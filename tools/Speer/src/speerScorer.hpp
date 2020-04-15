/*  $Id: speerScorer.hpp,v 1.8 2010/10/13 21:30:51 lanczyck Exp $  */

#ifndef SPEER_SCORER__HPP
#define SPEER_SCORER__HPP

#include <map>
#include <string>
#include <vector>
#include <corelib/ncbistd.hpp>
#include <corelib/ncbiexec.hpp>
#include <corelib/tempstr.hpp>

#include "speerAlignment.hpp"

USING_NCBI_SCOPE;

enum ESpeerScorerType {
    eRelativeEntropyScorer = 0,
    eEvolutionaryDistanceScorer,
    eEvolutionaryRateScorer
};

class CSpeerScorer_Base {

public:

    static const TSpeerScore UNCOMPUTED;

    CSpeerScorer_Base() : m_score(UNCOMPUTED), m_gapFraction(0.0) {};
    virtual ~CSpeerScorer_Base() {};

    TSpeerScore GetScore() {return m_score;}

    //  If 'forceRecompute' is false, simply return m_score if it has already been computed.
    //  If 'forceRecompute' is true, recompute the score from scratch and reset m_score.
    virtual TSpeerScore ComputeScore(bool forceRecompute) = 0;

    //  Compute the Z and P scores for the data in m_colScores.
    //  As a sanity check, returns the average of the Z-scores for all columns in m_colScores, 
    //  which should vanish modulo rounding errors.  Calls ComputeScore first only if necessary.
    double ComputeZPScores();

    const TColumnScoreMap& GetColumnScores() const { return m_colScores; }
    const TColumnScoreMap& GetColumnZScores() const { return m_colZScores; }
    const TColumnScoreMap& GetColumnPScores() const { return m_colPScores; }

    //  These return 'UNCOMPUTED' if the column is invalid.
    double GetColumnScore(unsigned int column) const; 
    double GetColumnZScore(unsigned int column) const;
    double GetColumnPScore(unsigned int column) const;

    //  If a score is not available for the given column, it is reported 'UNCOMPUTED'
    void GetColumnScores(unsigned int column, double& score, double& zscore, double& pscore) const;

    unsigned int NumAlignments() const { return m_alignments.size(); }

    void SetGapFraction(double gapFraction);

    string GetScorerName() const {return m_scorerName;}

protected:

    //  a map of column index to the score/z-score computed for that column.
    TColumnScoreMap m_colScores;
    TColumnScoreMap m_colZScores;
    TColumnScoreMap m_colPScores;

    TSpeerScore m_score;

    vector< CSpeerAlignment* > m_alignments;
    double m_gapFraction;  //  the maximum gap fraction allowed for a column to be scored

    string m_scorerName;

    CSpeerScorer_Base(vector< CSpeerAlignment* > alignments, double gapFraction) : m_score(UNCOMPUTED), m_alignments(alignments) {
        SetGapFraction(gapFraction);
    }

    //  For the set of alignments in m_alignments, return the set of columns for
    //  which none of the alignments exceeds the maximum gap fraction.
    void FindNonGappyColumns(set<unsigned int>& nonGappyColumns);

};


ostream& operator <<(ostream& os, const CSpeerScorer_Base& scorer);



//  This subclass exists to generate a combination of the Z-scores generated 
//  by other subclass instances at each alignment column.
class CSpeerCompositeZScorer : public CSpeerScorer_Base {

public:

    //  Generate a composite score from the individual scorers provided.
    CSpeerCompositeZScorer(const map< CSpeerScorer_Base*, double>& scorerMap);
    virtual ~CSpeerCompositeZScorer() {};

    virtual TSpeerScore ComputeScore(bool forceRecompute);

private:

    map< CSpeerScorer_Base*, double> m_scorerMap;
};



class CSpeerRelativeEntropyScorer : public CSpeerScorer_Base {

    typedef map<unsigned int, TAADataMap> TColFreqsMap;
    typedef map<unsigned int, TAADataMap>::iterator TColFreqsMapIt;

public:

    //  If any of the CSpeerAlignment pointers are null, m_alignments will remain empty.
    //  NumAlignments() can be used to test if construction was successful.
    CSpeerRelativeEntropyScorer(vector< CSpeerAlignment*> alignments, double gapFraction);
    virtual ~CSpeerRelativeEntropyScorer() {};

    virtual TSpeerScore ComputeScore(bool forceRecompute);

private:

    void GetColumnFrequencies(const CSpeerColumn& column, TAADataMap& aaFreqs);

    //  Returns sum(aa)[(f1 - f2)*log_2(f1/f2)], where log_2(x) = log(x)/log(2).
    //  Only sums over residues represented as *uppercase* characters.
    double ComputeAlignmentPairContributionForColumn(const TAADataMap& column1Freqs, const TAADataMap& column2Freqs);

    void ComputeAlignmentColumnFrequencies(map< unsigned int, TColFreqsMap >& alignmentColFreqs);
};



class CSpeerEvolutionaryDistanceScorer : public CSpeerScorer_Base {

public:

    CSpeerEvolutionaryDistanceScorer(vector< CSpeerAlignment*> alignments, double gapFraction);
    virtual ~CSpeerEvolutionaryDistanceScorer() {};

    virtual TSpeerScore ComputeScore(bool forceRecompute);

private:

    //  For each alignment, save the contribution to SED/GED.  The GED contribution
    //  has index equal to m-alignments.size().
    map< unsigned int, TColumnScoreMap > m_colScoresForAlignment;

    //  If the alignment is a subfamily, this returns a contribution to SED [Eq. 3],
    //  otherwise, this returns the value of GED [Eq. 4] for each column defined in colScores.
    void ComputeFamilyED(const CSpeerAlignment& alignment, const TWeightMap& weights, TColumnScoreMap& colScores);

};


class CSpeerEvolutionRateScorer : public CSpeerScorer_Base {

    static const string m_rate4SiteExecutableName;

public:

    //  Non-negative jobid will be used to constructed non-default intermediate filenames, otherwise
    //  default filenames will be used.  If running in parallel, provide job id to avoid file contention
    //  issues of multiple processes writing to the same default filename.
    CSpeerEvolutionRateScorer(vector< CSpeerAlignment*> alignments, double gapFraction, int jobid = -1);
    virtual ~CSpeerEvolutionRateScorer() {};

    virtual TSpeerScore ComputeScore(bool forceRecompute);

    void SetExecutablePath(const string& path) { m_executablePath = path; }
    string GetExecutablePath() const { return m_executablePath; }

private:

    string m_executablePath;
    int m_jobid;

    TExitCode RunRate4Site(const string& fileIn, const string& fileOut, string& err);
    void ParseRate4SiteOutput(unsigned int alignmentIndex, CNcbiIstream& is);
    bool ParseRate4SiteOutputLine(CTempString& s, unsigned int& num, string::value_type& c, double& rate);
};

#endif // SPEER_SCORER__HPP
