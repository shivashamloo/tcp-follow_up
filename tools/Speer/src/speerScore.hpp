/*  $Id: speerScore.hpp,v 1.5 2010/10/13 21:30:51 lanczyck Exp $  */

#ifndef SPEER_SCORE__HPP
#define SPEER_SCORE__HPP

#include <map>
#include <string>
#include <vector>
#include <corelib/ncbistd.hpp>
#include "speerData.hpp"
#include "speerScorer.hpp"
#include "speerColumn.hpp"
#include "speerAlignment.hpp"

USING_NCBI_SCOPE;
//USING_SCOPE(objects);
//USING_SCOPE(cd_utils);

class CSpeerScore {

public:

    CSpeerScore() : m_jobid (-1) {};
    ~CSpeerScore();

    //  set 'jobid' to a non-negative integer when running SPEER on multiple processors to ensure
    //  intermediate filenames are unique across all active processes; 
    //  a negative job id will generate default filenames instead of ones constructed from the jobid.
    void SetJobid(int jobid) { m_jobid = jobid; }
    int GetJobid() const { return m_jobid; }

    //  Negative (or zero) values for weight are allowed.  Null scorers will be ignored.
    void AddScorer(CSpeerScorer_Base* scorer, double weight = 1.0);

    //  Creates the indicated scorer type and adds it to the set of scorers.
    //  Negative (or zero) values for weight are allowed.  
    void AddScorer(ESpeerScorerType type, vector< CSpeerAlignment* >& alignments, double gapFractionThreshold, double weight = 1.0, string pathToExecutable = kEmptyStr);

    //  Score for each column is computed as sum_on_scorers(score*weight).
    //  The results are stored in the CSpeerScorer_Base subclass returned.
    //  Client assumes ownership of the returned object.
    CSpeerScorer_Base* Compute();

private:

    //  Map between the scorer and the weight used when computing the total score.
    typedef map< CSpeerScorer_Base*, double> TScorers;
    typedef TScorers::iterator TScorersIt;
    typedef TScorers::const_iterator TScorersCit;
    typedef TScorers::value_type TScorersVT;

    TScorers m_scorers;
    int m_jobid;
};

#endif // SPEER_SCORE__HPP
