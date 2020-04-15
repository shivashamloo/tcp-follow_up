/*  $Id: speerScorer.cpp,v 1.10 2010/10/13 21:30:51 lanczyck Exp $  */

#include <ncbi_pch.hpp>
#include <math.h>
#include <corelib/ncbi_limits.h>
#include "speerData.hpp"
#include "speerUtils.hpp"
#include "speerScorer.hpp"
#include "zScoreGenerator.hpp"

const TSpeerScore CSpeerScorer_Base::UNCOMPUTED = kMax_Double;

ostream& operator <<(ostream& os, const CSpeerScorer_Base& scorer)
{
    unsigned int col;
    double score, zscore, pscore;
    const TColumnScoreMap& colScores = scorer.GetColumnScores();
    const TColumnScoreMap& colZScores = scorer.GetColumnZScores();
    const TColumnScoreMap& colPScores = scorer.GetColumnPScores();

    TColumnScoreMapCit scoreCit, zscoreCit, pscoreCit;
    TColumnScoreMapCit scoreEnd = colScores.end(), zscoreEnd = colZScores.end(), pscoreEnd = colPScores.end();

    for (scoreCit = colScores.begin(); scoreCit != scoreEnd; ++scoreCit) {
        col = scoreCit->first;

        zscoreCit = colZScores.find(col);
        pscoreCit = colPScores.find(col);

        score = scoreCit->second;
        zscore = (zscoreCit != zscoreEnd) ? zscoreCit->second : CSpeerScorer_Base::UNCOMPUTED;
        pscore = (pscoreCit != pscoreEnd) ? pscoreCit->second : CSpeerScorer_Base::UNCOMPUTED;

        os << col << "    " << score << ",  Z = " << zscore << ", P = " << pscore << NcbiEndl;
    }

    return os;
    
}

void CSpeerScorer_Base::GetColumnScores(unsigned int column, double& score, double& zscore, double& pscore) const
{
    TColumnScoreMapCit scoreCit = m_colScores.find(column);
    TColumnScoreMapCit zscoreCit = m_colZScores.find(column); 
    TColumnScoreMapCit pscoreCit = m_colPScores.find(column);

    score = UNCOMPUTED;
    zscore = UNCOMPUTED;
    pscore = UNCOMPUTED;

    if (scoreCit != m_colScores.end()) {
        score = scoreCit->second;
    }
    if (zscoreCit != m_colZScores.end()) {
        zscore = zscoreCit->second;
    }
    if (pscoreCit != m_colPScores.end()) {
        pscore = pscoreCit->second;
    }
}

double CSpeerScorer_Base::GetColumnScore(unsigned int column) const
{
    double result = UNCOMPUTED;
    TColumnScoreMapCit cit = m_colScores.find(column);
    if (cit != m_colScores.end()) {
        result = cit->second;
    }
    return result;
}
double CSpeerScorer_Base::GetColumnZScore(unsigned int column) const
{
    double result = UNCOMPUTED;
    TColumnScoreMapCit cit = m_colZScores.find(column);
    if (cit != m_colZScores.end()) {
        result = cit->second;
    }
    return result;
}
double CSpeerScorer_Base::GetColumnPScore(unsigned int column) const
{
    double result = UNCOMPUTED;
    TColumnScoreMapCit cit = m_colPScores.find(column);
    if (cit != m_colPScores.end()) {
        result = cit->second;
    }
    return result;
}

double CSpeerScorer_Base::ComputeZPScores() 
{

    double result = 0.0;
    vector<double> scores, zscores, pscores;
    unsigned int i, nScores;
    TColumnScoreMapIt scoreIt, scoreEnd, zIt, zEnd;
    CZScoreGenerator zscorer;

    m_colZScores.clear();
    m_colPScores.clear();

    if (m_colScores.size() == 0 || m_score == UNCOMPUTED) {
        ComputeScore(true);
    }

    nScores = m_colScores.size();
    if (nScores > 0) {

        scoreEnd = m_colScores.end();
        for (scoreIt = m_colScores.begin(); scoreIt != scoreEnd; ++scoreIt) {
            scores.push_back(scoreIt->second);
        }

        zscorer.GetP(scores, pscores, &zscores);

        for (scoreIt = m_colScores.begin(), i = 0; scoreIt != scoreEnd; ++scoreIt, ++i) {
            m_colZScores[scoreIt->first] = zscores[i];
            m_colPScores[scoreIt->first] = pscores[i];
        }
    }


    //  Sum up the z-scores, as a sanity check (it should vanish).
    zEnd = m_colZScores.end();
    for (zIt = m_colZScores.begin(); zIt != zEnd; ++zIt) {
        result += zIt->second;
    }
    result /= m_colZScores.size();

    return result;
}

void CSpeerScorer_Base::SetGapFraction(double gapFraction) 
{
    if (gapFraction < 0.0) {
        m_gapFraction = 0.0;
    } else if (gapFraction > 1.0) {
        m_gapFraction = 1.0;
    } else {
        m_gapFraction = gapFraction;
    }
}

void CSpeerScorer_Base::FindNonGappyColumns(set<unsigned int>& nonGappyColumns)
{
    double gapFraction;
    unsigned int s, col, maxCol = kMax_UInt;
    unsigned int nAlignments = m_alignments.size();
    CSpeerAlignment::TColumnMapCit colCitS, colEndS;
    set<unsigned int> notFullyConservedColumnSet;
    set<unsigned int> gappyColumnSet;
    set<unsigned int>::iterator setIt, setEnd, fullyConsColEnd;

    nonGappyColumns.clear();

    //  Flag columns in which the gap percentage threshold is exceeded in any subfamily alignment.
    for (s = 0; s < nAlignments; ++s) {
        const CSpeerAlignment::TColumnMap& columnMapS = m_alignments[s]->GetColumnMap();
        if (columnMapS.size() < maxCol) maxCol = columnMapS.size();
        colEndS = columnMapS.end();

        for (colCitS = columnMapS.begin(); colCitS != colEndS; ++colCitS) {
            gapFraction = colCitS->second.GetGapFraction();
            if (gapFraction > m_gapFraction) {
                gappyColumnSet.insert(colCitS->first);
//                LOG_POST("    alignment " << s << ", col " << colCitS->first << ", high gapFraction = " << gapFraction << ":  " << colCitS->second.GetResidues());
            }
            if (! colCitS->second.IsSingleResidueColumn()) {
                notFullyConservedColumnSet.insert(colCitS->first);
            }
        }
    }

    //  Invert the set of gappy columns to record the set of non-gappy columns.
    setEnd = gappyColumnSet.end();
    fullyConsColEnd = notFullyConservedColumnSet.end();
    for (col = 0; col < maxCol; ++col) {
        if (gappyColumnSet.find(col) == setEnd) {

            //  do not include columns that are 100% conserved
            if (notFullyConservedColumnSet.find(col) != fullyConsColEnd) {
                nonGappyColumns.insert(col);
//            } else {
//                LOG_POST("    excluding fully conserved column " << col << " from results");
            }
        }
    }

}

////////////////////////////////////////////////////////////
//
//   Composite-Z scorer
//
////////////////////////////////////////////////////////////


CSpeerCompositeZScorer::CSpeerCompositeZScorer(const map< CSpeerScorer_Base*, double>& scorerMap) : CSpeerScorer_Base(), m_scorerMap(scorerMap)
{
    bool hasNullScorer = false;
    map< CSpeerScorer_Base*, double >::const_iterator cit = scorerMap.begin(), cend = scorerMap.end();
    for (; cit != cend; ++cit) {
        if (cit->first == NULL) {
            hasNullScorer = true;
            break;
        }
    }
    if (hasNullScorer) {
        m_scorerMap.clear();
    }

    m_scorerName = "Composite Z-score";
}

TSpeerScore CSpeerCompositeZScorer::ComputeScore(bool forceRecompute) 
{
    //  Stop if there was a problem with the input alignments.
    if (m_scorerMap.size() == 0) {
        return UNCOMPUTED;
    }

    //  Don't recompute if it's unnecessary.
    if (!forceRecompute && m_score != UNCOMPUTED) {
        return m_score;
    }

    TSpeerScore result = 0.0;
    double weight;
    TColumnScoreMapCit scoreMapCit, scoreMapCend;
    map< CSpeerScorer_Base*, double >::iterator scorerIt = m_scorerMap.begin(), scorerEnd = m_scorerMap.end();

    //  The composite score S is a linear combination of the Z-scores of the
    //  individual scorers held in m_scorerMap.  They are weighted by the
    //  value associated with the scorer in this map.

    //  Clean out preexisting data.
    m_colScores.clear();
    m_colZScores.clear();
    m_colPScores.clear();

    for (; scorerIt != scorerEnd; ++scorerIt) {
        const TColumnScoreMap& zScores = scorerIt->first->GetColumnZScores();
        weight = scorerIt->second;
        if (weight == 0.0) continue;

        scoreMapCend = zScores.end();
        for (scoreMapCit = zScores.begin(); scoreMapCit != scoreMapCend; ++scoreMapCit) {

            //  They should be, but ensure map elements are properly initialized to zero.
            if (m_colScores.find(scoreMapCit->first) == m_colScores.end()) {
                m_colScores[scoreMapCit->first] = 0;
            }

            m_colScores[scoreMapCit->first] += (scoreMapCit->second * weight);
        }
    }

    if (m_colScores.size() > 0) {
        //  Compute the average over the number of alignment combinations, 
        //  and then average over all columns used.
        scoreMapCend = m_colScores.end();
        for (scoreMapCit  = m_colScores.begin(); scoreMapCit != scoreMapCend; ++scoreMapCit) {
            result += scoreMapCit->second;
        }
        result /= m_colScores.size();
    }

    m_score = result;
    return result;
}



////////////////////////////////////////////////////////////
//
//   Relative-entropy scorer
//
////////////////////////////////////////////////////////////


CSpeerRelativeEntropyScorer::CSpeerRelativeEntropyScorer(vector< CSpeerAlignment*> alignments, double gapFraction) : CSpeerScorer_Base(alignments, gapFraction)
{
    bool hasNullAlignment = false;
    for (unsigned int i = 0; i < m_alignments.size(); ++i) {
        if (m_alignments[i] == NULL) {
            hasNullAlignment = true;
            break;
        }
    }
    if (hasNullAlignment) {
        m_alignments.clear();
    }

    m_scorerName = "Relative entropy";
}

TSpeerScore CSpeerRelativeEntropyScorer::ComputeScore(bool forceRecompute) 
{
    //  Stop if there was a problem with the input alignments.
    if (m_alignments.size() == 0) {
        return UNCOMPUTED;
    }

    //  Don't recompute if it's unnecessary.
    if (!forceRecompute && m_score != UNCOMPUTED) {
        return m_score;
    }

    TSpeerScore result = 0.0;
    unsigned int s, t;
    unsigned int col;
    unsigned int nAlignments = m_alignments.size();
    unsigned int nCombinations = nAlignments*(nAlignments - 1)/2;
    double contribution;

    TColumnScoreMapIt scoreMapIt, scoreMapEnd;
    CSpeerAlignment::TColumnMapCit colCitS, colCitT, colEndS, colEndT;
//    CSpeerColumn column;

    //  Map to hold pre-computed frequency maps for each alignment.
    map< unsigned int, TColFreqsMap > alignmentColFreqs;
    TColFreqsMapIt colFreqsItS, colFreqsEndS, colFreqsItT, colFreqsEndT;

    //  If computing, need to clean out preexisting data.
    m_colScores.clear();
    m_colZScores.clear();
    m_colPScores.clear();

/*
    //  Determine which columns exceed the gap percentage threshold in one
    //  or more of the input alignments.
    nCols = kMax_UInt;
    for (s = 0; s < nAlignments; ++s) {
        const CSpeerAlignment::TColumnMap& columnMapS = m_alignments[s]->GetColumnMap();
        colEndS = columnMapS.end();
        for (colCitS = columnMapS.begin(); colCitS != colEndS; ++colCitS) {
            gapFraction = colCitS->second.GetGapFraction();
            if (gapFraction > m_gapFraction) {
                nonGappyColumns.insert(colCitS->first);
                cerr << "    alignment " << s << ", col " << colCitS->first << ", high gapFraction = " << gapFraction << ":  " << colCitS->second.GetResidues() << NcbiEndl;
            }
        }

        //  In case alignments have different length (they shouldn't but...),
        //  the *smallest* alignment length is the maximum possible length of
        //  the intersection of all alignments.
        if (columnMapS.size() < nCols) nCols = columnMapS.size();
    }
*/
/*
    //  Initialize m_colScores to only have entries for columns in which each   
    //  alignment satisfies the gap percentage threshold.
    FindNonGappyColumns(nonGappyColumns);
    nonGappyColumnsEnd = nonGappyColumns.end();
    for (nonGappyColumnsIt = nonGappyColumns.begin(); nonGappyColumnsIt != nonGappyColumnsEnd; ++nonGappyColumnsIt) {
//    for (col = 0; col < nCols; ++col) {
//        if (nonGappyColumns.find(col) == nonGappyColumnsEnd) {
        m_colScores[col] = 0;
//            cerr << "    column " << col << " satisfies gap contraint\n";
//        }
    }
    scoreMapEnd = m_colScores.end();


    //  Pre-compute column frequencies for non-gappy columns.
    for (s = 0; s < nAlignments; ++s) {
        TColFreqsMap& colFreqsMap = m_alignmentColFreqs[s];
        for (scoreMapIt = m_colScores.begin(); scoreMapIt != scoreMapEnd; ++scoreMapIt) {
            col = scoreMapIt->first;
            if (m_alignments[s]->GetColumn(col, column)) {
                TAADataMap& colFreqs = colFreqsMap[col];
//                cerr << "s = " << s << "  col = " << col << "  residues = " << column.GetResidues() << NcbiEndl;
                GetColumnFrequencies(column, colFreqs);
            }
        }
    }
*/

    ComputeAlignmentColumnFrequencies(alignmentColFreqs);

    //  Loop over unique subgroup combinations (s,t).
    scoreMapEnd = m_colScores.end();
    for (s = 0; s < nAlignments; ++s) {
        TColFreqsMap& colFreqsMapS = alignmentColFreqs[s];
        colFreqsEndS = colFreqsMapS.end();

        for (t = s + 1; t < nAlignments; ++t) {
            TColFreqsMap& colFreqsMapT = alignmentColFreqs[t];
            colFreqsEndT = colFreqsMapT.end();

            //  For each validated column, compute the contribution from all amino acid
            //  types for the (s, t) pair of alignments.
            for (scoreMapIt  = m_colScores.begin(); scoreMapIt != scoreMapEnd; ++scoreMapIt) {
                col = scoreMapIt->first;
                colFreqsItS = colFreqsMapS.find(col);
                colFreqsItT = colFreqsMapT.find(col);
                if (colFreqsItS == colFreqsEndS || colFreqsItT == colFreqsEndT) continue;

                contribution = ComputeAlignmentPairContributionForColumn(colFreqsItS->second, colFreqsItT->second);
                scoreMapIt->second += contribution;
//                LOG_POST("   s " << s << "  t " << t << "  col " << col << " contrib " << contribution);
            }
        }

    }

    //  Compute the average over the number of alignment combinations, 
    //  and then average over all columns used.
    for (scoreMapIt  = m_colScores.begin(); scoreMapIt != scoreMapEnd; ++scoreMapIt) {
        scoreMapIt->second /= nCombinations;
        result += scoreMapIt->second;
    }
    result /= m_colScores.size();

    m_score = result;
    return result;
}

void CSpeerRelativeEntropyScorer::GetColumnFrequencies(const CSpeerColumn& column, TAADataMap& aaFreqs)
{
    unsigned int nGaps, nIndeterminate;
    unsigned int nRows = column.GetColumnLength();
    unsigned int nMaxGaps = (unsigned int) (nRows * m_gapFraction);
    double freq;

    nGaps = column.GetNumGaps();
    if (nGaps > nMaxGaps) {
        cout << "Skipping column " << column.GetColumnIndex() << " due to " << nGaps << " gaps (max allowed = " << nMaxGaps << ")\n";
        return;
    }

    TCountsMap aaCounts;
    TCountsMapIt aaIt, aaEnd;
    column.Counts(aaCounts, eSpeerConvertToUpper);

    nIndeterminate = aaCounts[AminoAcidChar[eAA_IndeterminateB]] + aaCounts[AminoAcidChar[eAA_IndeterminateZ]] + aaCounts[AminoAcidChar[eAA_IndeterminateX]];
    if (nIndeterminate > 0) {
        cout << "Warning:  skipped " << nIndeterminate << " ambiguous residues (B, Z, X) in column " << column.GetColumnIndex() << " in the RE computation.\n";
    }

    aaEnd = aaCounts.end();
    for (aaIt = aaCounts.begin(); aaIt != aaEnd; ++aaIt) {

        //  Completely skip the indeterminate residues.
        if (aaIt->first == AminoAcidChar[eAA_IndeterminateB]
            || aaIt->first == AminoAcidChar[eAA_IndeterminateZ]
            || aaIt->first == AminoAcidChar[eAA_IndeterminateX]) {
            continue;
        }

        //  Compensate for any indeterminate rows.
//        term = (double)aaIt->second + refFreqMap[aaIt->first] * 20.0;
//        freq = term/((double)(nRows - nIndeterminate + 20));
//        cerr << "    " << aaIt->second << ", " << refFreqMap[aaIt->first] << "  term, freq:  " << term << ", " << freq <<  NcbiEndl;
        freq = ((double)aaIt->second + refFreqMap[aaIt->first] * 20.0)/((double)(nRows - nIndeterminate + 20));
        aaFreqs[aaIt->first] = freq;
    }

}

void CSpeerRelativeEntropyScorer::ComputeAlignmentColumnFrequencies(map< unsigned int, TColFreqsMap >& alignmentColFreqs)
{
    unsigned int s, col;
    unsigned int nAlignments = m_alignments.size();
    TColumnScoreMapIt scoreMapIt, scoreMapEnd;
    TColFreqsMapIt colFreqsItS, colFreqsEndS, colFreqsItT, colFreqsEndT;
    set<unsigned int> nonGappyColumns;
    set<unsigned int>::iterator nonGappyColumnsIt, nonGappyColumnsEnd;
    CSpeerColumn column;

    alignmentColFreqs.clear();

    //  Initialize m_colScores to only have entries for columns in which each   
    //  alignment satisfies the gap percentage threshold.
    FindNonGappyColumns(nonGappyColumns);
    nonGappyColumnsEnd = nonGappyColumns.end();
    for (nonGappyColumnsIt = nonGappyColumns.begin(); nonGappyColumnsIt != nonGappyColumnsEnd; ++nonGappyColumnsIt) {
//    for (col = 0; col < nCols; ++col) {
//        if (nonGappyColumns.find(col) == nonGappyColumnsEnd) {
        col = *nonGappyColumnsIt;
        m_colScores[col] = 0;
//        LOG_POST("    column " << col << " satisfies gap contraint");
//        }
    }
    scoreMapEnd = m_colScores.end();


    //  Pre-compute column frequencies for non-gappy columns.
    for (s = 0; s < nAlignments; ++s) {
        TColFreqsMap& colFreqsMap = alignmentColFreqs[s];
        for (scoreMapIt = m_colScores.begin(); scoreMapIt != scoreMapEnd; ++scoreMapIt) {
            col = scoreMapIt->first;
            if (m_alignments[s]->GetColumn(col, column)) {
                TAADataMap& colFreqs = colFreqsMap[col];
//                LOG_POST("s = " << s << "  col = " << col << "  residues = " << column.GetResidues());
                GetColumnFrequencies(column, colFreqs);
            }
        }
    }
}



double CSpeerRelativeEntropyScorer::ComputeAlignmentPairContributionForColumn(const TAADataMap& column1Freqs, const TAADataMap& column2Freqs)
{
    static const double invLogBaseEOf2 = 1.0/log(2.0);
    static const string allowedCharacters = "ACDEFGHIKLMNPQRSTUVWY";  //  NOTE:  -, B, Z, X excluded!
    string::const_iterator allowedIt = allowedCharacters.begin(), allowedEnd = allowedCharacters.end();

    double result = 0;
    TAADataMapCit col1It, col2It;
    TAADataMapCit col1End = column1Freqs.end(), col2End = column1Freqs.end();

/*
    cerr << "col 1 freqs:\n";
    for (col1It = column1Freqs.begin(); col1It != col1End; ++col1It) {
        cerr << "    " << col1It->first << "  " << col1It->second << NcbiEndl;
    }
    cerr << "col 2 freqs:\n";
    for (col2It = column2Freqs.begin(); col2It != col2End; ++col2It) {
        cerr << "    " << col2It->first << "  " << col2It->second << NcbiEndl;
    }
*/
    for (; allowedIt != allowedEnd; ++allowedIt) {
        col1It = column1Freqs.find(*allowedIt);
        col2It = column2Freqs.find(*allowedIt);
        if (col1It != col1End && col2It != col2End && col2It->second != 0.0) {
//            term = (col1It->second - col2It->second)*log(col1It->second/col2It->second);
//            result += term;
            result += (col1It->second - col2It->second)*log(col1It->second/col2It->second);
        }
    }
    result *= invLogBaseEOf2;
    return result;
}

////////////////////////////////////////////////////////////
//
//   Evolutionary distance scorer
//
////////////////////////////////////////////////////////////

CSpeerEvolutionaryDistanceScorer::CSpeerEvolutionaryDistanceScorer(vector< CSpeerAlignment*> alignments, double gapFraction) : CSpeerScorer_Base(alignments, gapFraction)
{
    bool hasNullAlignment = false;
    for (unsigned int i = 0; i < m_alignments.size(); ++i) {
        if (m_alignments[i] == NULL) {
            hasNullAlignment = true;
            break;
        }
    }
    if (hasNullAlignment) {
        m_alignments.clear();
    }

    m_scorerName = "Evolutionary distance";
}

TSpeerScore CSpeerEvolutionaryDistanceScorer::ComputeScore(bool forceRecompute)  
{
    //  Stop if there was a problem with the input alignments.
    if (m_alignments.size() == 0) {
        return UNCOMPUTED;
    }

    //  Don't recompute if it's unnecessary.
    if (!forceRecompute && m_score != UNCOMPUTED) {
        return m_score;
    }

    TSpeerScore result = 0.0;
    TSpeerScore sum, denom;

    bool hasParent;
    unsigned int i, col;
    unsigned int nAlignments = m_alignments.size();

    set<unsigned int> nonGappyColumns;
    set<unsigned int>::iterator nonGappyColumnsIt, nonGappyColumnsEnd;
    map< unsigned int, TWeightMap > alignmentWeightMap;
    map< unsigned int, TColumnScoreMap > alIt, alEnd;

    TColumnScoreMapIt scoreMapIt, scoreMapEnd, colIt;
//    CSpeerColumn column;

    //  If computing, need to clean out preexisting data.
    m_colScores.clear();
    m_colZScores.clear();
    m_colPScores.clear();

    //  Initialize m_colScores to only have entries for columns in which each   
    //  alignment satisfies the gap percentage threshold.
    FindNonGappyColumns(nonGappyColumns);
    nonGappyColumnsEnd = nonGappyColumns.end();
    for (nonGappyColumnsIt = nonGappyColumns.begin(); nonGappyColumnsIt != nonGappyColumnsEnd; ++nonGappyColumnsIt) {
        col = *nonGappyColumnsIt;
        m_colScores[col] = 0;
//        LOG_POST("    column " << col << " satisfies gap contraint");
    }
    scoreMapEnd = m_colScores.end();


    //  Initialize the row weights for each subfamily alignment, along with the parent,
    //  given index = nAlignments.  Then compute the contributions to SED/GED for 
    //  each column of that alignment.  Assumes all subfamilies have the same parent
    //  alignment.
    TWeightMap baseWeights;
    hasParent = (m_alignments[0]->GetBaseAlignment() != NULL);
    if (!hasParent) {
        cerr << "Error:  no parent alignment detected.  Can't compute weights!!\n";
        return UNCOMPUTED;
    } else {
        //  Deal with the parent alignment
        m_colScoresForAlignment[nAlignments] = m_colScores;
        ComputeSequenceWeights(*(m_alignments[0]->GetBaseAlignment()), baseWeights, m_gapFraction);
        ComputeFamilyED(*(m_alignments[0]->GetBaseAlignment()), baseWeights, m_colScoresForAlignment[nAlignments]);
    }


    for (i = 0; i < nAlignments; ++i) {

        m_colScoresForAlignment[i] = m_colScores;
//        TWeightMap& weights = alignmentWeightMap[i];
//        ComputeSequenceWeights(*(m_alignments[i]), weights, m_gapFraction);
//        ComputeFamilyED(*(m_alignments[i]), weights, m_colScoresForAlignment[i]);
        ComputeFamilyED(*(m_alignments[i]), baseWeights, m_colScoresForAlignment[i]);

    }


    //  Sum over families and normalize by GED to compute the final ED score in each column.
    scoreMapIt = m_colScores.begin(), scoreMapEnd = m_colScores.end();
    for (; scoreMapIt != scoreMapEnd; ++scoreMapIt) {
        sum = 0.0;
        denom = 1.0;
        col = scoreMapIt->first;
        if (hasParent) {
            colIt = m_colScoresForAlignment[nAlignments].find(col);
            if (colIt != m_colScoresForAlignment[nAlignments].end() && colIt->second != 0) {
                denom = colIt->second;
            }
        }
        for (i = 0; i < nAlignments; ++i) {
            colIt = m_colScoresForAlignment[i].find(col);
            if (colIt != m_colScoresForAlignment[i].end()) {
                sum += colIt->second;
            }
        }    
        //LOG_POST("summing for " << col << "; SED = " << sum << ";  GED = " << denom << ";  SED/GED = " << sum/denom);
        sum /= denom;
        result += sum;
        scoreMapIt->second = sum;
    }

    m_score = result;
    return result;
}

void CSpeerEvolutionaryDistanceScorer::ComputeFamilyED(const CSpeerAlignment& alignment, const TWeightMap& weights, TColumnScoreMap& colScores)
{
    unsigned int i, j, col, len, nPairs;
    unsigned int fromRow = (alignment.GetBaseAlignment() != NULL) ? alignment.GetFromRow() : 0;
//    unsigned int toRow = (alignment.GetBaseAlignment() != NULL) ? alignment.GetToRow() : alignment.NumRows() - 1;
    string columnStr;
    double wi, wj;
    TSpeerScore sum, dij, contribution;
    TResiduePair resPair;
    const CSpeerColumn* column;
    TColumnScoreMapIt scoreMapIt = colScores.begin(), scoreMapEnd = colScores.end();
    TWeightMapCit weightCit, weightEnd = weights.end();
    string::iterator strIt1, strIt2, strEnd;

//    cerr << "from row " << fromRow << "; toRow " << toRow << endl;

    for (; scoreMapIt != scoreMapEnd; ++scoreMapIt) {
        sum = 0.0;
        col = scoreMapIt->first;
        column = alignment.GetColumn(col);
        if (column) {
            columnStr = column->GetResidues();
            len = columnStr.length();


            strEnd = columnStr.end();
            nPairs = 0;
            for (i = fromRow, strIt1 = columnStr.begin(); strIt1 != strEnd; ++i, ++strIt1) {
                weightCit = weights.find(i);
                wi = (weightCit != weightEnd) ? weightCit->second : 0.0;
                resPair.first = *strIt1;
                strIt2 = strIt1;
                ++strIt2;
                for (j = i+1; strIt2 != strEnd; ++j, ++strIt2) {
                    ++nPairs;
                    resPair.second = *strIt2;

                    //  No contribution from identical residue pairs
                    if (*strIt1 != *strIt2) {
                        weightCit = weights.find(j);
                        wj = (weightCit != weightEnd) ? weightCit->second : 0.0;
                        if (wi != 0 && wj != 0) {
                            dij = speerPairDistanceMap[resPair];
                            contribution = ((wi*wj)/(wi + wj)) * dij;
//                            if (col == 429) {
//                                LOG_POST("   " << col << "  (" << resPair.first << ", " << resPair.second << ") " << i << "  " << j << "  " << dij << "  " << contribution);
//                            }
                            sum += contribution;
                        }
//                    } else {
//                        if (col == 429) {
//                            LOG_POST("   " << col << " (" << resPair.first << ", " << resPair.second << ") " << i << "  " << j << "  0  0");
//                        }
                    }

                }
            }
            if (nPairs > 0) {
                sum /= (double) nPairs;
            } else {
                sum = 0;
            }
//            if (col == 429) {
//                LOG_POST("   " << col << ", " << columnStr << ":  nPairs = " << nPairs);
//            }
        }
        scoreMapIt->second = sum;
    }
}


////////////////////////////////////////////////////////////
//
//   Evolutionary rate scorer
//
////////////////////////////////////////////////////////////

const string CSpeerEvolutionRateScorer::m_rate4SiteExecutableName = "rate4site.exe";

CSpeerEvolutionRateScorer::CSpeerEvolutionRateScorer(vector< CSpeerAlignment*> alignments, double gapFraction, int jobid) : CSpeerScorer_Base(alignments, gapFraction), m_jobid(-1) 
{
    bool hasNullAlignment = false;
    for (unsigned int i = 0; i < m_alignments.size(); ++i) {
        if (m_alignments[i] == NULL) {
            hasNullAlignment = true;
            break;
        }
    }
    if (hasNullAlignment) {
        m_alignments.clear();
    }
    if (jobid >= 0) {
        m_jobid = jobid;
    }


    m_executablePath = ".";

    m_scorerName = "Evolutionary rate";
}


TSpeerScore CSpeerEvolutionRateScorer::ComputeScore(bool forceRecompute)  
{
    //  Stop if there was a problem with the input alignments.
    if (m_alignments.size() == 0) {
        return UNCOMPUTED;
    }

    //  Don't recompute if it's unnecessary.
    if (!forceRecompute && m_score != UNCOMPUTED) {
        return m_score;
    }

    TExitCode exitCode;
    TSpeerScore result = 0.0;
    unsigned int nAlignments = m_alignments.size();
    unsigned int i, col;
    string rate4SiteIn, rate4SiteOut, rate4SiteErr, fastaErr;

    set<unsigned int> nonGappyColumns;
    set<unsigned int>::iterator nonGappyColumnsIt, nonGappyColumnsEnd;
    TColumnScoreMapIt scoreMapIt, scoreMapEnd, colIt;
//    CSpeerColumn column;

    //  If computing, need to clean out preexisting data.
    m_colScores.clear();
    m_colZScores.clear();
    m_colPScores.clear();

    //  Initialize m_colScores to only have entries for columns in which each   
    //  alignment satisfies the gap percentage threshold.
    FindNonGappyColumns(nonGappyColumns);
    nonGappyColumnsEnd = nonGappyColumns.end();
    for (nonGappyColumnsIt = nonGappyColumns.begin(); nonGappyColumnsIt != nonGappyColumnsEnd; ++nonGappyColumnsIt) {
        col = *nonGappyColumnsIt;
        m_colScores[col] = 0;
//        LOG_POST("    column " << col << " satisfies gap contraint");
    }
    scoreMapEnd = m_colScores.end();

    for (i = 0; i < nAlignments; ++i) {

        //  make sure in/out filenames differ;
        //  implementation of GetTmpName on the PC allows for duplicate name
        //  when not enough time elapses between two calls.
        rate4SiteIn = CFile::GetTmpName() + "_in";
        if (m_jobid >= 0) {
            rate4SiteOut = NStr::IntToString(m_jobid) + ".res";
        } else {
            rate4SiteOut = CFile::GetTmpName() + "_out";  
        }

        //  Create tmp file containing alignment i as input, scoped to ensure flushing/closing of file.
        {{
            CNcbiOfstream ofs(rate4SiteIn.c_str());
            if (!m_alignments[i]->WriteAsFasta(ofs, fastaErr)) {
                ERR_POST("Could not write Fasta for subfamily " << i+1 << " to temporary file " << rate4SiteIn << ".\nStopping.\n");
                return UNCOMPUTED;
            } else {
                LOG_POST("Fasta for subfamily " << i+1 << " written to temporary file " << rate4SiteIn << " (it will be removed when SPEER terminates normally).");
            }
        }}
        
        LOG_POST("Running 'rate4Site' for subfamily " << i + i << " (please wait; this may take a while...)");
        exitCode = RunRate4Site(rate4SiteIn, rate4SiteOut, rate4SiteErr);
        if (exitCode == 0) {
            LOG_POST("    ... 'rate4Site' has finished for subfamily " << i + i << ".\n");
        } else {
            LOG_POST("    ... 'rate4Site' has returned with error code " << exitCode << " for subfamily " << i + i << ".\n    ***  There may have been a problem -> inspect any results carefully!  ***");
        }

        if (rate4SiteErr.length() > 0) {
            LOG_POST("    'rate4Site' message:\n    " << rate4SiteErr);
        }

        if (!CFile(rate4SiteOut).Exists()) {
            ERR_POST("No output file '" << rate4SiteOut << "' for subfamily " << i+1 << " exists.\nStopping.\n");
            return UNCOMPUTED;
        }

        CNcbiIfstream ifs(rate4SiteOut.c_str());
        ParseRate4SiteOutput(i, ifs);

        //  Clean up temporary files.
        CFile(rate4SiteIn).Remove();
        CFile(rate4SiteOut).Remove();

    }

    //  Clean up other auxiliary files created by rate4Site.
    if (m_jobid >= 0) {
        string jobIdStr = NStr::IntToString(m_jobid);
        CFile(jobIdStr + string(".tree")).Remove();
        CFile(jobIdStr + string(".res")).Remove();
        CFile(jobIdStr + string("_Org.res")).Remove();
    } else {
        CFile("TheTree.txt").Remove();
        CFile("r4sOrig.res").Remove();
    }
    CFile("r4s.res").Remove();

    //  Compute an average rate as the return value.
    if (m_colScores.size() > 0) {
        scoreMapIt = m_colScores.begin();
        for (; scoreMapIt != scoreMapEnd; ++scoreMapIt) {
            result += scoreMapIt->second;
        }
        result /= m_colScores.size();
    }
    m_score = result;
    return result;
}


TExitCode CSpeerEvolutionRateScorer::RunRate4Site(const string& inputFile, const string& outputFile, string& err)
{
    TExitCode result = -1;
    string command = m_executablePath + CDirEntry::GetPathSeparator() + m_rate4SiteExecutableName;
    string jobIdStr = NStr::IntToString(m_jobid);

    command += " -s " + inputFile + " -o " + outputFile;
    if (m_jobid >= 0) {
        command += " -x " + jobIdStr + ".tree";
        command += " -y " + jobIdStr + "_Org.res";
    }

    err.erase();
    try {
        NcbiCerr << "Running rate4site:  this may take a while..." << NcbiEndl;
        NcbiCerr << command << NcbiEndl;

        result = CExec::System(command.c_str());

        NcbiCerr << "    rate4site has completed." << NcbiEndl;


    } catch (CException& e) {
        err = "NCBI exception in rate4site:  " + e.GetMsg();
    } catch (std::exception& e) {
        err = "Standard exception in rate4site:  " + string(e.what());
    } catch (...) {
        err = "Unknown exception in rate4site.\n";
    }
    return result;
}

void CSpeerEvolutionRateScorer::ParseRate4SiteOutput(unsigned int alignmentIndex, CNcbiIstream& is)
{
    CSpeerAlignment::SRowData rowData;
    CStreamLineReader lineReader(is);
    unsigned int num;
    double rate;
    string::size_type charPos, lastCharPos;
    string::value_type c;
    string referenceSequence = kEmptyStr;
    CTempString line;
    TColumnScoreMapIt colScoreIt, colScoreEnd;

    if (!m_alignments[alignmentIndex]->GetRowData(0, rowData)) {
        ERR_POST("Could not get reference sequence from alignment " << alignmentIndex + 1 << " in ParseRate4SiteOutput.");
        return;
    }
    referenceSequence = rowData.seq;

    //  rate4Site returns data for each column in the alignment for which the reference sequence 
    //  is not a gap.  Some of those columns exceed the gap-percent threshold for Speer and need
    //  to be ignored.  Some columns with a gap on the reference sequence *do* satisfy the gap-percent
    //  threshold cutoff but we will not be able to obtain data in those cases.  
    //  [Although we could reorder the alignment rows to avoid this in the future.]

    charPos = referenceSequence.find_first_not_of('-', 0);
    lastCharPos = referenceSequence.find_last_not_of('-');

/*
    if (charPos != string::npos) {
        cerr << "first charPos = " << charPos << ", " << referenceSequence[charPos] << endl;
    } else {
        cerr << "no characters in ref seq!  ??" << endl;
    }

    if (lastCharPos != string::npos) {
        cerr << "last charPos = " << lastCharPos << ", " << referenceSequence[lastCharPos] << endl;
    }
*/

    colScoreEnd = m_colScores.end();
    while (!lineReader.AtEOF() && charPos <= lastCharPos) {

        //  Skip to next non-gap position in the sequence.
        if (referenceSequence[charPos]  == '-') {
            ++charPos;
            continue;
        }

        //  Skip header/footer info in the output file.
        line = NStr::TruncateSpaces(*++lineReader);
        if (line.empty() || line[0] == '#') {
            continue;
        }

        //  Find and increment column score if this character position satisfied gap fraction threshold.
        colScoreIt = m_colScores.find(charPos);
        if (colScoreIt != colScoreEnd) {

            //  If can't parse a non-header/footer, don't increment character position.
            if (!ParseRate4SiteOutputLine(line, num, c, rate)) {
                ERR_POST("Unexpected line format; line skipped:  '" << line << "'");
                continue;
            }

            colScoreIt->second += rate;
//            cerr << "compare to ref seq:  " << num << "/" << charPos << ";  " << c << " " << referenceSequence[charPos] << ";  " << rate << endl;
        }
        ++charPos;
    }
}

bool CSpeerEvolutionRateScorer::ParseRate4SiteOutputLine(CTempString& s, unsigned int& num, string::value_type& c, double& rate) 
{
    static const string delims = " \t\n";

    bool result = false;
    list<string> tokens;
    list<string>::iterator lit;

    try {
        NStr::Split(s, delims, tokens);
        if (tokens.size() > 2) {
            lit = tokens.begin();
            num = NStr::StringToUInt(*lit);
            c = (*++lit)[0];
            rate = NStr::StringToDouble(*++lit);
            result = true;
//            cerr << "line results:  " << s << "\n" << num << ";  " << c << "; " << rate << endl;
        }
    } catch (...) {
    }
    return result;
}
