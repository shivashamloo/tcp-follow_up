/*  $Id: speerScore.cpp,v 1.10 2010/10/13 21:30:51 lanczyck Exp $  */

#include <ncbi_pch.hpp>
#include "speerScore.hpp"

CSpeerScore::~CSpeerScore() 
{
    TScorersIt it = m_scorers.begin(), itEnd = m_scorers.end();
    for (; it != itEnd; ++it) {
        delete it->first;
    }
    m_scorers.clear();
}

void CSpeerScore::AddScorer(CSpeerScorer_Base* scorer, double weight) 
{
    if (scorer) {
        m_scorers.insert(TScorersVT(scorer, weight));
    }
}


void CSpeerScore::AddScorer(ESpeerScorerType type, vector< CSpeerAlignment* >& alignments, double gapFractionThreshold, double weight, string pathToExecutable)
{
    CSpeerScorer_Base* scorer = NULL;
    CSpeerEvolutionRateScorer* erScorer = NULL;

    if (alignments.size() > 0) {
        switch (type) {
        case eRelativeEntropyScorer:
            scorer = new CSpeerRelativeEntropyScorer(alignments, gapFractionThreshold);
            break;
        case eEvolutionaryDistanceScorer:
            scorer = new CSpeerEvolutionaryDistanceScorer(alignments, gapFractionThreshold);
            break;
        case eEvolutionaryRateScorer:
            erScorer = new CSpeerEvolutionRateScorer(alignments, gapFractionThreshold, m_jobid);
            if (erScorer && pathToExecutable.length() > 0) {
                erScorer->SetExecutablePath(pathToExecutable);
            }
            scorer = erScorer;                
            break;
        default:
            break;
        }
    }

    if (scorer) {
        m_scorers.insert(TScorersVT(scorer, weight));
    }
}


CSpeerScorer_Base* CSpeerScore::Compute() 
{
    CSpeerScorer_Base* result = NULL;
    TSpeerScore combinedScalarResult = 0.0, scalarResult = 0.0;
    TScorersIt it = m_scorers.begin(), itEnd = m_scorers.end();

    if (m_scorers.size() > 0) {

        result = new CSpeerCompositeZScorer(m_scorers);
        if (!result) return NULL;

        //  Compute results for each individual scorer.
        for (; it != itEnd; ++it) {
            if (it->first) {
                combinedScalarResult += (it->first->ComputeZPScores()) * it->second;
            }
//            LOG_POST(it->first->GetScorerName() << " scorer results:\n" << *(it->first));
        }

        //  Trigger computation of the final, composite scores.
        scalarResult = result->ComputeZPScores();
//        LOG_POST(result->GetScorerName() << " scorer results:\n" << *result);

//        LOG_POST("Aggregate result (sanity check:  should be ~0.0):  " << scalarResult);
//        LOG_POST("   combined individual results:  " << combinedScalarResult);
    }

    return result;
}
