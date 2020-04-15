#include <ncbi_pch.hpp>
#include <math.h>
#include "zScoreGenerator.hpp"
#include "zToPTable.hpp"

USING_NCBI_SCOPE;

void CZScoreGenerator::Initialize(const vector<double>& scores)
{
    unsigned int j;
    double sum = 0.0;

    m_mean = 0.0;
    m_stdDev = 1;
    m_numScoresInput = scores.size();
    if (m_numScoresInput == 0) return;

    //  Get mean for z-score calculation.
    for (j = 0; j < m_numScoresInput; ++j) {
        m_mean += scores[j];
    }    
    m_mean /= m_numScoresInput;

    //  Get standard deviation for z-score calculation.
    for (j = 0; j < m_numScoresInput; ++j) {
        sum += (scores[j] - m_mean)*(scores[j] - m_mean);
    }    

    if (m_numScoresInput > 1) {
        m_stdDev = sqrt(sum/(m_numScoresInput - 1));
    } else {
        m_stdDev = 0.0;  //  only one score in this case; one sample -> std dev vanishes/ill-defined
    }

//    cerr << "num scores = " << m_numScoresInput << ":  mean = " << m_mean << "; stdDev = " << m_stdDev << endl;
}

double CZScoreGenerator::GetZ(double score) const
{
    double result = 0.0;
    if (m_numScoresInput > 0 && m_stdDev != 0) {
        result = (score - m_mean)/m_stdDev;
    }
    return result;
}

void CZScoreGenerator::GetZ(const vector<double>& scores, vector<double>& zScores)
{
    unsigned int i, n = scores.size();
    double zscore;

    Initialize(scores);

    zScores.clear();
    for (i = 0; i < n; ++i) {
        zscore = GetZ(scores[i]);
        zScores.push_back(zscore);
    }
}

double CZScoreGenerator::GetP(double score, double* zScore) const
{
    double z = GetZ(score);

    if (zScore) {
        *zScore = z;
    }

    return CZ_To_P_Table::GetP(z);
}

void CZScoreGenerator::GetP(const vector<double>& scores, vector<double>& pScores, vector<double>* zScores)
{
    unsigned int i, n = scores.size();
    double pscore, zscore;

    Initialize(scores);

    pScores.clear();
    if (zScores) zScores->clear();

    for (i = 0; i < n; ++i) {
        pscore = GetP(scores[i], &zscore);
        pScores.push_back(pscore);
        if (zScores) zScores->push_back(zscore);
    }
}

