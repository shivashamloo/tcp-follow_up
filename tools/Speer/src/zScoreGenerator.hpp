#ifndef Z_SCORE_GENERATOR__HPP
#define Z_SCORE_GENERATOR__HPP

#include <corelib/ncbistd.hpp>

USING_NCBI_SCOPE;

class CZScoreGenerator {

public:

    //  Initialize with a dummy vector of scores; requires scores to be passed the first time used
    //  to generate a Z or P value.
    CZScoreGenerator() {
        vector<double> scores;
        Initialize(scores);
    }

    CZScoreGenerator(const vector<double>& scores) {
        Initialize(scores);
    }

    //  Uses existing mean, deviation to compute Z for score.
    //  If object was not properly initialized, return 0.
    double GetZ(double score) const;

    //  Uses existing mean, deviation to compute p-value for score.  Optionally returns the Z-score.
    //  If object was not properly initialized, return 0.
    double GetP(double score, double* zScore = NULL) const;

    //  Initializes object with the first argument, then computes Z/p for each of those values.
    void GetZ(const vector<double>& scores, vector<double>& zScores);
    void GetP(const vector<double>& scores, vector<double>& pScores, vector<double>* zScores = NULL);

    double GetMean() const {return m_mean;}
    double GetStdDev() const {return m_stdDev;}

private:

    unsigned int m_numScoresInput;
    double m_mean;
    double m_stdDev;

    void Initialize(const vector<double>& scores);

};

#endif
