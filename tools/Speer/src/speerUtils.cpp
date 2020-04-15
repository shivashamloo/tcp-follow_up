/*  $Id: speerUtils.cpp,v 1.7 2008/11/04 00:12:47 lanczyck Exp $  */

#include <ncbi_pch.hpp>
#include "speerUtils.hpp"
#include "speerColumn.hpp"

USING_NCBI_SCOPE;

unsigned int CountResidueTypes(TCountsMap& counts)
{
    static const string allowedCharactersL = "acdefghiklmnpqrstuvwy";
    static const string allowedCharactersU = "ACDEFGHIKLMNPQRSTUVWY";

    unsigned int result = 0;
    unsigned int i = 0, nLetters = allowedCharactersU.length();

    for (; i < nLetters; ++i) {
        if (counts[allowedCharactersL[i]] > 0 || counts[allowedCharactersU[i]]) {
            ++result;
        }
    }
    return result;
}

void ComputeSequenceWeights(const CSpeerAlignment& alignment, TWeightMap& weights, double maxGapFraction)
{
    unsigned int col, row;
    unsigned int nGaps, nGoodColumns, nPairs, nResidueTypes, nMaxGaps;
    unsigned int nRows = alignment.NumRows();
    unsigned int nCols = alignment.NumColumns();
    unsigned char ucharLower, ucharUpper;
    ESpeerCaseSensitivity caseSensitivity = alignment.GetCaseSensitivity();
    double invResiduesTypes, gapFraction;
    string residues;
    string::iterator resIt, resEnd;
    TWeightMapIt weightIt, weightEnd;
    TCountsMap countMap;
    const CSpeerColumn* column;

    if (maxGapFraction > 1.0) {
        ERR_POST("Invalid maxGapFraction = " << maxGapFraction << "; reset it to 1.0 (any number allowed).");
        maxGapFraction = 1.0;
    } else if (maxGapFraction < 0.0) {
        ERR_POST("Invalid maxGapFraction = " << maxGapFraction << "; reset it to 0.0 (no gaps allowed).");
        maxGapFraction = 0.0;
    }

    nMaxGaps = (unsigned int) (nRows * maxGapFraction);
    //LOG_POST("For " << nRows << " rows, skip columns with more than " << nMaxGaps << " gaps.");

    //  Initialize the weights map.
    for (row = 0; row < nRows; ++row) {
        weights[row] = 0.0;
    }
    
    //  Loop over alignment columns.
    for (col = 0, nGoodColumns = 0; col < nCols; ++col) {

        column = alignment.GetColumn(col);
        if (!column) continue;

        nGaps = column->GetNumGaps();
        gapFraction = column->GetGapFraction();
        if (gapFraction > maxGapFraction) {
//            cout << "Skipping column " << col << ":  " << nGaps << " gaps exceeds the maximum allowed (" << nMaxGaps << ")\n";
            continue;
        }
        ++nGoodColumns;

        column->Counts(countMap, caseSensitivity);

        //  Count the different types of *real* residues in the column.
        //  (i.e., skipping -, B, Z, X)
        nResidueTypes = CountResidueTypes(countMap);
        invResiduesTypes = (nResidueTypes > 0) ? 1.0/(double)(nResidueTypes) : 1.0;
//        cout << "Including column " << col << " with " << nGaps << " gaps in the weight   # types = " << nResidueTypes << NcbiEndl;

        //  Loop over residues in the column, in row order
        residues = column->GetResidues();
        if (caseSensitivity == eSpeerConvertToUpper) {
            NStr::ToUpper(residues);
        } else if (caseSensitivity == eSpeerConvertToLower) {
            NStr::ToLower(residues);
        }
        resEnd = residues.end();
        for (row = 0, resIt = residues.begin(); resIt != resEnd; ++row, ++resIt) {
            if (*resIt == '-') {
                nPairs = countMap[*resIt];
            } else {
                //  For the weights, treat 'a' and 'A' together even if the alignment is
                //  treated case sensitively otherwise.
                if (caseSensitivity == eSpeerCaseSensitive) {
                    ucharLower = tolower((unsigned char) *resIt);
                    ucharUpper = toupper((unsigned char) *resIt);
                    nPairs = countMap[ucharLower] + countMap[ucharUpper];
                } else {
                    nPairs = countMap[*resIt];
                }
            }

            //  nPairs should always >= 1 due to self-pairing, but just in case...
            if (nPairs > 0) {
                weights[row] += invResiduesTypes/(double) nPairs;
            }
        }
    }

    //  Normalize the weights by the number of columns that passed the gap % screen.
    if (nGoodColumns > 0) {
        weightEnd = weights.end();
        for (weightIt = weights.begin(); weightIt != weightEnd; ++weightIt) {
//            cout << "Row " << weightIt->first << " " << weightIt->second << " " << nGoodColumns << " ";
            weightIt->second /= nGoodColumns;
//            cout << weightIt->second << endl;
        }
    } else {
        weights.clear();
        ERR_POST("Error:  no columns found that have less than " << 100.0*maxGapFraction << "% gap characters!");
    }
}

