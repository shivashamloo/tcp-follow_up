#include <ncbi_pch.hpp>
#include <corelib/ncbistd.hpp>
#include <corelib/ncbi_limits.hpp>
#include <corelib/ncbistr.hpp>
#include "speerColumn.hpp"

ostream& operator<<(ostream& os, const CSpeerColumn& speerColumn)
{
    os << speerColumn.ToString() << NcbiEndl;
    return os;
}

CSpeerColumn::CSpeerColumn(unsigned int columnIndex, const string& columnResidues) : m_columnIndex(columnIndex), m_columnResidues(columnResidues), m_gapFraction(0.0), m_numGaps(0) 
{
    SetGapFraction();
}


CSpeerColumn::CSpeerColumn(const CSpeerColumn& rhs, unsigned int from, unsigned int to)
{
    m_columnIndex = rhs.m_columnIndex;
    m_columnResidues = rhs.GetResidues(from, to);
    SetGapFraction();
}

void CSpeerColumn::SetGapFraction()
{
    m_numGaps = 0;
    m_gapFraction = 0.0;
    if (m_columnResidues.length() > 0) {
        m_numGaps = Count('-');
        m_gapFraction = ((double) m_numGaps)/((double) m_columnResidues.length());
    }
}

string CSpeerColumn::ToString() const
{
    static const string spacer(":  ");
    string s = NStr::UIntToString(m_columnIndex) + spacer + m_columnResidues;
    double gapPerc = 100.0*m_gapFraction;
    if (m_numGaps > 0) {
        s += "\ncolumn has " + NStr::UIntToString(m_numGaps) + " gaps (" +  NStr::DoubleToString(gapPerc) + "%)";
//    } else {
//        s += "  (no gaps)";
    }
    return s;
}

string CSpeerColumn::GetResidues(unsigned int from, unsigned int to) const
{
    string result = kEmptyStr;
    if (from <= to && to < m_columnResidues.length()) {
        result = m_columnResidues.substr(from, to - from + 1);
    }
    return result;
}

unsigned int CSpeerColumn::Count(char c) const
{
    unsigned int n = count(m_columnResidues.begin(), m_columnResidues.end(), c);
    return n;
}

void CSpeerColumn::Counts(TCountsMap& counts, ESpeerCaseSensitivity caseSensitivity) const
{
    static const string allCharacters = "-abcdefghiklmnpqrstuvwxyzABCDEFGHIKLMNPQRSTUVWXYZ";
    static const string lowerCharacters = "-abcdefghiklmnpqrstuvwxyz";
    static const string upperCharacters = "-ABCDEFGHIKLMNPQRSTUVWXYZ";

    string allowedCharacters;
    string columnResidues = m_columnResidues;
    switch (caseSensitivity) {
    case eSpeerConvertToUpper:
        allowedCharacters = upperCharacters;
        NStr::ToUpper(columnResidues);
        break;
    case eSpeerConvertToLower:
        allowedCharacters = lowerCharacters;
        NStr::ToLower(columnResidues);
        break;
    case eSpeerCaseSensitive:
    default:
        allowedCharacters = allCharacters;
        break;
    }

    string::value_type c;
    string::const_iterator allowedIt = allowedCharacters.begin(), allowedEnd = allowedCharacters.end();
    string::const_iterator cit = columnResidues.begin(), cend = columnResidues.end();

    //  Initialize the map for each letter.
    counts.clear();
    for (; allowedIt != allowedEnd; ++allowedIt) counts[*allowedIt] = 0;
    
    for (; cit != cend; ++cit) {
        c = *cit;
        ++counts[c];
    }
}

double CSpeerColumn::GetPercentIdentity() const
{
    double result = 0.0;
    double nCombinations;
    double nMatches = 0.0;
    string::const_iterator strCit2;
    string::const_iterator strCit1 = m_columnResidues.begin();
    string::const_iterator strEnd = m_columnResidues.end();

    nCombinations = m_columnResidues.size() * (m_columnResidues.size() - 1.0)/2.0;

    for (; strCit1 != strEnd; ++strCit1) {
        strCit2 = strCit1;
        ++strCit2;
        for (; strCit2 != strEnd; ++strCit2) {
            if (*strCit1 == *strCit2) {
                nMatches += 1.0;
            }
        }
    }

    if (nCombinations > 0) {
        result = 100.0*nMatches/nCombinations;
    }
    return result;
}

bool CSpeerColumn::IsSingleResidueColumn() const
{
    bool result = true;
    unsigned int len = m_columnResidues.size();
    string::const_iterator strIt;
    string::const_iterator strFirst = m_columnResidues.begin();
    string::const_iterator strEnd = m_columnResidues.end();

    if (len > 0) {
        if (len > 1) {
            strIt = strFirst;
            ++strIt;
            while (result && strIt != strEnd) {
                result = (*strIt == *strFirst);
                ++strIt;
            }
        }
    } else {
        result = false;
    }
    return result;
}
