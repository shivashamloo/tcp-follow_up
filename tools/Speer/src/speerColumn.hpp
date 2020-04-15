/*  $Id: speerColumn.hpp,v 1.4 2010/10/13 21:30:51 lanczyck Exp $  */

#ifndef SPEERCOLUMN__HPP
#define SPEERCOLUMN__HPP

#include <map>
#include <string>
#include <vector>
#include <corelib/ncbistd.hpp>
#include <corelib/ncbi_limits.h>
#include "speerTypes.hpp"

USING_NCBI_SCOPE;

// class for representing the column of an alignment
class CSpeerColumn 
{

    friend ostream& operator <<(ostream& os, const CSpeerColumn& speerColumn);

public:

    CSpeerColumn() : m_columnIndex(kMax_UInt), m_columnResidues(kEmptyStr), m_gapFraction(0.0), m_numGaps(0) {};
    CSpeerColumn(unsigned int columnIndex, const string& columnResidues);

    CSpeerColumn(const CSpeerColumn& rhs) {
        m_columnIndex = rhs.m_columnIndex;
        m_columnResidues = rhs.m_columnResidues;
        m_gapFraction = rhs.m_gapFraction;
        m_numGaps = rhs.m_numGaps;
    }
    
    CSpeerColumn& operator=(const CSpeerColumn& rhs) {
        if (this != &rhs) {
            m_columnIndex = rhs.m_columnIndex;
            m_columnResidues = rhs.m_columnResidues;
            m_gapFraction = rhs.m_gapFraction;
            m_numGaps = rhs.m_numGaps;
        }
        return *this;
    }

    //  Create a column using a sub-string of the passed string of residues.
    CSpeerColumn(const CSpeerColumn& rhs, unsigned int from, unsigned int to);

    ~CSpeerColumn() {};

    unsigned int GetColumnIndex() const {return m_columnIndex;}
    unsigned int GetColumnLength() const {return m_columnResidues.length();}

    //  Return the number of times the character 'c' appears in m_columnResidues.
    //  Note:  this is case-sensitive!!
    unsigned int Count(char c) const;

    //  For convenience, gap content is pre-computed.
    unsigned int GetNumGaps() const { return m_numGaps; }
    double GetGapFraction() const { return m_gapFraction;}

    //  Return true if every residue in the column is the same
    // [faster than testing if GetPercentIdentity == 100].
    bool IsSingleResidueColumn() const;

    //  Compute a percent identity for the residues in the column:  100*#matches/#unique pairs.
    //  gap-gap pairs are considered an identical pair here.
    //  Return value is between 0.0 and 100.0.
    double GetPercentIdentity() const;

    //  Count the number of times the characters in m_columnResidues appear.
    //  Initializes counts to zero for all allowed characters based on
    //  the value of caseSensitivity given.  Note that the gap character
    //  and ambiguity residues B, Z and X are also counted individually.
    void Counts(TCountsMap& counts, ESpeerCaseSensitivity caseSensitivity) const;


    string GetResidues() {return m_columnResidues;}
    const string& GetResidues() const {return m_columnResidues;}

    //  Return an empty string if the from/to coordinates are out of range.
    //  Note:  from/to coordinates are zero-based.
    string GetResidues(unsigned int from, unsigned int to) const;

    //  Return the most common residue, along w/ the number of times it appeared.
    //  If a tie, return 'smallest' residue alphabetically.  If an error, a '-' is returned.
//    static char MostCommonResidue(const vector<char>& residues, unsigned int& appearances);

    string ToString() const;

private:

    unsigned int m_columnIndex;
    string m_columnResidues;
    double m_gapFraction;
    unsigned int m_numGaps;

    void SetGapFraction();
};

ostream& operator <<(ostream& os, const CSpeerColumn& speerColumn);

#endif // SPEER__HPP 
