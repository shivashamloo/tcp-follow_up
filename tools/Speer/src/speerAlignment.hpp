/*  $Id: speerAlignment.hpp,v 1.5 2008/10/10 16:21:44 lanczyck Exp $  */

#ifndef SPEERALIGNMENT__HPP
#define SPEERALIGNMENT__HPP

#include <map>
#include <string>
#include <vector>
#include <corelib/ncbistd.hpp>
#include <corelib/ncbi_limits.hpp>
#include <algo/structure/cd_utils/cuReadFastaWrapper.hpp>
#include <algo/structure/cd_utils/cuSeqAnnotFromFasta.hpp>
#include <algo/structure/cd_utils/cuCdFromFasta.hpp>
#include "speerColumn.hpp"

USING_NCBI_SCOPE;
USING_SCOPE(objects);
USING_SCOPE(cd_utils);

// class for representing an alignment for SPEER calculations
class CSpeerAlignment
{
    friend ostream& operator <<(ostream& os, const CSpeerAlignment& speerAlignment);

public:

    typedef map<unsigned int, CSpeerColumn> TColumnMap;
    typedef TColumnMap::iterator TColumnMapIt;
    typedef TColumnMap::const_iterator TColumnMapCit;
    typedef TColumnMap::value_type TColumnMapVT;

    struct SRowData {
        unsigned int row;
        string defline;
        string seq;
    };


    CSpeerAlignment(CCdFromFasta* cd, bool gaplessColumnsOnly, ESpeerCaseSensitivity caseSensitivity = eSpeerConvertToUpper);

    //  Create a subfamily using a sub-string of the passed string of residues.
    CSpeerAlignment(const CSpeerAlignment& rhs, unsigned int rowFrom, unsigned int rowTo, ESpeerCaseSensitivity caseSensitivity = eSpeerConvertToUpper);

    ~CSpeerAlignment() {};

    //  Create a subfamily using a sub-string of the passed string of residues.
    CSpeerAlignment* MakeSubfamily(unsigned int rowFrom, unsigned int rowTo, ESpeerCaseSensitivity caseSensitivity = eSpeerConvertToUpper) const;

    //  Lazy-initializes m_rowMap when necessary.
    bool GetRowData(unsigned int rowIndex, SRowData& rowData);

    bool GetColumn(unsigned int index, CSpeerColumn& column) const;
    const CSpeerColumn* GetColumn(unsigned int index) const;

    const TColumnMap& GetColumnMap() const { return m_columns; }

    bool IsSubfamily() const { return (m_baseAlignment != NULL);}

    unsigned int NumColumns() const { return m_columns.size(); }
    unsigned int NumRows() const { return m_rowTo - m_rowFrom + 1; }
    unsigned int GetFromRow() const { return m_rowFrom; }
    unsigned int GetToRow() const { return m_rowTo; }

    ESpeerCaseSensitivity GetCaseSensitivity() const { return m_caseSensitivity; }

    bool WriteAsClustal(ostream& s, string& err);
    bool WriteAsFasta(ostream& s, string& err);

    //  Return the most common residue, along w/ the number of times it appeared.
    //  If a tie, return 'smallest' residue alphabetically.  If an error, a '-' is returned.
//    static char MostCommonResidue(const vector<char>& residues, unsigned int& appearances);

    const CSpeerAlignment* GetBaseAlignment() const {return m_baseAlignment;}

private:

    const CSpeerAlignment* m_baseAlignment;   //  this is NULL if this object is not a subfamily
    ESpeerCaseSensitivity m_caseSensitivity;  //  how to treat the characters in the alignment
    unsigned int m_rowFrom;                   //  used for subfamilies; row m_rowFrom in baseAlignment is row 0 in this alignment
    unsigned int m_rowTo;                     //  used for subfamilies; row m_rowTo in baseAlignment is last row in this alignment
    TColumnMap m_columns;

    map<unsigned int, SRowData> m_rowMap;     //  Sequence data is lazy-initialized to the full sequence for each alignment row.
                                              //  If alignment was created with gapless columns only, only those columns appear
                                              //  in the sequence strings for each row.  Hence, these sequences may differ from
                                              //  that in the original CD.  Deflines are unmodified from the original FASTA.


    bool PopulateRowStrings();                //  only call when necessary; clears existing sequences in m_rowMap on each invokation
};


ostream& operator <<(ostream& os, const CSpeerAlignment& speerAlignment);


#endif // SPEERALIGNMENT__HPP 
