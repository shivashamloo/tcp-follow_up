
#include <ncbi_pch.hpp>
#include <corelib/ncbi_limits.hpp>
#include "speerAlignment.hpp"

//  Pad the string returned with whitespace out to 'padToLength', where there must be at least
//  one space at the end.  Hence, the actual defline has a max of padToLength - 1 characters.
string SpeerGetIDFromDefline(const string& defline, unsigned int rowId, unsigned int padToLength = 24)
{
    string result(padToLength, ' ');
    string defaultStr = "UnknownID-Row" + NStr::UIntToString(rowId + 1);
    string tmp = defaultStr;
    unsigned int len = defline.length();
    string::size_type firstPos = defline.find('>', 0);
    string::size_type spacePos, tabPos, pipePos;
// 
//    string::size_type firstSpace = defline.find(' ', 0);
//    string::size_type secondPipe = string::npos;

    if (firstPos == string::npos) {  //  no '>' character; start immediately
        firstPos = 0;
    } else if (len - firstPos > 1) { //  start after the first '>' character (when > 1 character)
        ++firstPos;
    } else {                         //  weird defline (e.g., ">"); start at the beginning
        firstPos = 0;
    }

    //  Truncate the defline immediately, assuming the relevant identifier is up front.
    if (defline.size() > 0) {

        //  If the defline is ">", use the default string.
        if (defline != ">") {
            tmp = defline.substr(firstPos, len - firstPos);
        }

        NStr::TruncateSpacesInPlace(tmp, NStr::eTrunc_Begin);
        if (tmp.length() >= padToLength) {
            tmp = tmp.substr(0, padToLength - 1);
        }
        NStr::TruncateSpacesInPlace(tmp, NStr::eTrunc_End);


        //  Chop off everything after the first whitespace character.
/*
        spacePos = tmp.find(' ', 0);
        if (spacePos != string::npos) {
            tmp = tmp.substr(0, spacePos);
        }
*/

        //  Replace all space characters with underscores.
        spacePos = tmp.find(' ', 0);
        while (spacePos != string::npos) {
            tmp.replace(spacePos, 1, 1, '_');
            spacePos = tmp.find(' ', spacePos + 1);
        }

        //  Replace all tab characters with underscores.
        tabPos = tmp.find('\t', 0);
        while (tabPos != string::npos) {
            tmp.replace(tabPos, 1, 1, '_');
            tabPos = tmp.find('\t', tabPos + 1);
        }

        //  Replace all pipe characters with underscores.
        pipePos = tmp.find('|', 0);
        while (pipePos != string::npos) {
            tmp.replace(pipePos, 1, 1, '_');
            pipePos = tmp.find('|', pipePos + 1);
        }

        if (tmp.length() == 0) {
            tmp = defaultStr;
        }

        result.replace(0, tmp.length(), tmp);
    }

    return result;
}

bool CSpeerAlignment::WriteAsClustal(ostream& os, string& err) 
{
    bool result = (os.good());
    err.erase();

    if (!result) {
        err = "WriteAsClustal:  stream is not writable.\n";
        return result;
    }

    if ((m_rowMap.size() == 0 || m_rowMap.begin()->second.seq.size() == 0) && !PopulateRowStrings()) {
        err = "Unexpected problem extracting sequence data for alignment.\n";
        return result;
    }

    SRowData rowData;
    map<unsigned int, SRowData>::iterator rowIt = m_rowMap.begin(), rowEnd = m_rowMap.end();


    //  Necessary header to be parsed as Clustal.
    //  Then, print identifier, followed by the sequence, all on one line vs. broken up in chunks
    os << "CLUSTAL V\n\n";
    for ( ; rowIt != rowEnd; ++rowIt) {

        os << SpeerGetIDFromDefline(rowIt->second.defline, rowIt->second.row);
        os << rowIt->second.seq << endl;
    }
    os << endl << endl;

    return result;
}

bool CSpeerAlignment::WriteAsFasta(ostream& os, string& err) 
{
    bool result = (os.good());
    err.erase();

    if (!result) {
        err = "WriteAsFasta:  stream is not writable.\n";
        return result;
    }

    if ((m_rowMap.size() == 0 || m_rowMap.begin()->second.seq.size() == 0) && !PopulateRowStrings()) {
        err = "Unexpected problem extracting sequence data for alignment.\n";
        return result;
    }

    SRowData rowData;
    map<unsigned int, SRowData>::iterator rowIt = m_rowMap.begin(), rowEnd = m_rowMap.end();

    for ( ; rowIt != rowEnd; ++rowIt) {

        os << rowIt->second.defline << endl;
        os << rowIt->second.seq << endl;
    }
    return result;
}

CSpeerAlignment::CSpeerAlignment(CCdFromFasta* cd, bool gaplessColumnsOnly, ESpeerCaseSensitivity caseSensitivity) : m_baseAlignment(NULL), m_caseSensitivity(caseSensitivity), m_rowFrom(0), m_rowTo(0)
{
    if (!cd) return;

    map<unsigned int, string> columns;
    map<unsigned int, string>::iterator colIt, colEnd;

    if (gaplessColumnsOnly) {
        cd->GetGaplessColumnsReadFromFile(columns);
    } else {
        cd->GetAllColumnsReadFromFile(columns);
    }

    colEnd = columns.end();
    for (colIt = columns.begin(); colIt != colEnd; ++colIt) {
        switch (caseSensitivity) {
        case eSpeerConvertToUpper:
            NStr::ToUpper(colIt->second);
            break;
        case eSpeerConvertToLower:
            NStr::ToLower(colIt->second);
            break;
        case eSpeerCaseSensitive:
        default:
            break;
        };

        m_columns.insert(TColumnMapVT(colIt->first, CSpeerColumn(colIt->first, colIt->second)));
    }

    m_rowTo = cd->GetNumRows() - 1;

    m_rowMap.clear();
    for (unsigned int row = m_rowFrom; row <= m_rowTo; ++row) {
        m_rowMap[row].row = row;
        m_rowMap[row].defline = cd->GetDeflineReadFromFile(row);
    }
}

CSpeerAlignment::CSpeerAlignment(const CSpeerAlignment& rhs, unsigned int rowFrom, unsigned int rowTo, ESpeerCaseSensitivity caseSensitivity)
{
    string subfamilyColumn, defline;
    TColumnMapCit cit, citEnd;
    map<unsigned int, SRowData>::const_iterator deflineCit, deflineEnd;

    m_rowFrom = rowFrom;
    m_rowTo = rowTo;
    m_baseAlignment = &rhs;

    // Iterate through rhs.m_columns and extract indicated substrings
    // to initialize this->m_columns.
    citEnd = rhs.m_columns.end();
    for (cit = rhs.m_columns.begin(); cit != citEnd; ++cit) {
        subfamilyColumn = cit->second.GetResidues(rowFrom, rowTo);
        switch (caseSensitivity) {
        case eSpeerConvertToUpper:
            NStr::ToUpper(subfamilyColumn);
            break;
        case eSpeerConvertToLower:
            NStr::ToLower(subfamilyColumn);
            break;
        case eSpeerCaseSensitive:
        default:
            break;
        };
        m_columns.insert(TColumnMapVT(cit->first, CSpeerColumn(cit->first, subfamilyColumn)));
    }

    m_rowMap.clear();
    deflineEnd = rhs.m_rowMap.end();
    for (unsigned int i = 0, row = m_rowFrom; row <= m_rowTo; ++i, ++row) {
        m_rowMap[i].row = i;
        deflineCit = rhs.m_rowMap.find(row);
        if (deflineCit != deflineEnd) {
            m_rowMap[i].defline = deflineCit->second.defline;
        } else {
            m_rowMap[i].defline = "UnknownID-Row" + NStr::UIntToString(i);
        }
    }

}

CSpeerAlignment* CSpeerAlignment::MakeSubfamily(unsigned int rowFrom, unsigned int rowTo, ESpeerCaseSensitivity caseSensitivity) const
{
    CSpeerAlignment* result = NULL;
    
    if (rowFrom <= rowTo && rowTo < m_rowTo - m_rowFrom + 1) {
        result = new CSpeerAlignment(*this, rowFrom, rowTo, caseSensitivity);
    }
    return result;
}

bool CSpeerAlignment::GetRowData(unsigned int rowIndex, SRowData& rowData) 
{
    bool result = false;
    if ((m_rowMap.size() == 0 || m_rowMap.begin()->second.seq.size() == 0) && !PopulateRowStrings()) {
        return result;
    }

    map<unsigned int, SRowData>::const_iterator rowIt = m_rowMap.find(rowIndex);
    if (rowIt != m_rowMap.end()) {
        rowData = rowIt->second;
        result = true;
    }
    return result;
}


const CSpeerColumn* CSpeerAlignment::GetColumn(unsigned int index) const
{
    TColumnMapCit cit = m_columns.find(index);
    if (cit != m_columns.end()) {
        return &(cit->second);
    } else {
        return NULL;
    }
}

bool CSpeerAlignment::GetColumn(unsigned int index, CSpeerColumn& column) const
{
    bool result = false;
    TColumnMapCit cit = m_columns.find(index);
    if (cit != m_columns.end()) {
        result = true;
        column = cit->second;
    }
    return result;
}

bool CSpeerAlignment::PopulateRowStrings()
{
    unsigned int row, col, nRows, nCols;    
    string colString;
    string::iterator strIt, strEnd;
    map<unsigned int, SRowData>::iterator rowIt, rowEnd;
    TColumnMapCit cit, citEnd = m_columns.end();

    //  Set aside enough space for each of the row strings.
    nCols = NumColumns();
    nRows = NumRows();
    for (row = 0; row < nRows; ++row) {
        m_rowMap[row].seq.erase();
        m_rowMap[row].seq.reserve(nCols);
    }
    rowEnd = m_rowMap.end();

//    cerr << nRows << " " << nCols << endl;

    // Iterate through m_columns and build up the strings for each row.
    col = 0;
    for (cit = m_columns.begin(); cit != citEnd; ++cit, ++col) {
        colString = cit->second.GetResidues();
        strEnd = colString.end();
        if (colString.length() != nRows) {
            ERR_POST("Unexpected column length " << colString.length() << " found for column " << cit->first << ".  Expected " << nRows << ".\nColumn string:  " << colString);
            m_rowMap.clear();
            return false;
        }
        for (row = 0, rowIt = m_rowMap.begin(), strIt = colString.begin(); row < nRows; ++row, ++rowIt, ++strIt) {
            rowIt->second.seq += *strIt;
//            if (row == 0) {
//                cerr << "    " << colString[row] << "; " << m_rowMap[row].seq << endl;
//            }
        }
    }

    return true;
}


ostream& operator<<(ostream& os, const CSpeerAlignment& speerAlignment)
{
    CSpeerAlignment::TColumnMapCit cit = speerAlignment.m_columns.begin(), cend = speerAlignment.m_columns.end();
    os << "Number of columns:  " << speerAlignment.m_columns.size() << NcbiEndl;
    for (; cit != cend; ++cit) {
        os << cit->second << NcbiEndl;
    }
    os << "========================================" << NcbiEndl;
    return os;
}

