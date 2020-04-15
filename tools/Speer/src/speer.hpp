/*  $Id: speer.hpp,v 1.5 2010/10/13 21:30:51 lanczyck Exp $  */

#ifndef SPEER__HPP
#define SPEER__HPP

#include <map>
#include <list>
#include <vector>
#include <corelib/ncbiexpt.hpp>
#include <corelib/ncbiapp.hpp>
#include <objects/cdd/Cdd.hpp>
#include "speerAlignment.hpp"
#include "speerScorer.hpp"

USING_NCBI_SCOPE;
USING_SCOPE(objects);

// class for standalone application
class CSpeerApp : public ncbi::CNcbiApplication
{

    static const string normalModeStr;
    static const string testingModeStr;

public:

    CSpeerApp();

    ~CSpeerApp();

    virtual void Init(void);
    virtual int  Run(void);
    virtual void Exit(void);

    //  Return the most common residue, along w/ the number of times it appeared.
    //  If a tie, return 'smallest' residue alphabetically.  If an error, a '-' is returned.
    static char MostCommonResidue(const vector<char>& residues, unsigned int& appearances);

    //  If there is no partition file given, looks for the 'extra' arguments on the command-line.
    //  Returns false if can't set up the partitioning - a fatal error.
    bool SetSubfamilyPartitioning(string& errStr);

private:

    ostream* m_outStream;
    vector<unsigned int> m_partitions;  //  number of rows in each subfamily alignment
    CSpeerAlignment* m_inputAlignment;
    vector<CSpeerAlignment*> m_subfamilies;

    string m_mode;  //  analyze a CD ("cd") or another database ("db")

    string m_rate4SiteDir;  //  name of the directory containing the 'rate4site' executable

    unsigned int m_outputOrder;  // 0 == column order; 1 == decreasing score; 2 == decreasing %ID

    int m_jobid;

    void OutputResults(CSpeerScorer_Base& compositeScorer);

    //  Check that the binary and the path to the binary exist when required.
    //  Return value is unaffected by whether or not user turns rate4site scoring term on.
    bool ValidateRate4Site(void) const;

};

#endif // SPEER__HPP 
