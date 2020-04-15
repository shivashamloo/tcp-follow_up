/*  $Id: speer.cpp,v 1.15 2010/10/13 21:30:51 lanczyck Exp $
 * ===========================================================================
 *
 *                            PUBLIC DOMAIN NOTICE
 *               National Center for Biotechnology Information
 *
 *  This software/database is a "United States Government Work" under the
 *  terms of the United States Copyright Act.  It was written as part of
 *  the author's official duties as a United States Government employee and
 *  thus cannot be copyrighted.  This software/database is freely available
 *  to the public for use. The National Library of Medicine and the U.S.
 *  Government have not placed any restriction on its use or reproduction.
 *
 *  Although all reasonable efforts have been taken to ensure the accuracy
 *  and reliability of the software and data, the NLM and the U.S.
 *  Government do not and cannot warrant the performance or results that
 *  may be obtained by using this software or data. The NLM and the U.S.
 *  Government disclaim all warranties, express or implied, including
 *  warranties of performance, merchantability or fitness for any particular
 *  purpose.
 *
 *  Please cite the author in any work or product based on this material.
 *
 * ===========================================================================
 *
 * Authors:  Chris Lanczycki
 *
 * File Description:
 *           
 *
 */

#include <ncbi_pch.hpp>
#include <corelib/ncbistd.hpp>
#include <corelib/ncbifile.hpp>

#include <objects/cdd/Cdd_id.hpp>
#include <objects/cdd/Cdd_id_set.hpp>
#include <objects/cdd/Global_id.hpp>

//#include <algo/structure/struct_util/su_sequence_set.hpp>

#include <serial/objistr.hpp>
#include <serial/objostr.hpp>
#include <algo/structure/bma_refine/BMAUtils.hpp>
#include <algo/structure/cd_utils/cuCdReadWriteASN.hpp>
#include <algo/structure/cd_utils/cuReadFastaWrapper.hpp>
#include <algo/structure/cd_utils/cuSeqAnnotFromFasta.hpp>
#include <algo/structure/cd_utils/cuCdFromFasta.hpp>

#include <math.h>
#include "speer.hpp"
#include "speerScore.hpp"
#include "speerUtils.hpp"

USING_NCBI_SCOPE;
USING_SCOPE(align_refine);
USING_SCOPE(objects);
USING_SCOPE(cd_utils);

template < class ASNClass >
static bool ReadASN(const char* filename, ASNClass *asnObject, bool isBinary, std::string *err)
{
    bool result = true;
    err->erase();
    try {
        ifstream fileStream(filename);
        if (!fileStream) throw runtime_error("Failed to open file!" + string(filename));
        if (isBinary) {
            fileStream >> MSerial_AsnBinary >> *asnObject;
        } else {
            fileStream >> MSerial_AsnText >> *asnObject;
        }
    } catch (runtime_error re) {
        *err = re.what();
        result = false;
    } catch (...) {
        *err = "Could not open or extract ASN object from file " + string(filename);
        result = false;
    }
    return result;
}

static bool ReadFasta(const char* filename, CCdFromFasta*& fastaCD, std::string *err)
{
    CCdFromFasta::Fasta2CdParams params;
    CDirEntry dirEntry(filename);
    params.cdName = dirEntry.GetBase();
    params.cdAcc = params.cdName;
    params.useLocalIds = true;
    params.useAsIs = true;
    params.masterIndex = 0;
    params.masterMethod = CSeqAnnotFromFasta::eFirstSequence;

    if (fastaCD) {
        delete fastaCD;
        fastaCD = NULL;
    }

    fastaCD = new CCdFromFasta(filename, params);
    if (fastaCD) {
        string source = "FASTA source file:  " + string(filename);
        fastaCD->AddCreateDate();
        fastaCD->AddSource(source);
        if (fastaCD->WasInputError() && err) {
            *err = fastaCD->GetFastaInputErrorMsg();
        }
    } else if (err) {
        *err = "Unknown parsing error during ReadFasta!\n";
    }

    return (fastaCD != NULL && !fastaCD->WasInputError());
}



const string CSpeerApp::normalModeStr = "prod";
const string CSpeerApp::testingModeStr = "test";


char CSpeerApp::MostCommonResidue(const vector<char>& residues, unsigned int& appearances) {

    char c, mostCommon = '-';
    map<char, unsigned int> charMap;
    map<char, unsigned int>::iterator it, itEnd;
    vector<char>::const_iterator cit = residues.begin(), citEnd = residues.end();
    unsigned int i, n = residues.size();

    appearances = 0;
    for (i = 0; i < n; ++i) {
        c = residues[i];
        if (charMap.count(c) == 0) 
            charMap.insert(map<char, unsigned int>::value_type(c, 1));
        else
            ++charMap[c];
    }

    itEnd = charMap.end();
    for (it = charMap.begin(); it != itEnd; ++it) {
        if (it->second > appearances) {
            mostCommon  = it->first;
            appearances = it->second;
        }
    }
    
    return mostCommon;

}


CSpeerApp::CSpeerApp() : m_outStream(NULL), m_inputAlignment(NULL), m_outputOrder(0), m_jobid(-1)
{
    SetVersion(CVersionInfo(1, 0, 0, "SPEER:  A program to predict subfamily specific sites in a protein multiple alignment"));
}

CSpeerApp::~CSpeerApp()
{
    delete m_inputAlignment;
    for (unsigned int i = 0; i < m_subfamilies.size(); ++i) delete m_subfamilies[i];
}

/////////////////////////////////////////////////////////////////////////////
//  Init test for all different types of arguments

void CSpeerApp::Init(void)
{

    // Create command-line argument descriptions class
    auto_ptr<CArgDescriptions> argDescr(new CArgDescriptions);

    // Specify USAGE context
    argDescr->SetUsageContext(GetArguments().GetProgramBasename(),
                              "SPEER(Specificity prediction using amino acids' Properties, Entropy and Evolution Rate)");

    //  Hide some of the standard arguments...
    HideStdArgs(fHideLogfile | fHideConffile | fHideVersion | fHideFullVersion | fHideDryRun | fHideXmlHelp);

    argDescr->AddKey("i", "FilenameIn", "full filename of input CD (ascii or add -b flag for binary) or a FASTA-formatted alignment (ascii only)", argDescr->eString);
    argDescr->AddOptionalKey("o", "OutputFileName", "filename for program output\n\n", argDescr->eString);
    argDescr->AddOptionalKey("pf", "FilenameIn", "Subfamily sizes:  file containing the partitioning of the input alignment into subfamilies, one integer per line (optionally followed by a descriptive string) giving the number of rows in the subfamily.  If there are leftover rows, they will constitute a subfamily.\n[Overrides any subfamily sizes given at the end of the command-line (see '....' below).]\n\n", argDescr->eString);
    argDescr->AddOptionalKey("Jobid", "job_id", "Job id used for overriding the names of files generated by rate4site.\n\n", argDescr->eInteger);
    argDescr->SetConstraint("Jobid", new CArgAllow_Integers(1, kMax_Int));


    argDescr->AddDefaultKey("wRE", "RelativeEntropyWeight", "weight assigned to the relative entropy term in the total score of a column", argDescr->eDouble, "1");
    argDescr->AddDefaultKey("wEDist", "EvolutionaryDistanceWeight", "weight assigned to the evolutionary distance term in the total score of a column", argDescr->eDouble, "1");
    argDescr->AddDefaultKey("wERate", "EvolutionaryRateWeight", "weight assigned to the evolutionary rate term in the total score of a column", argDescr->eDouble, "1");
    argDescr->SetConstraint("wRE", new CArgAllow_Doubles(0, 1));
    argDescr->SetConstraint("wEDist", new CArgAllow_Doubles(0, 1));
    argDescr->SetConstraint("wERate", new CArgAllow_Doubles(0, 1));

    argDescr->AddDefaultKey("f", "OutputFormat", "Specify output column order:  0 = column number order [default]; 1 = order of decreasing score; 2 = order of decreasing column %ID.", argDescr->eInteger, "0");

    argDescr->AddDefaultKey("g", "double", "maximum fraction of gaps allowed in a column (all subfamilies must exceed this fraction in each column scored)", argDescr->eDouble, "0.2");
    argDescr->SetConstraint("g", new CArgAllow_Doubles(0, 1));
    argDescr->AddOptionalKey("r4sDir", "DirectoryFolder", "The directory or folder containing the 'rate4site.exe' executable.\nIgnored if the evolutionary rate weight (-wERate) is zero.\n[Default is the same directory that contains the current working directory.]\n\n", argDescr->eString);

    argDescr->AddExtra(0, kMax_UInt, "Subfamily sizes:  partitioning of the input alignment into subfamilies using space-delimited positive integers to give the number of rows in each subfamily.  If there are leftover rows, they will constitute a subfamily.\n[These values are ignored when the -pf option (described above) is also specified.]\n** When present, these integers must be the final arguments on the command line.**\n", argDescr->eInteger);


//    argDescr->AddFlag("fa2cd", "for FASTA input only:  write out the CD created from the FASTA file.\nName of output will be the input file name with the file extension replaced by '.cn3'\nSuch files can be used in Cn3D or CDTree.");

//  Cmdline options for internal use only.
//#ifndef _PUBLIC_
//    argDescr->SetCurrentGroup("Options for internal use");
//    argDescr->SetCurrentGroup("");
//#endif

    // Setup arg.descriptions for this application
    SetupArgDescriptions(argDescr.release());
}


bool CSpeerApp::SetSubfamilyPartitioning(string& errStr)
{
    bool result = true;
    const CArgs& args = GetArgs();
    size_t nExtra = args.GetNExtra();
    string partitionFile = (args["pf"]) ? args["pf"].AsString() : kEmptyStr;
    string s;
    int intN;
    unsigned int i, n;

    m_partitions.clear();
    if (partitionFile.length() > 0) {
        CNcbiIfstream pf(partitionFile.c_str());
        if (pf.good()) {
            try {
                while (pf >> s) {

                    //  Will throw an exception if the conversion fails.
                    n = NStr::StringToUInt(s, NStr::fAllowLeadingSpaces | NStr::fAllowTrailingSpaces);
                    LOG_POST("partition of size " << n << " to be created.");
                    m_partitions.push_back(n);
                }

            } catch (...) {
                errStr = "While reading " + partitionFile + ", error converting '" + s + "' to a positive integer.\nStopping.";
                result = false;
            }
        } else {
            errStr = "Could not open file " + partitionFile + " to read partioning data.";
            result = false;
        }
        pf.close();

    } else if (nExtra > 0) {
        for (i = 1; i <= nExtra; ++i) {
            intN = args[i].AsInteger();
            if (intN > 0) {
                m_partitions.push_back((unsigned int) intN);
            } else {
                errStr = "Invalid partition size '" + NStr::IntToString(intN) + "' found:  must be a positive integer.";
                result = false;
                break;
            }
        }
    } else {
        errStr = "No partitioning information provided.\nGive subfamily sizes at the end of the command-line,\nor specify a file containing the subfamily sizes using the \"-pf\" option.";
        result = false;
    }
    return result;
}

/////////////////////////////////////////////////////////////////////////////
//  Main Run method; invoke requested tests (printout arguments obtained from command-line)
int CSpeerApp::Run(void)
{

    bool dumpFastaAsCd;
    string r4sDirFromCmdLine;
    string errStr, hitStr, dumpStr, columnResidues;
    string ranDbFile, cdInFile, outFile, cdAcc, compareMethodStr;
    unsigned int i;
    unsigned int nColsFasta, nSubfamilies;
    unsigned int nRows = 0, rowSum = 0;
    unsigned int rowFrom = 0, rowTo = 0;
    double gapFractionThreshold;  //  maximum fraction of gaps allowed in an alignment column
    double reWeight, eDistWeight, eRateWeight;
    vector<int> scores;
    vector<char> residues;

    CSpeerAlignment* subfamilyAlignment;

    //  Important:  sets up the reference frequencies!!
    SpeerInitializeDefaultGlobalData();

    // Get arguments
    CArgs args = GetArgs();
    cdInFile = args["i"].AsString();
    outFile = (args["o"]) ? args["o"].AsString() : "";
    m_outStream = (outFile.size() > 0) ? new ofstream(outFile.c_str(), ios::out) : &cout;
    SetDiagPostFlag(eDPF_Log);      // show only the provided message
    SetDiagStream(m_outStream);

    m_outputOrder = (unsigned int) args["f"].AsInteger();

    gapFractionThreshold = args["g"].AsDouble();
    reWeight = args["wRE"].AsDouble();
    eDistWeight = args["wEDist"].AsDouble();
    eRateWeight = args["wERate"].AsDouble();

    if (reWeight == 0.0 && eDistWeight == 0.0 && eRateWeight == 0.0) {
        ERR_POST("Invalid command:  At least one non-zero weight must be supplied via the -wRE, -wEDist and/or -wERate options.\nTerminating.");
        return -4;
    }

    //  If rate4site is to be used, make sure everything is OK.
    r4sDirFromCmdLine = (args["r4sDir"]) ? args["r4sDir"].AsString() : kEmptyStr;
    m_rate4SiteDir = CDirEntry::ConvertToOSPath(r4sDirFromCmdLine);
//    if (m_rate4SiteDir.length() > 0) {
//        LOG_POST("rate4Site directory = '"<< m_rate4SiteDir << "' (speer argument = '" << args["r4sDir"].AsString() << "')\n");
//    }

    m_jobid = (args["Jobid"]) ? args["Jobid"].AsInteger() : 0;
    if (eRateWeight > 0.0 && ! ValidateRate4Site()) {
        return -5;
    }

    //  Put this here to avoid the warning about failing to open the logfile.
    SetDiagPostLevel(eDiag_Info);   

    // read alignment input file
    string asnErr, fastaErr, clustalErr;
    CCdFromFasta* cdFromFasta = NULL;

    if (ReadFasta(cdInFile.c_str(), cdFromFasta, &fastaErr) && cdFromFasta) {

        dumpFastaAsCd = false;  //  (args["fa2cd"]);
        if (fastaErr.length() > 0) {
            ERR_POST("Error while parsing file " << cdInFile << " as a FASTA-formatted alignment.\nStopping.\n\n" << fastaErr);
            delete cdFromFasta;
            return -1;
        }

        //  Do this so I can use the rest of the code without modification.
        if (dumpFastaAsCd) {
            CCdd cd;
            CDirEntry inFilename(cdInFile);
            string outFilename = CDirEntry::MakePath(inFilename.GetDir(), inFilename.GetBase(), "cn3");

            cd.Assign(*((CCdd*) cdFromFasta));
            if (!WriteASNToFile(outFilename.c_str(), cd, false, &fastaErr)) {
                ERR_POST("Error writing FASTA to a CD-formatted file " << outFilename << ":\n" << fastaErr << "\n\n");
            }
        }
 
        string firstFastaSequence = cdFromFasta->GetSequenceReadFromFile(0);
        nColsFasta = firstFastaSequence.length();
        nRows = cdFromFasta->GetNumRows();
        LOG_POST("\nInput FASTA from file " << cdInFile);
        LOG_POST("Number of columns in the alignment (== length of first sequence including gaps) = " << nColsFasta);
        LOG_POST("Number of rows in the alignment = " << nRows << "\n");

    } else {
        ERR_POST("Input error reading from file " << cdInFile);
        if (fastaErr.length() > 0) {
            ERR_POST(fastaErr << "\n\n");
        }
        return -1;
    }

    if (!SetSubfamilyPartitioning(errStr)) {
        ERR_POST(" " << errStr);
        return -2;
    }

    nSubfamilies = m_partitions.size();
    LOG_POST("Partitioning the input alignment into " << nSubfamilies << " subfamilies.");

    //  Sanity check the partitioning.
    for (i = 0; i < nSubfamilies; ++i) {
        rowSum += m_partitions[i];
        if (gapFractionThreshold > 0.0) {
            LOG_POST("    Subfamily " << i + 1 << " consists of " << m_partitions[i] << " rows; only columns with at most " << 100.0*gapFractionThreshold << "% gaps (<= " << (unsigned int) (m_partitions[i]*gapFractionThreshold) << ") are analyzed.");
        } else {
            LOG_POST("    Subfamily " << i + 1 << " consists of " << m_partitions[i] << " rows; only gapless columns are analyzed.");
        }
    }


    if (rowSum != nRows) {
        if (rowSum > nRows) {
            ERR_POST("Partitioning error:  Too many rows (" << rowSum << ") specified; only " << nRows << " rows present.\n(Make sure all fasta deflines start with the '>' character.)\n");
            return -3;
        } else {
            LOG_POST("The first " << rowSum << " of the " << nRows << " rows in the input alignment\nhave been specified for inclusion in a subfamily.\nThe remaining rows are being ignored.");
        }
    }
    
    LOG_POST("\nSpeer scoring term weights:\n  Relative entropy term weight (-wRE)         = " << reWeight << "\n  Evolutionary distance term weight (-wEDist) = " << eDistWeight << "\n  Evolutionary rate term weight (-wERate)     = " << eRateWeight << "\n");

    LOG_POST("\n\n");


    //  Create the base alignment, and the specified subfamilies.
    m_inputAlignment = new CSpeerAlignment(cdFromFasta, false);

//    cout << "Input alignment:\n" << *m_inputAlignment;

    for (i = 0; i < nSubfamilies; ++i) {
        rowTo = rowFrom + m_partitions[i] - 1;
        subfamilyAlignment = m_inputAlignment->MakeSubfamily(rowFrom,rowTo);
        if (!subfamilyAlignment) {
            ERR_POST("Subfamily creation error:  Unable to make subfamily alignment " << i << " for rows " << rowFrom << " to " << rowTo << ".\nSkip and continue with other subfamilies.\n");
//        } else {
//            LOG_POST("Subfamily " << i << " (" << rowFrom << ", " << rowTo << "):\n" << *subfamilyAlignment);
        }
        m_subfamilies.push_back(subfamilyAlignment);

        rowFrom = rowTo + 1;
    }

    delete cdFromFasta;


    CSpeerScore speerScore;
    if (m_jobid > 0) {
        speerScore.SetJobid(m_jobid);
    }
    if (reWeight != 0.0) {
        speerScore.AddScorer(eRelativeEntropyScorer, m_subfamilies, gapFractionThreshold, reWeight);
    }
    if (eDistWeight != 0.0) {
        //  When combining evolutionary distance term, flip the Z-score since lower
        //  distances reflect more likely subfamily specificity.
        speerScore.AddScorer(eEvolutionaryDistanceScorer, m_subfamilies, gapFractionThreshold, -1.0*eDistWeight);
    }
    if (eRateWeight != 0.0) {
        //  When combining evolutionary rate term, flip the Z-score since lower
        //  rates reflect more likely subfamily specificity.
        speerScore.AddScorer(eEvolutionaryRateScorer, m_subfamilies, gapFractionThreshold, -1.0*eRateWeight, m_rate4SiteDir);
    }

    CSpeerScorer_Base* compositeScorer = speerScore.Compute();


    if (compositeScorer) {
        OutputResults(*compositeScorer);
        delete compositeScorer;

    } else {
        ERR_POST("Could not combine the scores while processing output!\nTerminating abnormally.");
        return -6;
    }

    return 0;

}

bool CSpeerApp::ValidateRate4Site() const
{
    static const char pathSep = CDirEntry::GetPathSeparator();
    static const string r4sExe = "rate4site.exe";

    bool result = false;
    string dotSlash = ".";
    string r4sPath = kEmptyStr;
    dotSlash += pathSep;

    //  Need to check in the current directory
    if (m_rate4SiteDir.length() == 0 || m_rate4SiteDir == "." || m_rate4SiteDir == dotSlash) {
        r4sPath = dotSlash + r4sExe;
    //  Need to check in other specified directory
    } else {
        r4sPath = m_rate4SiteDir;  //CDirEntry::ConvertToOSPath(m_rate4SiteDir);
        CDir r4sDir(r4sPath);
        if (!r4sDir.Exists()) {
            ERR_POST("Directory '" << m_rate4SiteDir << "' does not exist.\n");
            r4sPath = kEmptyStr;
        } else if (!r4sDir.IsDir()) {
            ERR_POST("You have specified a location '" << m_rate4SiteDir << "' that is not a directory.\n");
            r4sPath = kEmptyStr;
        } else if (!r4sDir.CheckAccess(CDirEntry::fDefaultDirOther)) {
            ERR_POST("You do not appear to have read & execute access for '" << m_rate4SiteDir << "'.\nPlease check mode flags for this directory, or ask an administrator to help you do so.\n");
            r4sPath = kEmptyStr;
        } else {
            r4sPath = CDirEntry::AddTrailingPathSeparator(r4sPath) + r4sExe;
        }
    }

    if (r4sPath.length() > 0) {
        CFile r4sFile(r4sPath);
        if (!r4sFile.Exists()) {
            ERR_POST("'rate4site.exe' not found in " << m_rate4SiteDir << ".\nCheck the directory given with the -r4s option, or move rate4site.exe to speer's directory.\n");
        } else if (!r4sFile.CheckAccess(CDirEntry::fExecute)) {
            ERR_POST("'rate4site.exe' does not have execute permissions.\nMake sure the proper mode flags are set, or ask an administrator to help you do so.\n");
        } else {
            LOG_POST("Rate4Site to be run using the executable at '" << r4sPath << "'.\n");
            result = true;
        }
    }

    if (m_jobid < 0) {
        ERR_POST("Jobid = " << m_jobid << "; it must not be negative!\n");
    }

    return result;
}


void CSpeerApp::OutputResults(CSpeerScorer_Base& compositeScorer)
{
    unsigned int col;
    double score, zscore, pscore, percId;
    CSpeerColumn column;
    map<unsigned int, double> percIDs;
    const TColumnScoreMap& compScoreMap = compositeScorer.GetColumnScores();
    TColumnScoreMapCit scoreMapCit, scoreMapCend;
    scoreMapCend = compScoreMap.end();

    typedef multimap<TSpeerScore, unsigned int>::value_type MMVT;
    multimap<TSpeerScore, unsigned int> inverseScoreMap;
    multimap<TSpeerScore, unsigned int>::reverse_iterator invScoreMapRit, invScoreMapRend;

    LOG_POST("#                  SPEER Results\n#====================================================");
    LOG_POST("# Column     %ID       Score      Z-score     P-value\n#====================================================\n");

    //  Collect formatted strings in a strstream
    CNcbiOstrstream strStrm;
    strStrm.setf(IOS_BASE::fixed, IOS_BASE::floatfield);


    if (m_outputOrder > 0) {
        //  Invert the map to order by appropriate value.
        for (scoreMapCit = compScoreMap.begin(); scoreMapCit != scoreMapCend; ++scoreMapCit) {
            percId = 0.0;
            col = scoreMapCit->first;
            if (m_inputAlignment->GetColumn(col, column)) {
                percId = column.GetPercentIdentity();
                percIDs[col] = percId;
            }
            if (m_outputOrder == 2) {
                inverseScoreMap.insert(MMVT(percId, col));
            } else {
                inverseScoreMap.insert(MMVT(scoreMapCit->second, col));
            }
        }

        invScoreMapRend = inverseScoreMap.rend();
        //  Compute %identity for each column of the input alignment in the results.
        for (invScoreMapRit = inverseScoreMap.rbegin(); invScoreMapRit != invScoreMapRend; ++invScoreMapRit) {
            col = invScoreMapRit->second;
            zscore = compositeScorer.GetColumnZScore(col);
            pscore = compositeScorer.GetColumnPScore(col);
            if (m_outputOrder == 2) {
                percId = invScoreMapRit->first;
                scoreMapCit = compScoreMap.find(col);
                score = (scoreMapCit != scoreMapCend) ? scoreMapCit->second : CSpeerScorer_Base::UNCOMPUTED;
            } else {
                score = invScoreMapRit->first;
                percId = percIDs[col];
            }
            strStrm << left << NcbiSetw(6) << col << right << NcbiSetprecision(1) << NcbiSetw(10) << percId << NcbiSetprecision(3) << NcbiSetw(12) << score << NcbiSetw(12) << zscore << NcbiSetw(12) << pscore << NcbiEndl;
        }

    } else {

        //  Compute %identity for each column of the input alignment in the results.
        for (scoreMapCit = compScoreMap.begin(); scoreMapCit != scoreMapCend; ++scoreMapCit) {
            percId = 0.0;
            col = scoreMapCit->first;
            score = scoreMapCit->second;
            zscore = compositeScorer.GetColumnZScore(col);
            pscore = compositeScorer.GetColumnPScore(col);
            if (m_inputAlignment->GetColumn(col, column)) {
                percId = column.GetPercentIdentity();
                percIDs[col] = percId;
            }
        
            strStrm << left << NcbiSetw(6) << col << right << NcbiSetprecision(1) << NcbiSetw(10) << percId << NcbiSetprecision(3) << NcbiSetw(12) << score << NcbiSetw(12) << zscore << NcbiSetw(12) << pscore << NcbiEndl;
        }
    }

    string s = CNcbiOstrstreamToString(strStrm);
    LOG_POST(s << "\n#====================================================\n");
}


/////////////////////////////////////////////////////////////////////////////
//  Cleanup


void CSpeerApp::Exit(void)
{
    SetDiagStream(0);
}

/////////////////////////////////////////////////////////////////////////////
//  MAIN
/////////////////////////////////////////////////////////////////////////////

int main(int argc, const char* argv[])
{

    //  Diagnostics code from Paul's struct_util demo.
    SetDiagStream(&NcbiCout); // send all diagnostic messages to cout
    SetDiagPostLevel(eDiag_Warning);   
    //    SetupCToolkitErrPost(); // reroute C-toolkit err messages to C++ err streams

    SetDiagTrace(eDT_Default);      // trace messages only when DIAG_TRACE env. var. is set
    SetDiagPostFlag(eDPF_Log);      // show only the provided message
    //#ifdef _DEBUG
    //    SetDiagPostFlag(eDPF_File);
    //    SetDiagPostFlag(eDPF_Line);
    //#else
    //UnsetDiagTraceFlag(eDPF_File);
    //UnsetDiagTraceFlag(eDPF_Line);
    //#endif

    // C++ object verification
    CSerialObject::SetVerifyDataGlobal(eSerialVerifyData_Always);
    CObjectIStream::SetVerifyDataGlobal(eSerialVerifyData_Always);
    CObjectOStream::SetVerifyDataGlobal(eSerialVerifyData_Always);

    // Execute main application function
    CSpeerApp speer;
    int result = speer.AppMain(argc, argv, 0, eDS_Default, 0);
    return result;
}

