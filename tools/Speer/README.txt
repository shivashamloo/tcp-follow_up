This file summarizes the usage of the 'SPEER' program for the prediction of
specificity determining sites in proteins.  This software is a C++
implementation of the original SPEER package described in
"Functional Specificity Lies within the Properties and Evolutionary Changes of
Amino Acids", Chakrabarti S, Bryant SH, Panchenko AR (2007).  J Mol Biol. 373(3):801-10.
PubMed ID:  17868687

The input file must be a FASTA-formatted multiple alignment representing a
protein family that contains two or more subfamilies.  Within the FASTA file,
each subfamily is defined by a subset of consecutive sequences, with each
subfamily's size specified by the user at the command-line.

Contact Chris Lanczycki (lanczyck@ncbi.nlm.nih.gov) or Saikat Chakrabarti
(chakraba@ncbi.nlm.nih.gov) for comments or reports of problems with this
program.  

Binaries for the PC, Mac and Linux are available; the source code, based on
the NCBI C++ Toolkit, is also available (instructions are in the file
BUILDING_SPEER.txt).  SPEER makes use of the Rate4Site program, which is
distributed with SPEER by consent of its authors (Mayrose, I., Graur, D.,
Ben-Tal, N., and Pupko, T. 2004. Comparison of site-specific rate-inference
methods: Bayesian methods are superior. Mol Biol Evol 21: 1781-1791). 


Software developed at NCBI is in the public domain, as per the notice at the end of this file.  

Last modified:  $Date: 2010/10/13 21:30:51 $

**  10/30/2008:  
    Initial version
**  11/17/2008:
    Minor modifications to text and command-line parameters
**  10/13/2010:
    Exclude 100% identical columns from output; minor modification to command-line parameters

============================================================================
Contents:
1)  Overview
2)  Scoring function
3)  Input formats
4)  Program output
5)  Command-line options
6)  Example
============================================================================

=============================
1)  Overview:

SPEER (Chakrabarti et al., 2007) combines Euclidean distances based on amino
acids' physico-chemical properties, evolutionary rate and combined relative
entropy to predict specificity determining sites. All three terms account for
the variability of sites within the subfamilies in terms of their
physico-chemical properties, evolutionary rates and amino acid types. The first
and the third terms also approximate the variability of physico-chemical
properties and amino acid types between the subfamilies.  


=============================
2)  Scoring Function

The SPEER scoring function is described in detail in the references at the top
of this file.  The following is a qualitative description.  Most simply, it is a
linear combination of three terms each with a user-definable weight between 0
and 1.

The SPEER score S_i of a query alignment column 'i' is computed for each column
whose gap fraction does not exceed a user-specified percentage (20% by default).

SPEER score == S_i = w_RE*(Z_RE)i + w_ED*(Z_ED)i + w_ER*(Z_ER)i,

where Z_RE, Z_ED, and Z_ER are the Z-scores computed for that column's relative
entropy (RE) score, Euclidean distance (ED) score, and evolutionary rate (ER)
score, respectively.   Z-scores are calculated by subtracting the mean value and
dividing by the standard deviation of the score distribution obtained for all
columns in a given alignment (to account for different background conservation
levels in different families).  Although the user is able to experiment with
variable weights for each term, leaving each weight equal to one appears to work
well in general.   

Consult the published references for details of the scoring scheme employed.

IMPORTANT:  We note that whenever w_ER != 0, computing Z_ER can take several
orders-of-magnitude longer than the other two terms of S_i.  As such, if
run-time is a major consideration, the ER scoring can be turned off by setting
w_ER to zero at the command-line.  Z_ER is computed based on the unmodified set
of column scores produced by Rate4Site. 

Rate4Site NOTE:
From the Rate4Site web site (http://www.tau.ac.il/~itaymay/cp/rate4site.html):

    "The Rate4Site program may suffer from underflow problems when a large
    number of sequences are used (typically more than 200 [in a SPEER
    subfamily]). To bypass this problem Rate4Site can be compiled using the
    Makefile_slow file. The executable produced when compiling with this option
    can handle a larger number of sequences but is slower."  

If you encounter such problems, see the file 'BUILDING_SPEER.txt' for
instructions on how to build Rate4Site to circumvent this issue.

=============================
3)  Input Formats:

The input files must be a FASTA-formatted multiple alignment containing
sequences from two or more subfamilies.  The user must organize the input
alignment such that all sequences representing subfamily 1 appear first,
followed by those for subfamily 2, etc.  No extra subfamily annotation of the
alignment is required, but you must provide the number of aligned sequences in
each subfamily at the command line or in a file specified with the -pf command
line option (see below).

While the Rate4Site executable, rate4site.exe, is expected to be in the same
directory as the SPEER executable itself, the -r4sDir command line option allows
you to provide an alternate location.  

IMPORTANT:  Be sure to keep the executable's name as 'rate4site.exe', wherever
you choose place it. 

=============================
4)  Program output:

SPEER provides output for each query column where the fraction of gaps is below
a user-defined value.  In a tabular format, the program provides the column
number, SPEER score, the corresponding Z-score of the SPEER score, and the
probability P of obtaining that SPEER score or higher purely by chance from a
normal distribution.  For reference, a percent identity for each column is
provided in the second column.  This is useful information because it is
difficult to make predictions of subfamily specificity when there has been no
apparent evolutionary change at a site.  Note that all columns with 100%
identity have been excluded from the output since they are unlikely to be
involved in conveying specificity to that site.


The output can be presented in one of three orders, specified by the -f option.
The default (-f 0) is alignment column order, where the column number refers to
the position on the first sequence of the input alignment.  To order the columns
from highest-to-lowest SPEER score, use -f 1.  Finally, -f 2 print the output in
descending order by column percent identity.


=============================
5)  Command line options:

There are two mandatory arguments:

      i)  -i, the file containing the alignment in FASTA format,  


      ii) the number of sequences in each subfamily in the input alignment as
          a)  integers at the end of the command line, following all other options
              -- OR -- 
          b)  integers, one per line, in the file given by the "-pf" option.


E.g., to specify that the alignment in myAlignment.FSA has subfamilies of size
25, 5 and 8:

./SPEER -i myAlignment.FSA  25 5 8

-- OR --

./SPEER -i myAlignment.FSA -pf subfamSizes.txt

where the file 'subfamSizes.txt' contains '25', '5' and '8' each on a separate
line. 


We reiterate that by default, all weights (-wRE, -wEDist, -wERate) are set to 1.
If you want to turn any term off in the SPEER score you must explicitly set it
to zero at the command line.  At least one weight must have a non-zero value.


You can view a summary of command line options by typing:

./SPEER -h

More verbose help is also available:

./SPEER -help


For reference, here is the verbose help output.

USAGE
  SPEER [-h] [-help] -i FilenameIn [-o OutputFileName] [-pf FilenameIn]
    [-Jobid job_id] [-wRE RelativeEntropyWeight]
    [-wEDist EvolutionaryDistanceWeight] [-wERate EvolutionaryRateWeight]
    [-f OutputFormat] [-g double] [-r4sDir DirectoryFolder] [....]

DESCRIPTION
   SPEER(Specificity prediction using amino acids' Properties, Entropy and
   Evolution Rate)

REQUIRED ARGUMENTS
 -i <String>
   full filename of input CD (ascii or add -b flag for binary) or a
   FASTA-formatted alignment (ascii only)

OPTIONAL ARGUMENTS
 -h
   Print USAGE and DESCRIPTION;  ignore other arguments
 -help
   Print USAGE, DESCRIPTION and ARGUMENTS description;  ignore other arguments
 -o <String>
   filename for program output
   
 -pf <String>
   Subfamily sizes:  file containing the partitioning of the input alignment
   into subfamilies, one integer per line (optionally followed by a
   descriptive string) giving the number of rows in the subfamily.  If there
   are leftover rows, they will constitute a subfamily.
   [Overrides any subfamily sizes given at the end of the command-line (see
   '....' below).]
   
 -Jobid <Integer, 1..2147483647>
   Job id used for overriding the names of files generated by rate4site.
   
 -wRE <Real, 0..1>
   weight assigned to the relative entropy term in the total score of a column
   Default = `1'
 -wEDist <Real, 0..1>
   weight assigned to the evolutionary distance term in the total score of a
   column
   Default = `1'
 -wERate <Real, 0..1>
   weight assigned to the evolutionary rate term in the total score of a
   column
   Default = `1'
 -f <Integer>
   Specify output column order:  0 = column number order [default]; 1 = order
   of decreasing score; 2 = order of decreasing column %ID.
   Default = `0'
 -g <Real, 0..1>
   maximum fraction of gaps allowed in a column (all subfamilies must exceed
   this fraction in each column scored)
   Default = `0.2'
 -r4sDir <String>
   The directory or folder containing the 'rate4site.exe' executable.
   Ignored if the evolutionary rate weight (-wERate) is zero.
   [Default is the same directory that contains the current working
   directory.]
   
 .... <Integer>
   Subfamily sizes:  partitioning of the input alignment into subfamilies
   using space-delimited positive integers to give the number of rows in each
   subfamily.  If there are leftover rows, they will constitute a subfamily.
   [These values are ignored when the -pf option (described above) is also
   specified.]
   ** When present, these integers must be the final arguments on the command
   line.**


=============================
6)  Example:

In the 'example' directory of the SPEER distribution is a FASTA-formatted
alignment, 'cd00423.FSA' (other examples of protein family alignments,
subfamilies and experimentally annotated sites can be found at
ftp://ftp.ncbi.nih.gov/pub/SPEER/Alignments_and_Sites).  The first 25 sequences
of 'cd00423.FSA' comprise one subfamily and the last 8 sequences comprise a
second subfamily.  Subfamily-specific sites for this alignment as predicted by
SPEER have the highest SPEER score, and are at the top of the list of columns in
the sample output file 'cd00423_bySpeerScore.out'. 


./SPEER -i cd00423.FSA -o cd00423_bySpeerScore.out -f 1 25 8


Importantly, the size of the two subfamilies are given as the last arguments on
the command line.  As stated above, you may choose to put these values in a file
and use the -pf option instead.

Note that the default scoring term weights (-wRE 1  -wEDist 1  -wERate 1) are
implied here, along with the '-f 1' option needed to specify ordering the output
by SPEER score rather than by column number. 

We also have included the sample output file 'cd00423_byColumn.out'.  It was
produced when running SPEER for the same alignment but ordered by column index,
simply by dropping '-f 1' from the command given above.  

Finally, whenever -wERate is non-zero, rate4site.exe will run.  Because
this is by far the most computationally intensive step, we permit rate4site.exe 
to print messages to the terminal as a sort of progress report.  For 
reference, an example of the rate4site.exe output for the command above is in
the file 'rate4site.out'.  To avoid seeing rate4site.exe output, at the command
line redirect it to a file (or /dev/null): 

./SPEER -i cd00423.FSA -o cd00423_bySpeerScore.out -f 1 25 8 >& rate4site.out



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
