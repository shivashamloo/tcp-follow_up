// $Id: maseFormat.h,v 1.1 2008/11/17 16:45:32 lanczyck Exp $

#ifndef ___MASE_FORMAT
#define ___MASE_FORMAT

#include "sequenceContainer.h"

class maseFormat{
public:
	static sequenceContainer read(istream &infile, const alphabet* alph);
	static void write(ostream &out, const sequenceContainer& sd);
	//readUnAligned: the input sequences do not need to be aligned (not all sequences are the same length).
	static sequenceContainer readUnAligned(istream &infile, const alphabet* alph);
};

#endif

/* EXAMPLE OF THE FORMAT:

;;this is the place for general remarks.
;here we put sequence specific remark. 
Langur
KIFERCELARTLKKLGLDGYKGVSLANWVCLAKWESGYNTEATNYNPGDESTDYGIFQINSRYWCNNGKPGAVDACHISCSALLQNNIADAVACAKRVVSDQGIRAWVAWRNHCQNKDVSQYVKGCGV
;
Baboon
KIFERCELARTLKRLGLDGYRGISLANWVCLAKWESDYNTQATNYNPGDQSTDYGIFQINSHYWCNDGKPGAVNACHISCNALLQDNITDAVACAKRVVSDQGIRAWVAWRNHCQNRDVSQYVQGCGV
;
Human
KVFERCELARTLKRLGMDGYRGISLANWMCLAKWESGYNTRATNYNAGDRSTDYGIFQINSRYWCNDGKPGAVNACHLSCSALLQDNIADAVACAKRVVRDQGIRAWVAWRNRCQNRDVRQYVQGCGV
;
Rat
KTYERCEFARTLKRNGMSGYYGVSLADWVCLAQHESNYNTQARNYDPGDQSTDYGIFQINSRYWCNDGKPRAKNACGIPCSALLQDDITQAIQCAKRVVRDQGIRAWVAWQRHCKNRDLSGYIRNCGV
;
Cow
KVFERCELARTLKKLGLDGYKGVSLANWLCLTKWESSYNTKATNYNPSSESTDYGIFQINSKWWCNDGKPNAVDGCHVSCSELMENDIAKAVACAKKIVSEQGITAWVAWKSHCRDHDVSSYVEGCTL
;
Horse
KVFSKCELAHKLKAQEMDGFGGYSLANWVCMAEYESNFNTRAFNGKNANGSSDYGLFQLNNKWWCKDNKRSSSNACNIMCSKLLDENIDDDISCAKRVVRDKGMSAWKAWVKHCKDKDLSEYLASCNL

*/


