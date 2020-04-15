# $Date: 2010/10/13 21:30:51 $

APP = SPEER.exe

REQUIRES = objects

SRC = speerData \
      speerUtils \
      speerScore \
      speerScorer \
      speerColumn \
      speerAlignment \
      speer \
      zScoreGenerator \
      zToPTable


LIB =   xcd_utils \
        id1cli id1 \
        $(BLAST_LIBS) \
        cdd \
        ncbimime \
        scoremat \
        cn3d \
        seqset $(SEQ_LIBS) \
        mmdb \
        taxon1 \
        pub \
        medline \
        biblio \
        general \
        xser sequtil xobjread \
        xregexp $(PCRE_LIB) \
        xutil xconnect \
        xncbi $(SOBJMGR_LIBS)


CPPFLAGS = $(ORIG_CPPFLAGS) -D_PUBLIC_

LDFLAGS = $(FAST_LDFLAGS)

LIBS = $(NETWORK_LIBS) $(DL_LIBS) $(PCRE_LIBS) $(ORIG_LIBS) 

