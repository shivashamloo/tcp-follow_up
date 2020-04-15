all: secator

OPT = -O3
COMP = gcc

secator: main.o lecfi.o ecrfi.o secatorbionj.o dissimilarities_clustering.o weighting.o treecut.o \
	alignment.o bionj.o treeward.o wardbionj.o treebionj.o tools.o
	$(COMP) -o secator main.o lecfi.o ecrfi.o secatorbionj.o dissimilarities_clustering.o \
	weighting.o treecut.o alignment.o bionj.o treeward.o wardbionj.o treebionj.o tools.o -lm $(OPT)

clean:
	rm -f *.o

main.o: main.c main.h lecfi.h ecrfi.h divide.h alignment.h treebionj.h secatorbionj.h treeward.h treecut.h
	$(COMP) -o main.o -c main.c $(OPT)

lecfi.o: lecfi.c main.h lecfi.h  
	$(COMP) -o lecfi.o -c lecfi.c $(OPT)

ecrfi.o: ecrfi.c main.h ecrfi.h  
	$(COMP) -o ecrfi.o -c ecrfi.c $(OPT)

secatorbionj.o: secatorbionj.c dissimilarities_clustering.h weighting.h treecut.h
	$(COMP) -o secatorbionj.o -c secatorbionj.c $(OPT)

dissimilarities_clustering.o: dissimilarities_clustering.c main.h dissimilarities_clustering.h tools.h
	$(COMP) -o dissimilarities_clustering.o -c dissimilarities_clustering.c $(OPT)

weighting.o: weighting.c main.h weighting.h
	$(COMP) -o weighting.o -c weighting.c $(OPT)

treecut.o: treecut.c main.h treecut.h dissimilarities_clustering.h tools.h
	$(COMP) -o treecut.o -c treecut.c $(OPT)

alignment.o: alignment.c main.h 
	$(COMP) -o alignment.o -c alignment.c $(OPT)

bionj.o: bionj.c bionj.h 
	$(COMP) -o bionj.o -c bionj.c $(OPT)

treeward.o: treeward.c main.h treeward.h
	$(COMP) -o treeward.o -c treeward.c $(OPT)

wardbionj.o: wardbionj.c main.h wardbionj.h
	$(COMP) -o wardbionj.o -c wardbionj.c $(OPT)

treebionj.o: treebionj.c main.h treebionj.h bionj.h
	$(COMP) -o treebionj.o -c treebionj.c $(OPT)

tools.o: tools.c main.h tools.h 
	$(COMP) -o tools.o -c tools.c $(OPT)

