CC = g++ -std=c++11
CFLAGS  = -Wall -O3 -std=c++11
LDFLAGS = -lpthread -LBMEAN/spoa/build/lib -lspoa

all: CONSENT-correction CONSENT-polishing explode merge reformatPAF

reformatPAF: reformatPAF.o
	$(CC) -o bin/reformatPAF src/reformatPAF.o

reformatPAF.o: src/reformatPAF.cpp
	$(CC) -o src/reformatPAF.o -c src/reformatPAF.cpp $(CFLAGS)

explode: explode.o
	$(CC) -o bin/explode src/explode.o

explode.o: src/explode.cpp
	$(CC) -o src/explode.o -c src/explode.cpp $(CFLAGS)

merge: merge.o
	$(CC) -o bin/merge src/merge.o

merge.o: src/merge.cpp
	$(CC) -o src/merge.o -c src/merge.cpp $(CFLAGS)

CONSENT-correction: alignmentPiles.o alignmentWindows.o reverseComplement.o utils.o correctionAlignment.o correctionDBG.o correctionMSA.o CONSENT-correction.o DBG.o main.o
	$(CC) -o bin/CONSENT-correction src/main.o src/utils.o src/correctionDBG.o src/correctionAlignment.o src/correctionMSA.o src/reverseComplement.o src/alignmentPiles.o src/alignmentWindows.o src/CONSENT-correction.o src/DBG.o BMEAN/bmean.o BMEAN/utils.o BMEAN/Complete-Striped-Smith-Waterman-Library/src/*.o $(LDFLAGS)

CONSENT-correction.o: src/CONSENT-correction.cpp src/CONSENT-correction.h src/alignmentPiles.h src/alignmentWindows.h src/correctionAlignment.h src/correctionDBG.h src/correctionMSA.h src/DBG.h
	$(CC) -o src/CONSENT-correction.o -c src/CONSENT-correction.cpp $(CFLAGS)

CONSENT-polishing: alignmentPiles.o alignmentWindows.o reverseComplement.o utils.o correctionAlignment.o correctionDBG.o correctionMSA.o CONSENT-polishing.o DBG.o main.o
	$(CC) -o bin/CONSENT-polishing src/main.o src/utils.o src/correctionDBG.o src/correctionAlignment.o src/correctionMSA.o src/reverseComplement.o src/alignmentPiles.o src/alignmentWindows.o src/CONSENT-polishing.o src/DBG.o BMEAN/bmean.o BMEAN/utils.o BMEAN/Complete-Striped-Smith-Waterman-Library/src/*.o $(LDFLAGS)

CONSENT-polishing.o: src/CONSENT-polishing.cpp src/CONSENT-polishing.h src/alignmentPiles.o src/alignmentWindows.h src/correctionAlignment.h src/correctionDBG.h src/correctionMSA.h src/DBG.h
	$(CC) -o src/CONSENT-polishing.o -c src/CONSENT-polishing.cpp $(CFLAGS)

reverseComplement.o: src/reverseComplement.cpp
	$(CC) -o src/reverseComplement.o -c src/reverseComplement.cpp $(CFLAGS)

alignmentPiles.o: src/alignmentPiles.cpp src/Overlap.h src/utils.h
	$(CC) -o src/alignmentPiles.o -c src/alignmentPiles.cpp $(CFLAGS)

alignmentWindows.o: src/alignmentWindows.cpp src/Overlap.h src/reverseComplement.h
	$(CC) -o src/alignmentWindows.o -c src/alignmentWindows.cpp $(CFLAGS)

correctionAlignment.o: src/correctionAlignment.cpp src/utils.h
	$(CC) -o src/correctionAlignment.o -c src/correctionAlignment.cpp $(CFLAGS)

correctionDBG.o: src/correctionDBG.cpp src/utils.h
	$(CC) -o src/correctionDBG.o -c src/correctionDBG.cpp $(CFLAGS)

correctionMSA.o: src/correctionMSA.cpp src/utils.h src/correctionDBG.h
	$(CC) -o src/correctionMSA.o -c src/correctionMSA.cpp $(CFLAGS)

DBG.o: src/DBG.cpp src/reverseComplement.h
	$(CC) -o src/DBG.o -c src/DBG.cpp $(CFLAGS)

utils.o: src/utils.cpp
	$(CC) -o src/utils.o -c src/utils.cpp $(CFLAGS)

main.o: src/main.cpp src/CONSENT-correction.h src/CONSENT-polishing.h
	$(CC) -o src/main.o -c src/main.cpp $(CFLAGS)

clean:
	rm src/*.o bin/CONSENT-correction bin/CONSENT-polishing bin/explode
