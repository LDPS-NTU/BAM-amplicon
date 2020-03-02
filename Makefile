#CC=gcc -Wall
CC=gcc
SRC=./src
BIN=./bin


all:
#	$(CC) -lz -lm $(SRC)/BamPoly.c $(SRC)/BamPattern.c $(SRC)/BamAmpSummary.c $(SRC)/BamDepthDist.c $(SRC)/BamStat.c $(SRC)/BamRegionTxt.c  $(SRC)/BamQual.c $(SRC)/BamCovTxt.c $(SRC)/BamCovTxtDel.c $(SRC)/BamCov.c $(SRC)/BamSam.c $(SRC)/BamCommonLibrary.c $(SRC)/Bamtools.c $(SRC)/BamDelTxt.c $(SRC)/BamGapTxt.c $(SRC)/BamPos.c $(SRC)/BamSinglePoint.c $(SRC)/BamSinglePointQuality.c $(SRC)/BamReadQual.c $(SRC)/BamInsTxt.c $(SRC)/BamDetect.c -o $(BIN)/bamtools 
	$(CC) -O2 -lz -lm -o $(BIN)/bam-utility \
		$(SRC)/BamCreateRegionTable.c \
		$(SRC)/BamMappingLength.c \
		$(SRC)/BamDelTxt.c	\
		$(SRC)/BamSinglePointQuality.c \
		$(SRC)/BamPoly.c \
		$(SRC)/BamPattern.c \
		$(SRC)/BamAmpSummary.c \
		$(SRC)/BamDepthDist.c \
		$(SRC)/BamStat.c \
		$(SRC)/BamRegionTxt.c \
		$(SRC)/BamQual.c \
		$(SRC)/BamCommonLibrary.c \
		$(SRC)/BamTest.c \
		$(SRC)/BamTrim.c \
		$(SRC)/Bamtools.c
