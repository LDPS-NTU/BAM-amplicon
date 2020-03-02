#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include "zlib.h"
#include "BamStruct.h"






char	CIGAR(int A);
int	tag(void *stream);

int	BinToEncodedSeq_Top(uint8_t A,int MSB);
int	BinToEncodedSeq(uint8_t A);

void 	PrintAlignmenetHeader(alignmentHeader *AlignmentHeader);
void 	PrintSequence(uint8_t *seq,int32_t seq_length);
void 	PrintCIGAR(uint32_t *cigar,uint16_t n_cigar_op);
void 	PrintPosInformation(posInformation *PosInformation,int position);
int 	EndPosition(int32_t pos, uint16_t n_cigar_op, uint32_t *cigar);
void 	MemoryCopy(uint8_t **stream,void **ptr,int length,int size);
void 	PrintBlank(int size);

void	decompressBlock(z_stream *infstream, uint8_t *stream_i, uint8_t *stream_o, int len_data, FILE *file_bam_i);
int	BGFZBlock(FILE *file_ptr);
int	refInformation(void *stream, char *chr_name, int *chr_length);

//void	alignment_Stat(uint8_t *stream, alignmentHeader *AlignmentHeader, int *pair, int *cross, int *only, int *unmap, char **chr_name);
void 	alignment_List(uint8_t *stream, alignmentHeader *AlignmentHeader, FILE *file_o, char **chr_name);
void	alignment_Coverage(uint8_t *stream, alignmentHeader *AlignmentHeader, posCoverage *PosCoverage, uint32_t ref_length);
void 	alignment_ListAndCoverage(uint8_t *stream, alignmentHeader *AlignmentHeader, FILE *file_o, char **chr_name);
void	alignment_Range(uint8_t *stream, alignmentHeader *AlignmentHeader, uint8_t *in_region, uint64_t *map_in, uint64_t *map_out, uint32_t ref_length);
void	strcat_triple(char	*main, char *first, char *second, char *third, int length);

void	ResetBaseDist(covBaseDistribution *CovBaseDistribution);
void	ResetRegionDist(covRegionDistribution	*CovRegionDistribution);

void	PrintRegionDist (FILE *file_stat_o, FILE *file_ratio_o, bedTable *BedTable, posCoverage *PosCoverage, bamHeader *BamHeader, int ref_ID, int threshold, covRegionDistribution *RegionDist_All, int operation);

void	PrintBaseDist (FILE *file_base, posCoverage *PosCoverage, uint8_t *in_region,char *chr_name, int chr_length, covBaseDistribution *BaseDist_All, int operation);

void	alignment_CoverageTxt(uint8_t *stream, alignmentHeader *AlignmentHeader, posCoverage *PosCoverage, uint32_t *Coverage, uint32_t ref_length);
void	alignment_DepthDist(uint8_t *stream, alignmentHeader *AlignmentHeader, posCoverage *PosCoverage, uint32_t ref_length);

void BuildBaiTable (FILE *file_bai_i, baiTable *BaiTable );
uint64_t ReturnOffset (baiTable *BaiTable, int chr_ID, uint32_t start );
uint64_t CatchOffset (FILE *file_bai_i, int chr_ID, uint32_t start);
uint8_t *CatchBamHeader (FILE *file_bam_i, bamHeader *BamHeader, uint8_t *stream_i, uint8_t *stream_o, uint8_t *buffer, uint8_t *top);

int SelectReads (uint8_t *stream, alignmentHeader *AlignmentHeader, int QualCut, int NumQualCut);
int	decodeBedFile(FILE *file_bed_i, bedTable *BedTable, bamHeader *BamHeader);
int	concernTargetRegion(bedTable *TargetTable, int ref_ID, uint32_t first_pos, uint32_t last_pos);
void	createBedTable(FILE *file_bed_i, bedTable *BedTable, bamHeader *BamHeader, toolsFlags *ToolsFlags);

char Bin2SeqTop(uint8_t A,int MSB);
int32_t	alignment_EndPosition(uint8_t *stream, alignmentHeader *AlignmentHeader, uint32_t ref_length);

void	createRegionTable(FILE *file_bed_i, bedTable *BedTable, bamHeader *BamHeader, toolsFlags *ToolsFlags);

uint8_t *Catch_BamHeader (FILE *file_bam_i, bamHeader *BamHeader, uint8_t *stream_i, uint8_t *stream_o, uint8_t *buffer, toolsFlags *ToolsFlags);
void	alignment_Quality(uint8_t *stream, alignmentHeader *AlignmentHeader, posQuality *PosQuality, uint32_t ref_length);
void	alignment_DepthDist_TR(uint8_t *stream, alignmentHeader *AlignmentHeader, posCoverage *PosCoverage, uint32_t start, uint32_t end);
void	alignment_DepthTxt(uint8_t *stream, alignmentHeader *AlignmentHeader, posCoverage *PosCoverage, uint32_t ref_length);



int	BamDepthDist(FILE *file_bam_i, FILE *file_bai_i, FILE *file_bed_i, char *chromosome, uint32_t start, uint32_t end, toolsFlags *ToolsFlags);
int	BamAmp(FILE *file_bam_i, FILE *file_bai_i, FILE *file_bed_i, char *chr, uint32_t start, uint32_t end, toolsFlags *ToolsFlags);
int	BamSinglePointQuality(FILE *file_bam_i, FILE *file_bai_i, char *chromosome, uint32_t position, toolsFlags *ToolsFlags , char ref);
int	BamStat(FILE *file_bam_i, toolsFlags *ToolsFlags);
int	BamPattern(FILE *file_bam_i, FILE *file_bai_i,  FILE *file_bed_i, char *chromosome, uint32_t start, uint32_t end, toolsFlags *ToolsFlags);
int	BamPoly(FILE *file_bam_i, FILE *file_bai_i,  FILE *file_bed_i, char *chromosome, uint32_t start, uint32_t end, toolsFlags *ToolsFlags);
int	BamMappingLength(FILE *file_bam_i, FILE *file_bai_i, FILE *file_bed_i, char *chr, uint32_t start, uint32_t end, toolsFlags *ToolsFlags);
int	BamTest(FILE *file_bam_i, FILE *file_bai_i, FILE *file_bed_i, FILE *file_fasta_i, FILE *file_fai_i, toolsFlags *ToolsFlags);
int	BamDelTxt(FILE *file_bam_i, FILE *file_bai_i, FILE *file_bed_i, toolsFlags *ToolsFlags);
int	BamTrim(FILE *file_bam_i, FILE *file_bai_i, FILE *file_bed_i, toolsFlags *ToolsFlags);
int	BamDetect(FILE *file_bam_i, FILE *file_bai_i, char *chromosome, toolsFlags *ToolsFlags);
int	BamReadQual(FILE *file_bam_i, int flag_hide);
int	BamPos(FILE *file_bam_i, FILE *file_length_o, int flag_hide , char *dirname_region, int flag_bin);
int	BamQual(FILE *file_bam_i, FILE *file_length_o, int flag_hide , char *dirname_region, int flag_bed);
int	BamCovTxtDel(FILE *file_bam_i, FILE *file_length_o, int flag_hide, char *chr_dir);
int	BamCov(FILE *file_bam_i, FILE *file_length_o, int flag_hide );
int	BamSam(FILE *file_bam_i, FILE *file_length_o, int flag_hide );
int	BamInsTxt(FILE *file_bam_i, FILE *file_length_o, int flag_hide);
int	BamGapTxt(FILE *file_bam_i, FILE *file_length_o, int flag_hide);
