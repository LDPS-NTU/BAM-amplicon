#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <getopt.h> 
#include "zlib.h"




#include "BamCommonLibrary.h"
//#include "Function.h"




void 	alignment_List(uint8_t *stream, alignmentHeader *AlignmentHeader, FILE *file_o, char **chr_name);
void	alignment_Coverage(uint8_t *stream, alignmentHeader *AlignmentHeader, posCoverage *PosCoverage, uint32_t ref_length);
void 	alignment_ListAndCoverage(uint8_t *stream, alignmentHeader *AlignmentHeader, FILE *file_o, char **chr_name);
void	strcat_triple(char	*main, char *first, char *second, char *third, int length);

void	usage(){
	printf("Program: bam-utility (for BAM file)\n");
	printf("Version: 1.0\n");
	printf("Usage:   bam-utility -m <Mode> -b <bam_filename> [Options]\n");
	printf("\n");
	//printf("Mode:	depthdist	about depth's distribution\n");
	printf("Mode:	depthdist	show depths\n");
	//printf("	ampsummary	about amplicon's depth and cover rate\n");
	printf("	ampsummary	show depths and coverage rates within specified amplicons\n");
	printf("	stat		show BAM statistics\n");
	printf("	pattern		show frequency of patterns within a specified region\n");
	printf("	quality		show distribution of quality scores for a specified position\n");
	printf("	poly		show read sequences with homopolymers\n");
	printf("	length		show distribution of mapped read lenghts\n");
	printf("	del		show frequency of deletions with position info\n");
	printf("	trim		trim the BAM file using amplicons provided by the BED file\n");
//	printf("	poly		\n");
	printf("\n");

	printf("Options:\n");
	printf("	-b, --bam    [FILE]	input BAM file\n");
	printf("	-r, --bed    [FILE]	the range of specified regions from BED file\n");
	printf("	-c, --chr    [STR]	chromosome name\n");
	printf("	-s, --start  [INT]	start position (1-based, i.e. first base is 1)\n");
	printf("	-e, --end    [INT]	end position (1-based, i.e. first base is 1)\n");
	printf("	-g, --fasta  [FILE]	fasta file (reference file)\n");
	printf("	-t, --target [FILE]	the ragne of amplicon regions from BED file\n");
	printf("	-f, --filter [TYPE]	filter criteria\n");
	printf("				[readqual, mapqual]\n");
	printf("	-l, --compact		show only average depth and coverage rate (depthdist mode) \n");
	printf("	-a, --origin		show original sequence (include insertion)\n");
	printf("	-d, --duplicate		show duplicate\n");
	printf("	-v, --verbose		show the processed status \n");
	printf("	-n, --column_name	remove the first row (header)\n");
	printf("	-u, --threshold	[INT]	coverage above the threshold\n");
//	printf("	/* For Ion Torrent Data */\n");
//	printf("	-w, --flow		modify the sequence by flow signals (FZ)\n");
//	printf("	-z, --zm		change the FZ values to ZM values\n");

	printf("	-p, --pattern	[STR]	show the percentage of ALT pattern specified (pattern mode only)\n");
	printf("	-h, --help		show the manual page\n");
	printf("	\n");
}


int main(int argc, char *argv[]){

	FILE	*file_length_o;
	FILE	**file_sam_array_o;
	FILE	**file_cov_array_o;
	FILE	*file_region;

	uint8_t *stream_i;//[65536];
	uint8_t *stream_o;//[65536];
	uint8_t	*buffer;
	uint8_t	*top;

	int	i;
	int	len_data;
	int	counter;
	int	ref_ID	= 0;

	uint32_t CRC32;
	uint32_t I_size;

	uint16_t unknown;


//New Variable
//
	FILE	*file_bam_i = NULL;
	FILE	*file_bai_i = NULL;
	FILE	*file_bed_i = NULL;
	FILE	*file_fasta_i = NULL;
	FILE	*file_fai_i = NULL;

	char	*mode;
	char	*bam;
	char	*bai_1;
	char	*bai_2;
	char	*fasta;
	char	*fai;
	char	*bed;
	char	*target;
	char	*dir;
	char	*chromosome;
	uint32_t	start = 0;
	uint32_t	end = 0;


	int	flag_bam;
	int	flag_fasta = 0;
	int	flag_bed;
	int	flag_dir;
	int	flag_start	= 0;
	int	flag_end	= 0;

//END



	z_stream infstream;
	bgfzTail BgfzTail;
	
	alignmentHeader	AlignmentHeader;
	posCoverage	*PosCoverage;
	toolsFlags	ToolsFlags;



	char	magic[4];
	char	**chr_name;
	int	*chr_length;
	int32_t	l_text;
	char	*text;
	int32_t	n_ref;
	uint8_t *address;
	char	filename[name_len];
	char	filename_bai[name_len];
	char	dirname_region[name_len];
	char	dirname_ref[name_len];

	int	flag_length=0;
	int	flag_header=0;
	int	flag_cov = 0;
	int	flag_covtxt = 0;
	int	flag_deltxt = 0;
	int	flag_instxt = 0;
	int	flag_gaptxt = 0;
	int	flag_detect = 0;
	int	flag_covtxtdel = 0;
	int	flag_sam = 0;
	int	flag_region = 0;
	int	flag_regiontxt = 0;
	int	flag_singlepoint = 0;
	int	flag_singlequality = 0;
	int	flag_readqual = 0;
	int	flag_pos = 0;
	int	flag_stat = 0;
	int	flag_qual = 0;
	int	flag_hide = 0;
	int	flag_op = 0;
	int	flag_dirname_region= 0;
	int	flag_dirname_ref= 0;

	uint8_t	*in_region;
	int	valid_region;
	int	valid_all_region = 0;
	int	have_region = 0;
	int	map_in = 0;
	int	map_out = 0;
	int	map_all_in = 0;
	int	map_all_out = 0;
	int	pair = 0;

	int	cross = 0;
	int	only = 0;
	int	unmap = 0;
	char	chr[1000];
	char	ref;

	const char* const short_options = "m:b:r:c:s:e:g:t:o:f:p:u:vndhlaMzw";
	struct option long_options[] = {
		{"mode"	,	1,	NULL,	'm'},
		{"bam"	,	1,	NULL,	'b'},
		{"bed"	,	1,	NULL,	'r'},
		{"chr"	,	1,	NULL,	'c'},
		{"start",	1,	NULL,	's'},
		{"end",		1,	NULL,	'e'},
		{"target",	1,	NULL,	't'},
		{"fasta",	1,	NULL,	'g'},
		{"output",	1,	NULL,	'o'},
		{"filter",	1,	NULL,	'f'},
		{"pattern",	1,	NULL,	'p'},
		{"threshold",	1,	NULL,	'u'},
		{"verbose",	0,	NULL,	'v'},
		{"column_name",	0,	NULL,	'n'},
		{"compact",	0,	NULL,	'l'},
		{"origin",	0,	NULL,	'a'},
		{"duplicate",	0,	NULL,	'd'},
		{"mapq",	0,	NULL,	'M'},
		{"zm",		0,	NULL,	'z'},
		{"flow",	0,	NULL,	'w'},
		{"help",	0,	NULL,	'h'}
	};

	int	option_index=0;
	int	c;

	if (argc == 1){
		usage();
		return -1;
	}

	flag_bam	= 0;
	flag_hide 	= 0;
	flag_bed	= 0;


	//ALL Flag is Zero initially

	memset(&ToolsFlags,0,sizeof(toolsFlags));


	while ((c = getopt_long (argc, argv, short_options, long_options, NULL)) != -1){
		switch (c){
			case 'm':	
				mode = strdup(optarg);	
				//printf("mode:\t%s\n",optarg);	
				break;
			case 'b':
				//Check Bam & Bai
				bam = strdup(optarg);
				bai_1 = calloc(strlen(bam)+6,sizeof(char));
				strcpy(bai_1, optarg);	strcat( bai_1, ".bai");
				bai_2 = strdup(optarg);	bai_2[strlen(bai_2)-1] = 'i';

				if ( (access(bam, R_OK) != -1) && ((access(bai_1, R_OK) != -1) || (access(bai_2, R_OK) != -1))){	
					flag_bam = 1;
				}else {
					if (access(bam, R_OK) == -1){
						printf("bam:\t%s\n",bam);	
						printf("Please Check Bam File.\n");
					}else {
						printf("bai_1:\t%s\n",bai_1);	
						printf("bai_2:\t%s\n",bai_2);
						printf("Please Check Bai File.\n");
					}
					return -1;
				}
				break;
			case 'g':
				fasta	= strdup(optarg);
				fai	= calloc(strlen(fasta)+6,sizeof(char));
				strcpy(fai, optarg);	strcat(fai, ".fai");
				if ( (access(fasta, R_OK) != -1) && (access(fai, R_OK) != -1) ){	
					flag_fasta = 1;
				}else {
					if (access(fasta, R_OK) == -1){
						printf("fasta:\t%s\n", fasta);	
						printf("Please Check Fasta File.\n");
					}else {
						printf("fai:\t%s\n", fai);	
						printf("Please Check Fai File.\n");
					}
				}
				break;
			case 'r':	
				//Check Bed File
				bed = strdup(optarg);
				if ( access(bed, R_OK) != -1){
					//printf("bed:\t%s\n",bed);
					flag_bed = 1;
					/*Pattern*/
					ToolsFlags.flag_pattern = 1;
				}else {
					printf("Please Check Bed File.\n");
					return -1;
				}
				break;
			case 'c':
				chromosome = strdup(optarg);
				ToolsFlags.chromosome = strdup(optarg);
				//printf("chr:\t%s\n",optarg);
				break;
			case 's':	
				start = atoi(optarg);
				ToolsFlags.start = atoi(optarg);
				//printf("start:\t%s\n",optarg);	
				flag_start = 1;
				if (start < 1){
					printf("\n[Error!]\tStart Point:[%d] is smaller than 1.\n\n",start);
					usage();
					return -1;
				}
				break;
			case 'e':
				end = atoi(optarg);	
				ToolsFlags.end	= atoi(optarg);
				flag_end = 1;
				//printf("end:\t%s\n",optarg);	
				break;
			case 't':	
				target = strdup(optarg);
				if ( access(target, R_OK) != -1){
					ToolsFlags.flag_target = 1;
				}else {
					printf("Please Check Bed File.\n");
					return -1;
				}
				break;
			case 'o':	printf("output:\t%s\n",optarg);	break;
			case 'f':	
				ToolsFlags.flag_filter = 1;
				//printf("filter:\t%s\n",optarg);	
				break;
			case 'p':
				ToolsFlags.flag_pattern = 1;
				ToolsFlags.pattern = strdup(optarg);
				break;
			case 'v':	
				flag_hide = 1;
				ToolsFlags.flag_hide = 1;
				break;
			case 'n':	
				ToolsFlags.flag_columnName = 1;
				break;
			case 'l':	
				ToolsFlags.flag_simple = 1;
				break;
			case 'a':	
				ToolsFlags.flag_origin = 1;
				break;
			case 'd':	
				ToolsFlags.flag_dup = 1;
				break;
			case 'M':	
				ToolsFlags.flag_mapq = 1;
				break;
			case 'u':
				ToolsFlags.flag_coverage = 1;
				ToolsFlags.coverage_threshold = strdup(optarg);
				break;
			case 'w':
				ToolsFlags.flag_flow	= 1;
				break;
			case 'z':
				ToolsFlags.flag_zm	= 1;
				break;
			case 'h':	
				usage();
				return 0;
			default: 	
				printf("Error: Parameter\n");	
				return -1;
		}
	}
	//Piroity 
	//Bed File > Region > All
	
	if (flag_bam == 1){	
		file_bam_i = fopen(bam,"rb");	
		if (access(bai_1, R_OK) != -1){	
			file_bai_i = fopen(bai_1,"rb");	
			//printf("bai_1:\t%s\n",bai_1);
		} else if (access(bai_2, R_OK) != -1){			
			file_bai_i = fopen(bai_2,"rb");	
			//printf("bai_2:\t%s\n",bai_2);
		}else {
			usage();
			return -1;
		}
	} else{	
		usage();
		return -1;	
	}
	if (flag_fasta == 1){
		file_fasta_i	= fopen(fasta, "r");	
		file_fai_i	= fopen(fai, "r");	
	}



	if (flag_bed == 1){	
		file_bed_i = fopen(bed,"r");
	}

	if (ToolsFlags.flag_target == 1){
		ToolsFlags.file_target = fopen(target, "r");
	}

	if (flag_start == 1){
		if (ToolsFlags.chromosome == NULL){
			printf("\n[Error!]\tChromosome is missed.\n");
			return -1;
		}
	}

	if (flag_end == 1){
		if (ToolsFlags.chromosome == NULL){
			printf("\n[Error!]\tChromosome is missed.\n");
			return -1;
		}
 		if (flag_start == 0){
			printf("\n[Error!]\tStart position is missed.\n");
			return -1;
		}
	}


	//start--;
	if (strcmp ( mode, "depthdist")==0){		BamDepthDist(file_bam_i, file_bai_i, file_bed_i, chromosome, start, end, &ToolsFlags);
	}else if (strcmp ( mode, "indel")==0){
	}else if (strcmp ( mode, "ampsummary")==0){	BamAmp(file_bam_i, file_bai_i, file_bed_i, chromosome, start, end, &ToolsFlags);
	}else if (strcmp ( mode, "quality")==0){	BamSinglePointQuality(file_bam_i, file_bai_i, chromosome, start, &ToolsFlags, 'N');
	}else if (strcmp ( mode, "stat")==0){		BamStat(file_bam_i, &ToolsFlags);
	}else if (strcmp ( mode, "pattern")==0){	BamPattern(file_bam_i, file_bai_i, file_bed_i, chromosome, start, end, &ToolsFlags);
	}else if (strcmp ( mode, "poly")==0){		BamPoly(file_bam_i, file_bai_i, file_bed_i, chromosome, start, end, &ToolsFlags);
	}else if (strcmp ( mode, "length")==0){		BamMappingLength(file_bam_i, file_bai_i, file_bed_i, chromosome, start, end, &ToolsFlags);
	}else if (strcmp ( mode, "del")==0){		BamDelTxt(file_bam_i, file_bai_i, file_bed_i, &ToolsFlags);
	}else if (strcmp ( mode, "trim")==0){		
		if (ToolsFlags.flag_target == 0){	usage(); return -1;}
		BamTrim(file_bam_i, file_bai_i, file_bed_i, &ToolsFlags);
	}else if (strcmp ( mode, "iontorrent")==0){	
		if (flag_fasta == 0){	usage(); return -1;}
		if (ToolsFlags.flag_zm == 0 && ToolsFlags.flag_flow == 0){	usage(); return -1;}
		ToolsFlags.flag_header	= 1;	BamTest(file_bam_i, file_bai_i, file_bed_i, file_fasta_i, file_fai_i, &ToolsFlags);
	}else {
		printf("[Warning]No mode or mode is error\n");
		usage();
		return -1;
	}


	if (flag_bam == 1){	fclose(file_bam_i);	fclose(file_bai_i);}
	if (flag_fasta == 1){	fclose(file_fasta_i);	fclose(file_fai_i);}
	if (flag_bed == 1){	fclose(file_bed_i);	}

	return 0;

	////////////////////////////////////////////////////////////////////////////
/*
	memset(dirname_region,0,name_len);

	for (i=1;i < argc;i++){
		if (strcmp(argv[i],"-bam")==0){
			i++;
		}else if (strcmp(argv[i],"-length")==0){
			if ((file_length_o = fopen(argv[i+1],"r"))!=NULL){
				printf("Error: file [%s] does exist\n",argv[i+1]);
				fclose(file_length_o);
				return -1;
			}else {
				file_length_o = fopen(argv[i+1],"w");
				flag_length=1;
				i++;
			}
		}else if (strcmp(argv[i],"-ref_dir")==0){
			if (access(argv[i+1],0) < 0){
				printf("Reference Dir doesn't exist\n");
				return -1;
			}else {
				strcpy(dirname_ref, argv[i+1]);
				flag_dirname_ref = 1;
			}
			i++;
		}else if (strcmp(argv[i],"-region_dir")==0){
			if (access(argv[i+1],0) < 0){
				printf("Region Dir doesn't exist\n");
				return -1;
			}else {
				strcpy(dirname_region, argv[i+1]);
				flag_dirname_region = 1;
			}
			i++;
		}else if (strcmp(argv[i],"-header")==0){
			flag_header=1;
		}else if (strcmp(argv[i],"-rmdup")==0){
			ToolsFlags.flag_dup = 1;
		}else if (strcmp(argv[i],"-op")==0){
			if (flag_op == 0){
				flag_op = 1;	
			}else {
				printf("Error: Multi operation\n",argv[i+1]);
				return -1;		
			}
			
			if (strcmp(argv[i+1],"cov")==0)	{		flag_cov = 1;
			}else if(strcmp(argv[i+1],"covtxt")==0){	flag_covtxt = 1;
			}else if(strcmp(argv[i+1],"deltxt")==0){	flag_deltxt = 1;
			}else if(strcmp(argv[i+1],"instxt")==0){	flag_instxt = 1;
			}else if(strcmp(argv[i+1],"gaptxt")==0){	flag_gaptxt = 1;
			}else if(strcmp(argv[i+1],"detect")==0){	flag_detect = 1;
			}else if(strcmp(argv[i+1],"covtxtdel")==0){	flag_covtxtdel = 1;
			}else if(strcmp(argv[i+1],"sam")==0){		flag_sam = 1;
			}else if(strcmp(argv[i+1],"region")==0){	flag_region = 1;
			}else if(strcmp(argv[i+1],"regiontxt")==0){	flag_regiontxt = 1;
			}else if(strcmp(argv[i+1],"singlepoint")==0){	flag_singlepoint = 1;
			}else if(strcmp(argv[i+1],"singlequality")==0){	flag_singlequality = 1;
			}else if(strcmp(argv[i+1],"readqual")==0){	flag_readqual = 1;
			}else if(strcmp(argv[i+1],"pos")==0){		flag_pos = 1;
			}else if(strcmp(argv[i+1],"stat")==0){		flag_stat = 1;
			}else if(strcmp(argv[i+1],"qual")==0){		flag_qual = 1;
			}else{
				printf("-op: %s doesn't exist\n",argv[i+1]);
				return -1;	
			}
			i++;
		}else if (strcmp(argv[i],"-alt")==0){
			ToolsFlags.alt = argv[i+1][0];
			i++;
		}else if (strcmp(argv[i],"-type")==0){
			if (strcmp(argv[i+1],"snp")==0){	ToolsFlags.type = 0;	
			}else if (strcmp(argv[i+1],"del")==0){	ToolsFlags.type = 1;	
			}else if (strcmp(argv[i+1],"ins")==0){	ToolsFlags.type = 2;
			}else {	printf("Error\n");	return -1;}
			i++;
		}else {
			printf("-bam:\tinput bam file\n");
			printf("-length:\toutput chromosome length file\n");
			printf("-header:\tprint header information\n");
			printf("-hide:\tprint process information\n");
			printf("-op:\toperation [ cov, sam ,region ]\n");
			printf("-region_dir: dir of region\n");
			return -1;
		}
	}

	if (flag_bam == 0){
		printf("Error: No bam file\n");
		return -1;
	}else if (flag_region == 1 && flag_dirname_region == 0){
		printf("Error No Dirname Region\n");
		return -1;
	}else if (flag_covtxtdel == 1 && flag_dirname_ref == 0){
		printf("Error No Dirname of Reference\n");
		return -1;
	}else if (flag_op == 0){
		printf("Error: No opertation\n");
		return -1;	
	}


	if (flag_regiontxt == 1){
		return 0;
	}else if (flag_detect == 1){
		if ((file_bai_i = fopen(filename_bai, "rb"))==NULL){
			printf("Error: file [%s] doesn't exist\n", filename_bai);
			return -1;
		}
		ToolsFlags.length = ToolsFlags.end - ToolsFlags.start;
		BamDetect(file_bam_i, file_bai_i, chr, &ToolsFlags);	
		return 0;
	}else if (flag_singlequality == 1){
		if ((file_bai_i = fopen(filename_bai, "rb"))==NULL){
			printf("Error: file [%s] doesn't exist\n", filename_bai);
			return -1;
		}
		start--;
		BamSinglePointQuality(file_bam_i, file_bai_i, chr, start, &ToolsFlags, ref);	
		return 0;
	}else if (flag_readqual == 1){		BamReadQual(file_bam_i, flag_hide);	return 0;
	}else if (flag_pos == 1){	BamPos(file_bam_i, file_length_o, flag_hide, dirname_region, 0);	return 0;
	}else if (flag_qual == 1){	BamQual(file_bam_i, file_length_o, flag_hide, dirname_region, 0);	return 0;
	}else if (flag_covtxtdel == 1){	BamCovTxtDel(file_bam_i, file_length_o, flag_hide, dirname_ref);	return 0;
	}else if (flag_cov == 1){	BamCov(file_bam_i, file_length_o, flag_hide);	return 0;
	}else if (flag_sam == 1){	BamSam(file_bam_i, file_length_o, flag_hide);	return 0;
	}else if (flag_deltxt == 1){	BamDelTxt(file_bam_i, file_length_o, flag_hide);	return 0;
	}else if (flag_instxt == 1){	BamInsTxt(file_bam_i, file_length_o, flag_hide);	return 0;
	}else if (flag_gaptxt == 1){	BamGapTxt(file_bam_i, file_length_o, flag_hide);	return 0;
	}else {
		printf("Error\n");
		return -1;
	}
*/


	return 0;
}
