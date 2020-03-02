#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include "zlib.h"

#include "BamCommonLibrary.h"
//#include "Function.h"




void 	alignment_List(uint8_t *stream, alignmentHeader *AlignmentHeader, FILE *file_o, char **chr_name);
void	alignment_Coverage(uint8_t *stream, alignmentHeader *AlignmentHeader, posCoverage *PosCoverage, uint32_t ref_length);
void 	alignment_ListAndCoverage(uint8_t *stream, alignmentHeader *AlignmentHeader, FILE *file_o, char **chr_name);
void	strcat_triple(char	*main, char *first, char *second, char *third, int length);


int main(int argc, char *argv[]){

	FILE	*file_bam_i;
	FILE	*file_bai_i;
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

	int	flag_bam=0;
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
	int	start = 0;
	int	end = 0;
	char	chr[1000];
	char	ref;



	if (argc == 1){
		printf("-bam:\tinput bam file\n");
		printf("-length:\toutput chromosome length file\n");
		printf("-header:\tprint header information\n");
		printf("-op:\toperation [ cov, sam , region ]\n");
		printf("-region_dir: dir of region\n");
		return -1;
	}
	memset(dirname_region,0,name_len);
	memset(&ToolsFlags,0,sizeof(toolsFlags));



	for (i=1;i < argc;i++){
		if (strcmp(argv[i],"-bam")==0){
			if ((file_bam_i = fopen(argv[i+1],"rb"))==NULL){
				printf("Error: file [%s] doesn't exist\n",argv[i+1]);
				return -1;
			}
			strcpy(filename_bai, argv[i+1]);
			strcat(filename_bai, ".bai");
			flag_bam=1;
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
		}else if (strcmp(argv[i],"-hide")==0){
			flag_hide=1;
			ToolsFlags.flag_hide = 1;
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
		}else if (strcmp(argv[i],"-chr")==0){
			strcpy(chr, argv[i+1]);
			i++;
		}else if (strcmp(argv[i],"-start")==0){
			start = atoi(argv[i+1]);
			ToolsFlags.start = start;
			i++;
		}else if (strcmp(argv[i],"-end")==0){	
			end = atoi(argv[i+1]);
			ToolsFlags.end = end;
			i++;
		}else if (strcmp(argv[i],"-ref")==0){
			ref = argv[i+1][0];
			ToolsFlags.ref = ref;
		//	strcpy(&ref, argv[i+1]);
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


	if(flag_stat == 1){	BamStat(file_bam_i, file_length_o, flag_hide);	return 0;
	}else if (flag_region == 1){		BamRegion(file_bam_i, flag_hide, dirname_region, 0);	return 0;
	}else if (flag_regiontxt == 1){
		if ((file_bai_i = fopen(filename_bai, "rb"))==NULL){
			printf("Error: file [%s] doesn't exist\n", filename_bai);
			return -1;
		}
		BamRegionTxt(file_bam_i, file_bai_i, file_length_o, flag_hide, chr, start, end);	
		return 0;
	}else if (flag_detect == 1){
		if ((file_bai_i = fopen(filename_bai, "rb"))==NULL){
			printf("Error: file [%s] doesn't exist\n", filename_bai);
			return -1;
		}
		ToolsFlags.length = ToolsFlags.end - ToolsFlags.start;
		BamDetect(file_bam_i, file_bai_i, chr, &ToolsFlags);	
		return 0;
	}else if (flag_singlepoint == 1){
		if ((file_bai_i = fopen(filename_bai, "rb"))==NULL){
			printf("Error: file [%s] doesn't exist\n", filename_bai);
			return -1;
		}
		start--;
		BamSinglePoint(file_bam_i, file_bai_i, chr, start, &ToolsFlags);	
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
	}else if (flag_covtxt == 1){	BamCovTxt(file_bam_i, file_length_o, flag_hide);	return 0;
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



	return 0;
}
