#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include "zlib.h"

#include "BamCommonLibrary.h"




void	PrintBaseDist	(FILE *file_base, posCoverage *PosCoverage, uint8_t *in_region, char *chr_name, int chr_length, covBaseDistribution *BaseDist, int operation);

int	BamRegion(FILE *file_bam_i, int flag_hide , char *dirname_region){

	FILE	*file_count_o;
	FILE	*file_region_i;
	FILE	*file_region_o;
	FILE	*file_base_stat_o;

	FILE	*file_region_stat_0_o;
	FILE	*file_region_stat_10_o;
	FILE	*file_region_stat_20_o;
	FILE	*file_region_stat_30_o;
	FILE	*file_region_stat_40_o;
	FILE	*file_region_stat_50_o;
	FILE	*file_region_stat_60_o;
	FILE	*file_region_stat_70_o;
	FILE	*file_region_stat_80_o;
	FILE	*file_region_stat_90_o;
	FILE	*file_region_stat_100_o;
	
	FILE	*file_region_ratio_0_o;
	FILE	*file_region_ratio_10_o;
	FILE	*file_region_ratio_20_o;
	FILE	*file_region_ratio_30_o;
	FILE	*file_region_ratio_40_o;
	FILE	*file_region_ratio_50_o;
	FILE	*file_region_ratio_60_o;
	FILE	*file_region_ratio_70_o;
	FILE	*file_region_ratio_80_o;
	FILE	*file_region_ratio_90_o;
	FILE	*file_region_ratio_100_o;

	uint8_t *stream_i;//[65536];
	uint8_t *stream_o;//[65536];
	uint8_t	*buffer;
	uint8_t	*top;

	int	i;
	int	len_data;
	int	counter;
	int	ref_ID	= -1;

	uint32_t CRC32;
	uint32_t I_size;

	uint16_t unknown;

	z_stream infstream;
	bgfzTail BgfzTail;
	
	alignmentHeader	AlignmentHeader;
	posCoverage	*PosCoverage;

	char	magic[4];
	char	**chr_name;
	uint32_t	*chr_length;
	int32_t	l_text;
	char	*text;
	int32_t	n_ref;
	uint8_t *address;
	char	filename[name_len];
	char	line[10000];


	uint8_t	*in_region;
	uint64_t	valid_region = 0;
	uint64_t	valid_all_region = 0;
	int	have_region = 0;
	int	start,end;
	uint64_t	map_in;
	uint64_t	map_out;
	uint64_t	map_all_in;
	uint64_t	map_all_out;

	bamHeader	BamHeader;

	covBaseDistribution BaseDist_All;
	covRegionDistribution RegionDist_0_All;
	covRegionDistribution RegionDist_10_All;
	covRegionDistribution RegionDist_20_All;
	covRegionDistribution RegionDist_30_All;
	covRegionDistribution RegionDist_40_All;
	covRegionDistribution RegionDist_50_All;
	covRegionDistribution RegionDist_60_All;
	covRegionDistribution RegionDist_70_All;
	covRegionDistribution RegionDist_80_All;
	covRegionDistribution RegionDist_90_All;
	covRegionDistribution RegionDist_100_All;
	
	
	memset(&BaseDist_All, 0, sizeof(covBaseDistribution));
	memset(&RegionDist_0_All, 0, sizeof(covRegionDistribution));
	memset(&RegionDist_10_All, 0, sizeof(covRegionDistribution));
	memset(&RegionDist_20_All, 0, sizeof(covRegionDistribution));
	memset(&RegionDist_30_All, 0, sizeof(covRegionDistribution));
	memset(&RegionDist_40_All, 0, sizeof(covRegionDistribution));
	memset(&RegionDist_50_All, 0, sizeof(covRegionDistribution));
	memset(&RegionDist_60_All, 0, sizeof(covRegionDistribution));
	memset(&RegionDist_70_All, 0, sizeof(covRegionDistribution));
	memset(&RegionDist_80_All, 0, sizeof(covRegionDistribution));
	memset(&RegionDist_90_All, 0, sizeof(covRegionDistribution));
	memset(&RegionDist_100_All, 0, sizeof(covRegionDistribution));
	

	stream_i= calloc(65536,sizeof(uint8_t));
	stream_o= calloc(65536,sizeof(uint8_t));
	buffer	= calloc(65536*2,sizeof(uint8_t));	

	/* Initial*/
	stream_i[0] = 120;
	stream_i[1] = 156;
	/*---------------*/


	CatchBamHeader (file_bam_i, &BamHeader, stream_i, stream_o, buffer, address, top);
	counter	= top - address;

	if (flag_hide != 1){
		for (i = 0;i < n_ref;i++){
			printf("%s\t%d\n",BamHeader.chr_name[i], BamHeader.chr_length[i]);
		}
	}


	file_region_o 	= fopen("Region_InOut.txt","w");
	file_base_stat_o 	= fopen("Base_stat.txt","w");




	file_region_stat_0_o 	= fopen("Region_stat_0.txt","w");
	file_region_stat_10_o 	= fopen("Region_stat_10.txt","w");
	file_region_stat_20_o 	= fopen("Region_stat_20.txt","w");
	file_region_stat_30_o 	= fopen("Region_stat_30.txt","w");
	file_region_stat_40_o 	= fopen("Region_stat_40.txt","w");
	file_region_stat_50_o 	= fopen("Region_stat_50.txt","w");
	file_region_stat_60_o 	= fopen("Region_stat_60.txt","w");
	file_region_stat_70_o 	= fopen("Region_stat_70.txt","w");
	file_region_stat_80_o 	= fopen("Region_stat_80.txt","w");
	file_region_stat_90_o 	= fopen("Region_stat_90.txt","w");
	file_region_stat_100_o 	= fopen("Region_stat_100.txt","w");
	
	file_region_ratio_0_o 	= fopen("Region_ratio_0.txt","w");
	file_region_ratio_10_o 	= fopen("Region_ratio_10.txt","w");
	file_region_ratio_20_o 	= fopen("Region_ratio_20.txt","w");
	file_region_ratio_30_o 	= fopen("Region_ratio_30.txt","w");
	file_region_ratio_40_o 	= fopen("Region_ratio_40.txt","w");
	file_region_ratio_50_o 	= fopen("Region_ratio_50.txt","w");
	file_region_ratio_60_o 	= fopen("Region_ratio_60.txt","w");
	file_region_ratio_70_o 	= fopen("Region_ratio_70.txt","w");
	file_region_ratio_80_o 	= fopen("Region_ratio_80.txt","w");
	file_region_ratio_90_o 	= fopen("Region_ratio_90.txt","w");
	file_region_ratio_100_o 	= fopen("Region_ratio_100.txt","w");
	


	//Start Unzip and Process Bam File
	while ( (len_data = BGFZBlock(file_bam_i)) > 0){
		decompressBlock(&infstream, stream_i, stream_o, len_data, file_bam_i);
		fread(&BgfzTail	,sizeof(bgfzTail),1,file_bam_i);

		memcpy(buffer,address,sizeof(uint8_t)*counter);
		memcpy(&buffer[counter],stream_o,BgfzTail.I_size);
		top	= buffer + counter + BgfzTail.I_size;
		address = buffer;
		
		if (top - address > sizeof(alignmentHeader)){
			memcpy(&AlignmentHeader, address, sizeof(alignmentHeader));
		}
	
		while ( top - address >= (AlignmentHeader.block_size + 4) ){

			//Chromosome Change
			if (AlignmentHeader.refID != ref_ID){
				if (ref_ID != -1){
					if (flag_hide == 0 ){
						printf("[Bam File Unzip %d / %d ] %s done\n",ref_ID+1, BamHeader.n_ref, BamHeader.chr_name[ref_ID]);
					}
					if (have_region == 1){
		PrintBaseDist (file_base_stat_o, PosCoverage, in_region, chr_name[ref_ID], chr_length[ref_ID], &BaseDist_All , 1);
		PrintRegionDist (file_region_stat_0_o  , file_region_ratio_0_o   ,file_region_i, PosCoverage, &BamHeader, ref_ID, 0  , &RegionDist_0_All  , 1);
		PrintRegionDist (file_region_stat_10_o , file_region_ratio_10_o  ,file_region_i, PosCoverage, &BamHeader, ref_ID, 9  , &RegionDist_10_All , 1);
		PrintRegionDist (file_region_stat_20_o , file_region_ratio_20_o  ,file_region_i, PosCoverage, &BamHeader, ref_ID, 19 , &RegionDist_20_All , 1);
		PrintRegionDist (file_region_stat_30_o , file_region_ratio_30_o  ,file_region_i, PosCoverage, &BamHeader, ref_ID, 29 , &RegionDist_30_All , 1);
		PrintRegionDist (file_region_stat_40_o , file_region_ratio_40_o  ,file_region_i, PosCoverage, &BamHeader, ref_ID, 39 , &RegionDist_40_All , 1);
		PrintRegionDist (file_region_stat_50_o , file_region_ratio_50_o  ,file_region_i, PosCoverage, &BamHeader, ref_ID, 49 , &RegionDist_50_All , 1);
		PrintRegionDist (file_region_stat_60_o , file_region_ratio_60_o  ,file_region_i, PosCoverage, &BamHeader, ref_ID, 59 , &RegionDist_60_All , 1);
		PrintRegionDist (file_region_stat_70_o , file_region_ratio_70_o  ,file_region_i, PosCoverage, &BamHeader, ref_ID, 69 , &RegionDist_70_All , 1);
		PrintRegionDist (file_region_stat_80_o , file_region_ratio_80_o  ,file_region_i, PosCoverage, &BamHeader, ref_ID, 79 , &RegionDist_80_All , 1);
		PrintRegionDist (file_region_stat_90_o , file_region_ratio_90_o  ,file_region_i, PosCoverage, &BamHeader, ref_ID, 89 , &RegionDist_90_All , 1);
		PrintRegionDist (file_region_stat_100_o, file_region_ratio_100_o  ,file_region_i, PosCoverage, &BamHeader, ref_ID, 99, &RegionDist_100_All, 1);
						fclose(file_region_i);
						have_region = 0;
					}
					free(in_region);
					free(PosCoverage);
					fprintf(file_region_o,"%s\t%lu\t%lu\t%lu\n",chr_name[ref_ID],map_in,map_out, valid_region);
					map_all_in += map_in;
					map_all_out += map_out;
					valid_all_region += valid_region;
					valid_region = 0;
					map_in = 0;
					map_out = 0;
					
				}	
				if ( AlignmentHeader.refID != -1){	
					PosCoverage = calloc(chr_length[AlignmentHeader.refID],sizeof(posCoverage));
					in_region = calloc(chr_length[AlignmentHeader.refID],sizeof(uint8_t));
					strcat_triple(filename, dirname_region, chr_name[AlignmentHeader.refID], "_region.txt", name_len);
					if ((file_region_i = fopen(filename,"r")) != NULL){
						have_region = 1;
						while (fscanf(file_region_i,"%*s %d %d",&start,&end)!=EOF){
							fgets( line, 10000, file_region_i);
				//			printf("%s\t%d\t%d\n", chr_name[AlignmentHeader.refID], start, end);
							for (i = start;i < end;i++){
								if (in_region[i] == 0){
									in_region[i] = 1;
									valid_region++;
								}
							}	
						}							
					}
				}
				ref_ID = AlignmentHeader.refID;	
			}
			if (ref_ID == -1){
				break;	
			}
			address += sizeof(alignmentHeader);	
			if ((AlignmentHeader.FLAG&4) == 0 && (AlignmentHeader.FLAG&1024) == 0){
				alignment_DepthTxt(address, &AlignmentHeader, PosCoverage, chr_length[ref_ID]);
				if (have_region == 1){
					alignment_Range(address, &AlignmentHeader, in_region, &map_in, &map_out,chr_length[ref_ID]);
				}else {
					map_out++;	
				}
			}
			address += (AlignmentHeader.block_size + 4) - sizeof(alignmentHeader);

			if (top - address > sizeof(alignmentHeader)){
				memcpy(&AlignmentHeader, address, sizeof(alignmentHeader));
			}else {
				break;	
			}
		}
		counter = top - address;
		if (ref_ID == -1){
			break;	
		}
	}
	
	if (ref_ID != -1){
		if (flag_hide == 0){	printf("[Bam File Unzip %d / %d ] %s done\n",ref_ID+1,n_ref,chr_name[ref_ID]);	}

		if (have_region == 1){
			PrintBaseDist (file_base_stat_o, PosCoverage, in_region, chr_name[ref_ID], chr_length[ref_ID], &BaseDist_All , 1);
			PrintRegionDist (file_region_stat_0_o  , file_region_ratio_0_o   ,file_region_i, PosCoverage, &BamHeader, ref_ID, 0  , &RegionDist_0_All  , 1);
			PrintRegionDist (file_region_stat_10_o , file_region_ratio_10_o  ,file_region_i, PosCoverage, &BamHeader, ref_ID, 9  , &RegionDist_10_All , 1);
			PrintRegionDist (file_region_stat_20_o , file_region_ratio_20_o  ,file_region_i, PosCoverage, &BamHeader, ref_ID, 19 , &RegionDist_20_All , 1);
			PrintRegionDist (file_region_stat_30_o , file_region_ratio_30_o  ,file_region_i, PosCoverage, &BamHeader, ref_ID, 29 , &RegionDist_30_All , 1);
			PrintRegionDist (file_region_stat_40_o , file_region_ratio_40_o  ,file_region_i, PosCoverage, &BamHeader, ref_ID, 39 , &RegionDist_40_All , 1);
			PrintRegionDist (file_region_stat_50_o , file_region_ratio_50_o  ,file_region_i, PosCoverage, &BamHeader, ref_ID, 49 , &RegionDist_50_All , 1);
			PrintRegionDist (file_region_stat_60_o , file_region_ratio_60_o  ,file_region_i, PosCoverage, &BamHeader, ref_ID, 59 , &RegionDist_60_All , 1);
			PrintRegionDist (file_region_stat_70_o , file_region_ratio_70_o  ,file_region_i, PosCoverage, &BamHeader, ref_ID, 69 , &RegionDist_70_All , 1);
			PrintRegionDist (file_region_stat_80_o , file_region_ratio_80_o  ,file_region_i, PosCoverage, &BamHeader, ref_ID, 79 , &RegionDist_80_All , 1);
			PrintRegionDist (file_region_stat_90_o , file_region_ratio_90_o  ,file_region_i, PosCoverage, &BamHeader, ref_ID, 89 , &RegionDist_90_All , 1);
			PrintRegionDist (file_region_stat_100_o, file_region_ratio_100_o  ,file_region_i, PosCoverage, &BamHeader, ref_ID, 99, &RegionDist_100_All, 1);
			have_region = 0;
		}
		free(PosCoverage);
		fprintf(file_region_o,"%s\t%lu\t%lu\t%lu\n",chr_name[ref_ID],map_in,map_out,valid_region);
		map_all_in += map_in;
		map_all_out += map_out;
		valid_all_region += valid_region;
	}
	
	PrintBaseDist (file_base_stat_o, PosCoverage, in_region, "", 0, &BaseDist_All , 2);
	PrintRegionDist (file_region_stat_0_o  , file_region_ratio_0_o   ,file_region_i, PosCoverage, &BamHeader, ref_ID, 0  , &RegionDist_0_All  , 2);
		PrintRegionDist (file_region_stat_10_o , file_region_ratio_10_o  ,file_region_i, PosCoverage, &BamHeader, ref_ID, 9  , &RegionDist_10_All , 2);
		PrintRegionDist (file_region_stat_20_o , file_region_ratio_20_o  ,file_region_i, PosCoverage, &BamHeader, ref_ID, 19 , &RegionDist_20_All , 2);
		PrintRegionDist (file_region_stat_30_o , file_region_ratio_30_o  ,file_region_i, PosCoverage, &BamHeader, ref_ID, 29 , &RegionDist_30_All , 2);
		PrintRegionDist (file_region_stat_40_o , file_region_ratio_40_o  ,file_region_i, PosCoverage, &BamHeader, ref_ID, 39 , &RegionDist_40_All , 2);
		PrintRegionDist (file_region_stat_50_o , file_region_ratio_50_o  ,file_region_i, PosCoverage, &BamHeader, ref_ID, 49 , &RegionDist_50_All , 2);
		PrintRegionDist (file_region_stat_60_o , file_region_ratio_60_o  ,file_region_i, PosCoverage, &BamHeader, ref_ID, 59 , &RegionDist_60_All , 2);
		PrintRegionDist (file_region_stat_70_o , file_region_ratio_70_o  ,file_region_i, PosCoverage, &BamHeader, ref_ID, 69 , &RegionDist_70_All , 2);
		PrintRegionDist (file_region_stat_80_o , file_region_ratio_80_o  ,file_region_i, PosCoverage, &BamHeader, ref_ID, 79 , &RegionDist_80_All , 2);
		PrintRegionDist (file_region_stat_90_o , file_region_ratio_90_o  ,file_region_i, PosCoverage, &BamHeader, ref_ID, 89 , &RegionDist_90_All , 2);
		PrintRegionDist (file_region_stat_100_o, file_region_ratio_100_o  ,file_region_i, PosCoverage, &BamHeader, ref_ID, 99, &RegionDist_100_All, 2);
	fprintf(file_region_o,"Total\t%lu\t%lu\t%lu\t%f\t%f\n",
		map_all_in,
		map_all_out,
		valid_all_region,
		(float)map_all_in /(map_all_in+map_all_out),
		(float)map_all_out /(map_all_in+map_all_out));

	free(stream_i);	
	free(stream_o);	
	free(buffer);	
	return 0;
}
