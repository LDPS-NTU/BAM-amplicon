#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include "zlib.h"

#include "BamCommonLibrary.h"




int	BamAmp(FILE *file_bam_i, FILE *file_bai_i, FILE *file_bed_i, char *chr, uint32_t start, uint32_t end, toolsFlags *ToolsFlags){

	FILE	*file_count_o;
	FILE	*file_region_o;
	FILE	*file_base_stat_o;

	FILE	*file_stat_0_o;
	FILE	*file_stat_10_o;
	FILE	*file_stat_20_o;
	FILE	*file_stat_30_o;
	FILE	*file_stat_40_o;
	FILE	*file_stat_50_o;
	FILE	*file_stat_60_o;
	FILE	*file_stat_70_o;
	FILE	*file_stat_80_o;
	FILE	*file_stat_90_o;
	FILE	*file_stat_100_o;
	FILE	*file_stat_200_o;
	FILE	*file_stat_500_o;
	
	FILE	*file_ratio_0_o;
	FILE	*file_ratio_10_o;
	FILE	*file_ratio_20_o;
	FILE	*file_ratio_30_o;
	FILE	*file_ratio_40_o;
	FILE	*file_ratio_50_o;
	FILE	*file_ratio_60_o;
	FILE	*file_ratio_70_o;
	FILE	*file_ratio_80_o;
	FILE	*file_ratio_90_o;
	FILE	*file_ratio_100_o;
	FILE	*file_ratio_200_o;
	FILE	*file_ratio_500_o;

	uint8_t *stream_i;//[65536];
	uint8_t *stream_o;//[65536];
	uint8_t	*buffer;
	uint8_t *address;
	uint8_t	*top;
	
	int	i,j,m;
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



	char	filename[name_len];
	char	line[LINE_MAX_LEN];


	uint8_t	*on_target;
	uint64_t	valid_region = 0;
	uint64_t	valid_all_region = 0;
	int	chr_in_target = 0;
	uint64_t	map_in = 0;
	uint64_t	map_out = 0;
	uint64_t	map_all_in = 0;
	uint64_t	map_all_out = 0;

	bamHeader	BamHeader;
	bedTable	BedTable;

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
	covRegionDistribution RegionDist_200_All;
	covRegionDistribution RegionDist_500_All;

	uint32_t	num_target = 0;	

	char	chromosome_name[name_len];
	uint64_t	offset_beg;
	uint64_t	offset_bgzf;
	uint64_t	offset_decomp;
	int	index;






/////////////////////////////////////////////////////
	//Start 
	////Bam Header
	////
	stream_i= calloc(65536,sizeof(uint8_t));
	stream_o= calloc(65536,sizeof(uint8_t));
	buffer	= calloc(65536*3,sizeof(uint8_t));	

	/* Initial*/
	stream_i[0] = 120;
	stream_i[1] = 156;
	/*---------------*/
	
	top = CatchBamHeader (file_bam_i, &BamHeader, stream_i, stream_o, buffer, top);
	counter	= top - buffer;
	address = buffer;

	/* comment on 2020/03/20 YC
	if (ToolsFlags->flag_hide == 1){
		for (i = 0;i < BamHeader.n_ref;i++){
			printf("%s\t%d\n",BamHeader.chr_name[i], BamHeader.chr_length[i]);
		}
	}
	*/

	//Bed File
	if (file_bed_i != NULL){
		while	(fgets(line, LINE_MAX_LEN, file_bed_i) != NULL){
			num_target++;
		}
		//printf("%d\n",num_target);

		BedTable.table_start = calloc (BamHeader.n_ref, sizeof(uint32_t));
		BedTable.table_end = calloc (BamHeader.n_ref, sizeof(uint32_t));
    	
		BedTable.start = calloc (num_target, sizeof(uint32_t));
		BedTable.end = calloc (num_target, sizeof(uint32_t));
		
		rewind(file_bed_i);
		
		for(i = 0;i < num_target;i++){
			fscanf(file_bed_i, "%s %d %d", chromosome_name, &BedTable.start[i], &BedTable.end[i]);
			fgets(line, LINE_MAX_LEN, file_bed_i);	
			
			if (i == 0){
				for(j = 0; j < BamHeader.n_ref; j++){
					if (strcmp(chromosome_name, BamHeader.chr_name[j])==0){
						BedTable.table_start[j] = i;
						break;
					}
				}
			}else if (strcmp(chromosome_name, BamHeader.chr_name[j])!=0){
				BedTable.table_end[j] = i;
				for(j = 0; j < BamHeader.n_ref; j++){
					if (strcmp(chromosome_name, BamHeader.chr_name[j])==0){
						BedTable.table_start[j] = i;
						break;
					}
				}
			}
			if (j == BamHeader.n_ref){
				printf("Error\n");
				return -1;
			}
			//printf("%s %d %d\n", chromosome_name, BedTable.start[i], BedTable.end[i]);
		}
		BedTable.table_end[j] = i;
		////////////////////////////////////////////
	}else {
		printf("Total\n");
		BedTable.table_start = calloc (BamHeader.n_ref, sizeof(uint32_t));
		BedTable.table_end = calloc (BamHeader.n_ref, sizeof(uint32_t));
		BedTable.table_max = calloc (BamHeader.n_ref, sizeof(uint32_t));
		
		BedTable.start = calloc (BamHeader.n_ref, sizeof(uint32_t));
		BedTable.end = calloc (BamHeader.n_ref, sizeof(uint32_t));
		num_target = BamHeader.n_ref;
			
		for(j = 0; j < BamHeader.n_ref; j++){
				BedTable.start[j] = 0;
				BedTable.end[j] = BamHeader.chr_length[j];
				BedTable.table_start[j] = j;
				BedTable.table_end[j] = j+1;
				BedTable.table_max[j] = BamHeader.chr_length[j];

		}				
	}
//////////////////////////////////////////////////////////////////////////////////////////












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
	memset(&RegionDist_200_All, 0, sizeof(covRegionDistribution));
	memset(&RegionDist_500_All, 0, sizeof(covRegionDistribution));
	





	file_region_o 	= fopen("Region_InOut.txt","w");
	file_base_stat_o 	= fopen("Base_stat.txt","w");




	file_stat_0_o 	= fopen("Region_stat_0.txt","w");
	file_stat_10_o 	= fopen("Region_stat_10.txt","w");
	file_stat_20_o 	= fopen("Region_stat_20.txt","w");
	file_stat_30_o 	= fopen("Region_stat_30.txt","w");
	file_stat_40_o 	= fopen("Region_stat_40.txt","w");
	file_stat_50_o 	= fopen("Region_stat_50.txt","w");
	file_stat_60_o 	= fopen("Region_stat_60.txt","w");
	file_stat_70_o 	= fopen("Region_stat_70.txt","w");
	file_stat_80_o 	= fopen("Region_stat_80.txt","w");
	file_stat_90_o 	= fopen("Region_stat_90.txt","w");
	file_stat_100_o 	= fopen("Region_stat_100.txt","w");
	file_stat_200_o 	= fopen("Region_stat_200.txt","w");
	file_stat_500_o 	= fopen("Region_stat_500.txt","w");
	
	file_ratio_0_o 	= fopen("Region_ratio_0.txt","w");
	file_ratio_10_o 	= fopen("Region_ratio_10.txt","w");
	file_ratio_20_o 	= fopen("Region_ratio_20.txt","w");
	file_ratio_30_o 	= fopen("Region_ratio_30.txt","w");
	file_ratio_40_o 	= fopen("Region_ratio_40.txt","w");
	file_ratio_50_o 	= fopen("Region_ratio_50.txt","w");
	file_ratio_60_o 	= fopen("Region_ratio_60.txt","w");
	file_ratio_70_o 	= fopen("Region_ratio_70.txt","w");
	file_ratio_80_o 	= fopen("Region_ratio_80.txt","w");
	file_ratio_90_o 	= fopen("Region_ratio_90.txt","w");
	file_ratio_100_o 	= fopen("Region_ratio_100.txt","w");
	file_ratio_200_o 	= fopen("Region_ratio_200.txt","w");
	file_ratio_500_o 	= fopen("Region_ratio_500.txt","w");
	


	//Start Unzip and Process Bam File
	while ( (len_data = BGFZBlock(file_bam_i)) > 0){
		decompressBlock(&infstream, stream_i, stream_o, len_data, file_bam_i);
		fread(&BgfzTail	,sizeof(bgfzTail),1,file_bam_i);
		

//		printf("Johnnash\t%u\t%u\t%u\t%d\n", buffer, address, top, counter );
//		getchar();	
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
					if (ToolsFlags->flag_hide == 1 ){
						printf("[Bam File Unzip %d / %d ] %s done\n",ref_ID+1, BamHeader.n_ref, BamHeader.chr_name[ref_ID]);
					}
					
					if (chr_in_target == 1){
						PrintBaseDist (file_base_stat_o, PosCoverage, on_target, BamHeader.chr_name[ref_ID], BamHeader.chr_length[ref_ID], &BaseDist_All , 1);
						PrintRegionDist (file_stat_0_o  , file_ratio_0_o   , &BedTable, PosCoverage, &BamHeader, ref_ID, 0  , &RegionDist_0_All  , 1);
						PrintRegionDist (file_stat_10_o , file_ratio_10_o  ,&BedTable, PosCoverage, &BamHeader, ref_ID, 9  , &RegionDist_10_All , 1);
						PrintRegionDist (file_stat_20_o , file_ratio_20_o  ,&BedTable, PosCoverage, &BamHeader, ref_ID, 19 , &RegionDist_20_All , 1);
						PrintRegionDist (file_stat_30_o , file_ratio_30_o  ,&BedTable, PosCoverage, &BamHeader, ref_ID, 29 , &RegionDist_30_All , 1);
						PrintRegionDist (file_stat_40_o , file_ratio_40_o  ,&BedTable, PosCoverage, &BamHeader, ref_ID, 39 , &RegionDist_40_All , 1);
						PrintRegionDist (file_stat_50_o , file_ratio_50_o  ,&BedTable, PosCoverage, &BamHeader, ref_ID, 49 , &RegionDist_50_All , 1);
						PrintRegionDist (file_stat_60_o , file_ratio_60_o  ,&BedTable, PosCoverage, &BamHeader, ref_ID, 59 , &RegionDist_60_All , 1);
						PrintRegionDist (file_stat_70_o , file_ratio_70_o  ,&BedTable, PosCoverage, &BamHeader, ref_ID, 69 , &RegionDist_70_All , 1);
						PrintRegionDist (file_stat_80_o , file_ratio_80_o  ,&BedTable, PosCoverage, &BamHeader, ref_ID, 79 , &RegionDist_80_All , 1);
						PrintRegionDist (file_stat_90_o , file_ratio_90_o  ,&BedTable, PosCoverage, &BamHeader, ref_ID, 89 , &RegionDist_90_All , 1);
						PrintRegionDist (file_stat_100_o, file_ratio_100_o  ,&BedTable, PosCoverage, &BamHeader, ref_ID, 99, &RegionDist_100_All, 1);
						PrintRegionDist (file_stat_200_o, file_ratio_200_o  ,&BedTable, PosCoverage, &BamHeader, ref_ID, 199, &RegionDist_200_All, 1);
						PrintRegionDist (file_stat_500_o, file_ratio_500_o  ,&BedTable, PosCoverage, &BamHeader, ref_ID, 499, &RegionDist_500_All, 1);
						chr_in_target = 0;
						free(on_target);
						free(PosCoverage);
						//printf("%s\t%lu\t%lu\t%lu\n",BamHeader.chr_name[ref_ID],map_in,map_out, valid_region);
					}
					fprintf(file_region_o,"%s\t%lu\t%lu\t%lu\n",BamHeader.chr_name[ref_ID],map_in,map_out, valid_region);

					//printf("%s\t%lu\t%lu\t%lu\n",BamHeader.chr_name[ref_ID],map_in,map_out, valid_region);
					map_all_in += map_in;
					map_all_out += map_out;
					valid_all_region += valid_region;
					valid_region = 0;
					map_in = 0;
					map_out = 0;
				}	
				if ( AlignmentHeader.refID != -1){	
					if (BedTable.table_start[AlignmentHeader.refID] != BedTable.table_end[AlignmentHeader.refID]){
						chr_in_target = 1;
						PosCoverage = calloc(BamHeader.chr_length[AlignmentHeader.refID],sizeof(posCoverage));
						on_target = calloc(BamHeader.chr_length[AlignmentHeader.refID],sizeof(uint8_t));
						
						for (i = BedTable.table_start[AlignmentHeader.refID]; i < BedTable.table_end[AlignmentHeader.refID];i++){
							for (m = BedTable.start[i]; m < BedTable.end[i]; m++){
								if (on_target[m] == 0){
									on_target[m] = 1;
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
			if ((AlignmentHeader.FLAG&4) == 0 && ((AlignmentHeader.FLAG&1024) == 0 || ToolsFlags->flag_dup==1)){
//			if ((AlignmentHeader.FLAG&4) == 0 && (AlignmentHeader.FLAG&256) == 0 && ((AlignmentHeader.FLAG&1024) == 0 || ToolsFlags->flag_dup==1)){
				
				if (chr_in_target == 1){
					alignment_DepthTxt(address, &AlignmentHeader, PosCoverage, BamHeader.chr_length[ref_ID]);
					alignment_Range(address, &AlignmentHeader, on_target, &map_in, &map_out, BamHeader.chr_length[ref_ID]);
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
		if (ToolsFlags->flag_hide == 1){	printf("[Bam File Unzip %d / %d ] %s done\n",ref_ID+1,BamHeader.n_ref,BamHeader.chr_name[ref_ID]);	}

		if (chr_in_target == 1){
			PrintBaseDist (file_base_stat_o, PosCoverage, on_target, BamHeader.chr_name[ref_ID], BamHeader.chr_length[ref_ID], &BaseDist_All , 1);
			PrintRegionDist (file_stat_0_o  , file_ratio_0_o   ,&BedTable, PosCoverage, &BamHeader, ref_ID, 0  , &RegionDist_0_All  , 1);
			PrintRegionDist (file_stat_10_o , file_ratio_10_o  ,&BedTable, PosCoverage, &BamHeader, ref_ID, 9  , &RegionDist_10_All , 1);
			PrintRegionDist (file_stat_20_o , file_ratio_20_o  ,&BedTable, PosCoverage, &BamHeader, ref_ID, 19 , &RegionDist_20_All , 1);
			PrintRegionDist (file_stat_30_o , file_ratio_30_o  ,&BedTable, PosCoverage, &BamHeader, ref_ID, 29 , &RegionDist_30_All , 1);
			PrintRegionDist (file_stat_40_o , file_ratio_40_o  ,&BedTable, PosCoverage, &BamHeader, ref_ID, 39 , &RegionDist_40_All , 1);
			PrintRegionDist (file_stat_50_o , file_ratio_50_o  ,&BedTable, PosCoverage, &BamHeader, ref_ID, 49 , &RegionDist_50_All , 1);
			PrintRegionDist (file_stat_60_o , file_ratio_60_o  ,&BedTable, PosCoverage, &BamHeader, ref_ID, 59 , &RegionDist_60_All , 1);
			PrintRegionDist (file_stat_70_o , file_ratio_70_o  ,&BedTable, PosCoverage, &BamHeader, ref_ID, 69 , &RegionDist_70_All , 1);
			PrintRegionDist (file_stat_80_o , file_ratio_80_o  ,&BedTable, PosCoverage, &BamHeader, ref_ID, 79 , &RegionDist_80_All , 1);
			PrintRegionDist (file_stat_90_o , file_ratio_90_o  ,&BedTable, PosCoverage, &BamHeader, ref_ID, 89 , &RegionDist_90_All , 1);
			PrintRegionDist (file_stat_100_o, file_ratio_100_o  ,&BedTable, PosCoverage, &BamHeader, ref_ID, 99, &RegionDist_100_All, 1);
			PrintRegionDist (file_stat_200_o, file_ratio_200_o  ,&BedTable, PosCoverage, &BamHeader, ref_ID, 199, &RegionDist_200_All, 1);
			PrintRegionDist (file_stat_500_o, file_ratio_500_o  ,&BedTable, PosCoverage, &BamHeader, ref_ID, 499, &RegionDist_500_All, 1);
			free(PosCoverage);
			chr_in_target = 0;
		}
		fprintf(file_region_o,"%s\t%lu\t%lu\t%lu\n",BamHeader.chr_name[ref_ID],map_in,map_out,valid_region);
		map_all_in += map_in;
		map_all_out += map_out;
		valid_all_region += valid_region;
	}
	
	PrintBaseDist (file_base_stat_o, PosCoverage, on_target, "", 0, &BaseDist_All , 2);
	PrintRegionDist (file_stat_0_o  , file_ratio_0_o   ,&BedTable, PosCoverage, &BamHeader, ref_ID, 0  , &RegionDist_0_All  , 2);
		PrintRegionDist (file_stat_10_o , file_ratio_10_o  ,&BedTable, PosCoverage, &BamHeader, ref_ID, 9  , &RegionDist_10_All , 2);
		PrintRegionDist (file_stat_20_o , file_ratio_20_o  ,&BedTable, PosCoverage, &BamHeader, ref_ID, 19 , &RegionDist_20_All , 2);
		PrintRegionDist (file_stat_30_o , file_ratio_30_o  ,&BedTable, PosCoverage, &BamHeader, ref_ID, 29 , &RegionDist_30_All , 2);
		PrintRegionDist (file_stat_40_o , file_ratio_40_o  ,&BedTable, PosCoverage, &BamHeader, ref_ID, 39 , &RegionDist_40_All , 2);
		PrintRegionDist (file_stat_50_o , file_ratio_50_o  ,&BedTable, PosCoverage, &BamHeader, ref_ID, 49 , &RegionDist_50_All , 2);
		PrintRegionDist (file_stat_60_o , file_ratio_60_o  ,&BedTable, PosCoverage, &BamHeader, ref_ID, 59 , &RegionDist_60_All , 2);
		PrintRegionDist (file_stat_70_o , file_ratio_70_o  ,&BedTable, PosCoverage, &BamHeader, ref_ID, 69 , &RegionDist_70_All , 2);
		PrintRegionDist (file_stat_80_o , file_ratio_80_o  ,&BedTable, PosCoverage, &BamHeader, ref_ID, 79 , &RegionDist_80_All , 2);
		PrintRegionDist (file_stat_90_o , file_ratio_90_o  ,&BedTable, PosCoverage, &BamHeader, ref_ID, 89 , &RegionDist_90_All , 2);
		PrintRegionDist (file_stat_100_o, file_ratio_100_o  ,&BedTable, PosCoverage, &BamHeader, ref_ID, 99, &RegionDist_100_All, 2);
		PrintRegionDist (file_stat_200_o, file_ratio_200_o  ,&BedTable, PosCoverage, &BamHeader, ref_ID, 199, &RegionDist_200_All, 2);
		PrintRegionDist (file_stat_500_o, file_ratio_500_o  ,&BedTable, PosCoverage, &BamHeader, ref_ID, 499, &RegionDist_500_All, 2);
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
