#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#include "zlib.h"

#include "BamCommonLibrary.h"

int	BamDepthDist(FILE *file_bam_i, FILE *file_bai_i, FILE *file_bed_i, char *chromosome, uint32_t start, uint32_t end, toolsFlags *ToolsFlags){

	FILE	*file_cov_o;

	uint8_t *stream_i;//[65536];
	uint8_t *stream_o;//[65536];
	uint8_t	*buffer;
	uint8_t	*top;

	int	i,j,m,n,k;
	int	len_data;
	int	counter;
	int	ref_ID	= -1;

	uint32_t CRC32;
	uint32_t I_size;

	uint16_t unknown;

	z_stream infstream;
	bgfzTail BgfzTail;
	
	alignmentHeader	AlignmentHeader;

	uint8_t *address;
	char	filename[name_len];
	char	dirname_region[name_len];

	int	flag_header=0;
	int	flag_cov = 0;
	int	flag_sam = 0;
	int	flag_region = 0;
	int	flag_region_file = 0;
	int	flag_dirname_region= 0;
	int	All_Coverage;
	int	chr_ID = -2;
	uint64_t	offset;

	int32_t	n_ref_bai;
	int32_t	n_bin;
	uint32_t	bin;
	
	uint64_t	offset_beg;
	uint64_t	offset_bgzf;
	uint64_t	offset_decomp;
	int64_t	diff;
	int	region2bin;
	posCoverage 	*PosCoverage;
	bamHeader	BamHeader;
	bedTable	BedTable;

	bedTable	TargetTable;
	
	char	line[LINE_MAX_LEN];
	int	num_target	= 0;
	char	chromosome_name[LINE_MAX_LEN];
	int32_t	end_position;


	uint32_t	max_position;
	uint64_t	Total_Depth = 0;
	uint64_t	temp_Depth = 0;
	uint64_t	Total_Length = 0;
	uint64_t	Cover_Length = 0;
	double	Ave_Depth;
	int32_t	max_cover_index;

	struct timeval start_clock, end_clock;
	long long  sum_clock = 0;
	int	num_reads;

	int	*threshold;
	uint64_t	*num_coverage_chr;
	uint64_t	*num_coverage_all;
	uint64_t	coverage;
	
	int	num_threshold = 0;
	int	index_threshold = 0;
	char	*token;
	char	*duplicate_string;



	if (ToolsFlags->flag_coverage == 1){
		duplicate_string = strdup(ToolsFlags->coverage_threshold);
		//token = strtok( ToolsFlags->coverage_threshold, ",");
		token = strtok( duplicate_string, ",");
		while( token != NULL ){
			num_threshold++;
			/* Get next token: */
			token = strtok( NULL, ",");
		}

		threshold = calloc(num_threshold, sizeof(int));
		num_coverage_all = calloc(num_threshold, sizeof(uint64_t));
		num_coverage_chr = calloc(num_threshold, sizeof(uint64_t));

		num_threshold	= 0;
		//printf("%s\n", ToolsFlags->coverage_threshold);
		token = strtok( ToolsFlags->coverage_threshold, ",");
		while( token != NULL ){
			/* While there are tokens in "string" */
			threshold[num_threshold] = atoi(token);
//			printf("%d\n", atoi(token));
			token = strtok( NULL, ",");
			if (threshold[num_threshold] == 0){
				printf("Error: value of threshold\n");
				return -1;
			}
			num_threshold++;
			/* Get next token: */
		}
//		printf("%d\n",num_threshold);
		printf("CHR");
		for ( index_threshold = 0; index_threshold < num_threshold; index_threshold++){
			printf("\t%dX", threshold[index_threshold]);
		}
		printf("\n");
		
	}

	

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

	if (ToolsFlags->flag_hide != 1){
		for (i = 0;i < BamHeader.n_ref;i++){
			printf("%s\t%d\n",BamHeader.chr_name[i], BamHeader.chr_length[i]);
		}
	}

	if (ToolsFlags->flag_target){
		decodeBedFile(ToolsFlags->file_target, &TargetTable, &BamHeader);
	}

	createRegionTable(file_bed_i, &BedTable, &BamHeader, ToolsFlags);

		
	//BAI File
	for (j = 0; j < BamHeader.n_ref; j++){
		ref_ID = j;
		if (BedTable.table_end[ref_ID] != BedTable.table_start[ref_ID]){
			offset_beg = CatchOffset (file_bai_i, ref_ID, BedTable.start[BedTable.table_start[ref_ID]]);
			offset_bgzf = offset_beg >> 16;
			offset_decomp = offset_beg&65535;

			PosCoverage = calloc(BamHeader.chr_length[ref_ID],sizeof(posCoverage));
	//		printf("%lu\t%lu\t%lu\n",offset_beg, offset_bgzf, offset_decomp);			

			if (offset_beg == 0){

			}else {
			//Offset Bam File
				fseek(file_bam_i,offset_bgzf,SEEK_SET);
				len_data = BGFZBlock(file_bam_i);
				decompressBlock(&infstream, stream_i, stream_o, len_data, file_bam_i);
				fread(&BgfzTail	,sizeof(bgfzTail),1,file_bam_i);
				memcpy(buffer,stream_o,BgfzTail.I_size);
				address = buffer + offset_decomp;
				counter = BgfzTail.I_size - offset_decomp;
			}
				//printf("%u\t%u\t%u\n", buffer, counter, BgfzTail.I_size);
	
			//Start Unzip and Process Bam File	
			while ( (len_data = BGFZBlock(file_bam_i)) > 0 ){
//gettimeofday(&start_clock, NULL);
				decompressBlock(&infstream, stream_i, stream_o, len_data, file_bam_i);
//gettimeofday(&end_clock, NULL);
//sum_clock	+= ((end_clock.tv_sec - start_clock.tv_sec)*1000000+(end_clock.tv_usec - start_clock.tv_usec));
//printf("%d\t%lu\t%f\t%d\t%lu\t%d\t%d\n", k, sum_clock, (float)sum_clock / (k+1), num_reads, ((end_clock.tv_sec - start_clock.tv_sec)*1000000+(end_clock.tv_usec - start_clock.tv_usec)), ToolsFlags->num_ref, ToolsFlags->num_alt);

				fread(&BgfzTail	,sizeof(bgfzTail),1,file_bam_i);
				if (address - buffer > counter){
					memcpy(buffer, address, sizeof(uint8_t)*counter);
				}else {
					memmove(buffer, address, sizeof(uint8_t)*counter);
				}
				memcpy(&buffer[counter], stream_o, BgfzTail.I_size);
				top	= buffer + counter + BgfzTail.I_size;
				address = buffer;
				//printf("%u\t%u\t%u\n", buffer, counter, BgfzTail.I_size);
				if (top - address > sizeof(alignmentHeader)){
					memcpy(&AlignmentHeader, address, sizeof(alignmentHeader));
				}
				//printf("%u\t%u\t%u\n", top, address, AlignmentHeader.refID);
				while ( top - address >= (AlignmentHeader.block_size + 4) ){
					if (AlignmentHeader.refID > ref_ID || AlignmentHeader.refID < (ref_ID - 1)){
//					if (AlignmentHeader.refID > ref_ID){
						break;
					}
					if (AlignmentHeader.refID == ref_ID && AlignmentHeader.pos > BedTable.table_max[ref_ID]){
						break;	
					} 

					address += sizeof(alignmentHeader);
				//	printf("%u\n",AlignmentHeader.pos);
					//Core_START
					if (AlignmentHeader.refID == ref_ID){
						if ((AlignmentHeader.FLAG&1024)==0 || ToolsFlags->flag_dup==1 ){
							if (!ToolsFlags->flag_filter || SelectReads(address, &AlignmentHeader, 30, 30)){
								if (ToolsFlags->flag_target){
									end_position	= alignment_EndPosition(address, &AlignmentHeader, BamHeader.chr_length[ref_ID]);
									max_cover_index	= concernTargetRegion( &TargetTable, ref_ID, AlignmentHeader.pos, end_position);
									//printf("%u\t%u\t%d\n", AlignmentHeader.pos, end_position, max_cover_index);
									if (max_cover_index >= 0){
										alignment_DepthDist_TR(address, &AlignmentHeader, PosCoverage, TargetTable.start[max_cover_index], TargetTable.end[max_cover_index]);	
									}
								}else {
									alignment_DepthDist_TR(address, &AlignmentHeader, PosCoverage, 0, BamHeader.chr_length[ref_ID]);
								}
								//alignment_DepthDist(address, &AlignmentHeader, PosCoverage, BamHeader.chr_length[ref_ID]);
							}
						}
					}
					//Core_END
			
					address += (AlignmentHeader.block_size + 4) - sizeof(alignmentHeader);

					if (top - address > sizeof(alignmentHeader)){
						memcpy(&AlignmentHeader, address, sizeof(alignmentHeader));
					}else {
						break;	
					}
				}
			//	printf("%u\t%u\t%u\n", BedTable.start[BedTable.table_start[ref_ID]], AlignmentHeader.pos, BedTable.table_max[ref_ID]);
				if (AlignmentHeader.refID > ref_ID || AlignmentHeader.refID < (ref_ID - 1)){
		//		if (AlignmentHeader.refID > ref_ID){
					break;
				}
				if (AlignmentHeader.refID == ref_ID && AlignmentHeader.pos > BedTable.table_max[ref_ID]){
					break;
				}
				counter = top - address;
			}
			if (ToolsFlags->flag_hide == 0){
				printf("[Bam File Unzip %d / %d ] %s done\n",ref_ID+1, BamHeader.n_ref, BamHeader.chr_name[ref_ID]);
			}
				
			max_position = 0;

			if (ToolsFlags->flag_coverage == 1){
				for (m = BedTable.table_start[j]; m < BedTable.table_end[j];m++){
					//printf("%d\t%d\n",BedTable.table_start[j],BedTable.table_end[j]);
					if ( BedTable.end[m] > max_position){
						if (BedTable.start[m] < max_position){
							BedTable.start[m] = max_position;
						}
					
						for (n = BedTable.start[m]; n < BedTable.end[m];n++){
								coverage = (PosCoverage[n].A +
								PosCoverage[n].C +
								PosCoverage[n].G +
								PosCoverage[n].T +
								PosCoverage[n].N +
								PosCoverage[n].DEL);
								Total_Length++;
								for ( index_threshold = 0; index_threshold < num_threshold; index_threshold++){
									if (coverage >= threshold[index_threshold]){	
										num_coverage_chr[index_threshold]++;	
										num_coverage_all[index_threshold]++;	
									}
								}
						}
						max_position = BedTable.end[m];	
					}
				}
				printf("%s",BamHeader.chr_name[ref_ID]);
				for ( index_threshold = 0; index_threshold < num_threshold; index_threshold++){
					printf("\t%u",num_coverage_chr[index_threshold]);
					num_coverage_chr[index_threshold] = 0;
				}
				printf("\n");
				
		
			}else if (ToolsFlags->flag_simple == 1){
				for (m = BedTable.table_start[j]; m < BedTable.table_end[j];m++){
					//printf("%d\t%d\n",BedTable.table_start[j],BedTable.table_end[j]);
					if ( BedTable.end[m] > max_position){
						if (BedTable.start[m] < max_position){
							BedTable.start[m] = max_position;
						}
					
						for (n = BedTable.start[m]; n < BedTable.end[m];n++){
								temp_Depth = (PosCoverage[n].A +
								PosCoverage[n].C +
								PosCoverage[n].G +
								PosCoverage[n].T +
								PosCoverage[n].N +
								PosCoverage[n].DEL);
								if (temp_Depth > 0){
									Cover_Length++;
									Total_Depth += temp_Depth;
								}
								Total_Length++;
						}
						max_position = BedTable.end[m];	
					}
				}
			}else {
				for (m = BedTable.table_start[j]; m < BedTable.table_end[j];m++){
					//printf("%d\t%d\n",BedTable.table_start[j],BedTable.table_end[j]);
					if ( BedTable.end[m] > max_position){
						if (BedTable.start[m] < max_position){
							BedTable.start[m] = max_position;
						}
				
						for (n = BedTable.start[m]; n < BedTable.end[m];n++){
							printf("%s\t",BamHeader.chr_name[j]);
							printf("%u\t",n + 1);
							printf("%u\t",PosCoverage[n].A);
							printf("%u\t",PosCoverage[n].C);
							printf("%u\t",PosCoverage[n].G);
							printf("%u\t",PosCoverage[n].T);
							printf("%u\t",PosCoverage[n].N);
							printf("%u\t",PosCoverage[n].DEL);
							printf("%u",	PosCoverage[n].A +
									PosCoverage[n].C +
									PosCoverage[n].G +
									PosCoverage[n].T +
									PosCoverage[n].N +
									PosCoverage[n].DEL);					
							printf("\n");
						}
						max_position = BedTable.end[m];	
					}
				}
			}
			free(PosCoverage);
		}
	}
	
	if (ToolsFlags->flag_coverage == 1){
		printf("Total");
		for ( index_threshold = 0; index_threshold < num_threshold; index_threshold++){
			printf("\t%u",num_coverage_all[index_threshold]);
		}
		printf("\n");
		printf("Rate");
		for ( index_threshold = 0; index_threshold < num_threshold; index_threshold++){
			printf("\t%.5lf", (double)num_coverage_all[index_threshold] /Total_Length);
		}
		printf("\n");

	}else if (ToolsFlags->flag_simple == 1){
		printf("%f\t%f\n", (double)Total_Depth / Cover_Length, (double)Cover_Length / Total_Length);
	}
	free(stream_i);	
	free(stream_o);	
	free(buffer);	
	return 0;
}

