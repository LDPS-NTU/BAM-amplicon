#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include "zlib.h"

#include "BamCommonLibrary.h"

int	BamMappingLength(FILE *file_bam_i, FILE *file_bai_i, FILE *file_bed_i, char *chr, uint32_t start, uint32_t end, toolsFlags *ToolsFlags){

	FILE	*file_count_o;
	FILE	*file_region_o;
	FILE	*file_base_stat_o;

	uint8_t *stream_i;//[65536];
	uint8_t *stream_o;//[65536];
	uint8_t	*buffer;
	uint8_t *address;
	uint8_t	*top;
	
	int	i,j,m;
	int	len_data;
	int	counter;
	int	ref_ID	= -1;
	int	length;
	uint32_t CRC32;
	uint32_t I_size;
	uint32_t end_position;
	uint16_t unknown;

	z_stream infstream;
	bgfzTail BgfzTail;
	
	alignmentHeader	AlignmentHeader;



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

	uint32_t	num_target = 0;	

	char	chromosome_name[name_len];
	uint64_t	offset_beg;
	uint64_t	offset_bgzf;
	uint64_t	offset_decomp;
	int	index;

	uint64_t	*length_array;
	uint32_t	max_mapping_length;
	max_mapping_length	= 300;
	length_array	= calloc(sizeof(uint64_t), max_mapping_length);

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

	if (ToolsFlags->flag_hide != 1){
		for (i = 0;i < BamHeader.n_ref;i++){
			printf("%s\t%d\n",BamHeader.chr_name[i], BamHeader.chr_length[i]);
		}
	}
	
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
	//	printf("Total\n");
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
					if (ToolsFlags->flag_hide == 0 ){
						printf("[Bam File Unzip %d / %d ] %s done\n",ref_ID+1, BamHeader.n_ref, BamHeader.chr_name[ref_ID]);
					}
					
					if (chr_in_target == 1){
						chr_in_target = 0;
						free(on_target);
						//printf("%s\t%lu\t%lu\t%lu\n",BamHeader.chr_name[ref_ID],map_in,map_out, valid_region);
					}

					//printf("%s\t%lu\t%lu\t%lu\n",BamHeader.chr_name[ref_ID],map_in,map_out, valid_region);
				}	
				if ( AlignmentHeader.refID != -1){	
					if (BedTable.table_start[AlignmentHeader.refID] != BedTable.table_end[AlignmentHeader.refID]){
						chr_in_target = 1;
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
				
				if (chr_in_target == 1){
					end_position	= alignment_EndPosition(address, &AlignmentHeader, BamHeader.chr_length[ref_ID]);
					length	= end_position - AlignmentHeader.pos;
					if (length >= max_mapping_length){
						
						length_array = realloc(length_array, sizeof(uint64_t)*(length+10));
						if (length_array !=NULL){
							for (i = max_mapping_length; i < max_mapping_length+10;i++){
								length_array[i] = 0;
							}
							max_mapping_length = max_mapping_length+10;
							length_array[length]++;
						}else {
							printf("Error\n");
							exit(1);
						}

					}else {
						length_array[length]++;
					}
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
		if (ToolsFlags->flag_hide == 0){	
			printf("[Bam File Unzip %d / %d ] %s done\n",ref_ID+1,BamHeader.n_ref,BamHeader.chr_name[ref_ID]);	
		}

	}
	for (i = 0;i < max_mapping_length;i++){
		printf("%d\t%ld\n",i, length_array[i]);
	}	
	free(stream_i);	
	free(stream_o);	
	free(buffer);	
	return 0;
}
