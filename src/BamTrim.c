#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include "zlib.h"

#include "BamCommonLibrary.h"

int	measureTargetRegion(bedTable *TargetTable, int ref_ID, char *ref_name, uint32_t first_pos, uint32_t last_pos);
void	readTrim(uint8_t *stream, alignmentHeader *AlignmentHeader, uint32_t ref_length);

int	BamTrim(FILE *file_bam_i, FILE *file_bai_i, FILE *file_bed_i, toolsFlags *ToolsFlags){

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

	char	line[LINE_MAX_LEN];

	uint8_t	*on_target;
	uint64_t	valid_region = 0;
	uint64_t	valid_all_region = 0;
	int	chr_in_target = 0;
	int32_t	end_position;
	int32_t	max_cover_index;

	bamHeader	BamHeader;
	bedTable	BedTable;
	bedTable	TargetTable;

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

	/* comment 2020/03/30 YC
	if (ToolsFlags->flag_hide == 1){
		for (i = 0;i < BamHeader.n_ref;i++){
			printf("%s\t%d\n",BamHeader.chr_name[i], BamHeader.chr_length[i]);
		}
	}
	*/

	// For Amplicon
	if (ToolsFlags->flag_target){
		decodeBedFile(ToolsFlags->file_target, &TargetTable, &BamHeader);
	}

	//
	createRegionTable(file_bed_i, &BedTable, &BamHeader, ToolsFlags);
	//Bed File
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
					if (ToolsFlags->flag_hide == 1 ){
						printf("[Bam File Unzip %d / %d ] %s done\n",ref_ID+1, BamHeader.n_ref, BamHeader.chr_name[ref_ID]);
					}
				}	
				ref_ID = AlignmentHeader.refID;	
			}
			if (ref_ID == -1){	break;	}

			address += sizeof(alignmentHeader);	
			if ((AlignmentHeader.FLAG&4) == 0 && ((AlignmentHeader.FLAG&1024) == 0 || ToolsFlags->flag_dup==1)){
				end_position	= alignment_EndPosition(address, &AlignmentHeader, BamHeader.chr_length[ref_ID]);
				max_cover_index	= measureTargetRegion( &TargetTable, ref_ID,BamHeader.chr_name[ref_ID], AlignmentHeader.pos, end_position);
				//printf("%u\t%u\t%d\n", AlignmentHeader.pos, end_position, max_cover_index);
				if (max_cover_index >= 0){
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

		if (ref_ID == -1){	break;	}
	}
	if (ref_ID == -1){
		if (ToolsFlags->flag_hide == 1){
			printf("[Bam File Unzip %d / %d ] %s done\n",ref_ID+1,BamHeader.n_ref,BamHeader.chr_name[ref_ID]);	
		}
	}
	free(stream_i);	
	free(stream_o);	
	free(buffer);	
	return 0;
}

void	readTrim(uint8_t *stream, alignmentHeader *AlignmentHeader, uint32_t ref_length){

	int	i,j;

//	char	*read_name;
	uint32_t *cigar;
//	uint8_t	*seq;
//	char	*qual;

	int	op;
	int	op_len;
	char	residue;
	int	length;
	int	index;
	int	flag = 0;

//	MemoryCopy(&stream,(void **)&read_name	, AlignmentHeader->l_read_name	, sizeof(char));
	stream	+= AlignmentHeader->l_read_name;

	MemoryCopy(&stream,(void **)&cigar	, AlignmentHeader->n_cigar_op	, sizeof(uint32_t));
//	MemoryCopy(&stream,(void **)&seq	, (AlignmentHeader->l_seq+1)/2	, sizeof(uint8_t));

	length = 0;
	index = AlignmentHeader->pos;
	for (i = 0;i < AlignmentHeader->n_cigar_op;i++){
		op	= cigar[i] & 7;
		op_len	= cigar[i] >> 4;
		if (CIGAR(op) == 'M'){
			if (index + op_len > ref_length){
				op_len = ref_length - index;
			}
		length+=op_len;
		index+=op_len;
		
		}else if (CIGAR(op) == 'I'){
			length += op_len;
		}else if (CIGAR(op) == 'D'){
			if (op_len > 3 ){
				if (flag==0){
					printf("%d_%d",index+1,op_len);
					flag = 1;
				}else {
					printf(",%d_%d",index+1,op_len);
				}
			}
			index+=op_len;
		}else if (CIGAR(op) == 'S'){
			length += op_len;
		}else if (CIGAR(op) == 'N'){
			index += op_len;
		}
		// CIGAR(op) == 'H' don't care about it.
		
	}
	if(flag==1){
		printf("\n");
	}
//	free(read_name);
	free(cigar);
//	free(seq);
		
}

int	measureTargetRegion(bedTable *TargetTable, int ref_ID, char *ref_name, uint32_t first_pos, uint32_t last_pos){
	int	max_cover_index = -1;
	float	max_cover = 0;
	int	i;
	int	min_distance	= 999999999;
	int	distance;
	int	start;
	int	end;
	int	min_start	= 999999999;
	int	min_end		= 999999999;
	int	region_index	= -1;

	for (i =TargetTable->table_start[ref_ID]; i < TargetTable->table_end[ref_ID]; i++){
		if (first_pos < TargetTable->end[i] && last_pos > TargetTable->start[i]){
			if (first_pos >=  TargetTable->start[i]){	distance = first_pos - TargetTable->start[i];}
			else {						distance = TargetTable->start[i] - first_pos;}

			if (last_pos >= TargetTable->end[i]){	distance += (last_pos - TargetTable->end[i]);}
			else {					distance += (TargetTable->end[i] - last_pos);}
			start	= TargetTable->start[i] - first_pos;
			end	= last_pos - TargetTable->end[i];

			if (distance < min_distance){
				region_index	= i;
				min_distance = distance;
				min_start	= start;
				min_end	= end;
				max_cover_index = i;
			}
			/*
			*/
		}
	}
	printf("%u\t%d\t%d\t%s:%u-%u\t%s:%u-%u\n",min_distance, min_start, min_end, ref_name, first_pos, last_pos, ref_name, TargetTable->start[region_index], TargetTable->end[region_index]);
	return max_cover_index;
}
 
