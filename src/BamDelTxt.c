#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include "zlib.h"

#include "BamCommonLibrary.h"
void	alignment_DeletionTxt(uint8_t *stream, alignmentHeader *AlignmentHeader, FILE *file_deletion_o, uint32_t ref_length);

int	BamDelTxt(FILE *file_bam_i, FILE *file_bai_i, FILE *file_bed_i, toolsFlags *ToolsFlags){

	FILE	*file_cov_o;

	uint8_t *stream_i;//[65536];
	uint8_t *stream_o;//[65536];
	uint8_t	*buffer;
	uint8_t	*top;

	int	i,j;
	int	counter;
	int	len_data;
	int	ref_ID	= -1;

	uint32_t CRC32;
	uint32_t I_size;

	uint16_t unknown;

	z_stream infstream;
	bgfzTail BgfzTail;
	
	bamHeader	BamHeader;
	bedTable	BedTable;

	uint64_t	offset_beg;
	uint64_t	offset_bgzf;
	uint64_t	offset_decomp;

	alignmentHeader	AlignmentHeader;

	char	magic[4];
	int32_t	l_text;
	char	*text;
	int32_t	n_ref;
	uint8_t *address;
	char	filename[name_len];
	char	dirname_region[name_len];

	int	flag_bam=0;
	int	flag_outlen=0;
	int	flag_header=0;
	int	flag_cov = 0;
	int	flag_sam = 0;
	int	flag_region = 0;
	int	flag_region_file = 0;
	int	flag_dirname_region= 0;
	int	flag_break	= 0;

	stream_i= calloc(65536,sizeof(uint8_t));
	stream_o= calloc(65536,sizeof(uint8_t));
	buffer	= calloc(65536*2,sizeof(uint8_t));	

	/* Initial*/
	stream_i[0] = 120;
	stream_i[1] = 156;
	/*---------------*/


	top = CatchBamHeader (file_bam_i, &BamHeader, stream_i, stream_o, buffer, top);
	counter	= top - buffer;
	address = buffer;

	
	/* comment 2020/03/30 YC
	if (ToolsFlags->flag_hide != 1){
		for (i = 0;i < BamHeader.n_ref;i++){
			printf("%s\t%d\n",BamHeader.chr_name[i], BamHeader.chr_length[i]);
		}
	}
	*/
	//Start Unzip and Process Bam File
	createRegionTable(file_bed_i, &BedTable, &BamHeader, ToolsFlags);
	
	for (j = 0; j < BamHeader.n_ref; j++){
		ref_ID	= j;
		if (BedTable.table_end[ref_ID] != BedTable.table_start[ref_ID]){
			offset_beg = CatchOffset (file_bai_i, ref_ID, BedTable.start[BedTable.table_start[ref_ID]]);
			offset_bgzf = offset_beg >> 16;
			offset_decomp = offset_beg&65535;
//			printf("%d\t%d\n", BedTable.start[BedTable.table_start[ref_ID]], BedTable.end[BedTable.table_start[ref_ID]]);
			//printf("%lu\t%lu\t%lu\n",offset_beg, offset_bgzf, offset_decomp);			
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

			strcat_triple(filename, "", BamHeader.chr_name[ref_ID], "_del.txt", name_len);
			file_cov_o = fopen(filename,"w");
			//Start Unzip and Process Bam File	
			while ( (len_data = BGFZBlock(file_bam_i)) > 0){
				decompressBlock(&infstream, stream_i, stream_o, len_data, file_bam_i);
				fread(&BgfzTail	,sizeof(bgfzTail),1,file_bam_i);

				memcpy(buffer,address,sizeof(uint8_t)*counter);
				
			//	printf("%u\t%u\t%u\t%u\n", buffer, &buffer[counter], stream_o, BgfzTail.I_size);

				memcpy(&buffer[counter], stream_o, BgfzTail.I_size);

			//	printf("Johnnash\n");	
				top	= buffer + counter + BgfzTail.I_size;
				address = buffer;
				if (top - address > sizeof(alignmentHeader)){
					memcpy(&AlignmentHeader, address, sizeof(alignmentHeader));
				}
	
				while ( top - address >= (AlignmentHeader.block_size + 4) ){
			//		printf("%d\t%d\t%d\n",j, AlignmentHeader.refID, AlignmentHeader.pos);
			//		printf("%d\t%d\t%d\n",flag_break, top - address, (AlignmentHeader.block_size + 4));	
					if (AlignmentHeader.refID != ref_ID){
						flag_break = 1;
						break;
					}
					address += sizeof(alignmentHeader);	
					//Core_START
					alignment_DeletionTxt(address, &AlignmentHeader, file_cov_o, BamHeader.chr_length[ref_ID]);
					//Core_END
					address += (AlignmentHeader.block_size + 4) - sizeof(alignmentHeader);

					if (top - address > sizeof(alignmentHeader)){
						memcpy(&AlignmentHeader, address, sizeof(alignmentHeader));
					}else {
						break;	
					}
				}

			//	printf("%d\t%d\t%d\tJohnnash - 2\n",flag_break, top - address, (AlignmentHeader.block_size + 4));	
				if (flag_break	== 1){	flag_break = 0;	break;
				}else {			counter = top - address;}
			}
		//	fclose(file_cov_o);
		}
		if (ToolsFlags->flag_hide == 1){
			printf("[Bam File Unzip %d / %d ] %s done\n",ref_ID+1,BamHeader.n_ref,BamHeader.chr_name[ref_ID]);
		}
	}
	
	free(stream_i);	
	free(stream_o);	
	free(buffer);	
	return 0;
}





void	alignment_DeletionTxt(uint8_t *stream, alignmentHeader *AlignmentHeader, FILE *file_deletion_o, uint32_t ref_length){

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
					fprintf(file_deletion_o,"%d_%d",index+1,op_len);
					flag = 1;
				}else {
					fprintf(file_deletion_o,",%d_%d",index+1,op_len);
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
		fprintf(file_deletion_o,"\n");
	}
//	free(read_name);
	free(cigar);
//	free(seq);
}

