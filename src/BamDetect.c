#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include "zlib.h"

#include "BamCommonLibrary.h"




void	alignment_Detect(uint8_t *stream, alignmentHeader *AlignmentHeader, uint32_t ref_length, toolsFlags *ToolsFlags);

int	BamDetect(FILE *file_bam_i, FILE *file_bai_i, char *chromosome, toolsFlags *ToolsFlags){

	FILE	*file_cov_o;

	uint8_t *stream_i;//[65536];
	uint8_t *stream_o;//[65536];
	uint8_t	*buffer;
	uint8_t	*top;

	int	i,j,m,n;
	int	len_data;
	int	counter;
	int	ref_ID	= -1;

	uint32_t CRC32;
	uint32_t I_size;

	uint16_t unknown;

	z_stream infstream;
	bgfzTail BgfzTail;
	
	alignmentHeader	AlignmentHeader;

	char	magic[4];
	char	**chr_name;
	int	*chr_length;
	int32_t	l_text;
	char	*text;
	int32_t	n_ref;
	uint8_t *address;
	char	filename[name_len];
	char	dirname_region[name_len];

	int	flag_header=0;
	int	flag_in = 0;
	int	All_Coverage;
	int	chr_ID = -2;
	uint64_t	offset;

	char	BAMI[4];
	int32_t	n_ref_bai;
	int32_t	n_bin;
	uint32_t	bin;
	int32_t	n_chunk;
	uint64_t	chunk_beg;
	uint64_t	chunk_end;
	int32_t	n_intv;
	uint64_t	*ioffset;
	
	uint64_t	offset_beg;
	uint64_t	offset_bgzf;
	uint64_t	offset_decomp;
	int64_t	diff;
	int	region2bin;

	stream_i= calloc(65536,sizeof(uint8_t));
	stream_o= calloc(65536,sizeof(uint8_t));
	buffer	= calloc(65536*2,sizeof(uint8_t));	


	/* Initial*/
	stream_i[0] = 120;
	stream_i[1] = 156;
	/*---------------*/


	len_data = BGFZBlock(file_bam_i);
	decompressBlock(&infstream, stream_i, stream_o, len_data, file_bam_i);
	fread(&BgfzTail	,sizeof(bgfzTail),1,file_bam_i);

	memcpy(buffer,stream_o,65536);
	address = buffer;	
	
	memcpy(magic,address,sizeof(char)*4);
	memcpy(&l_text,address+4,sizeof(int32_t));
	text = calloc(l_text+1,sizeof (char));
	memcpy(text,address+8,sizeof(char)*l_text);
	memcpy(&n_ref,address+8+l_text,sizeof(int32_t));

	if (flag_header){
		printf("%d\n%s%d\n",l_text,text,n_ref);	
	}
	address += (12+l_text);
	chr_name	= malloc(n_ref*sizeof(char *));	
	chr_length	= calloc(n_ref,sizeof(int));
	for (i = 0;i < n_ref;i++){
		chr_name[i]	= calloc(name_len,sizeof(char));
		address += refInformation(address,chr_name[i],&chr_length[i]);
		if (strcmp(chromosome, chr_name[i]) == 0){
			chr_ID = i;
			if( ToolsFlags->start > ToolsFlags->end || ToolsFlags->end > chr_length[i]){
				printf("Error\n");
				return -1;	
			}
		}
	}
	
	if (ToolsFlags->flag_hide != 1){
		for (i = 0;i < n_ref;i++){
			printf("%s\t%d\n",chr_name[i], chr_length[i]);
		}
	}
	top = buffer + BgfzTail.I_size;
	counter	= top - address;

	fread(BAMI,4, sizeof(char), file_bai_i);
	fread(&n_ref_bai,1,sizeof(int32_t),file_bai_i);
	for(i = 0;i < n_ref_bai;i++){
		fread(&n_bin,1,sizeof(int32_t),file_bai_i);
		for (j = 0;j < n_bin;j++){
			fread(&bin,1,sizeof(uint32_t),file_bai_i);
			fread(&n_chunk,1,sizeof(int32_t),file_bai_i);
			for(m = 0;m < n_chunk;m++){
				fread(&chunk_beg,1,sizeof(uint64_t),file_bai_i);
				fread(&chunk_end,1,sizeof(uint64_t),file_bai_i);
			}
		}
		fread(&n_intv, 1, sizeof(int32_t), file_bai_i);
		ioffset = calloc(sizeof(uint64_t), n_intv);
		fread(ioffset, n_intv, sizeof(uint64_t), file_bai_i);
		if (chr_ID == i){
			offset_bgzf = ioffset[ToolsFlags->start/16384] >> 16;
			offset_decomp = ioffset[ToolsFlags->start/16384]&65535;
		}
		if (i == chr_ID){
			break;	
		}
	}

	fseek(file_bam_i,offset_bgzf,SEEK_SET);
	len_data = BGFZBlock(file_bam_i);
	decompressBlock(&infstream, stream_i, stream_o, len_data, file_bam_i);
	fread(&BgfzTail	,sizeof(bgfzTail),1,file_bam_i);
	memcpy(buffer,stream_o,BgfzTail.I_size);
	address = buffer + offset_decomp;
	counter = BgfzTail.I_size - offset_decomp;
//	printf("%lu\t%lu\t%lu\n", offset_bgzf, offset_decomp, counter);

	
	//Start Unzip and Process Bam File	
	while ( (len_data = BGFZBlock(file_bam_i)) > 0 ){
		decompressBlock(&infstream, stream_i, stream_o, len_data, file_bam_i);
		fread(&BgfzTail	,sizeof(bgfzTail),1,file_bam_i);

		memcpy(buffer, address, sizeof(uint8_t)*counter);
		memcpy(&buffer[counter], stream_o, BgfzTail.I_size);
		top	= buffer + counter + BgfzTail.I_size;
		address = buffer;
		if (top - address > sizeof(alignmentHeader)){
			memcpy(&AlignmentHeader, address, sizeof(alignmentHeader));
		}
		while ( top - address >= (AlignmentHeader.block_size + 4) ){
		//	printf("Johnnash\n");
			if (AlignmentHeader.refID != ref_ID){
				if (ref_ID != -1){
					if (ToolsFlags->flag_hide == 0){
						printf("[Bam File Unzip %d / %d ] %s done\n",ref_ID+1,n_ref,chr_name[ref_ID]);
					}	
				}
				ref_ID = AlignmentHeader.refID;
			}
			if (ref_ID == -1){	break;	}
			address += sizeof(alignmentHeader);
		
			//Core_START
			if (strcmp(chromosome, chr_name[ref_ID]) == 0){
				alignment_Detect(address, &AlignmentHeader, chr_length[ref_ID], ToolsFlags);
				flag_in = 1;
			}else if (flag_in == 1){
				break;	
			}
			//Core_END
			
			address += (AlignmentHeader.block_size + 4) - sizeof(alignmentHeader);
				if (top - address > sizeof(alignmentHeader)){
				memcpy(&AlignmentHeader, address, sizeof(alignmentHeader));
			}else {
				break;	
			}
		}
		counter = top - address;
		if (ref_ID == -1){	break;	}
		if (strcmp(chromosome, chr_name[ref_ID]) != 0 && flag_in == 1){	break;	}
		if (AlignmentHeader.pos > ToolsFlags->end && flag_in == 1){	break;	}
	}
	if (ref_ID != -1){
		if (ToolsFlags->flag_hide == 0){
			printf("[Bam File Unzip %d / %d ] %s done\n",ref_ID+1,n_ref,chr_name[ref_ID]);
		}
	}
	printf("%u\t%f\n", ToolsFlags->num_ref+ToolsFlags->num_alt, (float)ToolsFlags->num_alt /(ToolsFlags->num_ref+ToolsFlags->num_alt));	
	
	free(stream_i);	
	free(stream_o);	
	free(buffer);	
	return 0;
}

void	alignment_Detect(uint8_t *stream, alignmentHeader *AlignmentHeader, uint32_t ref_length, toolsFlags *ToolsFlags){

	int	i,j,k;

	char	*read_name;
	uint32_t *cigar;
	uint8_t	*seq;
	char	*qual;

	int	op;
	int	op_len;
	char	residue;
	int	length;
	int	index;
	int	flag = 0;
	MemoryCopy(&stream,(void **)&read_name	, AlignmentHeader->l_read_name	, sizeof(char));
	MemoryCopy(&stream,(void **)&cigar	, AlignmentHeader->n_cigar_op	, sizeof(uint32_t));
	MemoryCopy(&stream,(void **)&seq	, (AlignmentHeader->l_seq+1)/2	, sizeof(uint8_t));

	length = 0;
	index = AlignmentHeader->pos;
	flag = 0;
	for (i = 0;i < AlignmentHeader->n_cigar_op;i++){
		op	= cigar[i] & 7;
		op_len	= cigar[i] >> 4;
		if (CIGAR(op) == 'M'){
			if (index + op_len > ref_length){
				op_len = ref_length - index;
			}
			if (ToolsFlags->type == 0 && index <= ToolsFlags->start && ToolsFlags->start < index+op_len){
				length += (ToolsFlags->start-index);
				residue = Bin2SeqTop(seq[length >> 1],length&1);
				if (residue == ToolsFlags->alt){
					ToolsFlags->num_alt = ToolsFlags->num_alt + 1;	
				}else {
					ToolsFlags->num_ref = ToolsFlags->num_ref + 1;	
				}
				break;
			}else if (ToolsFlags->type == 1){
				if (index <= ToolsFlags->start && op_len >= ToolsFlags->length){
					 ToolsFlags->num_ref = ToolsFlags->num_ref + 1;	
				}else if (index <= ToolsFlags->start){
					flag = 1;
				}else if (index+op_len > ToolsFlags->end && flag == 1){
					 ToolsFlags->num_ref = ToolsFlags->num_ref + 1;	
				}
			}
			index += op_len;
			length += op_len;
		}else if (CIGAR(op) == 'I'){
			length += op_len;
		}else if (CIGAR(op) == 'D'){
			if (ToolsFlags->type == 1 && ToolsFlags->start == index && ToolsFlags->length == op_len){
				ToolsFlags->num_alt = ToolsFlags->num_alt + 1;
			}
			index+=op_len;
		}else if (CIGAR(op) == 'S'){
			length += op_len;
		}else if (CIGAR(op) == 'N'){
			if (ToolsFlags->type == 0 && index <= ToolsFlags->start && ToolsFlags->start < index+op_len){
				ToolsFlags->num_ref = ToolsFlags->num_ref + 1;	
				break;
			}
			index += op_len;
		}
		// CIGAR(op) == 'H' don't care about it.
		
	}
	
	free(read_name);
	free(cigar);
	free(seq);
}
