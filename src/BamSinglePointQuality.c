#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include "zlib.h"

#include "BamCommonLibrary.h"


void	alignment_SinglePointQuality(uint8_t *stream, alignmentHeader *AlignmentHeader, uint32_t *PosQuality, uint32_t position , uint32_t ref_length, char ref, uint8_t flag_mapq);

int	BamSinglePointQuality(FILE *file_bam_i, FILE *file_bai_i, char *chromosome, uint32_t position, toolsFlags *ToolsFlags , char ref){

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
	int	flag_cov = 0;
	int	flag_sam = 0;
	int	flag_region = 0;
	int	flag_region_file = 0;
	int	flag_dirname_region= 0;
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
	uint32_t	PosQuality[128*8];
	uint32_t	SumQuality[8];
	uint32_t	sumQuality = 0;
	
//	PosQuality = calloc(256,sizeof(uint32_t));
	memset(PosQuality,0,sizeof(uint32_t)*128*8);
	
	stream_i= calloc(65536,sizeof(uint8_t));
	stream_o= calloc(65536,sizeof(uint8_t));
	buffer	= calloc(65536*2,sizeof(uint8_t));	


	/* Initial*/
	stream_i[0] = 120;
	stream_i[1] = 156;
	/*---------------*/




	//Modify
		position--;
	//





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
			if( position > chr_length[i]){
				printf("[Error] Position is bigger than Size of %s\n",chromosome);
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
//			printf("%u\t%u\n",ioffset[0], position/16384);
			offset_bgzf = ioffset[position/16384] >> 16;
			offset_decomp = ioffset[position/16384]&65535;
	//		printf("OFF:%lu\t%lu\t%u\n", ioffset[position/16384], ioffset[position/16384] >> 16, ioffset[position/16384]&65535);
		}
		if (i == chr_ID){
			break;	
		}
	}

	//Offset Bam File
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

			if (AlignmentHeader.refID != ref_ID){
				if (ref_ID != -1){
					if (ToolsFlags->flag_hide == 0){
						printf("[Bam File Unzip %d / %d ] %s done\n",ref_ID+1,n_ref,chr_name[ref_ID]);
					}	
				}
				ref_ID = AlignmentHeader.refID;
			}

			if ((ref_ID == -1) || (strcmp(chromosome, chr_name[ref_ID]) != 0) || (AlignmentHeader.pos > position )){
				break;	
			}


			
			address += sizeof(alignmentHeader);
		
			//Core_START
			if ((AlignmentHeader.FLAG&1024)==0 || ToolsFlags->flag_dup==0 ){
				if (ToolsFlags->flag_mapq == 1){
					alignment_SinglePointQuality(address, &AlignmentHeader, PosQuality, position, chr_length[ref_ID], ref, 1);
				}else {
					alignment_SinglePointQuality(address, &AlignmentHeader, PosQuality, position, chr_length[ref_ID], ref, 0);
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
		counter = top - address;
		if ((ref_ID == -1) || (strcmp(chromosome, chr_name[ref_ID]) != 0) || (AlignmentHeader.pos > position )){
			break;	
		}
	}
	if (ref_ID != -1){
		if (ToolsFlags->flag_hide == 0){
			printf("[Bam File Unzip %d / %d ] %s done\n",ref_ID+1,n_ref,chr_name[ref_ID]);
		}
	}
	

	for(i = 0;i < 64;i++){
		printf("%d\t",i);
		for (j = 0;j < 8;j++){
			printf("%u\t",PosQuality[i+j*128]);
			SumQuality[j] += i*PosQuality[i+j*128];
			sumQuality += i*PosQuality[i+j*128];
		}
		printf("\n");
	}
	printf("SUM\t");
	for (j = 0;j < 8;j++){
		printf("%u\t",SumQuality[j]);
	}
	printf("\n");




//	for(i = 128;i < 256;i++){
//		printf("%u\t",PosQuality[i]);
//		sumQuality += (i-128)*PosQuality[i];
//	}
//	printf("\n");
	printf("%d\n", sumQuality);
	
	free(stream_i);	
	free(stream_o);	
	free(buffer);	
	return 0;
}

void	alignment_SinglePointQuality(uint8_t *stream, alignmentHeader *AlignmentHeader, uint32_t *PosQuality, uint32_t position, uint32_t ref_length, char ref, uint8_t flag_mapq){

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
	uint8_t	strands;

	uint32_t	baseQuality;
	
	MemoryCopy(&stream,(void **)&read_name	, AlignmentHeader->l_read_name	, sizeof(char));
	MemoryCopy(&stream,(void **)&cigar	, AlignmentHeader->n_cigar_op	, sizeof(uint32_t));
	MemoryCopy(&stream,(void **)&seq	, (AlignmentHeader->l_seq+1)/2	, sizeof(uint8_t));
	MemoryCopy(&stream,(void **)&qual	, AlignmentHeader->l_seq	, sizeof(char));

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
			if (index <= position && position < index+op_len){
				length += (position-index);
				residue = Bin2SeqTop(seq[length >> 1],length&1);
			
				//Care Mapping Qualuty
				if (flag_mapq == 1 && AlignmentHeader->MAPQ < qual[length]){
					baseQuality = AlignmentHeader->MAPQ;
				}else {
					baseQuality = qual[length];
				}
				//Care Strand
				if ((AlignmentHeader->FLAG&16)>>4){
					baseQuality += 512;
				}


				
				if ( residue == 'A'){
					PosQuality[baseQuality] += 1;
				}else if ( residue == 'C'){
					PosQuality[baseQuality+128] += 1;
				}else if ( residue == 'G'){
					PosQuality[baseQuality+256] += 1;
				}else if ( residue == 'T'){
					PosQuality[baseQuality+384] += 1;
				}else {

				}

				break;
			}
			index += op_len;
			length += op_len;
		}else if (CIGAR(op) == 'I'){
			length += op_len;
		}else if (CIGAR(op) == 'D'){
			index += op_len;
		}else if (CIGAR(op) == 'S'){
			length += op_len;
		}else if (CIGAR(op) == 'N'){
			index += op_len;
		}
		// CIGAR(op) == 'H' don't care about it.
		
	}
	
	free(read_name);
	free(cigar);
	free(seq);
}
