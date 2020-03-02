#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include "zlib.h"

#include "BamCommonLibrary.h"




void	alignment_InsTxt(uint8_t *stream, alignmentHeader *AlignmentHeader, FILE *file_deletion_o, uint32_t ref_length);

int	BamInsTxt(FILE *file_bam_i, FILE *file_length_o, int flag_hide){

	FILE	*file_cov_o;

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
	uint32_t	*Coverage;

	char	magic[4];
	char	**chr_name;
	int	*chr_length;
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
	int	All_Coverage;

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
	}
	if (flag_hide != 1){
		for (i = 0;i < n_ref;i++){
			printf("%s\t%d\n",chr_name[i], chr_length[i]);
		}
	}
	top = buffer + BgfzTail.I_size;
	counter	= top - address;



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
			if (AlignmentHeader.refID != ref_ID){
				if (ref_ID != -1){
					if (flag_hide == 0){
						printf("[Bam File Unzip %d / %d ] %s done\n",ref_ID+1,n_ref,chr_name[ref_ID]);
					}	
					fclose(file_cov_o);
				}
				if (AlignmentHeader.refID != -1){
					strcat_triple(filename, "", chr_name[AlignmentHeader.refID], "_ins.txt", name_len);
					file_cov_o = fopen(filename,"w");
				}
				ref_ID = AlignmentHeader.refID;
			}
			if (ref_ID == -1){	break;	}
			address += sizeof(alignmentHeader);
		
			//Core_START
			alignment_InsTxt(address, &AlignmentHeader, file_cov_o, chr_length[ref_ID]);
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
	}
	if (ref_ID != -1){
		if (flag_hide == 0){
			printf("[Bam File Unzip %d / %d ] %s done\n",ref_ID+1,n_ref,chr_name[ref_ID]);
		}
		fclose(file_cov_o);
	}
	
	free(stream_i);	
	free(stream_o);	
	free(buffer);	
	return 0;
}

void	alignment_InsTxt(uint8_t *stream, alignmentHeader *AlignmentHeader, FILE *file_deletion_o, uint32_t ref_length){

	int	i,j;

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
			if (flag==0){
				fprintf(file_deletion_o,"%d_%d_",index+1,op_len);
				flag = 1;
			} else {
				fprintf(file_deletion_o,",%d_%d_",index+1,op_len);
			}
			for (j = 0;j < op_len;j++){
				residue = Bin2SeqTop(seq[length >> 1],length&1);
				fprintf(file_deletion_o,"%c", residue);
				length ++;
			}
		}else if (CIGAR(op) == 'D'){
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
	free(read_name);
	free(cigar);
	free(seq);
}




