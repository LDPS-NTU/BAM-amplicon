#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <inttypes.h>
#include <ctype.h>
#include "zlib.h"

#include "BamCommonLibrary.h"
uint32_t taglen(uint8_t *stream, char type);
uint32_t print_tag(uint8_t *stream, char type);
//void	banded_smith_waterman(char *ref, char *seq, uint8_t *qual);
void	banded_smith_waterman(char *ref, char *seq, char *qual, uint16_t *flow_base_forward, uint16_t *flow_base_reverse, uint16_t *flow_ori);
void	print_alignment(char *ref, char *alt, uint32_t *cigar, uint16_t n_cigar_op);
void	dna_str_reverse(char *string);
void	dna_flow_base_reverse(uint16_t *flow_base, int length);
void	strrev(char s[]);
uint32_t FZ_2_ZM(uint8_t *stream, char type);
faiTable *_FaiTable_construct(FILE *file_fai);
void	create_new_sequence(int	start_flow, char *flow_order, uint16_t *flow_signal, uint32_t array_length, char *alt, char *qual, char *alt_new, char *qual_new, int alt_index, int alt_length, uint16_t *flow_base_forward, uint16_t *flow_base_reverse, uint16_t *flow_ori);
void	extract_ref_sequence(FILE *file_fasta_i, faiTable *FaiTable, char *ref_name, uint32_t start, uint32_t end, char *ref);

void	alignment_Test(uint8_t *stream, alignmentHeader *AlignmentHeader, uint32_t ref_length, char *flow_order, char *ref_name, FILE *file_fasta_i, faiTable *FaiTable, toolsFlags *ToolsFlags);

int	BamTest(FILE *file_bam_i, FILE *file_bai_i, FILE *file_bed_i, FILE *file_fasta_i, FILE *file_fai_i, toolsFlags *ToolsFlags){

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
	faiTable	*FaiTable; 

	uint64_t	offset_beg;
	uint64_t	offset_bgzf;
	uint64_t	offset_decomp;

	alignmentHeader	AlignmentHeader;

	char	magic[4];
	int32_t	l_text;
	char	*text;
	int32_t	n_ref;
	uint8_t *address;
	char	line[LINE_MAX_LEN];
	char	dirname_region[name_len];

	char	*flow_order;

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


	FaiTable = calloc(sizeof(faiTable), 1);

	while (fgets(line, 1000, file_fai_i) != NULL){
		FaiTable->num_chr++;
	}
	
	FaiTable->chr_name	= calloc(sizeof(char *), FaiTable->num_chr);
	for(i = 0; i < FaiTable->num_chr; i++){
		FaiTable->chr_name[i]	= calloc(sizeof(char), 100);
	}
	FaiTable->chr_length	= calloc(sizeof(uint32_t), FaiTable->num_chr);
	FaiTable->index	= calloc(sizeof(uint32_t), FaiTable->num_chr);
	FaiTable->line_len	= calloc(sizeof(uint32_t), FaiTable->num_chr);
	FaiTable->line_len_real	= calloc(sizeof(uint32_t), FaiTable->num_chr);
			
	fseek(file_fai_i, 0, SEEK_SET);
	for(i = 0; i < FaiTable->num_chr; i++){
		fscanf(file_fai_i, "%s\t%u\t%u\t%u\t%u", FaiTable->chr_name[i], &FaiTable->chr_length[i], &FaiTable->index[i], &FaiTable->line_len[i], &FaiTable->line_len_real[i]);
	}


	top = Catch_BamHeader(file_bam_i, &BamHeader, stream_i, stream_o, buffer, ToolsFlags);
	//top = CatchBamHeader(file_bam_i, &BamHeader, stream_i, stream_o, buffer,top);
	counter	= top - buffer;
	address = buffer;

	printf("%s",BamHeader.text);

	while (strncmp("@RG", BamHeader.text, 3) != 0){
		BamHeader.text = strchr(BamHeader.text, '\n') + 1;
	}
	while (strncmp("FO:", BamHeader.text, 3) != 0){
		BamHeader.text = strchr(BamHeader.text, '\t') + 1;
	}

	flow_order	= strtok(BamHeader.text+3, "\t");
	//printf("%s\n", flow_order);
	
	if (ToolsFlags->flag_hide != 1){
		for (i = 0;i < BamHeader.n_ref;i++){
			printf("%s\t%d\n",BamHeader.chr_name[i], BamHeader.chr_length[i]);
		}
	}
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
					//printf("%d\t%d\t%d\n",j, AlignmentHeader.refID, AlignmentHeader.pos);
					if (AlignmentHeader.refID != ref_ID){
						flag_break = 1;
						break;
					}
					address += sizeof(alignmentHeader);	
					//Core_START
					alignment_Test(address, &AlignmentHeader, BamHeader.chr_length[ref_ID], flow_order, BamHeader.chr_name[ref_ID], file_fasta_i, FaiTable, ToolsFlags);
					//Core_END
					address += (AlignmentHeader.block_size + 4) - sizeof(alignmentHeader);

					if (top - address > sizeof(alignmentHeader)){
						memcpy(&AlignmentHeader, address, sizeof(alignmentHeader));
					}else {
						break;	
					}
				}

				if (flag_break	== 1){	flag_break = 0;	break;
				}else {			counter = top - address;}
			}
		}
		if (ToolsFlags->flag_hide == 0){
			printf("[Bam File Unzip %d / %d ] %s done\n",ref_ID+1,BamHeader.n_ref,BamHeader.chr_name[ref_ID]);
		}
	}
	
	free(stream_i);	
	free(stream_o);	
	free(buffer);	
	return 0;
}

void	alignment_Test(uint8_t *stream, alignmentHeader *AlignmentHeader, uint32_t ref_length, char *flow_order, char* ref_name, FILE *file_fasta_i, faiTable *FaiTable, toolsFlags *ToolsFlags){

	int	i,j;
	char	*read_name;
	uint32_t *cigar;
	uint8_t	*seq;
	char	*qual;
	char	*tok;

	int	op;
	int	op_len;
	char	residue;
	uint32_t	index;
	int	flag = 0;

	FILE	*file_p;

	int	alt_index	= 0;
	int	alt_length	= 0;

	char	*alt;
	char	*ref;
	char	*alt_new;
	char	*qual_new;
	uint16_t	*flow_base_forward;
	uint16_t	*flow_base_reverse;
	uint16_t	*flow_ori;

	uint32_t	start_position = 0;
	uint32_t	end_position;
	char	tag[3];
	char	type;
	char	sub_type;
	uint16_t	*flow_signal;
	uint32_t	array_length;
	uint32_t	block_size	= (AlignmentHeader->block_size - sizeof(alignmentHeader) + 4);
	uint8_t	*stream_temp	= stream;
	uint8_t	*stream_temp_2;
	uint8_t	flag_FZ;
	uint8_t	tag_NM	= -1;
	uint32_t start_flow;
	uint32_t	m, n, k;

	MemoryCopy(&stream,(void **)&read_name	, AlignmentHeader->l_read_name	, sizeof(char));
//	stream	+= AlignmentHeader->l_read_name;
	MemoryCopy(&stream,(void **)&cigar	, AlignmentHeader->n_cigar_op	, sizeof(uint32_t));
//	stream	+= AlignmentHeader->n_cigar_op*4;
	MemoryCopy(&stream,(void **)&seq	, (AlignmentHeader->l_seq+1)/2	, sizeof(uint8_t));
//	stream	+= (AlignmentHeader->l_seq+1)/2;
	MemoryCopy(&stream,(void **)&qual	, AlignmentHeader->l_seq	, sizeof(char));
//	stream	+= AlignmentHeader->l_seq;

	stream_temp_2	= stream;

	alt	= calloc(AlignmentHeader->l_seq + 1, sizeof(char));
	//Sequence Extract
	for (i = 0; i < AlignmentHeader->l_seq; i++){	alt[i]	= Bin2SeqTop(seq[i >> 1],i&1);	}
	for (i = 0; i < AlignmentHeader->l_seq; i++){	qual[i]	+= 33;	}



	printf("%s\t", read_name);
	printf("%d\t", AlignmentHeader->FLAG);
	printf("%s\t", ref_name);
	printf("%d\t", AlignmentHeader->pos+1);
	printf("%d\t", AlignmentHeader->MAPQ);

	if (ToolsFlags->flag_flow == 0){	/* Only ZM */
		for (i = 0;i < AlignmentHeader->n_cigar_op;i++){
			op	= cigar[i] & 7;
			op_len	= cigar[i] >> 4;
			printf("%d%c",op_len, CIGAR(op));
		}
		printf("\t*\t0\t0\t%s\t%s", alt, qual);
		stream	= stream_temp_2;
		while (block_size > (stream - stream_temp)){
			printf("\t");
			type = stream[2];
			if (strncmp(stream, "FZ", 2) == 0){
				stream	+= FZ_2_ZM(stream,type);
				printf("\tZP:B:f,0.007,0.007,0.000001");
			}else {
				stream += print_tag(stream, type);
			}
		}
		printf("\n");
		return;
	}else {
		//Check NM
		while (block_size > (stream - stream_temp)){
			type = stream[2];
			if (strncmp(stream, "NM", 2) == 0){
				tag_NM	= stream[3];
				break;
			}else {
				stream += taglen(stream, type);
			}
		}

		if (tag_NM == 0){
			for (i = 0;i < AlignmentHeader->n_cigar_op;i++){
				op	= cigar[i] & 7;
				op_len	= cigar[i] >> 4;
				printf("%d%c",op_len, CIGAR(op));
			}
			printf("\t*\t0\t0\t%s\t%s", alt, qual);
			stream	= stream_temp_2;
			while (block_size > (stream - stream_temp)){
				printf("\t");
				type = stream[2];
				if (strncmp(stream, "FZ", 2) == 0 && ToolsFlags->flag_zm == 1){
					stream	+= FZ_2_ZM(stream,type);
					printf("\tZP:B:f,0.007,0.007,0.000001");
				}else {
					stream += print_tag(stream, type);
				}
			}
			printf("\n");
		}else {
			flag_FZ	= 0;
			start_flow	= 0;
			stream	= stream_temp_2;
			while (block_size > (stream - stream_temp) && flag_FZ < 2){
				type = stream[2];
				if (strncmp(stream, "FZ", 2) == 0){
					memcpy(&array_length, stream+4, sizeof(uint32_t));
					flow_signal	= malloc (array_length * sizeof(uint16_t));
					memcpy(flow_signal, stream+8, array_length*sizeof(uint16_t));
					stream	+= (8+array_length*2);
					flag_FZ++;
				}else if (strncmp(stream, "ZF", 2) == 0){
					memcpy(&start_flow, stream+3, sizeof(uint8_t));
					stream	+= 4;
					flag_FZ++;
				}else {
					stream += taglen(stream, type);
				}
			}
			alt_new	= calloc(2*AlignmentHeader->l_seq, sizeof(char));
			qual_new	= calloc(2*AlignmentHeader->l_seq, sizeof(char));
			flow_base_forward	= calloc(2*AlignmentHeader->l_seq, sizeof(uint16_t));
			flow_base_reverse	= calloc(2*AlignmentHeader->l_seq, sizeof(uint16_t));
			flow_ori	= calloc(2*AlignmentHeader->l_seq, sizeof(uint16_t));
		
			//Alignment End and Sequence Start
			alt_index	= 0;
			alt_length	= 0;
			end_position	= AlignmentHeader->pos;
			start_position	= AlignmentHeader->pos;
			for (i = 0;i < AlignmentHeader->n_cigar_op;i++){
				op	= cigar[i] & 7;
				op_len	= cigar[i] >> 4;
				switch(CIGAR(op)){
				case 'M':	
					end_position +=op_len;	
					alt_length += op_len;
					break;
				case 'D':	end_position +=op_len;	break;
				case 'I':	alt_length += op_len;	break;
				case 'S':	
					if (end_position == AlignmentHeader->pos && (AlignmentHeader->FLAG&16)==0){ 
						alt_index = op_len;
						alt_length += op_len;
					}else if (end_position != AlignmentHeader->pos && (AlignmentHeader->FLAG&16)!=0){
						alt_index = op_len;
						alt_length += op_len;
					}	
					break;
				default:	break;	
				}
			}
		
			ref	= calloc(end_position - start_position + 2, sizeof(char));
			extract_ref_sequence(file_fasta_i, FaiTable, ref_name, start_position, end_position, ref);
		
			//Original Alignment	
			//print_alignment(ref, alt, cigar, AlignmentHeader->n_cigar_op);
		
		
			if ((AlignmentHeader->FLAG&16)!=0){
				strrev(qual);
				dna_str_reverse(alt);
			}
		
			create_new_sequence(start_flow, flow_order, flow_signal, array_length, alt, qual, alt_new, qual_new, alt_index, alt_length, flow_base_forward, flow_base_reverse, flow_ori);
		
			if ((AlignmentHeader->FLAG&16)!=0){
				strrev(qual_new);
				dna_str_reverse(alt_new);
				dna_flow_base_reverse(flow_base_forward, strlen(alt_new));
				dna_flow_base_reverse(flow_base_reverse, strlen(alt_new));
			}
		
			banded_smith_waterman(ref, alt_new, qual_new, flow_base_forward, flow_base_reverse, flow_ori);
		
			stream	= stream_temp_2;
			while (block_size > (stream - stream_temp)){
				printf("\t");
				type = stream[2];
				if (strncmp(stream, "FZ", 2) == 0 && ToolsFlags->flag_zm == 1){
					stream	+= FZ_2_ZM(stream,type);
					printf("\tZP:B:f,0.007,0.007,0.000001");
		
				}else {
					stream	+= print_tag(stream, type);
				}
			}
			printf("\n");

			free(flow_signal);
			free(alt_new);
			free(qual_new);
			free(flow_base_forward);
			free(flow_base_reverse);
			free(flow_ori);
		}
	}

	free(read_name);
	free(cigar);
	free(seq);
	free(alt);
	free(qual);
}

void	dna_flow_base_reverse(uint16_t *qual, int length){
	int	i,j;
	uint16_t	temp;
	j	= length-1;
	for (i = 0; i < (length+1)/2; i++){
		temp	= qual[i];
		qual[i]	= qual[j];
		qual[j]	= temp;
		j--;
	}
}
void	dna_str_reverse(char *string){
	int	i,j;
	char	*reverse;
	int	length	= strlen(string);
	reverse	= calloc(length+1, sizeof(char));

	j = 0;
	for (i = length - 1; i >=0 ;i--){
		switch (string[i]){
		case 'A':	reverse[j]	= 'T';	break;
		case 'C':	reverse[j]	= 'G';	break;
		case 'G':	reverse[j]	= 'C';	break;
		case 'T':	reverse[j]	= 'A';	break;
		default:	printf("Error\n");	exit(-1);
		}
		j++;
	}
	strcpy(string, reverse);
	free(reverse);
}




uint32_t taglen(uint8_t *stream, char type){

	char	sub_type;
	uint32_t	array_length;

	switch (type){
	case 'A':	break;
	case 'c':	return 4;
	case 'C':	return 4;
	case 's':	return 5;
	case 'S':	return 5;
	case 'i':	return 7;
	case 'I':	return 7;
	case 'f':	return 7;
	case 'Z':	return (strlen(strtok(stream,"\n"))+1);
	case 'H':	printf("Tag Error: H (Please contact d01943008@ntu.edu.tw)\n");      exit(-1);
	case 'B':	
		sub_type = stream[3];
		memcpy(&array_length, stream+4, sizeof(uint32_t));
		stream	+= 8;
		switch (sub_type){
		case 'c':	return (8+array_length);
		case 'C':	return (8+array_length);
		case 's':	return (8+2*array_length);
		case 'S':	return (8+2*array_length);
		case 'i':	return (8+4*array_length);
		case 'I':	return (8+4*array_length);
		case 'f':	return (8+4*array_length);
		default:	printf("Tag Error: B (Please contact d01943008@ntu.edu.tw)\n");	exit(-1);
		}
	default:	printf("Error\n");	exit(-1);
	}
}

uint32_t print_tag(uint8_t *stream, char type){
	int	i;
	char	sub_type;
	uint32_t	array_length;
	//printf("%.2s:%c:", stream, type);
	printf("%.2s:", stream, type);
	switch (type){
	case 'A':	break;
	case 'c':	printf("i:%d", ((int8_t *)(stream+3))[0]);	return 4;
	case 'C':	printf("i:%u", ((uint8_t *)(stream+3))[0]);	return 4;
	case 's':	printf("i:%d", ((int16_t *)(stream+3))[0]);	return 5;
	case 'S':	printf("i:%u", ((uint16_t *)(stream+3))[0]);	return 5;
	case 'i':	printf("i:%d", ((int32_t *)(stream+3))[0]);	return 7;
	case 'I':	printf("i:%u", ((uint32_t *)(stream+3))[0]);	return 7;
	case 'f':	printf("f:%f", ((float *) (stream+3))[0]);	return 7;
	case 'Z':	printf("Z:%s", stream+3);		return (strlen(strtok(stream,"\n"))+1);
	case 'H':	printf("Tag Error: H (Please contact d01943008@ntu.edu.tw)\n");      exit(-1);
	case 'B':	
		sub_type = stream[3];
		printf("B:%c",sub_type);
		memcpy(&array_length, stream+4, sizeof(uint32_t));
//		stream	+= 8;
		switch (sub_type){
		case 'c':	//printf(",%d",((int8_t *)(stream+8))[0]); 
				for(i=0;i<array_length;i++){printf(",%d",((int8_t *)(stream+8))[i]);} 	return (8+array_length);
		case 'C': 	//printf(",%u",((uint8_t *)(stream+8))[0]); 
				for(i=0;i<array_length;i++){printf(",%u",((uint8_t *)(stream+8))[i]);}	return (8+array_length);
		case 's':	//printf(",%d",((int16_t *)(stream+8))[0]); 
				for(i=0;i<array_length;i++){printf(",%d",((int16_t *)(stream+8))[i]);}	return (8+2*array_length);
		case 'S':	//printf(",%u",((uint16_t *)(stream+8))[0]); 
				for(i=0;i<array_length;i++){printf(",%u",((uint16_t *)(stream+8))[i]);}	return (8+2*array_length);
		case 'i':	//printf(",%d",((int32_t *)(stream+8))[0]); 
				for(i=0;i<array_length;i++){printf(",%d",((int32_t *)(stream+8))[i]);}	return (8+4*array_length);
		case 'I':	//printf(",%u",((uint32_t *)(stream+8))[0]);
				for(i=0;i<array_length;i++){printf(",%u",((uint32_t *)(stream+8))[i]);}	return (8+4*array_length);
		case 'f':	//printf(",%f",((float *)(stream+8))[0]); 
				for(i=0;i<array_length;i++){printf(",%f",((float *)(stream+8))[i]);}	return (8+4*array_length);
		default:	printf("Tag Error: B (Please contact d01943008@ntu.edu.tw)\n");	exit(-1);
		}
	default:	printf("Error\n");	exit(-1);
	}
}

uint32_t FZ_2_ZM(uint8_t *stream, char type){
	int	i;
	char	sub_type;
	uint32_t	array_length;
	//printf("%.2s:%c:", stream, type);
	printf("ZM:", type);
	switch (type){
	case 'A':	break;
	case 'H':	printf("Tag Error: H (Please contact d01943008@ntu.edu.tw)\n");      exit(-1);
	case 'B':	
		sub_type = stream[3];
		printf("B:s",sub_type);
		memcpy(&array_length, stream+4, sizeof(uint32_t));
//		stream	+= 8;
		switch (sub_type){
//		case 'c':	for(i=0;i<array_length;i++){printf(",%d", 2.56*((int8_t *)(stream+8))[i]);} 	return (8+array_length);
//		case 'C': 	for(i=0;i<array_length;i++){printf(",%u", 2.56*((uint8_t *)(stream+8))[i]);}	return (8+array_length);
		case 's':	for(i=0;i<array_length;i++){printf(",%d", (int16_t) (0.5+2.56*((int16_t *)(stream+8))[i]));}	return (8+2*array_length);
		case 'S':	for(i=0;i<array_length;i++){printf(",%u", (uint16_t) (0.5+2.56*((uint16_t *)(stream+8))[i]));}	return (8+2*array_length);
		case 'i':	for(i=0;i<array_length;i++){printf(",%d", (int32_t) (0.5+2.56*((int32_t *)(stream+8))[i]));}	return (8+4*array_length);
		case 'I':	for(i=0;i<array_length;i++){printf(",%u", (uint16_t) (0.5+2.56*((uint32_t *)(stream+8))[i]));}	return (8+4*array_length);
		case 'f':	for(i=0;i<array_length;i++){printf(",%f", 0.5+2.56*((float *)(stream+8))[i]);}	return (8+4*array_length);
		default:	printf("Tag Error: B (Please contact d01943008@ntu.edu.tw)\n");	exit(-1);
		}
	default:	printf("Error\n");	exit(-1);
	}
}


void	print_alignment(char *ref, char *alt, uint32_t *cigar, uint16_t n_cigar_op){

	int	i,j;
	int	op;
	int	op_len;
	char	*alt_align;
	char	*ref_align;
	int	length	= 0;
	int	length_2	= 0;

	alt_align = calloc(strlen(ref)*2, sizeof(char));
	ref_align = calloc(strlen(ref)*2, sizeof(char));

	for (i = 0;i < n_cigar_op;i++){
		op	= cigar[i] & 7;
		op_len	= cigar[i] >> 4;
//		printf("%d%c ", op_len, CIGAR(op));
		if (CIGAR(op) == 'M'){
			strncat(alt_align, alt+length, op_len);
			strncat(ref_align, ref+length_2, op_len);
			length+=op_len;
			length_2 += op_len;
		}else if (CIGAR(op) == 'I'){
			strncat(alt_align, alt+length, op_len);
			length+=op_len;
			for (j = 0;j < op_len;j++){
				strncat(ref_align, "*", 1);
			}
		}else if (CIGAR(op) == 'D'){
			strncat(ref_align, ref+length_2, op_len);
			length_2 += op_len;
			for (j = 0;j < op_len;j++){
				strncat(alt_align, "-", 1);
			}
		}else if (CIGAR(op) == 'S'){
			length +=op_len;
		}else if (CIGAR(op) == 'N'){
		}
		// CIGAR(op) == 'H' don't care about it.	
	}
	printf("\n%s\n",ref_align);
	printf("%s\n",alt_align);
	for (i = 0; i < strlen(ref_align); i++){
		if (ref_align[i] != alt_align[i]){	printf("#");}else { printf(" ");}
	}
	printf("\n");
	free(ref_align);
	free(alt_align);
}



void	banded_smith_waterman(char *ref, char *seq, char *qual, uint16_t *flow_base_forward, uint16_t *flow_base_reverse, uint16_t *flow_ori){ 
	
	int	i,j,k,w;
	char	*ref_align;
	char	*alt_align;
	char	*alt_new;
	char	*qual_new;


	int	match	= 100;
	int	mismatch	= 400;
	float	gap_o	= 600;
	float	gap_e	= 100;
	int	op;
	int	op_len;

	int	**score_d;
	int	**score_p;
	int	**score_q;
	uint8_t	**direct_d;
	uint8_t	**direct_p;
	uint8_t	**direct_q;
	uint8_t	**direct;

	int	ref_len;
	int	seq_len;

	ref_len	= strlen(ref);
	seq_len	= strlen(seq);


	ref_align	= calloc(ref_len+seq_len, sizeof(char));
	alt_align	= calloc(ref_len+seq_len, sizeof(char));

	alt_new	= calloc(seq_len+1, sizeof(char));
	qual_new = calloc(seq_len+1, sizeof(char));
	
	score_d	= malloc((seq_len+1)*sizeof(int *));
	score_p	= malloc((seq_len+1)*sizeof(int *));
	score_q	= malloc((seq_len+1)*sizeof(int *));
	direct_d	= malloc((seq_len+1)*sizeof(uint8_t *));
	direct_p	= malloc((seq_len+1)*sizeof(uint8_t *));
	direct_q	= malloc((seq_len+1)*sizeof(uint8_t *));

	for (i = 0; i <= seq_len; i++){
		score_d[i]	= malloc((1+ref_len)*sizeof(int));
		score_p[i]	= malloc((1+ref_len)*sizeof(int));
		score_q[i]	= malloc((1+ref_len)*sizeof(int));
		direct_d[i]	= malloc((1+ref_len)*sizeof(uint8_t));
		direct_p[i]	= malloc((1+ref_len)*sizeof(uint8_t));
		direct_q[i]	= malloc((1+ref_len)*sizeof(uint8_t));
	}
	score_d[0][0]	= 0;
	score_p[0][0]	= gap_o;
	score_q[0][0]	= gap_o;

	for (i = 1;i <= seq_len; i++){
		score_d[i][0]	= -(gap_o+gap_e*i);
		score_p[i][0]	= -(gap_o+gap_e*i);
		score_q[i][0]	= -(gap_o+gap_e*i)*10;
		//score_q[i][0]	= -(gap_o+gap_e*i);
		direct_d[i][0]	= 1;
		direct_p[i][0]	= 1;
		direct_q[i][0]	= 1;

	}
	
	for (j = 1;j <= ref_len; j++){
		score_d[0][j]	= -(gap_o+gap_e*j);
		score_p[0][j]	= -(gap_o+gap_e*j)*10;
		score_q[0][j]	= -(gap_o+gap_e*j);
		direct_d[0][j]	= 2;
		direct_p[0][j]	= 2;
		direct_q[0][j]	= 2;
	}

	for (i = 1; i <= seq_len;i++){
		for (j = 1; j <= ref_len;j++){
			if (score_d[i-1][j]-(gap_o*flow_base_forward[i-1]/100+gap_e) > score_p[i-1][j]-gap_e){
				score_p[i][j]	= score_d[i-1][j]-(gap_o*flow_base_forward[i-1]/100+gap_e);
				direct_p[i][j]	= 0;
			}else {
				score_p[i][j]	= score_p[i-1][j]-gap_e;
				direct_p[i][j]	= 1;
			}
			if (score_d[i][j-1]-(gap_o+gap_e) > score_q[i][j-1]-gap_e){
				score_q[i][j]	= score_d[i][j-1]-(gap_o+gap_e);
				direct_q[i][j]	= 0;
			}else {
				score_q[i][j]	= score_q[i][j-1]-gap_e;
				direct_q[i][j]	= 2;
			}
			if (seq[i-1] == ref[j-1]){
				//w = match;
				w = flow_base_forward[i-1];
			}else {
				w = -mismatch;
			}

			if (score_p[i][j] >= score_q[i][j]){
				if (score_d[i-1][j-1]+w > score_p[i][j]){
					score_d[i][j]	= score_d[i-1][j-1]+w;
					direct_d[i][j]	= 0;
				}else {
					score_d[i][j]	= score_p[i][j];
					direct_d[i][j]	= 1;

				}
			}else {
				if (score_d[i-1][j-1]+w > score_q[i][j]){
					score_d[i][j]	= score_d[i-1][j-1]+w;
					direct_d[i][j]	= 0;	
				}else {
					score_d[i][j]	= score_q[i][j];
					direct_d[i][j]	= 2;
				}
			}
		}			
	}
	i = seq_len;
	j = ref_len;
	direct	= direct_d;
	while (i != 0 || j != 0){
		if (direct == direct_p){
			if (direct[i][j] == 0)	{	
				direct = direct_d;
			}else if (direct[i][j] == 1){	
				direct = direct_p;
			}
				if (flow_base_forward[i-1] > 60 && flow_base_reverse[i-1] > 60){
//					printf("%d\t%d\t%d\n", i-1, flow_base_forward[i-1], flow_ori[i-1]);
					strncat(alt_align, &seq[i-1], 1);
					strcat(ref_align, "*");
					strncat(alt_new, &seq[i-1], 1);
					strncat(qual_new, &qual[i-1], 1);
//				printf("%d\t%c\t%c\t%d\t%d\t%d\n", flow_base_forward[i-1], seq[i-1], '*', score_d[i][j], score_p[i][j], score_q[i][j]);
				}
			//	printf("\t\t\t%d\t%d\t%d\n", score_d[i-1][j], score_p[i-1][j], score_q[i-1][j]);
				i--;
		}else if (direct == direct_q){
			if (direct[i][j] == 0){
				direct = direct_d;
			}else if (direct[i][j] == 2){	
				direct = direct_q;
			}
				strcat(alt_align, "-");
				strncat(ref_align, &ref[j-1], 1);
//				printf("\t%c\t%c\t%d\t%d\t%d\n", '-', ref[j-1], score_d[i][j], score_p[i][j], score_q[i][j]);
				//printf("\t\t\t%d\t%d\t%d\n", score_d[i][j-1], score_p[i][j-1], score_q[i][j-1]);
				j--;
		}else {
			if (direct[i][j] == 0){
				strncat(ref_align, &ref[j-1], 1);
				strncat(alt_align, &seq[i-1], 1);

				strncat(alt_new, &seq[i-1], 1);
				strncat(qual_new, &qual[i-1], 1);
//				printf("(%d,%d)\t%d\t%c\t%c\t%d\t%d\t%d\n",i-1, j-1, flow_base_forward[i-1], seq[i-1], ref[j-1], score_d[i][j], score_p[i][j], score_q[i][j]);
				//printf("\t\t\t%d\t%d\t%d\n", score_d[i-1][j-1], score_p[i-1][j], score_q[i][j-1]);
				i--;
				j--;
				direct	= direct_d;
			}else if (direct[i][j] == 1){
				direct	= direct_p;
			}else if (direct[i][j] == 2){
				direct	= direct_q;
			}
		}
	}
	strrev(alt_align);
	strrev(ref_align);
	strrev(alt_new);
	strrev(qual_new);

	op_len	= 0;
	op	= 'M';
	for (i = 0; i < strlen(alt_align); i++){
		if (ref_align[i] != '*' && alt_align[i] != '-'){
			if (op == 'M'){	op_len++;}
			else {
				printf("%d%c", op_len, op);
				op	= 'M';
				op_len	= 1;	
			}
		}else if (ref_align[i] == '*'){
			if (op == 'I'){	op_len++;}
			else {
				printf("%d%c", op_len, op);
				op	= 'I';
				op_len	= 1;	
			}
		}else if (alt_align[i] == '-'){
			if (op == 'D'){	op_len++;}
			else {
				printf("%d%c", op_len, op);
				op	= 'D';
				op_len	= 1;	
			}
		}
	}
	printf("%d%c\t", op_len, op);
	printf("*\t0\t0\t%s\t%s", alt_new, qual_new);

//	printf("\n%s\n", ref_align);
//	printf("%s\n", alt_align);
	for (i = 0; i <= seq_len; i++){
		free(score_d[i]);
		free(score_p[i]);
		free(score_q[i]);
		free(direct_d[i]);
		free(direct_p[i]);
		free(direct_q[i]);
	}
	free(score_d);
	free(score_p);
	free(score_q);
	free(direct_d);
	free(direct_p);
	free(direct_q);
}

void strrev(char s[])
{
	int length = strlen(s) ;
	int c, i, j;

	for (i = 0, j = length - 1; i < j; i++, j--){
		c = s[i];
		s[i] = s[j];
		s[j] = c;
	}
}

faiTable *_FaiTable_construct(FILE *file_fai){
	faiTable *FaiTable;
	int	num_chr;
	int	i;
	char	line[1000];	

	FaiTable	= calloc(sizeof(faiTable), 1);

//	while (EOF != (fscanf(file_fai,"%*[^\n]"), fscanf(file_fai,"%*c"))){ 
	while (fgets(line, 1000, file_fai) != NULL){
		FaiTable->num_chr++;
		num_chr++;	
	}

	FaiTable->chr_name	= calloc(sizeof(char *), num_chr);
	for(i = 0; i < num_chr; i++){
		FaiTable->chr_name[i]	= calloc(sizeof(char), 100);
	}
	FaiTable->chr_length	= calloc(sizeof(uint32_t), num_chr);
	FaiTable->index	= calloc(sizeof(uint32_t), num_chr);
	FaiTable->line_len	= calloc(sizeof(uint32_t), num_chr);
	FaiTable->line_len_real	= calloc(sizeof(uint32_t), num_chr);

	fseek(file_fai, 0, SEEK_SET);
	for(i = 0; i < num_chr; i++){
		fscanf(file_fai, "%s\t%u\t%u\t%u\t%u", FaiTable->chr_name[i], &FaiTable->chr_length[i], &FaiTable->index[i], &FaiTable->line_len[i], &FaiTable->line_len_real[i]);
//		printf("%s\t%u\t%u\t%u\t%u\n", FaiTable.chr_name[i], FaiTable.chr_length[i], FaiTable.index[i], FaiTable.line_len[i], FaiTable.line_len_real[i]);
	}
	return	FaiTable;
}


void	create_new_sequence(int	start_flow, char *flow_order, uint16_t *flow_signal, uint32_t array_length, char *alt, char *qual, char *alt_new, char *qual_new, int alt_index, int alt_length, uint16_t *flow_base_forward, uint16_t *flow_base_reverse, uint16_t *flow_ori){

	int	i, j;
	int	num_base;
	int	alt_index_new	= 0;

	for (i = start_flow; i < array_length; i++){
		num_base = flow_signal[i]/100;			
		for (j = 0;j < num_base; j++){
			if (flow_order[i] == alt[alt_index]){
				alt_new[alt_index_new]	= flow_order[i];
				qual_new[alt_index_new] = qual[alt_index];
				flow_base_forward[alt_index_new]	= 100;
				flow_base_reverse[alt_index_new]	= 100;
				flow_ori[alt_index_new]	= flow_signal[i];

				alt_index_new++;
				alt_index++;
			//printf("%c\t%d\t%d\t%s\n", flow_order[i], flow_signal[i] - num_base*100,alt_index_new, alt_new);	
			}
		}

		if (alt_index >= alt_length){	break;}

		if (i == start_flow && flow_order[i] != alt[alt_index]){
		}else if (flow_signal[i] - num_base*100 > 50 && flow_order[i] == alt[alt_index]){
		//	printf("%c\t%d\t%d\n", flow_order[i], flow_signal[i] - num_base*100,alt_index_new);	
			alt_new[alt_index_new]	= flow_order[i];
			qual_new[alt_index_new] = qual[alt_index];
			flow_ori[alt_index_new]	= flow_signal[i];
			flow_base_forward[alt_index_new]	= flow_signal[i] - num_base*100;
			flow_base_reverse[alt_index_new]	= 100;
			flow_base_reverse[alt_index_new-num_base] = flow_signal[i] - num_base*100;
			alt_index_new++;
			alt_index++;
		}else if (flow_signal[i] - num_base*100 > 25 && flow_signal[i] > 50){
		//	printf("%c\t%d\t%d\n", flow_order[i], flow_signal[i] - num_base*100, alt_index_new);	
			alt_new[alt_index_new]	= flow_order[i];
			qual_new[alt_index_new] = 63;
			flow_ori[alt_index_new]	= flow_signal[i];
			flow_base_forward[alt_index_new]	= flow_signal[i] - num_base*100;
			flow_base_reverse[alt_index_new]	= 100;
			flow_base_reverse[alt_index_new-num_base] = flow_signal[i] - num_base*100;
			alt_index_new++;
		}
		if (alt_index >= alt_length){	break;}
	}

}

void	extract_ref_sequence(FILE *file_fasta_i, faiTable *FaiTable, char *ref_name, uint32_t start, uint32_t end, char *ref){
	int	i, j;
	uint32_t	index;

	for (i = 0; i < FaiTable->num_chr;i++){
		if (strcmp (FaiTable->chr_name[i], ref_name) == 0){
			index	= FaiTable->index[i] + start/FaiTable->line_len[i]*FaiTable->line_len_real[i] + (start % FaiTable->line_len[i]);
	//		printf("%s\t%d\t%d\t\n", FaiTable->chr_name[i], FaiTable->index[i], index);
			break;
		}
	}
	fseek(file_fasta_i, index, SEEK_SET);
	for (j = 0; j < (end - start); j++){
		ref[j]	= fgetc(file_fasta_i);
		if (ref[j] == '\n'){
			j--;	
		}else {
			ref[j] = toupper(ref[j]);	
		}
	}
}
