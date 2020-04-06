#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include "zlib.h"
#include <sys/time.h>



#include "BamCommonLibrary.h"




void	alignment_Pattern(uint8_t *stream, alignmentHeader *AlignmentHeader, uint32_t start, uint32_t end, uint32_t valid_start, uint32_t valid_end, uint32_t ref_length);
void	alignment_Pattern_Origin(uint8_t *stream, alignmentHeader *AlignmentHeader, uint32_t start, uint32_t end, uint32_t valid_start, uint32_t valid_end, uint32_t ref_length, toolsFlags *ToolsFlags);

int	BamPattern(FILE *file_bam_i, FILE *file_bai_i,  FILE *file_bed_i, char *chromosome, uint32_t start, uint32_t end, toolsFlags *ToolsFlags){

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

	bamHeader	BamHeader;
	bedTable	BedTable;
	bedTable	TargetTable;
	baiTable	BaiTable;
	char	line[LINE_MAX_LEN];
	int	num_target;
	char	chromosome_name[LINE_MAX_LEN];
	char	pattern_temp[LINE_MAX_LEN];
		
	struct timeval start_clock, end_clock;
	long long  sum_clock = 0;
	int	num_reads;
	int	end_position;
	int	max_cover_index;
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

	if (ToolsFlags->flag_target){
		decodeBedFile(ToolsFlags->file_target, &TargetTable, &BamHeader);
	}

	if (file_bed_i != NULL){
		if (ToolsFlags->flag_columnName){
			fgets(line, LINE_MAX_LEN, file_bed_i);
		}
		while	(fgets(line, LINE_MAX_LEN, file_bed_i) != NULL){
			num_target++;
		}

		BedTable.table_start = calloc (BamHeader.n_ref, sizeof(uint32_t));
		BedTable.table_end = calloc (BamHeader.n_ref, sizeof(uint32_t));
		BedTable.table_max = calloc (BamHeader.n_ref, sizeof(uint32_t));
		
		BedTable.start = calloc (num_target, sizeof(uint32_t));
		BedTable.end = calloc (num_target, sizeof(uint32_t));
		BedTable.ref_id = calloc (num_target, sizeof(uint32_t));
		BedTable.pattern	= calloc (num_target, sizeof(char *));

		rewind(file_bed_i);
		if (ToolsFlags->flag_columnName){
			fgets(line, LINE_MAX_LEN, file_bed_i);
		}
	
		for(i = 0;i < num_target;i++){
			memset(pattern_temp,0,LINE_MAX_LEN);
			fscanf(file_bed_i, "%s %d %d %*s %s", chromosome_name, &BedTable.start[i], &BedTable.end[i], pattern_temp);
			fgets(line, LINE_MAX_LEN, file_bed_i);	
			BedTable.start[i] -= 1;
			BedTable.pattern[i]	= strdup(pattern_temp);//calloc (500, sizeof(char));
			//printf("Johnnash\t%d\t%u\t%u\t%u\t%u\n",num_target, BedTable.start, BedTable.end, BedTable.ref_id, BedTable.pattern);
				
			//printf("%s\t%d\t%d\t%s\n", chromosome_name, BedTable.start[i], BedTable.end[i], BedTable.pattern[i]);
//			getchar();
			if (i == 0){
				for(j = 0; j < BamHeader.n_ref; j++){
					if (strcmp(chromosome_name, BamHeader.chr_name[j])==0){
						BedTable.table_start[j] = i;
						BedTable.table_max[j] = BedTable.end[i];
						BedTable.ref_id[i] = j;
						break;
					}
				}
			}else if (strcmp(chromosome_name, BamHeader.chr_name[j])!=0){
				BedTable.table_end[j] = i;
				if (BedTable.end[i] > BedTable.table_max[j]){
					BedTable.table_max[j] = BedTable.end[i];
				}
				for(j = 0; j < BamHeader.n_ref; j++){
					if (strcmp(chromosome_name, BamHeader.chr_name[j])==0){
						BedTable.table_start[j] = i;
						BedTable.ref_id[i] = j;
						break;
					}
				}
			}else {
				BedTable.ref_id[i] = j;
				if (BedTable.end[i] > BedTable.table_max[j]){
					BedTable.table_max[j] = BedTable.end[i];
				}	
			}
			if (j == BamHeader.n_ref){
				printf("Error\n");
				return -1;
			}
//		printf("%s %d %d\n", chromosome_name, BedTable.start[i], BedTable.end[i]);
		}
		BedTable.table_end[j] = i;
		
	}else if (end > 0){
		//printf("Region\n");
		BedTable.table_start = calloc (BamHeader.n_ref, sizeof(uint32_t));
		BedTable.table_end = calloc (BamHeader.n_ref, sizeof(uint32_t));
		BedTable.table_max = calloc (BamHeader.n_ref, sizeof(uint32_t));
		
		BedTable.start = calloc (1, sizeof(uint32_t));
		BedTable.end = calloc (1, sizeof(uint32_t));
		BedTable.ref_id = calloc (1, sizeof(uint32_t));
		num_target = 1;
		BedTable.start[0] = start-1;
		BedTable.end[0] = end;
		
		for(j = 0; j < BamHeader.n_ref; j++){
			if (strcmp(chromosome, BamHeader.chr_name[j])==0){
				BedTable.table_start[j] = 0;
				BedTable.table_end[j] = 1;
				BedTable.ref_id[0] = j;
				BedTable.table_max[j] = end;
				if( start > BamHeader.chr_length[j] || end > BamHeader.chr_length[j] || start > end){
					printf("[Error] Position is bigger than Size of %s\n",chromosome);
					return -1;	
				}
				break;
			}
		}
	}else if (start > 0){
		//printf("SinglePoint\n");
		BedTable.table_start = calloc (BamHeader.n_ref, sizeof(uint32_t));
		BedTable.table_end = calloc (BamHeader.n_ref, sizeof(uint32_t));
		BedTable.table_max = calloc (BamHeader.n_ref, sizeof(uint32_t));
		BedTable.start = calloc (1, sizeof(uint32_t));
		BedTable.end = calloc (1, sizeof(uint32_t));
		BedTable.ref_id = calloc (1, sizeof(uint32_t));
		num_target = 1;
		BedTable.start[0] = start-1;
		BedTable.end[0] = start;
		
		for(j = 0; j < BamHeader.n_ref; j++){
			if (strcmp(chromosome, BamHeader.chr_name[j])==0){
				BedTable.table_start[j] = 0;
				BedTable.table_end[j] = 1;
				BedTable.ref_id[0] = j;
				BedTable.table_max[j] = start+1;
				if( start > BamHeader.chr_length[j]){
					printf("[Error] Position is bigger than Size of %s\n",chromosome);
					return -1;	
				}
				break;
			}
		}	
	}else {
		//printf("Total\n");
		BedTable.table_start = calloc (BamHeader.n_ref, sizeof(uint32_t));
		BedTable.table_end = calloc (BamHeader.n_ref, sizeof(uint32_t));
		BedTable.table_max = calloc (BamHeader.n_ref, sizeof(uint32_t));
		
		BedTable.start = calloc (BamHeader.n_ref, sizeof(uint32_t));
		BedTable.end = calloc (BamHeader.n_ref, sizeof(uint32_t));
		BedTable.ref_id = calloc (BamHeader.n_ref, sizeof(uint32_t));
		num_target = BamHeader.n_ref;
			
		for(j = 0; j < BamHeader.n_ref; j++){
				BedTable.start[j] = 0;
				BedTable.end[j] = BamHeader.chr_length[j];
				BedTable.ref_id[j] = j;
				BedTable.table_start[j] = j;
				BedTable.table_end[j] = j+1;
				BedTable.table_max[j] = BamHeader.chr_length[j];
		}		
	}


	BuildBaiTable (file_bai_i, &BaiTable);
	if (!ToolsFlags->flag_pattern){
		if (ToolsFlags->flag_origin){
			printf("PATTERN\n");
		}else{
			printf("ID\tPATTERN\n");
		}
	}

	//BAI File
/*
	for (j = 0; j < BamHeader.n_ref; j++){
		ref_ID = j;
		for (k = BedTable.table_start[ref_ID]; k < BedTable.table_end[ref_ID]; k++){
*/
		for (k = 0; k < num_target; k++){
			ref_ID	= BedTable.ref_id[k];

//		if (BedTable.table_end[ref_ID] != BedTable.table_start[ref_ID]){
			if (file_bed_i != NULL){
				ToolsFlags->pattern	= BedTable.pattern[k];
				ToolsFlags->num_ref	= 0;
				ToolsFlags->num_alt	= 0;
//				num_reads	= 0;
//				printf("%s\t%s\n", ToolsFlags->pattern, BedTable.pattern[k]);
			}
//			printf("%lu\t%lu\n", ReturnOffset(&BaiTable, ref_ID, BedTable.start[k]), CatchOffset(file_bai_i, ref_ID, BedTable.start[k]));	
			offset_beg = ReturnOffset ( &BaiTable, ref_ID, BedTable.start[k]);
//			offset_beg = CatchOffset (file_bai_i, ref_ID, BedTable.start[k]);
//			printf("%lu\n",offset_beg);
//			offset_beg = CatchOffset_Bin (file_bai_i, ref_ID, BedTable.start[BedTable.table_start[ref_ID]]);
			offset_bgzf = offset_beg >> 16;
			offset_decomp = offset_beg&65535;
			
//			printf("%lu\t%lu\n", offset_bgzf, offset_decomp);	
			//Offset Bam File
			if (offset_beg == 0){
				rewind(file_bam_i);
				top = CatchBamHeader (file_bam_i, &BamHeader, stream_i, stream_o, buffer, top);
				counter	= top - buffer;
				address = buffer;

			}else {
				fseek(file_bam_i,offset_bgzf,SEEK_SET);
				len_data = BGFZBlock(file_bam_i);
//				printf("%lu\n",len_data);
				decompressBlock(&infstream, stream_i, stream_o, len_data, file_bam_i);
				fread(&BgfzTail	,sizeof(bgfzTail),1,file_bam_i);
				memcpy(buffer,stream_o,BgfzTail.I_size);
//				printf("%u\n",BgfzTail.I_size);
				address = buffer + offset_decomp;
				counter = BgfzTail.I_size - offset_decomp;
			}
	
			//Start Unzip and Process Bam File	
			while ( (len_data = BGFZBlock(file_bam_i)) > 0 ){

				decompressBlock(&infstream, stream_i, stream_o, len_data, file_bam_i);
				fread(&BgfzTail	,sizeof(bgfzTail),1,file_bam_i);
				if (address - buffer > counter){
					memcpy(buffer, address, sizeof(uint8_t)*counter);
				}else {
					memmove(buffer, address, sizeof(uint8_t)*counter);
				}
				//memcpy(buffer, address, sizeof(uint8_t)*counter);
				memcpy(&buffer[counter], stream_o, BgfzTail.I_size);
				top	= buffer + counter + BgfzTail.I_size;
				address = buffer;


				if (top - address > sizeof(alignmentHeader)){
					memcpy(&AlignmentHeader, address, sizeof(alignmentHeader));
				}
//				printf("%u\n", AlignmentHeader.block_size);
//				printf("%u\n", AlignmentHeader.refID);
				while ( top - address >= (AlignmentHeader.block_size + 4) ){
					if (AlignmentHeader.refID != ref_ID || AlignmentHeader.pos > BedTable.end[k]){
						break;
					}
					address += sizeof(alignmentHeader);
					//Core_START
					if (((AlignmentHeader.FLAG&4)==0) && ((AlignmentHeader.FLAG&1024)==0 || ToolsFlags->flag_dup==1) ){
						if (!ToolsFlags->flag_filter || SelectReads(address, &AlignmentHeader, 30, 30)){
							end_position	= alignment_EndPosition(address, &AlignmentHeader, BamHeader.chr_length[ref_ID]);
							if (end_position > BedTable.start[k] && AlignmentHeader.pos < BedTable.end[k]){
								if (ToolsFlags->flag_target){
									max_cover_index	= concernTargetRegion( &TargetTable, ref_ID, AlignmentHeader.pos, end_position);
									if (max_cover_index >= 0){
										if (ToolsFlags->flag_origin){
											alignment_Pattern_Origin(address, &AlignmentHeader, BedTable.start[k], BedTable.end[k], TargetTable.start[max_cover_index], TargetTable.end[max_cover_index], BamHeader.chr_length[ref_ID], ToolsFlags);
										}else{
											alignment_Pattern(address, &AlignmentHeader, BedTable.start[k], BedTable.end[k], TargetTable.start[max_cover_index], TargetTable.end[max_cover_index], BamHeader.chr_length[ref_ID]);
										}
									}								
								}else {
									if (ToolsFlags->flag_origin){
										alignment_Pattern_Origin(address, &AlignmentHeader, BedTable.start[k], BedTable.end[k], BedTable.start[k], BedTable.end[k], BamHeader.chr_length[ref_ID], ToolsFlags);
									}else{
										alignment_Pattern(address, &AlignmentHeader, BedTable.start[k], BedTable.end[k], BedTable.start[k], BedTable.end[k], BamHeader.chr_length[ref_ID]);
									}
								}
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
				//printf("%u\t%u\t%u\n", BedTable.start[BedTable.table_start[ref_ID]], AlignmentHeader.pos, BedTable.table_max[ref_ID]);
				if (AlignmentHeader.refID != ref_ID || AlignmentHeader.pos > BedTable.end[k]){
					break;
				}
				counter = top - address;
			}
			if (ToolsFlags->flag_hide == 1){
				printf("[Bam File Unzip %d / %d ] %s done\n",ref_ID+1, BamHeader.n_ref, BamHeader.chr_name[ref_ID]);
			}
			if (ToolsFlags->flag_pattern){
				printf("TOT.DP\tPAT.FREQ **(TOT.DP: total depth; PAT.FREQ: pattern frequency)\n");
				if (ToolsFlags->num_ref != 0){
					printf("%u\t%7.5f\n", ToolsFlags->num_ref, (float)ToolsFlags->num_alt /ToolsFlags->num_ref);
				}else {
					printf("%u\t%7.5f\n", 0, 0);
				}
			}
		}
	free(stream_i);	
	free(stream_o);	
	free(buffer);	
	return 0;
}

void	alignment_Pattern(uint8_t *stream, alignmentHeader *AlignmentHeader, uint32_t start, uint32_t end, uint32_t valid_start, uint32_t valid_end, uint32_t ref_length){

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
	uint32_t	end_tmp;

	MemoryCopy(&stream,(void **)&read_name	, AlignmentHeader->l_read_name	, sizeof(char));
	MemoryCopy(&stream,(void **)&cigar	, AlignmentHeader->n_cigar_op	, sizeof(uint32_t));
	MemoryCopy(&stream,(void **)&seq	, (AlignmentHeader->l_seq+1)/2	, sizeof(uint8_t));

	length = 0;
	index = AlignmentHeader->pos;

	printf("%s\t",read_name);
	//Check Valid Start Position
	if (valid_start	> start){
		for (k=start;k <valid_start;k++){
			printf(">");	
		}
		start	= valid_start;
	}

	for (k=start;k <index;k++){
		printf(">");	
	}
	end_tmp	= end;
	if (valid_end < end){	
		end	= valid_end;
	}

	for (i = 0;i < AlignmentHeader->n_cigar_op;i++){
		op	= cigar[i] & 7;
		op_len	= cigar[i] >> 4;
		if (CIGAR(op) == 'M'){
			if (index + op_len > ref_length){
				op_len = ref_length - index;
			}
			for (j = 0;j < op_len;j++){
				residue = Bin2SeqTop(seq[length >> 1],length&1);
				if (index >= start && index < end){
					printf("%c",residue);
				}
				length++;
				index++;
			}
		}else if (CIGAR(op) == 'I'){
			length += op_len;
		}else if (CIGAR(op) == 'D'){
			for (j = 0;j < op_len;j++){
				if (index >= start && index < end){
					printf("-");
				}
				index++;
			}
		}else if (CIGAR(op) == 'S'){
			length += op_len;
		}else if (CIGAR(op) == 'N'){
			for (j = 0;j < op_len;j++){
				if (index >= start && index < end){
					printf("N");
				}
				index++;
			}
		//	index += op_len;
		}
		// CIGAR(op) == 'H' don't care about it.
		
	}
	//Ed. index = 10, end=6, end_tmp=8;
	for (k=index;k < end;k++){
		printf("<");	
	}
	for (k = end; k < end_tmp;k++){
		printf("<");	
	}
//	printf("\t%d\t%d\t%d", index, end, end_tmp);
	printf("\n");	
	
	free(read_name);
	free(cigar);
	free(seq);
}



void	alignment_Pattern_Origin(uint8_t *stream, alignmentHeader *AlignmentHeader, uint32_t start, uint32_t end, uint32_t valid_start, uint32_t valid_end, uint32_t ref_length, toolsFlags *ToolsFlags){

	int	i,j,k;

	char	*read_name;
	uint32_t *cigar;
	uint8_t	*seq;
	char	*qual;

	int	op;
	int	op_len;
	char	residue;
	int	length	= 0;
	int	index;
	int	flag	= 0;
	int	index_pattern = 0;
	MemoryCopy(&stream,(void **)&read_name	, AlignmentHeader->l_read_name	, sizeof(char));
//	stream += AlignmentHeader->l_read_name;
	MemoryCopy(&stream,(void **)&cigar	, AlignmentHeader->n_cigar_op	, sizeof(uint32_t));
	MemoryCopy(&stream,(void **)&seq	, (AlignmentHeader->l_seq+1)/2	, sizeof(uint8_t));

	index = AlignmentHeader->pos;
	if (start < valid_start){
		start	= valid_start;
	}

	if (end > valid_end){
		end	= valid_end;
	}

//	printf("%s\t",read_name);
	for (i = 0;i < AlignmentHeader->n_cigar_op;i++){
		op	= cigar[i] & 7;
		op_len	= cigar[i] >> 4;

		if (op == 0){
//		if (CIGAR(op) == 'M'){
			if (index + op_len > ref_length){
				op_len = ref_length - index;
			}
			if (index + op_len < start){
				length += op_len;
				index += op_len;
			}else{
				if (start > index){
					length	+= (start-index);
					j	= (start-index);
					index	= start;
				}else {
					j = 0;
				}
				for (;j < op_len;j++){
					residue = Bin2SeqTop(seq[length >> 1],length&1);
					if (index < end){
//					if (flag == 0){
//						printf("%s\t",read_name);	
//					}
						if (ToolsFlags->flag_pattern){
							if ( (index_pattern > strlen(ToolsFlags->pattern)) || (residue != ToolsFlags->pattern[index_pattern])){
								ToolsFlags->num_ref += 1;
								return ;
							}else {
								index_pattern++;
							}
						}else {
							printf("%c",residue);
						}
						flag = 1;
					}
					length++;
					index++;
				}
			}
		}else if (op == 1){
//		}else if (CIGAR(op) == 'I'){
			if (index > start && index <= end){
				for (j = 0;j < op_len;j++){
					residue = Bin2SeqTop(seq[length >> 1],length&1);
					if (ToolsFlags->flag_pattern){
						if ( (index_pattern > strlen(ToolsFlags->pattern)) || (residue != ToolsFlags->pattern[index_pattern])){
							ToolsFlags->num_ref += 1;
							return ;
						}else {
							index_pattern++;
						}
					}else {
						printf("%c",residue);
					}
					length++;
				}
			}else {
				length+=op_len;
			}
		}else if (op == 2){
//		}else if (CIGAR(op) == 'D'){
			for (j = 0;j < op_len;j++){
				if (index >= start && index < end){
//					if (flag == 0){
//						printf("%s\t",read_name);	
//					}
					flag = 1;
				}
				index++;
			}
		}else if (op == 4){
//		}else if (CIGAR(op) == 'S'){
			length += op_len;
		}else if (op == 3){
//		}else if (CIGAR(op) == 'N'){
			for (j = 0;j < op_len;j++){
				if (index >= start && index < end){
				//	if (flag == 0){
				//		printf("%s\t",read_name);	
				//	}
					if (ToolsFlags->flag_pattern){
						if ( (index_pattern > strlen(ToolsFlags->pattern)) || (residue != ToolsFlags->pattern[index_pattern])){
							ToolsFlags->num_ref += 1;
							return ;
						}else {
							index_pattern++;
						}
					}else {
						printf("N");
					}
					flag = 1;
				}
				index++;
			}
		}
		// CIGAR(op) == 'H' don't care about it.
		
		// because insertion so > not >=
		if (index > end){
			break;
		}
		
	}
	if (flag == 1){
		if (ToolsFlags->flag_pattern){
			ToolsFlags->num_ref += 1;
			if (index_pattern == strlen(ToolsFlags->pattern)){
				ToolsFlags->num_alt += 1;
			}
		}else {
			printf("\n");	
		}
	}
	
//	free(read_name);
	free(cigar);
	free(seq);
}
