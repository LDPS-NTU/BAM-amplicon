#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include "zlib.h"

#include "BamCommonLibrary.h"




typedef struct bamStatRead_st{
	uint64_t	Total;
	uint64_t	Mapped;
	uint64_t	Mapped_Chimeric;
	uint64_t	UnMapped;
	uint64_t	UnMapped_ID;

	//Duplicated
	uint64_t	Dup_R1;
	uint64_t	Dup_R2;
	
	//------R1 R2
	uint64_t	R1;
	uint64_t	R2;
	uint64_t	SingleEnd;
	uint64_t	R1_Mapped;
	uint64_t	R2_Mapped;
	uint64_t	SingleEnd_Mapped;
	uint64_t	R1_Mapped_Properly;
	uint64_t	R2_Mapped_Properly;
	uint64_t	R1_Mapped_Chimeric;
	uint64_t	R2_Mapped_Chimeric;
	uint64_t	R1_UnMapped_ID;
	uint64_t	R2_UnMapped_ID;
	uint64_t	SingleEnd_UnMapped_ID;
	uint64_t	R1_UnMapped_ID_Properly;
	uint64_t	R2_UnMapped_ID_Properly;
	uint64_t	R1_UnMapped_NoID;
	uint64_t	R2_UnMapped_NoID;
	uint64_t	SingleEnd_UnMapped_NoID;
	
	//------Same Cross Singleton
	uint64_t	Same;
	uint64_t	Cross;
	uint64_t	Singleton;
	uint64_t	UnMapped_Properly;
	uint64_t	Same_Properly;
	uint64_t	Cross_Properly;
	uint64_t	Singleton_Properly;
	uint64_t	Mapped_Properly;
	uint64_t	Error;
}  __attribute__((packed)) bamStatRead;

void	alignment_Stat(uint8_t *stream, alignmentHeader *AlignmentHeader, char **chr_name, bamStatRead *BamStatRead, FILE *file_o);

int	BamStat(FILE *file_bam_i, toolsFlags *ToolsFlags){

	uint8_t *stream_i;//[65536];
	uint8_t *stream_o;//[65536];
	uint8_t	*buffer;
	uint8_t	*top;



	FILE	*file_tlength_o;
	FILE	*file_confuse_o;
	int	i;
	int	len_data;
	int	counter;
	int	ref_ID	= 0;

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

	uint8_t	*in_region;
	int	valid_region;
	int	valid_all_region = 0;
	int	have_region = 0;
	int	start,end;
	int	pair = 0;
	int	cross = 0;
	int	only = 0;
	int	unmap = 0;
	int	template_length;
	int32_t	data_length;
	int32_t	text_length;
	int32_t	l_name;
	bamStatRead	BamStatRead;



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
	memcpy(buffer, stream_o, BgfzTail.I_size);

	address = buffer;	
	memcpy(magic, address,sizeof(char)*4);
	memcpy(&l_text, address+4,sizeof(int32_t));
	text = calloc(l_text+1, sizeof(char));
	address += 8;
	data_length = BgfzTail.I_size - 8;


	while(l_text - data_length > 0){
		len_data = BGFZBlock(file_bam_i);
		decompressBlock(&infstream, stream_i, stream_o, len_data, file_bam_i);
		fread(&BgfzTail	,sizeof(bgfzTail),1,file_bam_i);
		memcpy(buffer,stream_o, BgfzTail.I_size);
		address = buffer;	
		
		l_text = l_text - data_length;
		data_length = BgfzTail.I_size;
	}

	memcpy(&n_ref, address+l_text, sizeof(int32_t));
	address += (4+l_text);

	chr_name	= malloc(n_ref*sizeof(char *));	
	chr_length	= calloc(n_ref,sizeof(int));

	top = buffer + BgfzTail.I_size;
	counter	= top - address;

	//printf("%u\t%u\n", l_text, n_ref);
	for (i = 0;i < n_ref;i++){
		if (top - address < 4){
			counter = top - address;
			len_data = BGFZBlock(file_bam_i);
			decompressBlock(&infstream, stream_i, stream_o, len_data, file_bam_i);
			fread(&BgfzTail	,sizeof(bgfzTail),1,file_bam_i);
			memcpy(buffer, address,sizeof(uint8_t)*counter);
			memcpy(&buffer[counter], stream_o,BgfzTail.I_size);
			top	= buffer + counter + BgfzTail.I_size;
			address = buffer;
		}

		memcpy(&l_name, address, sizeof(int32_t));
		address += 4;
		chr_name[i]	= calloc(l_name, sizeof(char));

		if (top - address < l_name + 4){
			counter = top - address;
			len_data = BGFZBlock(file_bam_i);
			decompressBlock(&infstream, stream_i, stream_o, len_data, file_bam_i);
			fread(&BgfzTail	,sizeof(bgfzTail),1,file_bam_i);
			memcpy(buffer,address,sizeof(uint8_t)*counter);
			memcpy(buffer+counter, stream_o, BgfzTail.I_size);
			top	= buffer + counter + BgfzTail.I_size;
			address = buffer;
		}

		memcpy(chr_name[i], address, sizeof(char)*l_name);
		memcpy(&chr_length[i], address+l_name, sizeof(int32_t));
		address += (l_name+4);
		//printf("%s\t%d\n",chr_name[i], chr_length[i]);
	}

	/* comment 2020/03/30 YC
	if (ToolsFlags->flag_hide == 1){
		for (i = 0;i < n_ref;i++){
			printf("%s\t%d\n",chr_name[i], chr_length[i]);
		}
	}
	*/
	counter	= top - address;
	
	memset((char *)&BamStatRead,0,sizeof(bamStatRead));
	file_tlength_o = fopen("Template_length.txt","w");
	if (file_tlength_o == NULL){
		printf("[Error] Can't Build the Template_length.txt\n Please Check Permission\n");
		return 0;
	}
	file_confuse_o = fopen("Confuse_Reads.txt","w");
	if (file_confuse_o == NULL){
		printf("[Error] Can't Build the Confuse_Reads.txt\n Please Check Permission\n");
		return 0;
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
	
//		printf("1:%s\t%d\t%u\t%u\n",chr_name[ref_ID],AlignmentHeader.pos, top-address,sizeof(AlignmentHeader.block_size));
		while ( top - address >= (AlignmentHeader.block_size + 4) ){
//			printf("%d\n",AlignmentHeader.refID);
			//Chromosome Change
			if (AlignmentHeader.refID != ref_ID){
				if (ToolsFlags->flag_hide == 1 && n_ref != 0){
					printf("[Bam File Unzip %d / %d ] %s done\n",ref_ID+1,n_ref,chr_name[ref_ID]);
				}
				ref_ID = AlignmentHeader.refID;
			}
			address += sizeof(alignmentHeader);
			alignment_Stat(address, &AlignmentHeader, chr_name, &BamStatRead, file_confuse_o);
			if (AlignmentHeader.tlen > 0){
				fprintf(file_tlength_o,"%d\n",AlignmentHeader.tlen);
			}
			address += (AlignmentHeader.block_size + 4) - sizeof(alignmentHeader);

			if (top - address > sizeof(alignmentHeader)){
				memcpy(&AlignmentHeader, address, sizeof(alignmentHeader));
			}else {
				break;	
			}
//			printf("%s\t%d\t%u\t%u\n",chr_name[ref_ID],AlignmentHeader.pos, top-address, AlignmentHeader.block_size);
		}
//		printf("2:%s\t%d\t%u\t%u\n",chr_name[ref_ID],AlignmentHeader.pos, top-address,sizeof(alignmentHeader));
		counter = top - address;

	}
	
	BamStatRead.Mapped_Chimeric = BamStatRead.R1_Mapped_Chimeric + BamStatRead.R2_Mapped_Chimeric;
	
	BamStatRead.Mapped	= BamStatRead.R1_Mapped + BamStatRead.R2_Mapped + BamStatRead.SingleEnd_Mapped;
	BamStatRead.UnMapped_ID	= BamStatRead.R1_UnMapped_ID + BamStatRead.R2_UnMapped_ID + BamStatRead.SingleEnd_UnMapped_ID;
	BamStatRead.UnMapped	= BamStatRead.R1_UnMapped_NoID + BamStatRead.R2_UnMapped_NoID + BamStatRead.SingleEnd_UnMapped_NoID + BamStatRead.UnMapped_ID;
	BamStatRead.Total	= BamStatRead.Mapped+BamStatRead.UnMapped;

	printf("Item\tMapped\tUnMap\tTotal\t(Mapped %)\n");
	printf("===== ===== ===== ===== ===== ===== ===== =====\n");
	printf("Sample\t%lu\t%lu\t%lu\t(%f)\n",BamStatRead.Mapped, BamStatRead.UnMapped , BamStatRead.Total, (float)BamStatRead.Mapped /BamStatRead.Total);
	printf("R1\t%lu\t%lu\t%lu\t(%f)\n",BamStatRead.R1_Mapped, BamStatRead.R1_UnMapped_ID+BamStatRead.R1_UnMapped_NoID, BamStatRead.R1,
		(float)BamStatRead.R1_Mapped / BamStatRead.R1);
	printf("R2\t%lu\t%lu\t%lu\t(%f)\n",BamStatRead.R2_Mapped, BamStatRead.R2_UnMapped_ID+BamStatRead.R2_UnMapped_NoID, BamStatRead.R2, 
		(float)BamStatRead.R2_Mapped / BamStatRead.R2);
	printf("----- ----- ----- ----- ----- ----- ----- -----\n");
	
	printf("UnChimeric_Sample\t%lu\t%lu\t%lu\t(%f)\n"
		,BamStatRead.Mapped-BamStatRead.Mapped_Chimeric
		,BamStatRead.UnMapped
		,BamStatRead.Mapped-BamStatRead.Mapped_Chimeric+BamStatRead.UnMapped,
		(float)(BamStatRead.Mapped-BamStatRead.Mapped_Chimeric)/ (BamStatRead.Mapped-BamStatRead.Mapped_Chimeric+BamStatRead.UnMapped));
	
	printf("UnChimeric_R1\t%lu\t%lu\t%lu\t(%f)\n"
		,BamStatRead.R1_Mapped-BamStatRead.R1_Mapped_Chimeric
		,BamStatRead.R1_UnMapped_ID+BamStatRead.R1_UnMapped_NoID
		,BamStatRead.R1-BamStatRead.R1_Mapped_Chimeric,
		(float)(BamStatRead.R1_Mapped-BamStatRead.R1_Mapped_Chimeric)/(BamStatRead.R1-BamStatRead.R1_Mapped_Chimeric));
	printf("UnChimeric_R2\t%lu\t%lu\t%lu\t(%f)\n",
		BamStatRead.R2_Mapped-BamStatRead.R2_Mapped_Chimeric,
		BamStatRead.R2_UnMapped_ID+BamStatRead.R2_UnMapped_NoID,
		BamStatRead.R2-BamStatRead.R2_Mapped_Chimeric, 
		(float)(BamStatRead.R2_Mapped-BamStatRead.R2_Mapped_Chimeric)/(BamStatRead.R2-BamStatRead.R2_Mapped_Chimeric));
	printf("----- ----- ----- ----- ----- ----- ----- -----\n");
	
	printf("Proper\t%lu\t%lu\t%lu\n",
		BamStatRead.Same_Properly + BamStatRead.Cross_Properly + BamStatRead.Singleton_Properly,
		BamStatRead.R1_UnMapped_ID_Properly + BamStatRead.R2_UnMapped_ID_Properly,
		BamStatRead.Same_Properly + BamStatRead.Cross_Properly + BamStatRead.Singleton_Properly + BamStatRead.R1_UnMapped_ID_Properly + BamStatRead.R2_UnMapped_ID_Properly);
	printf("R1_PP\t%lu\t%lu\t%lu\t(%f)\n",BamStatRead.R1_Mapped_Properly, BamStatRead.R1_UnMapped_ID_Properly, 
		BamStatRead.R1_Mapped_Properly+BamStatRead.R1_UnMapped_ID_Properly,
		(float)BamStatRead.R1_Mapped_Properly / (BamStatRead.R1_Mapped_Properly+BamStatRead.R1_UnMapped_ID_Properly));
	printf("R2_PP\t%lu\t%lu\t%lu\t(%f)\n",BamStatRead.R2_Mapped_Properly, BamStatRead.R2_UnMapped_ID_Properly, 
		BamStatRead.R2_Mapped_Properly+BamStatRead.R2_UnMapped_ID_Properly,
		(float)BamStatRead.R2_Mapped_Properly / (BamStatRead.R2_Mapped_Properly+BamStatRead.R2_UnMapped_ID_Properly));
	printf("Dupli ----- ----- ----- ----- ----- ----- -----\n");
	printf("R1_Dup\t%lu\n",BamStatRead.Dup_R1);
	printf("R2_Dup\t%lu\n",BamStatRead.Dup_R2);
	printf("===== ===== ===== ===== ===== ===== ===== =====\n");
	printf("Item\tProper\tTotal\t(Properly %)\n");
	printf("===== ===== ===== ===== ===== ===== ===== =====\n");
	if (BamStatRead.Same == 0){
		printf("SameChr\t%lu\t%lu\t(%f)\n", 0, 0, 0.0);
	}else{
		printf("SameChr\t%lu\t%lu\t(%f)\n",BamStatRead.Same_Properly, BamStatRead.Same, (float)BamStatRead.Same_Properly / BamStatRead.Same);
	}
	if (BamStatRead.Cross == 0){
		printf("Cross\t%lu\t%lu\t(%f)\n", 0, 0, 0.0);
	}else{
		printf("Cross\t%lu\t%lu\t(%f)\n",BamStatRead.Cross_Properly, BamStatRead.Cross, (float)BamStatRead.Cross_Properly / BamStatRead.Cross);
	}
	if (BamStatRead.Singleton == 0){
		printf("Single\t%lu\t%lu\t(%f)\n", 0, 0, 0.0);
	}else{
		printf("Single\t%lu\t%lu\t(%f)\n",BamStatRead.Singleton_Properly, BamStatRead.Singleton, (float)BamStatRead.Singleton_Properly / BamStatRead.Singleton);
	}
	printf("===== ===== ===== ===== ===== ===== ===== =====\n");
//	printf("Error\t%lu\n",BamStatRead.Error);
	
	free(stream_i);	
	free(stream_o);	
	free(buffer);	
	return 0;
}
void	alignment_Stat(uint8_t *stream, alignmentHeader *AlignmentHeader, char **chr_name, bamStatRead *BamStatRead, FILE *file_o){

	char	*read_name;
	uint32_t *cigar;
	int	i;
	int	op;
	int	op_len;
	
	if (AlignmentHeader->refID >= 0){
		MemoryCopy(&stream,(void **)&read_name	, AlignmentHeader->l_read_name	, sizeof(char));
//		MemoryCopy(&stream,(void **)&cigar	, AlignmentHeader->n_cigar_op	, sizeof(uint32_t));
//		printf("%d\t%d\t%d\t",AlignmentHeader->refID, AlignmentHeader->l_seq, AlignmentHeader->FLAG);
//		for (i = 0;i < AlignmentHeader->n_cigar_op;i++){
//			op	= cigar[i] & 7;
//			op_len	= cigar[i] >> 4;
//			printf("%d%c", op_len,CIGAR(op));
//		}
//		printf("\t%s\n", read_name);
		if ((AlignmentHeader->FLAG&4) == 0){
			if ((AlignmentHeader->FLAG&64) == 64){
				BamStatRead->R1++;
				BamStatRead->R1_Mapped++;
				if ((AlignmentHeader->FLAG&2) == 2){	BamStatRead->R1_Mapped_Properly++;}
				if ((AlignmentHeader->FLAG&256) == 256){	BamStatRead->R1_Mapped_Chimeric++;}
			}else if ((AlignmentHeader->FLAG&128) == 128){
				BamStatRead->R2++;	
				BamStatRead->R2_Mapped++;	
				if ((AlignmentHeader->FLAG&2) == 2){	BamStatRead->R2_Mapped_Properly++;}
				if ((AlignmentHeader->FLAG&256) == 256){	BamStatRead->R2_Mapped_Chimeric++;}
			}else {
				BamStatRead->SingleEnd++;	
				BamStatRead->SingleEnd_Mapped++;	
			}
		}else {
			BamStatRead->UnMapped++;
			BamStatRead->UnMapped_ID++;
			
			if ((AlignmentHeader->FLAG&64) == 64){
				BamStatRead->R1++;
				BamStatRead->R1_UnMapped_ID++;
			}else if ((AlignmentHeader->FLAG&128) == 128){
				BamStatRead->R2++;	
				BamStatRead->R2_UnMapped_ID++;
			}else {
				BamStatRead->SingleEnd++;	
				BamStatRead->SingleEnd_UnMapped_ID++;
			}	
		}
		if ((AlignmentHeader->FLAG&1024)!=0){
			if ((AlignmentHeader->FLAG&64) == 64){
				BamStatRead->Dup_R1++;
			}else if ((AlignmentHeader->FLAG&128) == 128){
				BamStatRead->Dup_R2++;
			}
			
		}else {	
			//Both are mapped and they are on Same Chromosome
			if (AlignmentHeader->next_refID == AlignmentHeader->refID){
				if ((AlignmentHeader->FLAG&4) == 0 && (AlignmentHeader->FLAG&8) == 0){
					BamStatRead->Same++;
					if ((AlignmentHeader->FLAG&2) == 2){	BamStatRead->Same_Properly++;}
				}else if ((AlignmentHeader->FLAG&4) == 0 &&(AlignmentHeader->FLAG&8) == 8){
					BamStatRead->Singleton++;
					if ((AlignmentHeader->FLAG&2) == 2){	BamStatRead->Singleton_Properly++;}
				}else if ((AlignmentHeader->FLAG&4) == 4){
					if ((AlignmentHeader->FLAG&2) == 2){	
			//			printf("%d\t%s\n",AlignmentHeader->FLAG,read_name);
						fprintf(file_o,"%d\t%s\n",AlignmentHeader->FLAG,read_name);
						BamStatRead->UnMapped_Properly++;
						if ((AlignmentHeader->FLAG&64) == 64){
							BamStatRead->R1_UnMapped_ID_Properly++;
						}else if ((AlignmentHeader->FLAG&128) == 128){
							BamStatRead->R2_UnMapped_ID_Properly++;
						}
					}
				}else{
					BamStatRead->Error++;	
				}
			//Both are mapped ,but they are on Different Chromosome
			}else if (AlignmentHeader->next_refID != AlignmentHeader->refID && AlignmentHeader->next_refID >= 0){
				if ((AlignmentHeader->FLAG&4) == 0 && (AlignmentHeader->FLAG&8) == 0){
					BamStatRead->Cross++;
					if ((AlignmentHeader->FLAG&2) == 2){	
						BamStatRead->Cross_Properly++;	
			//			printf("%d\t%s\t(Cross_PP)\n",AlignmentHeader->FLAG,read_name);	
						fprintf(file_o,"%d\t%s\t(Cross_PP)\n",AlignmentHeader->FLAG,read_name);	
					}
				}else if ((AlignmentHeader->FLAG&4) == 0 &&(AlignmentHeader->FLAG&8) == 8){
					BamStatRead->Singleton++;
					if ((AlignmentHeader->FLAG&2) == 2){	BamStatRead->Singleton_Properly++;}
				}else if ((AlignmentHeader->FLAG&4) == 4){
					if ((AlignmentHeader->FLAG&2) == 2){	
			//			printf("%d\t%s\n",AlignmentHeader->FLAG,read_name);
						fprintf(file_o,"%d\t%s\n",AlignmentHeader->FLAG,read_name);
						BamStatRead->UnMapped_Properly++;
						if ((AlignmentHeader->FLAG&64) == 64){
							BamStatRead->R1_UnMapped_ID_Properly++;
						}else if ((AlignmentHeader->FLAG&128) == 128){
							BamStatRead->R2_UnMapped_ID_Properly++;
						}					
					}
				}else{
					BamStatRead->Error++;	
				}	
			}else {
				if ((AlignmentHeader->FLAG&4) == 0 &&(AlignmentHeader->FLAG&8) == 8){
					BamStatRead->Singleton++;
					if ((AlignmentHeader->FLAG&2) == 2){	BamStatRead->Singleton_Properly++;}
				}else {
//					printf("%d\t%s\n", AlignmentHeader->FLAG, read_name);
//					BamStatRead->Error++;	
					//Singleton
					BamStatRead->Singleton++;
				}
			}
		}
		free(read_name);
	}else {
		if ((AlignmentHeader->FLAG&64) == 64){
			BamStatRead->R1++;
			BamStatRead->R1_UnMapped_NoID++;
		}else if ((AlignmentHeader->FLAG&128) == 128){
			BamStatRead->R2++;
			BamStatRead->R2_UnMapped_NoID++;
		}else {
			BamStatRead->SingleEnd++;
			BamStatRead->SingleEnd_UnMapped_NoID++;
		}	
	}
}
