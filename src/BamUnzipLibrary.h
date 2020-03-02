#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include "zlib.h"


typedef	struct bgfzHeader_st {

	uint8_t	NUM_31;
	uint8_t	NUM_139;
	uint8_t	NUM_8;
	uint8_t	NUM_4;
	uint32_t MTIME;
	uint8_t	XFL;
	uint8_t	OS;
	uint16_t XLEN;
	uint8_t	NUM_66;
	uint8_t	NUM_67;
	uint16_t NUM_2;
	uint16_t BSIZE;

} __attribute__((packed)) bgfzHeader;


typedef struct	alignmentHeader_st {

	int32_t	block_size;
	int32_t	refID;
	int32_t	pos;
	//bin_mq_nl
	uint8_t l_read_name;
	uint8_t	MAPQ;
	uint16_t bin;
	//flag_nc
	uint16_t n_cigar_op;
	uint16_t FLAG;
	int32_t	l_seq;
	int32_t next_refID;
	int32_t next_pos;
	int32_t	tlen;
} alignmentHeader;

char	CIGAR(int A);
int tag(void *stream);

int BinToEncodedSeq_Top(uint8_t A,int MSB);
int BinToEncodedSeq(uint8_t A);

void PrintAlignmenetHeader(alignmentHeader *AlignmentHeader);
void PrintSequence(uint8_t *seq,int32_t seq_length);
void PrintCIGAR(uint32_t *cigar,uint16_t n_cigar_op);
void PrintPosInformation(posInformation *PosInformation,int position);
int EndPosition(int32_t pos, uint16_t n_cigar_op, uint32_t *cigar);
void MemoryCopy(char **stream,void **ptr,int length,int size);
void PrintBlank(int size);


int alignment(FILE *file,int query_position,int *begin);


int	alignment(FILE *file ,int query_position,int *begin){

	alignmentHeader	AlignmentHeader;
	char	*read_name;
	uint32_t *cigar;
	uint8_t	*seq;
	char	*qual;
	int	i,j;
	int	op;
	int	op_len;
	int	length;
	char	*stream;
	char	*temp_stream;
	int temp = 0;
	int	result;
	int	gap;
	int	number;
	int	end_position;
	int	index;
	stream	= calloc(500,sizeof(char));	
	temp_stream = stream;
	if ( fread(&AlignmentHeader,sizeof(alignmentHeader),1,file) != 1){
		return 0;
	};
	fread(stream,sizeof(char),AlignmentHeader.block_size-32,file);

	if (query_position - 60 < AlignmentHeader.pos && query_position > AlignmentHeader.pos){
		MemoryCopy(&stream,(void **)&read_name   ,AlignmentHeader.l_read_name    ,sizeof(char));
		MemoryCopy(&stream,(void **)&cigar       ,AlignmentHeader.n_cigar_op     ,sizeof(uint32_t));
		end_position = EndPosition(AlignmentHeader.pos,AlignmentHeader.n_cigar_op,cigar)+1;

		if (query_position <= end_position){
		if (*begin == 0){
			PrintBlank(22-1);
			*begin = AlignmentHeader.pos;
			number = AlignmentHeader.pos % 10;
			for(i = 0;i < 100;i++){
				printf("%d",(number+i)%10);
			}	
			printf("\r\n");
		}
		gap = AlignmentHeader.pos - *begin;

		MemoryCopy(&stream,(void **)&seq         ,(AlignmentHeader.l_seq+1)/2    ,sizeof(uint8_t));
		PrintAlignmenetHeader(&AlignmentHeader);
		PrintBlank(gap);

		//for (i = 0;i < AlignmentHeader.l_seq;i++){
		//	BinToEncodedSeq_Top(seq[i >> 1],i&1);
		//}

		length = 0;
//		printf("\t");

		for (i = 0;i < AlignmentHeader.n_cigar_op;i++){
			op	= cigar[i] & 7;
			op_len	= cigar[i] >> 4;
//			printf("%d%c",op_len,CIGAR(op));
			index = 0;
			if (CIGAR(op) == 'M'){
				for (j = 0;j < op_len;j++){
					BinToEncodedSeq_Top(seq[index >> 1],index&1);
					index++;
				}
			}else if (CIGAR(op) == 'I'){
				index = index+op_len;
			}else if (CIGAR(op) == 'D'){
				for (j = 0;j < op_len;j++){
					printf("*");
				}
			}else if (CIGAR(op) == 'S'){
			}else if (CIGAR(op) == 'N'){
				index = index+op_len;	
			}
		}
		printf("\r\n");
		free(read_name);
		free(cigar);
		free(seq);
	//	free(qual);
		}
	}else if (query_position + 60 <= AlignmentHeader.pos){
		return 0;
	}
	stream = temp_stream;
	free(stream);

	return	(AlignmentHeader.block_size + 4);
}

int BinToEncodedSeq(uint8_t A){
	
	switch(A){
	case 0:	printf("=");
		break;
	case 1:	printf("A");
		break;
	case 2:	printf("C");
		break;
	case 3:	printf("M");
		break;
	case 4:	printf("G");
		break;
	case 5:	printf("R");
		break;
	case 6:	printf("S");
		break;
	case 7:	printf("V");
		break;
	case 8:	printf("T");
		break;
	case 9:	printf("W");
		break;
	case 10:printf("Y");
		break;
	case 11:printf("H");
		break;
	case 12:printf("K");
		break;
	case 13:printf("D");
		break;
	case 14:printf("B");
		break;
	case 15:printf("N");
		break;
	default:printf("Error\n");
		exit(1);
	}
}

int BinToEncodedSeq_Top(uint8_t A,int MSB){
	uint8_t	low 	= A & 15;
	uint8_t	high	= (A & 240)>>4;
	if (MSB == 0){
		BinToEncodedSeq(high);
	} else {
		BinToEncodedSeq(low);
	}
}

void PrintSequence(uint8_t *seq,int32_t seq_length){
	int	i;
	for (i = 0;i < seq_length;i++){
		BinToEncodedSeq_Top(seq[i >> 1],i&1);
	}
}
void PrintCIGAR(uint32_t *cigar,uint16_t n_cigar_op){
	int	i;
	for (i = 0;i < n_cigar_op;i++){
		printf("%d%c",cigar[i] >> 4,CIGAR(cigar[i] & 7));
	}
}

int	tag(void *stream){
	char	tag[2];
	char	val_type;
	uint32_t value;

	memcpy(&tag,stream,2);
	memcpy(&val_type,stream+2,1);
	stream = stream + 3;
	
	switch(val_type){
	case 'A': printf("Unknown A\n");
		  exit(1);
		break;
	case 'c': memcpy(&value,stream,1);
		break;
	case 'C':  memcpy(&value,stream,1);
		break;
	case 's':  memcpy(&value,stream,2);
		break;
	case 'S':  memcpy(&value,stream,2);
		break;
	case 'i':  memcpy(&value,stream,4);
		break;
	case 'I':  memcpy(&value,stream,4);
		break;
	case 'f': printf("Unknown f\n");
		  exit(1);
		break;
	case 'Z': printf("Unknown Z\n");
		  exit(1);
		break;
	case 'H': printf("Unknown H\n");
		  exit(1);
		break;
	case 'B': printf("Unknown B\n");
		  exit(1);
		break;
	default: printf("Tag Error\n");
		exit(1);
	}
	
	printf("%c%c %c %d\n",tag[0],tag[1],val_type,value);
		
	return	(3+sizeof(value));
	
}

char	CIGAR(int A){
	switch(A){
	case 0:	return 'M';
	case 1:	return 'I';
	case 2:	return 'D';
	case 3:	return 'N';
	case 4:	return 'S';
	case 5:	return 'H';
	case 6:	return 'P';
	case 7:	return '=';
	case 8:	return 'X';
	default : printf("Error CIGAR\n");
		exit(1);
	}
}

void PrintAlignmenetHeader(alignmentHeader *AlignmentHeader){
//	printf("%d ",AlignmentHeader->block_size);
	printf("%2d ",AlignmentHeader->refID+1);
	printf("%10d ",AlignmentHeader->pos+1);
//	printf("%u ",AlignmentHeader->l_read_name);
//	printf("%u ",AlignmentHeader->MAPQ);
//	printf("%u ",AlignmentHeader->bin);
//	printf("%u ",AlignmentHeader->n_cigar_op);
	printf("%4u ",AlignmentHeader->FLAG);
	printf("%2d ",AlignmentHeader->l_seq);
//	printf("%d ",AlignmentHeader->next_refID+1);
//	printf("%d ",AlignmentHeader->nex_pos+1);
//	printf("%d ",AlignmentHeader->tlen);
}

void PrintPosInformation(posInformation *PosInformation,int position){
	printf("position:%10d ",position);
	printf("genotype:%c ",PosInformation->genotype);	
	printf("het:%c ",PosInformation->het);	
	printf("A:%5d ",PosInformation->A);	
	printf("T:%5d ",PosInformation->T);	
	printf("G:%5d ",PosInformation->G);	
	printf("C:%5d\n",PosInformation->C);	
}
void PrintBlank(int size){
	for (;size > 0;size--){
		printf(" ");
	}
}
int EndPosition(int32_t pos, uint16_t n_cigar_op, uint32_t *cigar){
	
	int	i;
	int	op;
	int	op_len;
	int	flag = 1;
	int	end;

	end	= pos;
	for (i = 0;i < n_cigar_op;i++){
		op	= cigar[i] & 7;
		op_len	= cigar[i] >> 4;
		if (CIGAR(op) == 'M'){
			end = end + op_len;	
		}else if (CIGAR(op) == 'I'){
		}else if (CIGAR(op) == 'D'){
			end = end + op_len;	
		}else if (CIGAR(op) == 'S'){
		}else if (CIGAR(op) == 'N'){
			end = end + op_len;	
		}
	}
	return (end-1);
}


void MemoryCopy(char **stream,void **ptr,int length,int size){
	*ptr	= calloc(length,size);
	memcpy(*ptr,*stream,length*size);
	*stream	+= length*size;
}

