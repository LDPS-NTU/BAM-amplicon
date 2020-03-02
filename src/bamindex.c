#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>


int main(int argc, char *argv[]){
	FILE	*file_i;

	int	i,j,m,n;
	char	BAMI[4];
	int32_t	n_ref;
	int32_t	n_bin;
	uint32_t	bin;
	int32_t	n_chunk;
	uint64_t	chunk_beg;
	uint64_t	chunk_end;
	int32_t	n_intv;
	uint64_t	ioffset;
	file_i=fopen(argv[1],"rb");

	fread(BAMI,4,sizeof(char),file_i);
	fread(&n_ref,1,sizeof(int32_t),file_i);
	printf("%s\n",BAMI);
	printf("REF:%d\n",n_ref);
	for(i = 0;i < n_ref;i++){
		fread(&n_bin,1,sizeof(int32_t),file_i);
		printf("N_BIN:%d\n",n_bin);
		for (j = 0;j < n_bin;j++){
			fread(&bin,1,sizeof(uint32_t),file_i);
			fread(&n_chunk,1,sizeof(int32_t),file_i);
			printf("BIN:%u\n",bin);
//			printf("\tN_Chunk:%d\n",n_chunk);
			for(m = 0;m < n_chunk;m++){
				fread(&chunk_beg,1,sizeof(uint64_t),file_i);
				fread(&chunk_end,1,sizeof(uint64_t),file_i);
				printf("\tBEG:%u\n",chunk_beg);
				printf("\tEND:%u\n",chunk_end);
			}
		}
		fread(&n_intv,1,sizeof(int32_t),file_i);
		printf("N_INTV:%d\n",n_intv);
		for (n = 0;n < n_intv;n++){
			fread(&ioffset,1,sizeof(uint64_t),file_i);
//			printf("OFF:%u\n",ioffset);
		}
		getchar();
	}
	return 0;

}
