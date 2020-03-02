#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include "zlib.h"

#include "BamCommonLibrary.h"


int	BamPos(FILE *file_bam_i, FILE *file_length_o, int flag_hide , char *dirname_region, int flag_bin){

	FILE	*file_o;
	FILE	*file_pos_i;

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

	char	magic[4];
	char	**chr_name;
	int	*chr_length;
	int32_t	l_text;
	char	*text;
	int32_t	n_ref;
	uint8_t *address;
	char	filename[name_len];



	uint8_t	*in_region;
	int	valid_region;
	int	valid_all_region = 0;
	int	have_region = 0;
	int	position;


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
	address += (12+l_text);


	chr_name	= malloc(n_ref*sizeof(char *));	
	chr_length	= calloc(n_ref,sizeof(int));
	for (i = 0;i < n_ref;i++){
		chr_name[i]	= calloc(name_len,sizeof(char));
		address += refInformation(address,chr_name[i],&chr_length[i]);
	}
	if (file_length_o != NULL){
		for (i = 0;i < n_ref;i++){
			fprintf(file_length_o,"%s\t%d\n",chr_name[i], chr_length[i]);
		}
	}
	if (flag_hide != 1){
		for (i = 0;i < n_ref;i++){
			printf("%s\t%d\n",chr_name[i], chr_length[i]);
		}
	}
	top = buffer + BgfzTail.I_size;
	counter	= top - address;


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

			//Chromosome Change
			if (AlignmentHeader.refID != ref_ID){
				if (ref_ID != -1){
					if (flag_hide == 0 ){
						printf("[Bam File Unzip %d / %d ] %s done\n",ref_ID+1,n_ref,chr_name[ref_ID]);
					}
					if (have_region == 1){
						fclose(file_pos_i);
						have_region = 0;
					}
					if (have_region == 1){
						rewind(file_pos_i);
						while (fscanf(file_pos_i,"%*s%d",&position)!=EOF){
							i = position -1;
							printf("%s\t%d\t%d\t%d\t%d\t%d\n"
								,chr_name[ref_ID]
								,i+1
								,PosCoverage[i].A
								,PosCoverage[i].C
								,PosCoverage[i].G
								,PosCoverage[i].T);
						}							
					}
//					for(i = 0;i < chr_length[ref_ID];i++){
//						if (in_region[i] == 1){
//							printf("%s\t%d\t%d\t%d\t%d\t%d\n"
//								,chr_name[ref_ID]
//								,i+1
//								,PosCoverage[i].A
//								,PosCoverage[i].C
//								,PosCoverage[i].G
//								,PosCoverage[i].T);
//						}	
//					}
					free(in_region);
					free(PosCoverage);
				}	
				if ( AlignmentHeader.refID != -1){	
					PosCoverage = calloc(chr_length[AlignmentHeader.refID],sizeof(posCoverage));
					in_region = calloc(chr_length[AlignmentHeader.refID],sizeof(uint8_t));
					strcat_triple(filename, dirname_region, chr_name[AlignmentHeader.refID], ".txt", name_len);
					if ((file_pos_i = fopen(filename,"r")) != NULL){
						have_region = 1;
						while (fscanf(file_pos_i,"%*s%d",&position)!=EOF){
							in_region[position-1] = 1;
						}							
					}
				}
				ref_ID = AlignmentHeader.refID;	
			}
			if (ref_ID == -1){
				break;	
			}
			address += sizeof(alignmentHeader);	
			if ((AlignmentHeader.FLAG&4) == 0 && have_region == 1){
				alignment_Coverage(address, &AlignmentHeader, PosCoverage, chr_length[ref_ID]);
			}
			address += (AlignmentHeader.block_size + 4) - sizeof(alignmentHeader);

			if (top - address > sizeof(alignmentHeader)){
				memcpy(&AlignmentHeader, address, sizeof(alignmentHeader));
			}else {
				break;	
			}
		}
		counter = top - address;
		if (ref_ID == -1){
			break;	
		}
	}
	
	if (ref_ID != -1){
		if (flag_hide == 0){
			printf("[Bam File Unzip %d / %d ] %s done\n",ref_ID+1,n_ref,chr_name[ref_ID]);
		}
		if (have_region == 1){
			have_region = 0;
		}
		free(PosCoverage);
	}
	
	free(stream_i);	
	free(stream_o);	
	free(buffer);	
	return 0;
}
