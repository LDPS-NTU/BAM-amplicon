#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include "zlib.h"

#include "BamCommonLibrary.h"



void	PrintPosQuality(FILE *file_qual_o, char *chr_name, int position, posQuality PosQuality);

int	BamQual(FILE *file_bam_i, FILE *file_length_o, int flag_hide , char *dirname_region, int flag_bed){

	FILE	*file_count_o;
	FILE	*file_region_i;
	FILE	*file_qual_o;

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
	posQuality	*PosQuality;

	char	magic[4];
	char	**chr_name;
	int	*chr_length;
	int32_t	l_text;
	char	*text;
	int32_t	n_ref;
	uint8_t *address;
	char	filename[name_len];



	uint8_t	*in_region;
//	int	have_region = 0;
	int	start,end;
	int	flag_inTarget = (strlen(dirname_region) == 0);

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
	/* comment 2020/03/30 YC
	if (flag_hide != 1){
		for (i = 0;i < n_ref;i++){
			printf("%s\t%d\n",chr_name[i], chr_length[i]);
		}
	}
	*/
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
					if (flag_hide == 1 ){
						printf("[Bam File Unzip %d / %d ] %s done\n",ref_ID+1,n_ref,chr_name[ref_ID]);
					}
					if(strlen(dirname_region) != 0){
						if (flag_inTarget == 1){
							flag_inTarget = 0;
							fclose(file_region_i);
							for(i = 0;i < chr_length[ref_ID];i++){
								if (in_region[i]==1 ){
									PrintPosQuality(file_qual_o, chr_name[ref_ID], i, PosQuality[i]);
								}
							}
						}
					}else {
						for(i = 0;i < chr_length[ref_ID];i++){
							PrintPosQuality(file_qual_o, chr_name[ref_ID], i, PosQuality[i]);
						}
					}
					free(PosQuality);
					free(in_region);
				}	
				if ( AlignmentHeader.refID != -1){
					PosQuality = calloc(chr_length[AlignmentHeader.refID],sizeof(posQuality));
					in_region = calloc(chr_length[AlignmentHeader.refID],sizeof(uint8_t));
					if(strlen(dirname_region) != 0){
						strcat_triple(filename, dirname_region, chr_name[AlignmentHeader.refID], "_region.txt", name_len);
						if ((file_region_i = fopen(filename,"r")) != NULL){
							strcat_triple(filename, "", chr_name[AlignmentHeader.refID], "_qual.txt", name_len);
							file_qual_o = fopen(filename,"w");
							flag_inTarget = 1;
							while (fscanf(file_region_i,"%d %d",&start,&end)!=EOF){
								for (i = start;i < end;i++){
									if (in_region[i] == 0){
										in_region[i] = 1;
									}
								}
							}
						}
					}
				}
				ref_ID = AlignmentHeader.refID;	
			}
			if (ref_ID == -1){
				break;	
			}
			address += sizeof(alignmentHeader);	
			if (flag_inTarget == 1){
				if ((AlignmentHeader.FLAG&4) == 0){
				//if (strlen(dirname_region) == 0 || have_region == 1){
					alignment_Quality(address, &AlignmentHeader, PosQuality, chr_length[ref_ID]);
				}
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

	printf("Q.SCORE\tA_FOR\tC_FOR\tG_FOR\tT_FOR\tA_REV\tC_REV\tG_REV\tT_REV\t**(Q.SCORE: quality score; FOR: forward; REV: reverse)\n");

		
	if (ref_ID != -1){
		if (flag_hide == 1){
			printf("[Bam File Unzip %d / %d ] %s done\n",ref_ID+1,n_ref,chr_name[ref_ID]);
		}
		if(strlen(dirname_region) != 0){
			if (flag_inTarget == 1){
				for(i = 0;i < chr_length[ref_ID];i++){
					if (in_region[i] == 1){
						PrintPosQuality(file_qual_o, chr_name[ref_ID], i, PosQuality[i]);
					}	
				}
			}
		}else {
			for(i = 0;i < chr_length[ref_ID];i++){
				PrintPosQuality(file_qual_o, chr_name[ref_ID], i, PosQuality[i]);
			}
		}
		free(PosQuality);
		free(in_region);
	}
	

	free(stream_i);	
	free(stream_o);	
	free(buffer);	
	return 0;
}

void	PrintPosQuality(FILE *file_qual_o, char *chr_name, int position, posQuality PosQuality){
	if (PosQuality.cov_sum != 0){
		fprintf(file_qual_o,"%s\t%d\t%d\t%lf\t%lf\t%lf\t%lf\n",
			chr_name,
			position,
			PosQuality.cov_sum,
			(double)PosQuality.base_qual_sum /PosQuality.cov_sum,
			sqrt((double)(PosQuality.base_qual_pow-pow(PosQuality.base_qual_sum,2)/PosQuality.cov_sum)),
			(double)PosQuality.map_qual_sum /PosQuality.cov_sum,
			sqrt((double)(PosQuality.map_qual_pow-pow(PosQuality.map_qual_sum,2)/PosQuality.cov_sum))
		);
	}else {
		fprintf(file_qual_o,"%s\t%d\t%d\t%lf\t%lf\t%lf\t%lf\n",chr_name,position,0,0.0,0.0,0.0,0.0);
	}
}
