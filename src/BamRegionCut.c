#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include "zlib.h"

#include "BamCommonLibrary.h"




void	PrintBaseDist	(FILE *file_base, posCoverage *PosCoverage, uint8_t *in_region, char *chr_name, int chr_length, covBaseDistribution *BaseDist, int operation);

int	BamRegionCut(FILE *file_bam_i, int flag_hide , char *dirname_region, int cutDepth){

	FILE	*file_region_i;
	FILE	*file_region_stat_o;	
	FILE	*file_region_ratio_o;

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
	uint64_t	valid_region;
	uint64_t	valid_all_region = 0;
	int	have_region = 0;
	int	start,end;
	uint64_t	map_in;
	uint64_t	map_out;
	uint64_t	map_all_in;
	uint64_t	map_all_out;

	covBaseDistribution BaseDist_All;
	covRegionDistribution RegionDist_All;
	
	
	memset(&RegionDist_All , 0, sizeof(covRegionDistribution));
	

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
	if (flag_hide != 1){
		for (i = 0;i < n_ref;i++){
			printf("%s\t%d\n",chr_name[i], chr_length[i]);
		}
	}
	top = buffer + BgfzTail.I_size;
	counter	= top - address;


	sprintf(filename,"Region_stat_%d.txt",cutDepth);
	file_region_stat_o 	= fopen(filename,"w");	
	memset(filename, 0, sizeof(char)*name_len);
	
	sprintf(filename,"Region_ratio_%d.txt",cutDepth);
	file_region_ratio_o 	= fopen(filename,"w");
	memset(filename, 0, sizeof(char)*name_len);
	
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

			//Chromosome Change
			if (AlignmentHeader.refID != ref_ID){
				if (ref_ID != -1){
					if (flag_hide == 0 ){
						printf("[Bam File Unzip %d / %d ] %s done\n",ref_ID+1,n_ref,chr_name[ref_ID]);
					}
					if (have_region == 1){
			PrintRegionDist (file_region_stat_o, file_region_ratio_o, file_region_i, PosCoverage, chr_name[ref_ID], cutDepth, &RegionDist_All  , 1, in_region);
						fclose(file_region_i);
						have_region = 0;
					}
					free(in_region);
					free(PosCoverage);
				}	
				if ( AlignmentHeader.refID != -1){	
					PosCoverage = calloc(chr_length[AlignmentHeader.refID],sizeof(posCoverage));
					in_region = calloc(chr_length[AlignmentHeader.refID],sizeof(uint8_t));
					strcat_triple(filename, dirname_region, chr_name[AlignmentHeader.refID], "_region.txt", name_len);
					if ((file_region_i = fopen(filename,"r")) != NULL){
						have_region = 1;
						while (fscanf(file_region_i,"%d %d",&start,&end)!=EOF){
							for (i = start;i < end;i++){
								if (in_region[i] == 0){
									in_region[i] = 1;
									valid_region++;
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
			if ((AlignmentHeader.FLAG&4) == 0 && (AlignmentHeader.FLAG&1024) == 0){
				alignment_DepthTxt(address, &AlignmentHeader, PosCoverage, chr_length[ref_ID]);
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
			PrintRegionDist(file_region_stat_o  , file_region_ratio_o  ,file_region_i, PosCoverage, chr_name[ref_ID], cutDepth, &RegionDist_All, 1, in_region);
			have_region = 0;
		}
		free(PosCoverage);
	}
	
	PrintRegionDist (file_region_stat_o, file_region_ratio_o, file_region_i, PosCoverage, "", cutDepth, &RegionDist_All, 2, in_region);

	free(stream_i);	
	free(stream_o);	
	free(buffer);	
	return 0;
}
