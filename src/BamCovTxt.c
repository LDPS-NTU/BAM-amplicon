#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include "zlib.h"

#include "BamCommonLibrary.h"





int	BamCovTxt(FILE *file_bam_i, FILE *file_length_o, int flag_hide){

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

	int32_t	data_length;
	int32_t	text_length;
	int32_t	l_name;

	int	SIZE_alignmentHeader = sizeof(alignmentHeader);




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

	data_length = 65536 - 8;
	text_length = 65536 - 8;


	while(l_text - data_length > 0){
//		memcpy(text, address+(65536-text_length), sizeof(char)*text_length);
//		printf("%s\n",address+(65536-text_length));
		
		len_data = BGFZBlock(file_bam_i);
		decompressBlock(&infstream, stream_i, stream_o, len_data, file_bam_i);
		fread(&BgfzTail	,sizeof(bgfzTail),1,file_bam_i);
		memcpy(buffer,stream_o,65536);
		address = buffer;	
		
		l_text = l_text - text_length;
		data_length = 65536;
		text_length = 65536;
	}

//	memcpy(text,address+8,sizeof(char)*l_text);
//	memcpy(&n_ref,address+8+l_text,sizeof(int32_t));

//	if (flag_header){
//		printf("%d\n%s%d\n",l_text,text,n_ref);	
//	}
//	address += (12+l_text);

	memcpy(&n_ref,address+(65536-text_length)+l_text,sizeof(int32_t));
	address += ((65536-text_length)+4+l_text);

	chr_name	= malloc(n_ref*sizeof(char *));	
	chr_length	= calloc(n_ref,sizeof(int));

	top = buffer + BgfzTail.I_size;
	counter	= top - address;

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
//	top = buffer + BgfzTail.I_size;
	counter	= top - address;


	//Start Unzip and Process Bam File
	
	while ( (len_data = BGFZBlock(file_bam_i)) > 0){
		decompressBlock(&infstream, stream_i, stream_o, len_data, file_bam_i);
		fread(&BgfzTail	,sizeof(bgfzTail),1,file_bam_i);

		memcpy(buffer,address,sizeof(uint8_t)*counter);
		memcpy(&buffer[counter],stream_o,BgfzTail.I_size);
		top	= buffer + counter + BgfzTail.I_size;
		address = buffer;
		
		if (top - address > SIZE_alignmentHeader){
			memcpy(&AlignmentHeader, address, SIZE_alignmentHeader);
		}
	
		while ( top - address >= (AlignmentHeader.block_size + 4) ){
			if (AlignmentHeader.refID != ref_ID){
				if (ref_ID != -1){
					if (flag_hide == 0){
						printf("[Bam File Unzip %d / %d ] %s done\n",ref_ID+1,n_ref,chr_name[ref_ID]);
					}	
					for (i = 0;i < chr_length[ref_ID];i++){		
						if (Coverage[i] != 0){
							All_Coverage = PosCoverage[i].A+ PosCoverage[i].C+ PosCoverage[i].G+ PosCoverage[i].T+ PosCoverage[i].N;
							fprintf(file_cov_o,"%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",
							i+1,
							PosCoverage[i].A,
							PosCoverage[i].C,
							PosCoverage[i].G,
							PosCoverage[i].T,
							PosCoverage[i].N,
							Coverage[i]-All_Coverage,
							Coverage[i]);	
						}	
					}
					free(PosCoverage);
					free(Coverage);
					fclose(file_cov_o);
				}
				if (AlignmentHeader.refID != -1){
					PosCoverage = calloc(chr_length[AlignmentHeader.refID],sizeof(posCoverage));
					Coverage = calloc(chr_length[AlignmentHeader.refID],sizeof(uint32_t));
					strcat_triple(filename, "", chr_name[AlignmentHeader.refID], "_cov.txt", name_len);
					if ((file_cov_o = fopen(filename,"w"))==NULL){
						printf("Error\n");	
					}
				}
				ref_ID = AlignmentHeader.refID;
			}
			if (ref_ID == -1){	break;	}
			address += SIZE_alignmentHeader;
		
			//Core_START
			if ((AlignmentHeader.FLAG&4) == 0 && (AlignmentHeader.FLAG&1024) == 0){
				alignment_CoverageTxt(address, &AlignmentHeader, PosCoverage, Coverage, chr_length[ref_ID]);
			}
			//Core_END
			
			address += (AlignmentHeader.block_size + 4) - SIZE_alignmentHeader;

			if (top - address > SIZE_alignmentHeader){
				memcpy(&AlignmentHeader, address, SIZE_alignmentHeader);
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
		for (i = 0;i < chr_length[ref_ID];i++){
			if (Coverage[i] != 0){
				All_Coverage = PosCoverage[i].A+ PosCoverage[i].C+ PosCoverage[i].G+ PosCoverage[i].T+ PosCoverage[i].N;
				fprintf(file_cov_o,"%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",
				i+1,
				PosCoverage[i].A,
				PosCoverage[i].C,
				PosCoverage[i].G,
				PosCoverage[i].T,
				PosCoverage[i].N,
				Coverage[i] - All_Coverage,
				Coverage[i]);	
			}	
		}
		free(PosCoverage);
		free(Coverage);
		fclose(file_cov_o);
	}
	
	free(stream_i);	
	free(stream_o);	
	free(buffer);	
	return 0;
}

