#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "BamStruct.h"

void	createRegionTable(FILE *file_bed_i, bedTable *BedTable, bamHeader *BamHeader, toolsFlags *ToolsFlags){

	char	line[LINE_MAX_LEN];
	char	chromosome_name[LINE_MAX_LEN];
	int	num_target	= 0;
	int	i,j;

	int	start	= ToolsFlags->start;
	int	end	= ToolsFlags->end;
	char	*chromosome	= ToolsFlags->chromosome;

	if (file_bed_i != NULL){
		//Bed File
		while	(fgets(line, LINE_MAX_LEN, file_bed_i) != NULL){
			num_target++;
		}
		//printf("%d\n",num_target);

		BedTable->table_start = calloc (BamHeader->n_ref, sizeof(uint32_t));
		BedTable->table_end = calloc (BamHeader->n_ref, sizeof(uint32_t));
		BedTable->table_max = calloc (BamHeader->n_ref, sizeof(uint32_t));
		
		
		BedTable->start = calloc (num_target, sizeof(uint32_t));
		BedTable->end = calloc (num_target, sizeof(uint32_t));
	
		rewind(file_bed_i);
	
		for(i = 0;i < num_target;i++){
			fscanf(file_bed_i, "%s %d %d", chromosome_name, &BedTable->start[i], &BedTable->end[i]);
			fgets(line, LINE_MAX_LEN, file_bed_i);	
			if (i == 0){
				for(j = 0; j < BamHeader->n_ref; j++){
					if (strcmp(chromosome_name, BamHeader->chr_name[j])==0){
						BedTable->table_start[j] = i;
						BedTable->table_max[j] = BedTable->end[i];
						break;
					}
				}
			}else if (strcmp(chromosome_name, BamHeader->chr_name[j])!=0){
				BedTable->table_end[j] = i;
				if (BedTable->end[i] > BedTable->table_max[j]){
					BedTable->table_max[j] = BedTable->end[i];
				}
				for(j = 0; j < BamHeader->n_ref; j++){
					if (strcmp(chromosome_name, BamHeader->chr_name[j])==0){
						BedTable->table_start[j] = i;
						BedTable->table_max[j] = BedTable->end[i];
						break;
					}
				}
			}else {
				if (BedTable->end[i] > BedTable->table_max[j]){
					BedTable->table_max[j] = BedTable->end[i];
				}	
			}
			if (j == BamHeader->n_ref){
				printf("Error: %d is over the number\n", j);
				exit(-1);
			}
		//printf("%s %d %d\n", chromosome_name, BedTable->start[i], BedTable->end[i]);
		}
		BedTable->table_end[j] = i;
		
	}else if (end > 0){
		//printf("Region\n");
		BedTable->table_start = calloc (BamHeader->n_ref, sizeof(uint32_t));
		BedTable->table_end = calloc (BamHeader->n_ref, sizeof(uint32_t));
		BedTable->table_max = calloc (BamHeader->n_ref, sizeof(uint32_t));
		
		BedTable->start = calloc (1, sizeof(uint32_t));
		BedTable->end = calloc (1, sizeof(uint32_t));
		num_target = 1;
		BedTable->start[0] = start-1;
		BedTable->end[0] = end;
		
		for(j = 0; j < BamHeader->n_ref; j++){
			if (strcmp(chromosome, BamHeader->chr_name[j])==0){
				BedTable->table_start[j] = 0;
				BedTable->table_end[j] = 1;
				BedTable->table_max[j] = end;
				if( start > BamHeader->chr_length[j] || end > BamHeader->chr_length[j] || start > end){
					printf("[Error] Position is bigger than Size of %s\n",chromosome);
					exit(-1);	
				}
				break;
			}
		}
	}else if (start > 0){
		//printf("SinglePoint\n");
		BedTable->table_start = calloc (BamHeader->n_ref, sizeof(uint32_t));
		BedTable->table_end = calloc (BamHeader->n_ref, sizeof(uint32_t));
		BedTable->table_max = calloc (BamHeader->n_ref, sizeof(uint32_t));
		BedTable->start = calloc (1, sizeof(uint32_t));
		BedTable->end = calloc (1, sizeof(uint32_t));
		num_target = 1;
		BedTable->start[0] = start-1;
		BedTable->end[0] = start;
		
		for(j = 0; j < BamHeader->n_ref; j++){
			if (strcmp(chromosome, BamHeader->chr_name[j])==0){
				BedTable->table_start[j] = 0;
				BedTable->table_end[j] = 1;
				BedTable->table_max[j] = start+1;
				if( start > BamHeader->chr_length[j]){
					printf("[Error] Position is bigger than Size of %s\n",chromosome);
					exit(-1);	
				}
				break;
			}
		}	
	}else {
		//printf("Total\n");
		BedTable->table_start = calloc (BamHeader->n_ref, sizeof(uint32_t));
		BedTable->table_end = calloc (BamHeader->n_ref, sizeof(uint32_t));
		BedTable->table_max = calloc (BamHeader->n_ref, sizeof(uint32_t));
		
		BedTable->start = calloc (BamHeader->n_ref, sizeof(uint32_t));
		BedTable->end = calloc (BamHeader->n_ref, sizeof(uint32_t));
		num_target = BamHeader->n_ref;
			
		for(j = 0; j < BamHeader->n_ref; j++){
			BedTable->start[j] = 0;
			BedTable->end[j] = BamHeader->chr_length[j];
			BedTable->table_start[j] = j;
			BedTable->table_end[j] = j+1;
			BedTable->table_max[j] = BamHeader->chr_length[j];
		}		
	}
}
