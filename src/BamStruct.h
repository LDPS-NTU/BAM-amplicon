#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#define	name_len	10000
#define	LINE_MAX_LEN	10000
#define	LINE_LEN	10000

typedef struct faiTabel_st {
	uint32_t	num_chr;
	char	**chr_name;
	uint32_t	*chr_length;
	uint32_t	*index;
	uint32_t	*line_len;
	uint32_t	*line_len_real;
} __attribute__((packed)) faiTable;

typedef struct bedTable_st {
	uint32_t	*table_start;	
	uint32_t	*table_end;
	uint32_t	*table_max;

	uint32_t	*ref_id;
	uint32_t	*start;
	uint32_t	*end;
	char		**pattern;

} __attribute__((packed)) bedTable;

typedef struct baiTable_st {

	char	BAMI[4];
	int32_t	n_ref_bai;
	int32_t		*n_intv;
	uint64_t	**ioffset;

} __attribute__((packed)) baiTable;

typedef struct bamHeader_st {
	char	magic[4];
	char	**chr_name;
	int	*chr_length;
	int32_t	l_text;
	char	*text;
	int32_t	n_ref;
	
} __attribute__((packed)) bamHeader;

typedef struct toolsFlags_st {
	uint8_t	flag_dup;
	uint8_t	flag_hide;
	uint8_t	flag_simple;
	uint8_t	flag_mapq;
	uint8_t	flag_header;

	uint8_t	flag_coverage;
	char	*coverage_threshold;



	/* Ion Torrent */
	uint8_t	flag_zm;
	uint8_t	flag_flow;
	/*---*/	

	uint8_t	flag_origin;
	uint8_t	flag_type;
	uint8_t	flag_pattern;
	char	*pattern;
	
	uint8_t	flag_target;
	uint8_t	flag_columnName;

	FILE	*file_target;
	/*
	0:	SNP
	1:	Deletion
	2:	Insertion
	*/
	char	ref;
	char	alt;
	uint8_t	type;
	
	uint8_t	flag_total;
	uint8_t	flag_bed;
	uint8_t	flag_region;
	uint8_t	flag_point;
	uint8_t	flag_filter;

	char	*chromosome;
	uint32_t	start;
	uint32_t	end;
	uint32_t	length;
	uint32_t	num_ref;
	uint32_t	num_alt;
	
} __attribute__((packed)) toolsFlags;

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


typedef	struct bgfztail_st {
	uint32_t CRC32;
	uint32_t I_size;
} __attribute__((packed)) bgfzTail;

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
} __attribute__((packed)) alignmentHeader;

typedef struct posInformation_st{
	char	genotype;
	char	het;
	uint16_t	A;
	uint16_t	T;
	uint16_t	G;
	uint16_t	C;
} __attribute__((packed)) posInformation;

typedef struct count_st{
	uint16_t	A;
	uint16_t	T;
	uint16_t	G;
	uint16_t	C;
} __attribute__((packed)) count;

typedef struct posCoverage_st{
	uint32_t	A;
	uint32_t	T;
	uint32_t	G;
	uint32_t	C;
	uint32_t	N;
	uint32_t	DEL;
	uint32_t	ALL;
} __attribute__((packed)) posCoverage;

typedef struct posData_st{
	char	genotype;
	char	het;
} __attribute__((packed)) posData;

typedef struct posQuality_st{
	uint64_t	base_qual_sum;
	uint64_t	base_qual_pow;
	uint64_t	map_qual_sum;
	uint64_t	map_qual_pow;
	uint64_t	cov_sum;

	uint64_t	A_F[128];
	uint64_t	C_F[128];
	uint64_t	G_F[128];
	uint64_t	T_F[128];

	uint64_t	A_R[128];
	uint64_t	C_R[128];
	uint64_t	G_R[128];
	uint64_t	T_R[128];

	uint64_t	sum_A_F;
	uint64_t	sum_C_F;
	uint64_t	sum_G_F;
	uint64_t	sum_T_F;

	uint64_t	sum_A_R;
	uint64_t	sum_C_R;
	uint64_t	sum_G_R;
	uint64_t	sum_T_R;


} __attribute__((packed)) posQuality;

typedef struct covRegionDistribution_st{
	int	per_0;
	int	per_50;
	int	per_90;
	int	per_99;
	int	per_100;
	double	per_all;
} __attribute__((packed)) covRegionDistribution;


typedef struct covBaseDistribution_st{
	uint64_t	cov_0;
	uint64_t	cov_9;
	uint64_t	cov_19;
	uint64_t	cov_29;
	uint64_t	cov_39;
	uint64_t	cov_49;
	uint64_t	cov_59;
	uint64_t	cov_69;
	uint64_t	cov_79;
	uint64_t	cov_89;
	uint64_t	cov_99;
	uint64_t	cov_100;
	uint64_t	cov_all;
	uint64_t	cov_del;
} __attribute__((packed)) covBaseDistribution;

