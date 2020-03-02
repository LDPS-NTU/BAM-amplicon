//======================
//Operation of Bam File
//1.BamStat:	Statistics of Reads
//2.BamRegion:	Statistics of InOrOut Target Reads
//3.BamQual:	Qual information per site
//4.BamCovTxt:	Print Mutation Site
//5.BamCov:	Create Coverage information Every Site (binary file)
//6.BamSam:	Create Sam-liked File
//======================


int	BamStat		(FILE *file_bam_i, FILE *file_length_o, int flag_hide);
int	BamRegion	(FILE *file_bam_i, int flag_hide, char *dirname_region, int flag_bin);
int	BamQual		(FILE *file_bam_i, FILE *file_length_o, int flag_hide, char *dirname_region, int flag_bin);
int	BamCovTxt	(FILE *file_bam_i, FILE *file_length_o, int flag_hide);
int	BamCovTxtDel	(FILE *file_bam_i, FILE *file_length_o, int flag_hide, char *dirname_chr);
int	BamCov		(FILE *file_bam_i, FILE *file_length_o, int flag_hide);
int	BamSam		(FILE *file_bam_i, FILE *file_length_o, int flag_hide);
