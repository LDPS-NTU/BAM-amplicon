# Formats of Output Files

## del
* Output File Names  
  * The output files of "del" function is named as `[ref-name]_del.txt` e.g. `chr1_del.txt`.
* The Formats
```
# The first tag, POS, is the starting position of a deletion.
# The second tag, LENGTH, is the length of that deletion.
# For each read, the output could look like: POS1_LENGTH1, POS2_LENGTH2, POS3_LENGTH3, ...
# Each read forms a seperate line.
# The following exmaple shows the case with two reads.
# A deletion (Length = 1) at the position, 39421, can be observed for both two reads. 

[user@local]$ bam-utility -m del -b ./example/foo_del.bam
[user@local]$ cat chr1_del.txt
39421_1
39421_1
```

## ampsummary
* Output File Names
  * After executing `ampsummary`, there are three kinds of output files, including `Region_InOut.txt`, `Base_stat.txt`,  `Region_ratio_[threshold].txt`, and `Region_stat_[threshold].txt`.  
  * The file, `Region_InOut.txt`, shows how many reads are inside/outside the regions of amplicons.  
  * The file, `Base_stat.txt`, shows the distribution of the depth of the targeted positions.  
  * The files, `Region_stat_[threshold].txt`, shows the coverage rates within the region of amplicons (for a given threshold in depth). 
  * The files, `Region_ratio_[threshold].txt`, shows the summary of reads within the region of amplicons (for a given threshold in depth).  
* Formats
```
## Region_InOut.txt
# Each reference forms a seperate line.
# The 1st column, CHROM, is the reference name.
# The 2nd column, NUM_OF_READS_I, is the number of reads within the region of amplicons on CHROM
# The 3rd column, NUM_OF_READS_O, is the number of reads outside the region of amplicons on CHROM
# The 4th column, NUM_OF_BASES_I, is the total length of amplicons on CHROM
# The output is look like: CHROM NUM_OF_READS_I NUM_OF_READS_O NUM_OF_BASES_I

# The 'Total' line:
# The 1st column indicates that it is the 'Total' line.
# The 2nd column, TOTAL_NUM_OF_READS_I, is the number of reads within the regions of amplicons for the reference(s) above
# The 3rd column, TOTAL_NUM_OF_READS_O, is the number of reads outside the regions of amplicons for the reference(s) above
# The 4th column, TOTAL_NUM_OF_BASES_I, is the total length of amplicons in the reference(s) above
# The last two columns are the rates of the reads inside/outside the regions of the amplicons
# 4/(4+2) = 0.666667
# 2/(4+2) = 0.333333

[user@local]$ bam-utility -m ampsummary -b ./example/foo_amp.bam -r example/foo_amp.bed
[user@local]$ cat Region_InOut.txt
chr1	2	0	12
chr2	2	2	12
Total	4	2	24	0.666667	0.333333

-----
## Base_stat.txt
# Each reference forms a seperate line.
# The 1st column, CHROM, is the reference name.
# The 2nd column, DP_0, is the number of targeted positions where their depth are 0 on CHROM
# The 3rd column, DP_9, is the number of targeted positions where their depth are >= 1 but <= 9 on CHROM
# The 4th column, DP_19, is the number of targeted positions where their depth are >= 10 but <= 19 on CHROM
# The 5th column, DP_29, is the number of targeted positions where their depth are >= 20 but <= 29 on CHROM
# The 6th column, DP_39, is the number of targeted positions where their depth are >= 30 but <= 39 on CHROM
# The 7th column, DP_49, is the number of targeted positions where their depth are >= 40 but <= 49 on CHROM
# The 8th column, DP_59, is the number of targeted positions where their depth are >= 50 but <= 59 on CHROM
# The 9th column, DP_69, is the number of targeted positions where their depth are >= 60 but <= 69 on CHROM
# The 10th column, DP_79, is the number of targeted positions where their depth are >= 70 but <= 79 on CHROM
# The 11th column, DP_89, is the number of targeted positions where their depth are >= 80 but <= 89 on CHROM
# The 12th column, DP_99, is the number of targeted positions where their depth are >= 90 but <= 99 on CHROM
# The 13th column, DP_100, is the number of targeted positions where their depth are >= 100 on CHROM

# The 'Total' line:
# The 1st column indicates that it is the 'Total' line.
# The 2nd column, DP_0, is the number of targeted positions where their depth are 0 for the reference(s) above
# The 3rd column, DP_9, is the number of targeted positions where their depth are >= 1 but <= 9 for the reference(s) above
# The 4th column, DP_19, is the number of targeted positions where their depth are >= 10 but <= 19 for the reference(s) above
# The 5th column, DP_29, is the number of targeted positions where their depth are >= 20 but <= 29 for the reference(s) above
# The 6th column, DP_39, is the number of targeted positions where their depth are >= 30 but <= 39 for the reference(s) above
# The 7th column, DP_49, is the number of targeted positions where their depth are >= 40 but <= 49 for the reference(s) above
# The 8th column, DP_59, is the number of targeted positions where their depth are >= 50 but <= 59 for the reference(s) above
# The 9th column, DP_69, is the number of targeted positions where their depth are >= 60 but <= 69 for the reference(s) above
# The 10th column, DP_79, is the number of targeted positions where their depth are >= 70 but <= 79 for the reference(s) above
# The 11th column, DP_89, is the number of targeted positions where their depth are >= 80 but <= 89 for the reference(s) above
# The 12th column, DP_99, is the number of targeted positions where their depth are >= 90 but <= 99 for the reference(s) above
# The 13th column, DP_100, is the number of targeted positions where their depth are >= 100 for the reference(s) above
# The last two columns are the number of bases and the number of deletions in the the regions of the amplicons.

[user@local]$ cat Base_stat.txt
chr1	0	12	0	0	0	0	0	0	0	0	0	0
chr2	1	11	0	0	0	0	0	0	0	0	0	0
Total	1	23	0	0	0	0	0	0	0	0	0	0	35	2

-----
## Region_stat_[threshold].txt

# Results of different depth thresholds are stored in different files.
# For example, Region_stat_0.txt contains the statistics where the depth threshold is set to 0.

# Each reference forms a seperate line.
# The 1st column, CHROM, is reference name.
# The 2nd column, CR_0, is the number of amplicons where their coverage rates are 0%
# The 3rd column, CR_50, is the number of amplicons where their coverage rates are > 0% but < 50%
# The 4th column, CR_90, is the number of amplicons where their coverage rates are >= 50% but < 90%
# The 5th column, CR_99, is the number of amplicons where their coverage rates are >= 90% but < 100%
# The 6th column, CR_100, is the number of amplicons where their coverage rates are 100%

# The 'Total' line
# The 1st column indicates that it is the 'Totel' line.
# The 2nd column, CR_0, is the total number of amplicons where their coverage rates are 0% for the reference(s) above
# The 3rd column, CR_50, is the total number of amplicons where their coverage rates are > 0% but < 50% for the reference(s) above
# The 4th column, CR_90, is the total number of amplicons where their coverage rates are >= 50% but < 90% for the reference(s) above
# The 5th column, CR_99, is the total number of amplicons where their coverage rates are >= 90% but < 100% for the reference(s) above
# The 6th column, CR_100, is the total number of amplicons where their coverage rates are 100% for the reference(s) above
# The last column is the sum of the coverage rates of all amplicons.

# For exmaple, if there are 3 amplicons, the coverage rates are 1.000000, 0.875000 and 0.909091, respectively.
# The value in the last column will be shown as '2.784091'.

[user@local]$ cat Region_stat_0.txt
chr1	0	0	0	0	1
chr2	0	0	1	1	0
Total	0	0	1	1	1	2.784091

-----
## Region_ratio_[threshold].txt

# Results of different depth thresholds are stored in different files.
# For example, Region_ratio_0.txt contains the statistics where the depth threshold is set to 0.

# Each amplicon forms a seperate line.
# The 1st column, CHROM, is reference name.
# The 2nd column, S_POINT, is the starting position of the amplicon (0-based, BED File format).
# The 3rd column, E_POINT, is the ending position of the amplicon (1-based, BED File format).
# The 4th column, LENGTH, is the length of the amplicon.
# The 5th column, LENGTH_NC, is the non-covered length of the amplicon.
# The 6th column, NUM_OF_BASES_I, is the total number of bases in reads falling into the regions of the amplicon.
# The 7th column, AVE_OF_DEPTH, is the average depth of the amplicon.
# The 8th column, NUM_OF_BASES_I_B, is the number of bases inside the regions of the amplicon, but is adjusted by the coverage rate when the amplicon is overlapped with other amplicons (*BETA*)
# The 9th column, AVE_OF_DEPTH_B, is the average depth of the amplicon, but is adjusted by the coverage rate when the amplicon is overlapped with other amplicons (*BETA*)
# The 10th column, CR, is the coverage rate under the given threshold, e.g. 0, 10, and 20.
# The 11th column, Q3_OF_DEPPTH, is the third quartile (Q3) depth of the amplicon.

## The formulas
# The 4th column = (The 3rd column) - (The 1st column)
# The 7th column = (The 6th column) / ((The 4th column)-(The 5th column))
# The 9th column = (The 8th column) / ((The 4th column)-(The 5th column))

# The following example shows that the length of the amplicon is 12. (Note that the starting position is 0-based, while the ending position is 1-based.)
# There are 2 two reads covering the amplicon, and the total read length within the amplicon region is 18.
# The average depth is 1.5000, and the coverage rate of the amplicon is 100%. (Every positions of the amplicon is covered by at least one read.)

[user@local]$ cat Region_ratio_0.txt
chr1	39414	39426	12	0	18	1.500000	18	1.500000	1.000000	2
chr2	19414	19422	8	1	11	1.571429	5	0.714286	0.875000	0
chr2	19415	19426	11	1	16	1.600000	12	1.200000	0.909091	2

```
