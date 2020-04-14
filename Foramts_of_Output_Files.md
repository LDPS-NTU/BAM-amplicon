# Formats of Output Files

## del
* Output File Names  
  * The naming convention of the output files is as `[ref-name]_del.txt`, e.g. `chr1_del.txt`.
* Formats  
```
# The first tag, POS, is the starting position of a deletion.
# The second tag, LENGTH, is the length of that deletion.
# POS_LENGTH

[user@local]$ bam-utility -m del -b ./example/foo_del.bam
[user@local]$ cat chr1_del.txt
39421_1
39421_1
```

## ampsummary
* File Name  
  * After executing `ampsummary`, there are three kinds of output files, including `Region_InOut.txt`, `Region_ratio_[threshold].txt`, and `Region_stat_[threshold].txt`.  
  * The file, `Region_InOut.txt`, shows how many reads are inside or outside the regions of amplicons.  
  * The files, `Region_stat_[threshold].txt`, shows the coverage rates within amplicons in the specified region (for a given threshold). 
  * The files, `Region_ratio_[threshold].txt`, shows the summary of reads within amplicons in the specified region (for a given threshold).  
* Formats
```
[user@local]$ bam-utility -m ampsummary -b ./example/foo.bam -r example/foo.bed

-----
## Region_InOut.txt
# The first line:
# The 1st column, CHROM, is the reference name.
# The 2nd column, NUM_OF_READS_I, is the number of reads within the amplicons on CHROM
# The 3rd column, NUM_OF_READS_O, is the number of reads outside the amplicons on CHROM
# The 4th column, LENGTH, is the total length of the amplicons on CHROM
# CHROM NUM_OF_READS_I NUM_OF_READS_O NUM_OF_BASES_I

# The 'Total' line:
# The 1st column indicates that it is the 'total' line.
# The 2nd to the 4th columns are the same as the first line.
# The last two columns are the rates of the reads inside/outside the regions of the amplicons
# 2/(2+0) = 1.00000
# 0/(2+0) = 0.00000

[user@local]$ cat Region_InOut.txt
chr1	2	0	12
Total	2	0	12	1.000000	0.000000

-----
## Region_stat_0.txt
# The first line
# The 1st column, CHROM, is reference name
# The 2nd column, CR_0   , is the number of amplicons where their coverage rates are 0%
# The 3rd column, CR_50  , is the number of amplicons where their coverage rates are > 0% but <= 50%
# The 4th column, CR_90  , is the number of amplicons where their coverage rates are > 50% but <= 90%
# The 5th column, CR_99  , is the number of amplicons where their coverage rates are > 90% but < 100%
# The 6th column, CR_100 , is the number of amplicons where their coverage rates are 100%

# The total line
# The 1st column indicates that it is the 'Totel' line.
# The 2nd to the 6th columns are the same as the first line.
# The last one column is the sum of the coverage rates of all amplicons.
## For exmaple, if there are 2 amplicons, the cover rates are 0.97 and 0.84, respectively.
## The value in the last column is '0.97+0.84'.

[user@local]$ cat Region_stat_0.txt
chr1	0	0	0	0	1
Total	0	0	0	0	1	1.000000

-----
## Region_ratio_0.txt
# The 1st column, CHROM, is reference name of the amplicon
# The 2nd column, S_POINT ,is the start position of the amplicon (0-based, BED File format)
# The 3rd column, E_POINT , is the end position of the amplicon (1-based, BED File format)
# The 4th column, LENGTH, is the length of the amplicon
# The 5th column, LENGTH_NC, is the non-covered length of the amplicon
# The 6th column, NUM_OF_BASES_I, is the number of bases inside the regions of the amplicon
# The 7th column, AVE_OF_DEPTH, is the average depth of the amplicon
# The 8th column, NUM_OF_BASES_I_B, is the number of bases inside the regions of the amplicon, but is adjusted by cover rate while amplicons are overlapped (*BETA*)
# The 9th column, AVE_OF_DEPTH_B, is the average depth of the amplicon, but is adjusted by cover rate while amplicons are overlapped (*BETA*)
# The 10th column, CR, is the coverage rate (under the given threshold, e.g. 0, 10, and 20)
# The 11th column, NUM_OF_READS_I, is the number of reads inside the regions of the amplicon

## The formulas
# The 4th column = (The 3rd column) - (The 1st column)
# The 7th column = (The 6th column) / (The 4th column)
# The 9th column = (The 8th column) / (The 4th column)
# The 11th column = 1 - (The 5th column) / (The 4th column)


[user@local]$ cat Region_ratio_0.txt
chr1	39414	39426	12	0	18	1.500000	18	1.500000	1.000000	2

* The information: 
* This amplicon's length is 12. 
* There are 2 reads inside this amplicon and there are 18 bases cover this amplicon. 
* The average depth is 1.5000
* The cover rate of amplicon is 1.00000, or 100%. (all positions of amplicon have more than one reads covering.)
```