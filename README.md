# BAM-utility

## Getting Started
```bash
git clone https://github.com/YuCheng-Group/BAM-utility.git
cd BAM-utility
make
export PATH=$PATH:"$PWD/bin/"
```


## Users' Guide
```
Program: bam-utility (for BAM file)
Usage:  bam-utility -m <Mode> -b <bam_filename> [Options]

Mode:	depthdist	show depths
	quality		show distribution of quality scores at a specified position
	
	pattern		show frequency of patterns within a specified region	
	del		show frequency of deletions with position information
	poly		show read sequences with homopolymers

	stat		show BAM statistics
	length		show distribution of mapped read lengths
	ampsummary	show depths and coverage rates within specified amplicons
	trim		trim the BAM file using amplicons provided by the BED file

Options:
	-b, --bam    [FILE]	input BAM file
	-r, --bed    [FILE]	the ranges of specified regions from BED file
	-c, --chr    [STR]	chromosome name
	-s, --start  [INT]	start position (1-based, i.e. first base is 1)
	-e, --end    [INT]	end position (1-based, i.e. first base is 1)
	-g, --fasta  [FILE]	fasta file (reference file)
	-t, --target [FILE]	the ranges of amplicon regions from BED file
	-f, --filter [TYPE]	filter criteria
				[ readqual, mapqual]
	-l, --compact		show only average depth and coverage rate (depthdist mode only)
	-a, --origin		show the original sequence (include insertion)
	-d, --duplicate		show duplicate
	-v, --verbose		don't show the processed status
	-n, --column_name	remove the first row (header)
	-u, --threshold	[INT]	coverage above the threshold
	/* For Ion Torrent Data */
	-w, --flow		modify the sequence according to flow signals (FZ)
	-z, --zm		change the FZ values to ZM values
	-p, --pattern	[STR]	show the percentage of ALT pattern specified (pattern mode)
	-h, --help		show the manual page    
```

## Examples
### stat
```
[user@local]$ bam-utility -m stat -b ./example/foo.bam
chr1	248956422
Item	Mapped	UnMap	Total	(Mapped)
===== ===== ===== ===== ===== ===== ===== =====
Sample	2	0	2	(1.000000)
R1	1	0	1	(1.000000)
R2	1	0	1	(1.000000)
----- ----- ----- ----- ----- ----- ----- -----
UnChimeric_Sample	2	0	2	(1.000000)
UnChimeric_R1	1	0	1	(1.000000)
UnChimeric_R2	1	0	1	(1.000000)
----- ----- ----- ----- ----- ----- ----- -----
Proper	2	0	2
R1_PP	1	0	1	(1.000000)
R2_PP	1	0	1	(1.000000)
----- ----- ----- ----- ----- ----- ----- -----
Dupli	0
R1_Dup	0
R2_Dup	0
===== ===== ===== ===== ===== ===== ===== =====
Item	Proper	Total	(Properly)
===== ===== ===== ===== ===== ===== ===== =====
SameChr	2	2	(1.000000)
Cross	0	0	(-nan)
Single	0	0	(-nan)
===== ===== ===== ===== ===== ===== ===== =====
```
### pattern

```
[user@local]$ bam-utility -m pattern -b ./example/foo.bam -c chr1 -s 39418 -e 39420 -v
ID	PATTERN
read_001_R1	TGT
read_001_R2	TGT

[user@local]$ bam-utility -m pattern -b ./example/foo.bam -c chr1 -s 39418 -e 39420 -a -v
PATTERN
TGT
TGT

[user@local]$ bam-utility -m pattern -b ./example/foo.bam -c chr1 -s 39418 -e 39420 -p TGT -a -v
TOT.DP	PAT.FREQ **(TOT.DP: total depth; PAT.FREQ: pattern frequency)
2	1.00000
```

### depthdist
```
[user@local]$ bam-utility -m depthdist -b ./example/foo.bam -c chr1 -s 39414 -e 39420 -v
CHR	POS	A	C	G	T	N	DEL	TOTAL
chr1	39414	0	0	0	0	0	0	0
chr1	39415	1	0	0	0	0	0	1
chr1	39416	0	0	0	1	0	0	1
chr1	39417	0	0	1	0	0	0	1
chr1	39418	0	0	0	2	0	0	2
chr1	39419	0	0	2	0	0	0	2
chr1	39420	0	0	0	2	0	0	2

[user@local]$ bam-utility -m depthdist -b ./example/foo.bam -c chr1 -s 39414 -e 39420 -v -l
AVE.DP	COV.RT	**(AVE.DP: average depth; COV.RT: coverage rate)
1.5000	0.857143

[user@local]$ bam-utility -m depthdist -b ./example/foo.bam -c chr1 -s 39414 -e 39420 -v -l -u 2
CHR	chr1
THRESH	2X	** (THRESH: threshold)
TOT.DP	3
COV.RT	0.42857
```

```
[user@local]$ bam-utility -m quality -b ./example/foo.bam -c chr1 -s 39420 -v
Q.SCORE	A_FOR	C_FOR	G_FOR	T_FOR	A_REV	C_REV	G_REV	T_REV	**(Q.SCORE: quality score)
0	0	0	0	0	0	0	0	0	**(FOR: forward; REV: reverse)
1	0	0	0	0	0	0	0	0
2	0	0	0	0	0	0	0	0
3	0	0	0	0	0	0	0	0
4	0	0	0	0	0	0	0	0
5	0	0	0	0	0	0	0	0
6	0	0	0	0	0	0	0	0
7	0	0	0	0	0	0	0	0
8	0	0	0	0	0	0	0	0
9	0	0	0	0	0	0	0	0
10	0	0	0	0	0	0	0	0
11	0	0	0	0	0	0	0	0
12	0	0	0	0	0	0	0	0
13	0	0	0	0	0	0	0	0
14	0	0	0	0	0	0	0	0
15	0	0	0	0	0	0	0	0
16	0	0	0	0	0	0	0	0
17	0	0	0	0	0	0	0	0
18	0	0	0	0	0	0	0	0
19	0	0	0	0	0	0	0	0
20	0	0	0	0	0	0	0	0
21	0	0	0	0	0	0	0	0
22	0	0	0	0	0	0	0	0
23	0	0	0	0	0	0	0	0
24	0	0	0	0	0	0	0	0
25	0	0	0	0	0	0	0	0
26	0	0	0	0	0	0	0	0
27	0	0	0	0	0	0	0	0
28	0	0	0	0	0	0	0	0
29	0	0	0	0	0	0	0	0
30	0	0	0	0	0	0	0	0
31	0	0	0	0	0	0	0	0
32	0	0	0	0	0	0	0	0
33	0	0	0	0	0	0	0	0
34	0	0	0	0	0	0	0	0
35	0	0	0	0	0	0	0	0
36	0	0	0	1	0	0	0	1
37	0	0	0	0	0	0	0	0
38	0	0	0	0	0	0	0	0
39	0	0	0	0	0	0	0	0
40	0	0	0	0	0	0	0	0
41	0	0	0	0	0	0	0	0
42	0	0	0	0	0	0	0	0
43	0	0	0	0	0	0	0	0
44	0	0	0	0	0	0	0	0
45	0	0	0	0	0	0	0	0
46	0	0	0	0	0	0	0	0
47	0	0	0	0	0	0	0	0
48	0	0	0	0	0	0	0	0
49	0	0	0	0	0	0	0	0
50	0	0	0	0	0	0	0	0
51	0	0	0	0	0	0	0	0
52	0	0	0	0	0	0	0	0
53	0	0	0	0	0	0	0	0
54	0	0	0	0	0	0	0	0
55	0	0	0	0	0	0	0	0
56	0	0	0	0	0	0	0	0
57	0	0	0	0	0	0	0	0
58	0	0	0	0	0	0	0	0
59	0	0	0	0	0	0	0	0
60	0	0	0	0	0	0	0	0
61	0	0	0	0	0	0	0	0
62	0	0	0	0	0	0	0	0
63	0	0	0	0	0	0	0	0
```
