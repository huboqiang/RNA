RNA seq pipeline.

First build-up a dir-system class, like:
```python

class DirSystem(dict):
   def __init__(self):
      home_dir       =  os.path.abspath('./')

      self['dir'] = {                                                       \
         'raw_data'              : "%s/00.raw_data"         % (home_dir),   \
         'clean_data'            : "%s/00.clean_data"       % (home_dir),   \
         'tophat_dir'            : "%s/01.Tophat"           % (home_dir),   \
         'HTSeq_result_dir'      : "%s/02.HTSeq_result"     % (home_dir),   \
         'HTSeq_known_dir'       : "%s/02.1.HTSeq_known"    % (home_dir),   \
         'HTSeq_unknown_dir'     : "%s/02.2.HTSeq_unknown"  % (home_dir),   \
         'cufflinks_unknown_dir' : "%s/03.cufflinks_unknown"% (home_dir),   \
         'cuffquant_dir'         : "%s/04.cuffquant"        % (home_dir),   \
         'cuffnorm_dir'          : "%s/05.cuffnorm"         % (home_dir),   \
         'cuffquant_ERCC_dir'    : "%s/04.1.cuffquant.ERCC" % (home_dir),   \
         'cuffnorm_ERCC_dir'     : "%s/05.1.cuffnorm.ERCC"  % (home_dir),   \
         'cuffquant_k_dir'       : "%s/04.cuffquant_known"        % (home_dir),   \
         'cuffnorm_k_dir'        : "%s/05.cuffnorm_known"         % (home_dir),   \
         'cuffquant_ERCC_k_dir'  : "%s/04.1.cuffquant.ERCC_known" % (home_dir),   \
         'cuffnorm_ERCC_k_dir'   : "%s/05.1.cuffnorm.ERCC_known"  % (home_dir),   \
         'repeat_counts_dir'     : "%s/06.repeat_counts"    % (home_dir),   \
         'repeat_mrg_dir'        : "%s/07.merge_repFPKM"    % (home_dir),   \
			
```
to put the results. Raw fq data should be put into 00.0.raw_data at the very beginning.

Then , build-up a software system class:

```python
### Using your own path for these softwares.
class UsedSoftware(object):
   def __init__(self):
      self.py        = "/usr/bin/python"
      self.pl        = "/usr/bin/perl"
      self.tophat    = "/software/tophat-2.0.12.Linux_x86_64/tophat"
      self.cflk_dir  = "/software/cufflinks-2.2.1.Linux_x86_64"
      self.samtools  = "/software/samtools-0.1.18/samtools"
      self.bedtools  = "/software/bedtools-2.17.0/bin/bedtools"
      self.deseq     = "/usr/lib/python2.7/site-packages/HTSeq/scripts/count.py"
      self.bgzip     = "/usr/bin/bgzip"
      self.tabix     = "/usr/bin/tabix"
```

Next, find the input files. You can download these files in UCSC or so on and then using own-scripts to merge the ERCC information, and generate files in this format.
``` bash
==> samp_info <==
sample_name		sample_brief_name		stage			sample_group	ERCC_times	RFP_polyA	GFP_polyA	CRE_polyA	data_type	rename
A1					EPI_A1					Epiblast		RNA            2.50E-07    1.50E+04    1.50E+03    1.50E+02	 	PE          Epiblast_Epiblast_A1
A2					EPI_A2					Epiblast		RNA            2.50E-07    1.50E+04    1.50E+03    1.50E+02	 	PE          Epiblast_Epiblast_A2


==> bt_index_base <==
hg19 genome fa + ERCC + RGC information. Then generate the index using bowtie2-build.


==> rmsk_bed <==
chr1	10001	10468	1504	1.3	0.4	1.3	(249240153)	+	(CCCTAA)n	Simple_repeat	1	463	(0)	1
chr1	10469	11447	3612	11.4	27.0	1.3	(249239174)	C	TAR1	Satellite/telo	(399)	1712	483	2
chr1	11504	11675	437	23.5	18.6	3.5	(249238946)	C	L1MC	LINE/L1	(2236)	5646	5449	3
chr1	11678	11780	239	29.4	1.9	1.0	(249238841)	C	MER5B	DNA/hAT-Charlie	(74)	104	1	4
chr1	15265	15355	318	23.0	3.8	0.0	(249235266)	C	MIR3	SINE/MIR	(119)	143	49	5
chr1	16713	16749	203	16.2	0.0	0.0	(249233872)	+	(TGG)n	Simple_repeat	1	37	(0)	6
chr1	18907	19048	239	33.8	14.8	0.0	(249231573)	+	L2a	LINE/L2	2942	3104	(322)	7
chr1	19948	20405	652	34.6	8.5	4.2	(249230216)	+	L3	LINE/CR1	3042	3519	(970)	8
chr1	20531	20679	270	33.1	0.7	2.7	(249229942)	+	Plat_L3	LINE/CR1	2802	2947	(639)	9
chr1	21949	22075	254	27.9	4.7	3.9	(249228546)	+	MLT1K	LTR/ERVL-MaLR	15	142	(453)	10

==> ercc_info <==
ERCC_ID	length	concentration_in_Mix1(attomoles/ul)
ERCC-00002	1061	15000
ERCC-00003	1023	937.5
ERCC-00004	523	7500
ERCC-00009	984	937.5
ERCC-00012	994	0.11444092
ERCC-00013	808	0.91552734
ERCC-00014	1957	3.66210938
ERCC-00016	844	0.22888184
ERCC-00017	1136	0.11444092


==> genome_gtf_wl(with lncRNA)<==
chr1	noncoding	exon	11874	12227	0.000000	+	.	gene_id "DDX11L1"; transcript_id "NR_046018"; exon_number "1"
chr1	noncoding	exon	12613	12721	0.000000	+	.	gene_id "DDX11L1"; transcript_id "NR_046018"; exon_number "2"
chr1	Cufflinks	exon	12975	13052	0	+	.	gene_id "NONHSAG000001"; transcript_id "NONHSAT000004"; FPKM "0"; exon_number 6;
chr1	noncoding	exon	13221	14409	0.000000	+	.	gene_id "DDX11L1"; transcript_id "NR_046018"; exon_number "3"
chr1	noncoding	exon	14362	14829	0.000000	-	.	gene_id "WASH7P"; transcript_id "NR_024540"; exon_number "1"
chr1	noncoding	exon	14970	15038	0.000000	-	.	gene_id "WASH7P"; transcript_id "NR_024540"; exon_number "2"
chr1	noncoding	exon	15796	15947	0.000000	-	.	gene_id "WASH7P"; transcript_id "NR_024540"; exon_number "3"
chr1	Cufflinks	exon	16035	16310	0	-	.	gene_id "NONHSAG000002"; transcript_id "NONHSAT000005"; FPKM "0"; exon_number 5;
chr1	noncoding	exon	16607	16765	0.000000	-	.	gene_id "WASH7P"; transcript_id "NR_024540"; exon_number "4"
chr1	noncoding	exon	16858	17055	0.000000	-	.	gene_id "WASH7P"; transcript_id "NR_024540"; exon_number "5"

==> genome_gtf <==
chr15	refGene	exon	68346572	68346688	.	+	.	gene_id "PIAS1"; transcript_id "NM_016166"; exon_number "1"; exon_id "NM_016166.1"; gene_name "PIAS1";
chr15	refGene	CDS	68346665	68346688	.	+	0	gene_id "PIAS1"; transcript_id "NM_016166"; exon_number "1"; exon_id "NM_016166.1"; gene_name "PIAS1";
chr15	refGene	exon	68378644	68379088	.	+	.	gene_id "PIAS1"; transcript_id "NM_016166"; exon_number "2"; exon_id "NM_016166.2"; gene_name "PIAS1";
chr15	refGene	CDS	68378644	68379088	.	+	0	gene_id "PIAS1"; transcript_id "NM_016166"; exon_number "2"; exon_id "NM_016166.2"; gene_name "PIAS1";
chr15	refGene	exon	68434284	68434368	.	+	.	gene_id "PIAS1"; transcript_id "NM_016166"; exon_number "3"; exon_id "NM_016166.3"; gene_name "PIAS1";
chr15	refGene	CDS	68434284	68434368	.	+	2	gene_id "PIAS1"; transcript_id "NM_016166"; exon_number "3"; exon_id "NM_016166.3"; gene_name "PIAS1";
chr15	refGene	exon	68434628	68434675	.	+	.	gene_id "PIAS1"; transcript_id "NM_016166"; exon_number "4"; exon_id "NM_016166.4"; gene_name "PIAS1";
chr15	refGene	CDS	68434628	68434675	.	+	1	gene_id "PIAS1"; transcript_id "NM_016166"; exon_number "4"; exon_id "NM_016166.4"; gene_name "PIAS1";
chr15	refGene	exon	68438154	68438244	.	+	.	gene_id "PIAS1"; transcript_id "NM_016166"; exon_number "5"; exon_id "NM_016166.5"; gene_name "PIAS1";
chr15	refGene	CDS	68438154	68438244	.	+	1	gene_id "PIAS1"; transcript_id "NM_016166"; exon_number "5"; exon_id "NM_016166.5"; gene_name "PIAS1";

==> refGene_hg19 <==
A3GALT2	NM_001080438	chr1	-	33772366	33786699	33772366	33786699	5	33772366,33777652,33778101,33778407,33786676,	33773054,33777790,33778191,33778491,33786699,
AADACL3	NM_001103169	chr1	+	12776117	12788726	12776343	12785963	3	12776117,12780884,12785188,	12776347,12780948,12788726,
AADACL3	NM_001103170	chr1	+	12776117	12788726	12776179	12785963	4	12776117,12779476,12780884,12785188,	12776347,12779693,12780948,12788726,
AADACL3	NR_111984	chr1	+	12776117	12788726	12788726	12788726	3	12776117,12780884,12785188,	12776347,12780948,12788726,
AADACL4	NM_001013630	chr1	+	12704565	12727097	12704565	12726746	4	12704565,12711141,12721801,12725971,	12704733,12711358,12721865,12727097,
ABCA4	NM_000350	chr1	-	94458393	94586705	94458792	94586601	50	94458393,94461664,94463416,94466391,94466557,94467413,94470996,94473189,94473790,94474306,94476355,94476817,94480098,94481294,94485137,94486795,94487195,94487401,94488941,94490509,94495000,94495983,94496551,94497333,94502295,94502700,94505598,94506764,94508316,94508891,94510168,94512474,94514423,94517188,94520666,94522156,94526092,94528132,94528667,94543245,94544145,94544877,94546033,94548907,94564349,94568570,94574132,94576993,94578528,94586535,	94458798,94461751,94463666,94466484,94466661,94467548,94471138,94473296,94473853,94474427,94476485,94476941,94480246,94481410,94485315,94486965,94487270,94487507,94488974,94490604,94495187,94496082,94496676,94497599,94502344,94502906,94505683,94506958,94508454,94509031,94510300,94512649,94514513,94517254,94520871,94522378,94526315,94528309,94528873,94543443,94544262,94545017,94546274,94548997,94564547,94568698,94574272,94577135,94578622,94586705,
ABCB10	NM_012089	chr1	-	229652328	229694442	229653925	229694399	13	229652328,229654587,229657338,229661682,229662975,229665945,229667382,229675202,229676352,229677983,229683245,229684980,229693882,	229654157,229654622,229657382,229661863,229663055,229666155,229667478,229675338,229676499,229678118,229683448,229685181,229694442,
ABCD3	NM_001122674	chr1	+	94883932	94944260	94884034	94944106	9	94883932,94924162,94930330,94933474,94939321,94940698,94941169,94943814,94944079,	94884144,94924199,94930429,94933563,94939391,94940796,94941293,94943871,94944260,
ABCD3	NM_002858	chr1	+	94883932	94984219	94884034	94982685	23	94883932,94924162,94930330,94933474,94939321,94940698,94941169,94943814,94946019,94948725,94953097,94953249,94953447,94955280,94955458,94956739,94964157,94964338,94964500,94965050,94972093,94980701,94982607,	94884144,94924199,94930429,94933563,94939391,94940796,94941293,94943871,94946162,94948795,94953167,94953347,94953539,94955372,94955531,94956803,94964235,94964404,94964590,94965170,94972198,94980758,94984219,
ABL2	NM_001136000	chr1	-	179068461	179112224	179076852	179112179	13	179068461,179078342,179079416,179081443,179084012,179086466,179087721,179089324,179090729,179095511,179100445,179102446,179112067,	179078033,179078576,179079590,179081533,179084165,179086651,179087899,179089409,179091002,179095807,179100616,179102509,179112224,

==> intragenic_bed <==
chr1	11874	12227	Intragenic
chr1	12228	12612	Intragenic
chr1	12613	12721	Intragenic
chr1	12722	13220	Intragenic
chr1	13221	14829	Intragenic
chr1	14830	14969	Intragenic
chr1	14970	15038	Intragenic
chr1	15039	15795	Intragenic
chr1	15796	15947	Intragenic
chr1	15948	16606	Intragenic

==> rmsk_bed <==
chr1	10001	10468	1504	1.3	0.4	1.3	(249240153)	+	(CCCTAA)n	Simple_repeat	1	463	(0)	1
chr1	10469	11447	3612	11.4	27.0	1.3	(249239174)	C	TAR1	Satellite/telo	(399)	1712	483	2
chr1	11504	11675	437	23.5	18.6	3.5	(249238946)	C	L1MC	LINE/L1	(2236)	5646	5449	3
chr1	11678	11780	239	29.4	1.9	1.0	(249238841)	C	MER5B	DNA/hAT-Charlie	(74)	104	1	4
chr1	15265	15355	318	23.0	3.8	0.0	(249235266)	C	MIR3	SINE/MIR	(119)	143	49	5
chr1	16713	16749	203	16.2	0.0	0.0	(249233872)	+	(TGG)n	Simple_repeat	1	37	(0)	6
chr1	18907	19048	239	33.8	14.8	0.0	(249231573)	+	L2a	LINE/L2	2942	3104	(322)	7
chr1	19948	20405	652	34.6	8.5	4.2	(249230216)	+	L3	LINE/CR1	3042	3519	(970)	8
chr1	20531	20679	270	33.1	0.7	2.7	(249229942)	+	Plat_L3	LINE/CR1	2802	2947	(639)	9
chr1	21949	22075	254	27.9	4.7	3.9	(249228546)	+	MLT1K	LTR/ERVL-MaLR	15	142	(453)	10

```

Now running the script below:

```python

from __future__ import division
import module01_mapping_from_cln    as m01
import module02_RNA_Quantification  as m02
import module04_Stat_Sample         as m04

dir_name = DirSystem()
sftw_name= UsedSoftware()

samp_process = m01.Map_From_cln( samp_info,bt_index_base,genome_gtf,dir_name['dir'],sftw_name )
samp_process.load_samp()
samp_process.run_QC()
samp_process.run_tophat()

samp_mRNAQ = m02.RnaQuantification( samp_info,bt_index_base,genome_gtf,genome_gtf_wl,intragenic_bed, rmsk_bed,  dir_name['dir'],sftw_name )
samp_mRNAQ.load_samp()
samp_mRNAQ.RNA_QuantPipe()

samp_Stat   = m04.SampStat( samp_info,ercc_info, genome_gtf ,dir_name['dir'] )
samp_Stat.load_samp()
samp_Stat.Basic_Stat()
samp_Stat.ERCC_Stat()
samp_Stat.plot_regression()
samp_Stat.Repeat_Stat()
```

