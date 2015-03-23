#-*- coding:utf-8 -*-
from __future__ import division
import re,sys,os
import subprocess
import numpy   as np
import cPickle as pickle
import scipy   as sp
import time
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from scipy import stats 
import module_running_jobs as m_jobs
import module_GTFFeature   as m_gtf
import module_CountInfo    as m_cnt

class QuantPipe(dict):
   def __init__( self,M_samp_info,l_samp,genome_file,anno_file_refERCC,anno_file,intragenic_bed,rmsk_gtf,rmsk_bed,  dir_name,sftw_name ):
      self['samp_info'] = M_samp_info
      self['samp']      = l_samp
      self['infile'] = { 'genome_file':genome_file,'anno_file_refERCC':anno_file_refERCC,'anno_file':anno_file,'intragenic_bed':intragenic_bed,'rmsk_gtf':rmsk_gtf,'rmsk_bed':rmsk_bed }
      self['dir_name']  = dir_name
      self['sftw_name'] = sftw_name
      self.__load_name()
      
   def __load_name(self):
      self.tophat       = self['dir_name']['tophat_dir']
      self.HTS          = self['dir_name']['HTSeq_result_dir']
      self.HTS_k        = self['dir_name']['HTSeq_known_dir']
      self.HTS_u        = self['dir_name']['HTSeq_unknown_dir']
      self.cufflink_u   = self['dir_name']['cufflinks_unknown_dir']
      self.cuffquant    = self['dir_name']['cuffquant_dir']
      self.cuffnorm     = self['dir_name']['cuffnorm_dir']
      self.cuffquant_ercc=self['dir_name']['cuffquant_ERCC_dir']
      self.cuffnorm_ercc =self['dir_name']['cuffnorm_ERCC_dir']
      
      self.cuffquant_k     = self['dir_name']['cuffquant_k_dir']
      self.cuffnorm_k      = self['dir_name']['cuffnorm_k_dir']
      self.cuffquant_ercc_k= self['dir_name']['cuffquant_ERCC_k_dir']
      self.cuffnorm_ercc_k = self['dir_name']['cuffnorm_ERCC_k_dir']
      
      self.repeat_mrg_dir = self['dir_name']['repeat_mrg_dir']
      self.repeatCount  = self['dir_name']['repeat_counts_dir']
      
      
      home_dir          = os.path.abspath('./')
      self.script_dir   = "%s/scripts" % (home_dir)
      self.bin_dir      = "%s/bin"     % (home_dir)
      self.data_dir     = "%s/Database"% (home_dir)

   def run_HTSeq_known(self):
      
      sh_file        =  "%s/s04.HTSeq_known.sh"       %  (self.script_dir)
      sh_work_file   =  "%s/s04.HTSeq_known_work.sh"  %  (self.script_dir)
      
      py_exe         = self['sftw_name'].py
      samtools_exe   = self['sftw_name'].samtools
      deseq_exe      = self['sftw_name'].deseq
      
      sh_info = """
py_exe=$1
samtools_exe=$2
deseq_exe=$3
tophat_dir=$4
samp_name=$5
HTS_k_dir=$6
known_GTF=$7

$samtools_exe view  -H              $tophat_dir/$samp_name/accepted_hits.bam  > $tophat_dir/$samp_name/accepted_hits.header.sam
$samtools_exe sort  -n -m 200000000 $tophat_dir/$samp_name/accepted_hits.bam    $tophat_dir/$samp_name/accepted_hits.sort_name
$samtools_exe view  -o    $tophat_dir/$samp_name/accepted_hits.sort_name.sam    $tophat_dir/$samp_name/accepted_hits.sort_name.bam 

[ ! -d $HTS_k_dir/$samp_name ] && mkdir -p $HTS_k_dir/$samp_name

$py_exe $deseq_exe                                                                                                                             \\
   -s no -f sam -a 10                                                                                                                        \\
   -o $tophat_dir/$samp_name/accepted_hits.sort_name.gene.sam                                                                                \\
   $tophat_dir/$samp_name/accepted_hits.sort_name.sam  $known_GTF            >$HTS_k_dir/$samp_name/$samp_name.dexseq.txt                 && \\
grep -v -P '^ERCC-|^RGC-|MIR|SNORD|Mir|Snord' $HTS_k_dir/$samp_name/$samp_name.dexseq.txt > $HTS_k_dir/$samp_name/$samp_name.dexseq_clean.txt			&& \\
grep    -P '^ERCC-|^RGC-'                     $HTS_k_dir/$samp_name/$samp_name.dexseq.txt > $HTS_k_dir/$samp_name/$samp_name.dexseq_ERCC_RGCPloyA.txt	&& \\
grep "__no_feature" $tophat_dir/$samp_name/accepted_hits.sort_name.gene.sam | grep -v chrM |                                                 \\
	cat           $tophat_dir/$samp_name/accepted_hits.header.sam /dev/stdin | 		                                                         \\
	$samtools_exe view -Sb /dev/stdin >$tophat_dir/$samp_name/accepted_hits.genome.bam                                                          && \\
$samtools_exe sort  -m 200000000 $tophat_dir/$samp_name/accepted_hits.genome.bam    $tophat_dir/$samp_name/accepted_hits.genome.sort          
rm             $tophat_dir/$samp_name/accepted_hits.sort_name.sam $tophat_dir/$samp_name/accepted_hits.sort_name.bam   $tophat_dir/$samp_name/accepted_hits.sort_name.gene.sam  $tophat_dir/$samp_name/accepted_hits.genome.bam
      """
      sh_work = ""
      for samp in self['samp']:
         tophat_dir  =  self.tophat
         samp_name   =  self['samp_info']['samp_brief'][samp]
         known_GTF   =  self['infile']['anno_file']

         sh_work += "sh %s  %s %s %s  %s %s %s %s\n" % ( sh_file,  py_exe,samtools_exe,deseq_exe,  self.tophat, samp_name, self.HTS_k, known_GTF)
      
      my_job = m_jobs.running_jobs(sh_file,sh_work_file)
      my_job.load_sh_file(      sh_info )
      my_job.load_sh_work_file( sh_work )
      my_job.running_multi( cpu=8 )
#      my_job.running_SGE( vf="1g",maxjob=100 )

   
   ### 
   def run_cufflinks_u(self):

      sh_file      = "%s/s05.cufflinks_GenomeMapped.sh"      % (self.script_dir)
      sh_work_file = "%s/s05.cufflinks_GenomeMapped_work.sh" % (self.script_dir)
      
      cflk_dir     = self['sftw_name'].cflk_dir
      
      sh_info = """
cflk_dir=$1
in_bam=$2
gtf_file=$3
out_dir=$4

$cflk_dir/cufflinks      \\
   -p 8  -u              \\
   -o $out_dir           \\
   $in_bam
      """
      sh_work = ""
      for samp in self['samp']:
         brief_name = self['samp_info']['samp_brief'][samp]
         in_bam      = "%s/%s/accepted_hits.genome.sort.bam"% ( self.tophat,    brief_name   ) 
         out_dir     = "%s/%s"                              % ( self.cufflink_u,brief_name   )
         sh_work += "sh %s  %s  %s %s %s \n" % ( sh_file,  cflk_dir, in_bam, self['infile']['anno_file'], out_dir)
      
      my_job = m_jobs.running_jobs(sh_file,sh_work_file)
      my_job.load_sh_file(      sh_info )
      my_job.load_sh_work_file( sh_work )
      my_job.running_multi( cpu=8 )
#      my_job.running_SGE( vf="500m",maxjob=100 )

   def run_cuffcomp_novo_trans(self):

      sh_file      = "%s/s06.1.cuffcompare_novo.sh"      % (self.script_dir)
      sh_work_file = "%s/s06.1.cuffcompare_novo_work.sh" % (self.script_dir)
      
      cflk_dir     = self['sftw_name'].cflk_dir
      
      sh_info = """
cflk_dir=$1
out_prefix=$2
shift
shift

$cflk_dir/cuffcompare    \\
   -o  $out_prefix       \\
   -T  $@                \\
      """
      sh_work = ""
      out_prefix  = "%s/novo_lnc_raw"           % ( self.data_dir )
      l_in_samp   = [  "%s/%s/transcripts.gtf"  % ( self.cufflink_u,self['samp_info']['samp_brief'][samp] ) for samp in self['samp'] ]
      sh_work = "sh %s  %s %s %s" % ( sh_file, cflk_dir, out_prefix, " ".join(l_in_samp) )
      
      my_job = m_jobs.running_jobs(sh_file,sh_work_file)
      my_job.load_sh_file(      sh_info )
      my_job.load_sh_work_file( sh_work )
      my_job.running_multi( cpu=8 )
#      my_job.running_SGE( vf="500m",maxjob=100 )


   def run_HTSeq_unknown(self):
      
      sh_file        =  "%s/s07.HTSeq_unknown.sh"       %  (self.script_dir)
      sh_work_file   =  "%s/s07.HTSeq_unknown_work.sh"  %  (self.script_dir)
      
      py_exe         = self['sftw_name'].py
      samtools_exe   = self['sftw_name'].samtools
      deseq_exe      = self['sftw_name'].deseq
      
      if not os.path.isdir( self.HTS_u ):
         os.mkdir( self.HTS_u )
      
      sh_info = """
py_exe=$1
samtools_exe=$2
deseq_exe=$3
tophat_dir=$4
samp_name=$5
HTS_u_dir=$6
unknown_GTF=$7

$samtools_exe view  -H              $tophat_dir/$samp_name/accepted_hits.genome.sort.bam          > $tophat_dir/$samp_name/accepted_hits.header.sam
$samtools_exe sort  -n -m 200000000 $tophat_dir/$samp_name/accepted_hits.genome.sort.bam            $tophat_dir/$samp_name/accepted_hits.genome.sort_name
$samtools_exe view  -o              $tophat_dir/$samp_name/accepted_hits.genome.sort_name.sam  $tophat_dir/$samp_name/accepted_hits.genome.sort_name.bam 

[ ! -d $HTS_u_dir/$samp_name ] && mkdir -p $HTS_u_dir/$samp_name

$py_exe $deseq_exe                                                                                                                           \\
   -s no -f sam -a 10                                                                                                                        \\
   $tophat_dir/$samp_name/accepted_hits.genome.sort_name.sam  $unknown_GTF            >$HTS_u_dir/$samp_name/$samp_name.dexseq_NeoRaw.txt
      """
      sh_work = ""
      for samp in self['samp']:
         tophat_dir  =  self.tophat
         samp_name   =  self['samp_info']['samp_brief'][samp]
         unknown_GTF =  "%s/novo_lnc_raw.combined.gtf" % ( self.data_dir )
         sh_work += "sh %s  %s %s %s  %s %s %s %s\n" % ( sh_file,  py_exe,samtools_exe,deseq_exe,  self.tophat, samp_name, self.HTS_u, unknown_GTF)
      
      my_job = m_jobs.running_jobs(sh_file,sh_work_file)
      my_job.load_sh_file(      sh_info )
      my_job.load_sh_work_file( sh_work )
      my_job.running_multi( cpu=8 )
#      my_job.running_SGE( vf="1g",maxjob=100 )
   
   def RPKM_novo_trans(self):
      l_brief_samp    = [ "%s" % ( self['samp_info']['samp_brief'][samp] ) for samp in self['samp'] ]
      unknown_GTF     = "%s/novo_lnc_raw.combined.gtf" % ( self.data_dir )
      
      Gtf_Info = m_gtf.GTFFeature( unknown_GTF  )
      Gtf_Info.gene_intergenic( self['infile']['intragenic_bed'] )
      Cnt_Info = m_cnt.CountInfo(  self.HTS_u, l_brief_samp, "dexseq_NeoRaw", self.HTS  )

      Cnt_Info.generate_mat()
      Cnt_Info.load_mat()
      Cnt_Info.cal_RPKM( Gtf_Info.gene,self.tophat )
      
      rpkm_file = "%s/merge.%s.RPKM.xls" % ( self.HTS,"dexseq_NeoRaw" )
      
      Gtf_Info.load_gene_RPKM( rpkm_file )
      Gtf_Info.output_GTF()
      Gtf_Info.get_gene_info()
      
   def merge_novo_known_GTF(self):
      sh_file        =  "%s/Generate_Transcriptome.sh"       %  (self.script_dir)
      sh_work_file   =  "%s/Generate_Transcriptome_work.sh"  %  (self.script_dir)
      
      sh_info = """
known_GTF=$1
unknown_GTF=$2
merge_GTF=$3
merge_ERCC_GTF=$4

sed 's/XLOC_/novoXLOC_/g' $unknown_GTF | sed 's/TCONS_/novoTCONS_/g' >$unknown_GTF.tmp
grep -P "^chr" $known_GTF | cat /dev/stdin $unknown_GTF.tmp | bedtools sort -i /dev/stdin | grep -P "^chr" >$merge_GTF
rm $unknown_GTF.tmp
grep -P "^ERCC|^RGC" $known_GTF | cat $merge_GTF /dev/stdin >$merge_ERCC_GTF
      """
      
      known_GTF     = self['infile']['anno_file']
      unknown_GTF   = "%s/novo_lnc_raw.combined.FPKM0.5_rep0.25.multiExon.gtf" % ( self.data_dir )
      merge_GTF     = "%s/all.exon.sort.gtf"         % ( self.data_dir )
      merge_ERCC_GTF= "%s/all.exon.sort.ERCC.gtf"    % ( self.data_dir )
      self['infile']['anno_file_merge']      = merge_GTF
      self['infile']['anno_file_merge_ERCC'] = merge_ERCC_GTF
      
      sh_work = "sh %s   %s %s %s %s" % (  sh_file,  known_GTF , unknown_GTF, merge_GTF, merge_ERCC_GTF )
      my_job = m_jobs.running_jobs(sh_file,sh_work_file)
      my_job.load_sh_file(      sh_info )
      my_job.load_sh_work_file( sh_work )
      my_job.running_multi( cpu=1 )
#      my_job.running_SGE( vf="1g",maxjob=100 )
      
      
      
      
      
      
   def run_cuffquant(self):
      sh_file      = "%s/s08.cuffquant.sh"      % (self.script_dir)
      sh_work_file = "%s/s08.cuffquant_work.sh" % (self.script_dir)
      
      if not os.path.isdir( self.cuffquant ):
         os.mkdir( self.cuffquant )
      
      cflk_dir     = self['sftw_name'].cflk_dir
      
      sh_info = """
cflk_dir=$1
in_bam=$2
gtf_file=$3
out_dir=$4

$cflk_dir/cuffquant      \\
   -p 8  -u              \\
   -o $out_dir           \\
   $gtf_file             \\
   $in_bam
      """
      sh_work = ""
      for samp in self['samp']:
         brief_name = self['samp_info']['samp_brief'][samp]
         in_bam      = "%s/%s/accepted_hits.bam"   % ( self.tophat,    brief_name   ) 
         out_dir     = "%s/%s"                     % ( self.cuffquant ,brief_name   )
         sh_work += "sh %s  %s %s %s %s \n" % ( sh_file,  cflk_dir, in_bam, self['infile']['anno_file_merge'], out_dir)
      
      my_job = m_jobs.running_jobs(sh_file,sh_work_file)
      my_job.load_sh_file(      sh_info )
      my_job.load_sh_work_file( sh_work )
      my_job.running_multi( cpu=8 )
#      my_job.running_SGE( vf="1000m",maxjob=100 )      


   def run_cuffnorm(self):
      sh_file      = "%s/s09.cuffnorm.sh"      % (self.script_dir)
      sh_work_file = "%s/s09.cuffnorm_work.sh" % (self.script_dir)
      
      cflk_dir     = self['sftw_name'].cflk_dir
      
      if not os.path.isdir( self.cuffnorm ):
         os.mkdir( self.cuffnorm )

      l_brief = []
      l_cxb   = []
      for samp in self['samp']:
         brief_name = self['samp_info']['samp_brief'][samp]
         l_brief.append( brief_name  )
         l_cxb.append(   "%s/%s/abundances.cxb" % (self.cuffquant,brief_name) )
      
      
      list_brief = ",".join( l_brief )
      list_cxb   = " ".join( l_cxb   )
      
      sh_info = """
cflk_dir=$1

$cflk_dir/cuffnorm        \\
   -p 8  -o %s  -L %s     \\
   %s                     \\
   %s
      """ % ( self.cuffnorm, list_brief, self['infile']['anno_file_merge'], list_cxb )
      
      sh_work = "sh %s %s" % (sh_file,cflk_dir)
      my_job = m_jobs.running_jobs(sh_file,sh_work_file)
      my_job.load_sh_file(      sh_info )
      my_job.load_sh_work_file( sh_work )
      my_job.running_multi( cpu=8 )
#      my_job.running_SGE( vf="1000m",maxjob=100 )   
      






   def run_cuffquant_ERCC(self):
      sh_file      = "%s/s08.1.cuffquant.ERCC.sh"      % (self.script_dir)
      sh_work_file = "%s/s08.1.cuffquant.ERCC_work.sh" % (self.script_dir)
      
      cflk_dir     = self['sftw_name'].cflk_dir
      
      if not os.path.isdir( self.cuffquant_ercc ):
         os.mkdir( self.cuffquant_ercc )
      
      sh_info = """
cflk_dir=$1
in_bam=$2
gtf_file=$3
out_dir=$4

$cflk_dir/cuffquant      \\
   -p 8  -u              \\
   -o $out_dir           \\
   $gtf_file             \\
   $in_bam
      """
      sh_work = ""
      for samp in self['samp']:
         brief_name = self['samp_info']['samp_brief'][samp]
         in_bam      = "%s/%s/accepted_hits.bam"   % ( self.tophat,         brief_name   ) 
         out_dir     = "%s/%s"                     % ( self.cuffquant_ercc ,brief_name   )
         sh_work += "sh %s   %s %s %s %s \n" % ( sh_file,  cflk_dir,in_bam, self['infile']['anno_file_merge_ERCC'], out_dir)
      
      my_job = m_jobs.running_jobs(sh_file,sh_work_file)
      my_job.load_sh_file(      sh_info )
      my_job.load_sh_work_file( sh_work )
      my_job.running_multi( cpu=8 )
#      my_job.running_SGE( vf="1000m",maxjob=100 )      


   def run_cuffnorm_ERCC(self):
      sh_file      = "%s/s09.1.cuffnorm.ERCC.sh"      % (self.script_dir)
      sh_work_file = "%s/s09.1.cuffnorm.ERCC_work.sh" % (self.script_dir)
      if not os.path.isdir( self.cuffnorm_ercc ):
         os.mkdir( self.cuffnorm_ercc )
      
      cflk_dir     = self['sftw_name'].cflk_dir

      l_brief = []
      l_cxb   = []
      for samp in self['samp']:
         brief_name = self['samp_info']['samp_brief'][samp]
         l_brief.append( brief_name  )
         l_cxb.append(   "%s/%s/abundances.cxb" % (self.cuffquant_ercc,brief_name) )
      
      
      list_brief = ",".join( l_brief )
      list_cxb   = " ".join( l_cxb   )
      
      sh_info = """
cflk_dir=$1

$cflk_dir/cuffnorm        \\
   -p 8  -o %s  -L %s     \\
   %s                     \\
   %s
      """ % ( self.cuffnorm_ercc, list_brief, self['infile']['anno_file_merge_ERCC'], list_cxb )
      
      sh_work = "sh %s  %s" % (sh_file, cflk_dir)
      my_job = m_jobs.running_jobs(sh_file,sh_work_file)
      my_job.load_sh_file(      sh_info )
      my_job.load_sh_work_file( sh_work )
      my_job.running_multi( cpu=8 )
#      my_job.running_SGE( vf="1000m",maxjob=100 )   


   def makeGTF_withoutERCC(self):
      sh_file        =  "%s/RemoveERCC.sh"       %  (self.script_dir)
      sh_work_file   =  "%s/RemoveERCC_work.sh"  %  (self.script_dir)
      
      sh_info = """
known_GTF=$1
remove_ERCC_GTF=$2

grep -P "^chr" $known_GTF >$remove_ERCC_GTF
      """
      
      known_GTF      = self['infile']['anno_file']
      remove_ERCC_GTF= "%s.sort.gtf"    % ( ".".join( self['infile']['anno_file'].split(".")[:-3] )  )
      self['infile']['anno_file_remove_ERCC'] = remove_ERCC_GTF
      
      sh_work = "sh %s   %s %s" % (  sh_file,  known_GTF , remove_ERCC_GTF )
      my_job = m_jobs.running_jobs(sh_file,sh_work_file)
      my_job.load_sh_file(      sh_info )
      my_job.load_sh_work_file( sh_work )
      my_job.running_multi( cpu=1 )



   def run_cuffquant_k(self):
      sh_file      = "%s/s08.cuffquant_k.sh"      % (self.script_dir)
      sh_work_file = "%s/s08.cuffquant_k_work.sh" % (self.script_dir)
      
      if not os.path.isdir( self.cuffquant_k ):
         os.mkdir( self.cuffquant_k )
      
      cflk_dir     = self['sftw_name'].cflk_dir
      
      sh_info = """
cflk_dir=$1
in_bam=$2
gtf_file=$3
out_dir=$4

$cflk_dir/cuffquant      \\
   -p 8  -u              \\
   -o $out_dir           \\
   $gtf_file             \\
   $in_bam
      """
      sh_work = ""
      for samp in self['samp']:
         brief_name = self['samp_info']['samp_brief'][samp]
         in_bam      = "%s/%s/accepted_hits.bam"   % ( self.tophat,      brief_name   ) 
         out_dir     = "%s/%s"                     % ( self.cuffquant_k, brief_name   )
         sh_work += "sh %s  %s %s %s %s \n" % ( sh_file,  cflk_dir, in_bam, self['infile']['anno_file_remove_ERCC'], out_dir)
      
      my_job = m_jobs.running_jobs(sh_file,sh_work_file)
      my_job.load_sh_file(      sh_info )
      my_job.load_sh_work_file( sh_work )
      my_job.running_multi( cpu=8 )
#      my_job.running_SGE( vf="1000m",maxjob=100 )      


   def run_cuffnorm_k(self):
      sh_file      = "%s/s09.cuffnorm_k.sh"      % (self.script_dir)
      sh_work_file = "%s/s09.cuffnorm_k_work.sh" % (self.script_dir)
      
      cflk_dir     = self['sftw_name'].cflk_dir
      
      if not os.path.isdir( self.cuffnorm_k ):
         os.mkdir( self.cuffnorm_k )

      l_brief = []
      l_cxb   = []
      for samp in self['samp']:
         brief_name = self['samp_info']['samp_brief'][samp]
         l_brief.append( brief_name  )
         l_cxb.append(   "%s/%s/abundances.cxb" % (self.cuffquant_k,brief_name) )
      
      
      list_brief = ",".join( l_brief )
      list_cxb   = " ".join( l_cxb   )
      
      sh_info = """
cflk_dir=$1

$cflk_dir/cuffnorm        \\
   -p 8  -o %s  -L %s     \\
   %s                     \\
   %s
      """ % ( self.cuffnorm_k, list_brief, self['infile']['anno_file_remove_ERCC'], list_cxb )
      
      sh_work = "sh %s %s" % (sh_file,cflk_dir)
      my_job = m_jobs.running_jobs(sh_file,sh_work_file)
      my_job.load_sh_file(      sh_info )
      my_job.load_sh_work_file( sh_work )
      my_job.running_multi( cpu=8 )
#      my_job.running_SGE( vf="1000m",maxjob=100 )   
      






   def run_cuffquant_ERCC_k(self):
      sh_file      = "%s/s08.1.cuffquant_k.ERCC.sh"      % (self.script_dir)
      sh_work_file = "%s/s08.1.cuffquant_k.ERCC_work.sh" % (self.script_dir)
      
      cflk_dir     = self['sftw_name'].cflk_dir
      
      if not os.path.isdir( self.cuffquant_ercc_k ):
         os.mkdir( self.cuffquant_ercc_k )
      
      sh_info = """
cflk_dir=$1
in_bam=$2
gtf_file=$3
out_dir=$4

$cflk_dir/cuffquant      \\
   -p 8  -u              \\
   -o $out_dir           \\
   $gtf_file             \\
   $in_bam
      """
      sh_work = ""
      for samp in self['samp']:
         brief_name = self['samp_info']['samp_brief'][samp]
         in_bam      = "%s/%s/accepted_hits.bam"   % ( self.tophat,           brief_name   ) 
         out_dir     = "%s/%s"                     % ( self.cuffquant_ercc_k ,brief_name   )
         sh_work += "sh %s   %s %s %s %s \n" % ( sh_file,  cflk_dir,in_bam, self['infile']['anno_file'], out_dir)
      
      my_job = m_jobs.running_jobs(sh_file,sh_work_file)
      my_job.load_sh_file(      sh_info )
      my_job.load_sh_work_file( sh_work )
      my_job.running_multi( cpu=8 )
#      my_job.running_SGE( vf="1000m",maxjob=100 )      


   def run_cuffnorm_ERCC_k(self):
      sh_file      = "%s/s09.1.cuffnorm.ERCC_k.sh"      % (self.script_dir)
      sh_work_file = "%s/s09.1.cuffnorm.ERCC_k_work.sh" % (self.script_dir)
      if not os.path.isdir( self.cuffnorm_ercc_k ):
         os.mkdir( self.cuffnorm_ercc_k )
      
      cflk_dir     = self['sftw_name'].cflk_dir

      l_brief = []
      l_cxb   = []
      for samp in self['samp']:
         brief_name = self['samp_info']['samp_brief'][samp]
         l_brief.append( brief_name  )
         l_cxb.append(   "%s/%s/abundances.cxb" % (self.cuffquant_ercc_k,brief_name) )
      
      
      list_brief = ",".join( l_brief )
      list_cxb   = " ".join( l_cxb   )
      
      sh_info = """
cflk_dir=$1

$cflk_dir/cuffnorm        \\
   -p 8  -o %s  -L %s     \\
   %s                     \\
   %s
      """ % ( self.cuffnorm_ercc_k, list_brief, self['infile']['anno_file'], list_cxb )
      
      sh_work = "sh %s  %s" % (sh_file, cflk_dir)
      my_job = m_jobs.running_jobs(sh_file,sh_work_file)
      my_job.load_sh_file(      sh_info )
      my_job.load_sh_work_file( sh_work )
      my_job.running_multi( cpu=8 )
#      my_job.running_SGE( vf="1000m",maxjob=100 )   







   def run_repeat_count(self):

      sh_file      = "%s/s13.Repeat_Count.sh"      % (self.script_dir)
      sh_work_file = "%s/s13.Repeat_Count_work.sh" % (self.script_dir)
      py_Repeat_Intersect2Count = "%s/Repeat_Intersect2Count.py" % (self.bin_dir)
      
      samtools_exe   = self['sftw_name'].samtools
      bedtools_exe   = self['sftw_name'].bedtools
      py_exe         = self['sftw_name'].py
      
      if not os.path.isdir( self.repeatCount ):
         os.mkdir( self.repeatCount )
      
      sh_info = """
samtools_exe=$1
bedtools_exe=$2
py_exe=$3
in_bam=$4
gtf_bed=$5
py_Repeat_Intersect2Count=$6
out_dir=$7

$samtools_exe view -F 0x0004 $in_bam |                     \\
   grep -v ERCC-00* | grep -v RGC-CRE|                \\
   grep -v RGC-GFP  | grep -v RGC-mRFP |grep NH:i:1 | \\
   awk '{OFS="\\t"; print $3,$4,$4+length($10),$1 }' >${out_dir}/repeat_result.bed

$bedtools_exe intersect -sorted -loj -a $gtf_bed -b ${out_dir}/repeat_result.bed | \\
   $py_exe $py_Repeat_Intersect2Count /dev/stdin >${out_dir}/repeat_count.bed
      """
      sh_work = ""
      for samp in self['samp']:
         brief_name = self['samp_info']['samp_brief'][samp]
         in_bam      = "%s/%s/accepted_hits.bam"         % ( self.tophat        ,brief_name   ) 
         out_dir     = "%s/%s"                           % ( self.repeatCount   ,brief_name   )
         sh_work += "sh %s  %s %s %s  %s %s %s %s\n" % ( sh_file,  samtools_exe,bedtools_exe,py_exe,  in_bam, self['infile']['rmsk_bed'], py_Repeat_Intersect2Count, out_dir)
         
         if not os.path.isdir( out_dir ):
            os.mkdir( out_dir )
      
      my_job = m_jobs.running_jobs(sh_file,sh_work_file)
      my_job.load_sh_file(      sh_info )
      my_job.load_sh_work_file( sh_work )
      my_job.running_multi( cpu=8 )
#      my_job.running_SGE( vf="500m",maxjob=100 )
      