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

class Map_From_cln(dict):
   
   def __init__(self,samp_info, genome_file,anno_file, dir_name,sftw_name ):
      self['sam_info'] = {}
      self['sample'] = []
      self['infile'] = {'info_file':samp_info, 'genome_file':genome_file,'anno_file':anno_file }
      self['stage'] = { 'name':[] }
      self['dir_name']  = dir_name
      self['sftw_name'] = sftw_name
      
   def load_samp(self):
      
      self['sam_info']['samp_brief'] = {}
      self['sam_info']['type']     = {}
      self['sam_info']['stage']    = {}
      self['sam_info']['dilute']   = {}
      self['sam_info']['stage_sam']= {}
      self['sam_info']['RFP_mols'] = {}
      self['sam_info']['GFP_mols'] = {}
      self['sam_info']['CRE_mols'] = {}
      self['sam_info']['stage_sam']= {}
      self['sam_info']['data_type']= {}
      
      info_file = self['infile']['info_file']
      file = open(info_file,"r")
      in_h = file.readline()
      for line in file:
         line = line.strip('\n')
         f = line.split()
         samp       = f[0]
         brief_name = f[1]
         ltype      = f[2]
         stage      = f[3]
         ERCC_dilute= float( f[4] )
         RFP_mols   = float( f[5] )
         GFP_mols   = float( f[6] )
         CRE_mols   = float( f[7] )
         data_type  = f[8]  # PE or SE

         self['sample'].append( samp )
         self['sam_info']['samp_brief'][ samp ] = brief_name
         self['sam_info']['type'][  samp ]      = ltype
         self['sam_info']['stage'][ samp ]      = stage
         self['sam_info']['dilute'][samp ]      = ERCC_dilute
         self['sam_info']['RFP_mols'][ samp ]   = RFP_mols
         self['sam_info']['GFP_mols'][ samp ]   = GFP_mols
         self['sam_info']['CRE_mols'][ samp ]   = CRE_mols
         self['sam_info']['data_type'][samp ]   = data_type
         
         if stage not in self['stage']['name']:
            self['stage']['name'].append( stage )
         if stage not in self['sam_info']['stage_sam']:
            self['sam_info']['stage_sam'][stage] = []
            
         self['sam_info']['stage_sam'][stage].append( samp )
      file.close()
 
   def run_QC(self):
      home_dir    = os.path.abspath('./')
      raw_dir     = self['dir_name']['raw_data']
      cln_dir     = self['dir_name']['clean_data']
      
      if not os.path.isdir( cln_dir ):
         os.mkdir( cln_dir )
      
      script_dir   = "%s/scripts"         % (home_dir)
      bin_dir      = "%s/bin"             % (home_dir)
      
      pl_exe      = self['sftw_name'].pl
      pl_QC       = "%s/bin/QC.pl" % ( home_dir )
      
      sh_file       = "%s/scripts/QC.sh"      % (home_dir)
      sh_work_file  = "%s/scripts/QC_work.sh" % (home_dir)
      
      sh_info = """
pl_exe=$1
pl_QC=$2
in_dir=$3
out_dir=$4
samp=$5
data_type=$6

$pl_exe $pl_QC --indir $in_dir --outdir $out_dir --sample $samp --end $data_type
      """
      
      sh_work = ""
      for samp in self['sample']:
         if not os.path.isdir( "%s/%s" % (cln_dir,samp) ):
            os.mkdir(          "%s/%s" % (cln_dir,samp) )
         in_dir    = raw_dir
         out_dir   = cln_dir
         data_type = 2
         if self['sam_info']['data_type'][samp ] == "SE":
            data_type = 1
         sh_work  += " sh %s  %s %s %s %s %s %d\n" % ( sh_file, pl_exe,  pl_QC, in_dir,out_dir,samp,data_type )
      
      my_job = m_jobs.running_jobs(sh_file,sh_work_file)
      my_job.load_sh_file(      sh_info )
      my_job.load_sh_work_file( sh_work )
      my_job.running_multi( cpu=8 )
#      my_job.running_SGE( vf="400m",maxjob=100 )
      
      
   
   def run_tophat(self):
      home_dir     = os.path.abspath('./')
      
      cln_dir      = self['dir_name']['clean_data']
      tophat_dir   = self['dir_name']['tophat_dir']
      
      if not os.path.isdir(  tophat_dir):
         os.mkdir( tophat_dir )
      
      script_dir   = "%s/scripts"         % (home_dir)
      bin_dir      = "%s/bin"             % (home_dir)
      
      tophat_py    = self['sftw_name'].tophat
      
      sh_file      = "%s/s02.tophat.sh"      % (script_dir)
      sh_work_file = "%s/s02.tophat_work.sh" % (script_dir)
      
      sh_info = """
tophat_py=$1
cln_dir=$2
samp_name=$3
brief_name=$4
tophat_dir=$5
genome=$6
gtf_file=$7
PE2=$8

$tophat_py  \\
   -p 8 -G $gtf_file                                                 \\
   --library-type fr-unstranded                                      \\
   --transcriptome-index /datc/huboqiang/cir_dyj_V2/Database/refseqGene.ERCC_RGCPloyA.exon.sort \\
   -o $tophat_dir/$brief_name                                        \\
   $genome                                                           \\
   $cln_dir/$samp_name/1.cln.fq.gz  $PE2
      """ 
      sh_work = ""
      for samp in self['sample']:
         brief_name = self['sam_info']['samp_brief'][samp]
         PE2 = ""
         if self['sam_info']['data_type'][samp ] == "PE":
            PE2 = "%s/%s/2.cln.fq.gz" % (cln_dir,samp)
         sh_work += "sh %s  %s %s %s %s  %s %s %s %s\n" % ( sh_file, tophat_py, cln_dir, samp, brief_name, tophat_dir, self['infile']['genome_file'],self['infile']['anno_file'],PE2  )
      
      my_job = m_jobs.running_jobs(sh_file,sh_work_file)
      my_job.load_sh_file(      sh_info )
      my_job.load_sh_work_file( sh_work )
      my_job.running_multi( cpu=6 )
#      my_job.running_SGE( vf="7g",maxjob=100 )

      
