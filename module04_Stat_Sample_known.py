from __future__ import division
import re,sys,os
import module_StatInfo as Stat
import scipy
import numpy as np
from scipy import optimize
import matplotlib.pyplot as plt
import module_CountInfo  as m_cnt
import module_GTFFeature   as m_gtf
import module_Matrix       as m_mat
import module_running_jobs as m_jobs
import module_RepeatCount   as m_repcnt

class SampStat(dict):
   def __init__(self,samp_info,ercc_info,genome_gtf,dir_name):
      self['sam_info'] = {}
      self['sample'] = []
      self['infile'] = {'info_file':samp_info,'ercc_info':ercc_info,'anno_file':genome_gtf}
      self['stage'] = { 'name':[] }
      self['dir_name'] = dir_name
      self.ERCC_info   = { 'len':{}, 'mol':{} }
      self.__load_name()
      self.__load_ERCC_info()
      
   def __load_ERCC_info(self):
      f_file = open( self['infile']['ercc_info'],"r")
      f_file.readline()
      for line in f_file:
         line = line.strip('\n')
         f    = line.split()
         ERCC_id = f[0]
         self.ERCC_info['len'][ERCC_id] = int(f[1])
         self.ERCC_info['mol'][ERCC_id] = float(f[2])
      f_file.close()
      
   def __load_name(self):
      self.cln          = self['dir_name']['clean_data']
      self.tophat       = self['dir_name']['tophat_dir']
      self.HTS          = self['dir_name']['HTSeq_result_dir']
      self.HTS_k        = self['dir_name']['HTSeq_known_dir']
      self.HTS_u        = self['dir_name']['HTSeq_unknown_dir']
      self.cufflink_u   = self['dir_name']['cufflinks_unknown_dir']
      self.cuffquant    = self['dir_name']['cuffquant_dir']
      self.cuffnorm     = self['dir_name']['cuffnorm_dir']
      self.cuffquant_ercc =self['dir_name']['cuffquant_ERCC_dir']
      self.cuffnorm_ercc  =self['dir_name']['cuffnorm_ERCC_dir']
      self.cuffnorm_ercc_k=self['dir_name']['cuffnorm_ERCC_k_dir']
      self.repeatCount    = self['dir_name']['repeat_counts_dir']
      self.repeat_mrg     = self['dir_name']['repeat_mrg_dir']
      
      home_dir          = os.path.abspath('./')
      self.script_dir   = "%s/scripts" % (home_dir)
      self.bin_dir      = "%s/bin"     % (home_dir)
      self.data_dir     = "%s/Database"% (home_dir)
      self.statInfo     = "%s/StatInfo"% (home_dir)
      if not os.path.isdir( self.statInfo ):
         os.mkdir( self.statInfo )

   def load_samp(self):      
      self['sam_info']['samp_brief'] = {}
      self['sam_info']['brief_samp'] = {}
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
         self['sam_info']['brief_samp'][ brief_name ] = samp
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
   
   def __get_HTS_reads(self, HTS_info, samp):
      print samp
      idx = self['sample'].index( samp )
      return HTS_info.sam_tot_reads[idx]
      
   def __get_HTS_clean_split(self):
      sh_file        =  "%s/p.HTSeq_split.sh"       %  (self.script_dir)
      sh_work_file   =  "%s/p.HTSeq_split_work.sh"  %  (self.script_dir)
      sh_info = """
infile=$1
out_Refseq=$2
out_NONCODE=$3
out_NSMB=$4

grep -v -P '^NONHSAG|XLOC_' $infile >$out_Refseq
head -n 1 $infile >$out_NONCODE && grep -P '^NONHSAG' $infile >>$out_NONCODE
head -n 1 $infile >$out_NSMB    && grep -P '^XLOC'    $infile >>$out_NSMB

      """ 
      
      infile      = "%s/merge.dexseq_clean.gene.xls"          % ( self.HTS )
      out_Refseq  = "%s/merge.dexseq_clean_refseq.gene.xls"   % ( self.HTS )
      out_NONCODE = "%s/merge.dexseq_clean_NONCODE.gene.xls"  % ( self.HTS )
      out_NSMB    = "%s/merge.dexseq_clean_NSMB.gene.xls"     % ( self.HTS )
      sh_work = "sh %s  %s %s %s %s " % ( sh_file,infile,out_Refseq,out_NONCODE,out_NSMB )
      my_job = m_jobs.running_jobs(sh_file,sh_work_file)
      my_job.load_sh_file(      sh_info )
      my_job.load_sh_work_file( sh_work )
      my_job.running_multi( cpu=1 )
      
   def Basic_Stat(self):
      """
      Stat for QC, Tophat mapping, ERCC RGC count
      """
      out_file   = "%s/01.BasicInfo_QC_map_SpikeIn.xls" % (self.statInfo)
      f_out_file = open( out_file,"w" )
      out_info   = "Sample\tBrief_samp\tRaw_Reads\tClean_Reads\t"    +\
                   "Pre_Map_Reads\tAligned_Reads\tHTSseq_Known_Reads\t"    +\
                   "HTSeq_Refseq_Reads\tHTSeq_NONCODE_V4_Reads\tHTSeq_Nsmb_Reads\t" +\
                   "RFP_Reads\tGFP_Reads\tCRE_Reads\tERCC_Reads\t"   +\
                   "RFP_Mols\tGFP_Mols\tCRE_Mols\tERCC_Mols"
                   
      print >>f_out_file, out_info

      l_breif_samp = [ self['sam_info']['samp_brief'][samp] for samp in self['sample'] ]

      '''
         Load refseq reads
      '''
      HTS_info = m_cnt.CountInfo( self.HTS_k,l_breif_samp,"dexseq_clean",self.HTS )
      HTS_info.load_mat()
      HTS_info.sam_tot_reads()
      
      self.__get_HTS_clean_split()
      
      Refseq_info = m_cnt.CountInfo( self.HTS_k,l_breif_samp,"dexseq_clean_refseq",self.HTS )
      Refseq_info.load_mat()
      Refseq_info.sam_tot_reads()

      NONCODE_info = m_cnt.CountInfo( self.HTS_k,l_breif_samp,"dexseq_clean_NONCODE",self.HTS )
      NONCODE_info.load_mat()
      NONCODE_info.sam_tot_reads()
      
      NSMB_info = m_cnt.CountInfo( self.HTS_k,l_breif_samp,"dexseq_clean_NSMB",self.HTS )
      NSMB_info.load_mat()
      NSMB_info.sam_tot_reads()      
      
      '''
         Load other information
      '''
      for idx,samp in enumerate(self['sample']):
         brief_name  = self['sam_info']['samp_brief'][samp]
         QC_log      =  "%s/%s/log"                         % ( self.cln   ,samp )
         Tophat_log  =  "%s/%s/align_summary.txt"           % ( self.tophat,brief_name )
         HTSeq_SpikeIn= "%s/%s/%s.dexseq_ERCC_RGCPloyA.txt" % ( self.HTS_k ,brief_name, brief_name )
         
         QcStat_info = Stat.QcStat( QC_log )
         MapStat_info= Stat.TophatStat( Tophat_log )
         SpikeIn_info= Stat.SpikeIn( HTSeq_SpikeIn, self['infile']['ercc_info'] )
         
      
         QcStat_info.read_infile()
         MapStat_info.read_infile()
         SpikeIn_info.load_HTS_file()
         
         pre_map_read = MapStat_info['statInfo']['totalRead']
         aligned_read = MapStat_info['statInfo']['mappedRead']
         
         if self['sam_info']['data_type'][samp ] == "PE":         
            HTSseq_read  = self.__get_HTS_reads( HTS_info    ,samp ) * 2
            Refseq_read  = self.__get_HTS_reads( Refseq_info ,samp ) * 2
            NONCODE_read = self.__get_HTS_reads( NONCODE_info,samp ) * 2
            NSMB_read    = self.__get_HTS_reads( NSMB_info   ,samp ) * 2
            read_RFP   =  SpikeIn_info.RGC_count['RGC-mRFP'] * 2
            read_GFP   =  SpikeIn_info.RGC_count['RGC-GFP' ] * 2
            read_CRE   =  SpikeIn_info.RGC_count['RGC-CRE' ] * 2
            read_ERCC  =  SpikeIn_info.ERCC_total            * 2
         else:
            HTSseq_read  = self.__get_HTS_reads( HTS_info    ,samp )
            Refseq_read  = self.__get_HTS_reads( Refseq_info ,samp )
            NONCODE_read = self.__get_HTS_reads( NONCODE_info,samp )
            NSMB_read    = self.__get_HTS_reads( NSMB_info   ,samp )         
            read_RFP   =  SpikeIn_info.RGC_count['RGC-mRFP']
            read_GFP   =  SpikeIn_info.RGC_count['RGC-GFP' ]
            read_CRE   =  SpikeIn_info.RGC_count['RGC-CRE' ]
            read_ERCC  =  SpikeIn_info.ERCC_total
            
         
         mol_RFP  = self['sam_info']['RFP_mols'][ samp ]
         mol_GFP  = self['sam_info']['GFP_mols'][ samp ]
         mol_CRE  = self['sam_info']['CRE_mols'][ samp ]
         mol_ERCC = self['sam_info']['dilute'][   samp ] * 6.023*10**10
         
         out_info =  "%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%1.2e\t%1.2e\t%1.2e\t%1.2e"   \
            % ( samp, brief_name, QcStat_info.raw_reads, QcStat_info.cln_reads,  \
                pre_map_read, aligned_read, HTSseq_read,                         \
                Refseq_read,NONCODE_read,NSMB_read,                              \
                read_RFP, read_GFP,read_CRE, read_ERCC,                          \
                mol_RFP , mol_GFP ,mol_CRE , mol_ERCC )
         print  >>f_out_file, out_info
      f_out_file.close()
   
   def ERCC_Stat(self):
      """
      Using ERCC mols to estimate the total amount of refGene RNA-mols in a cell. 
      (ERCC_MOLs) / (ERCC_FPKM) = (mRNA_MOLs) / (mRNA_FPKM)
      """
      self.l_ERCC_name = []
      self.l_RGCs_name = []
      self.l_mRNA_name = []
      self.l_ERCC_FPKM = {}
      self.l_RGCs_FPKM = {}
      self.l_mRNA_FPKM = {}
      self.l_cirRNA_FPKM={}

      self.l_ERCC_HTSname = []
      self.l_RGCs_HTSname = []
      self.l_mRNA_HTSname = []
      self.l_ERCC_RPKM = {}
      self.l_RGCs_RPKM = {}
      self.l_mRNA_RPKM = {}


      self.l_ERCC_MOLs = {}
      self.l_RGCs_MOLs = {}
      self.l_mRNA_MOLs = {}
      self.l_cirRNA_MOLs={}
      self.l_mRNA_MOLs_HTSname = {}

      
      self.regression = {}
      
      self.__load_FPKM()
      self.__load_MOLs()      # ERCC RGC mols
      self.__get_mRNA_MOLs()  # get mRNA mols using ERCC_FPKM, ERCC_MOLs and mRNA_FPKM
      self.__load_Count()
      
      out_file   = "%s/02.ERCC_Mols.xls" % (self.statInfo)
      f_out_file = open( out_file,"w" )
      
      out_info = "Sample\tBrief_samp\tERCC_MOLs\tRGC_MOLs\tmRNA_MOLs\tRefSeq_mRNA_MOLs\tRegression_R\tRegression_P\tRefSeq_mRNA_FPKM>0.1"
      print >>f_out_file, out_info
      for samp in self['sample']:
         brief_name = self['sam_info']['samp_brief'][samp]
         ERCC_MOLs  = sum(    self.l_ERCC_MOLs[  brief_name ] )
         RGC_MOLs   = sum(    self.l_RGCs_MOLs[  brief_name ] )
         mRNA_MOLs  = np.sum( self.l_mRNA_MOLs[  brief_name ] )
         RefSeq_mRNA_MOLs           = np.sum(  self.l_mRNA_MOLs[ brief_name ][ self.mRNA_refSeq_index ] )
         RefSeq_mRNA_lFPKM          = np.array(self.l_mRNA_FPKM[ brief_name ],dtype=float )
         RefSeq_mRNA_lFPKM          = RefSeq_mRNA_lFPKM[ self.mRNA_refSeq_index ] 
         RefSeq_mRNA_Exps           = np.shape(  RefSeq_mRNA_lFPKM[RefSeq_mRNA_lFPKM >= 0.1] )[0]
         regression_R               = self.regression[brief_name]['r_value']
         regression_P               = self.regression[brief_name]['p_value']
         
         out_info = "%s\t%s\t%6.3e\t%6.3e\t%6.3e\t%6.3e\t%1.3f\t%1.3e\t%d" % ( samp,brief_name, ERCC_MOLs,RGC_MOLs,mRNA_MOLs, RefSeq_mRNA_MOLs,regression_R,regression_P,RefSeq_mRNA_Exps )
         print >>f_out_file, out_info
      
      f_out_file.close()
      
   
   def Repeat_Stat(self):

      l_brief_samp = [ self['sam_info']['samp_brief'][samp] for samp in self['sample'] ]
     
      l_align_reads = []
      for brief_samp in l_brief_samp:
         Tophat_log  =  "%s/%s/align_summary.txt"           % ( self.tophat,brief_samp )
         MapStat_info= Stat.TophatStat( Tophat_log )
         MapStat_info.read_infile()
         l_align_reads.append( MapStat_info['statInfo']['mappedRead'] )
         
      np_align_reads = np.array( l_align_reads )
     
      RepCnt = m_repcnt.RepeatCount( self.repeatCount,l_brief_samp, self.repeat_mrg )
      RepCnt.generate_mat()
      RepCnt.element_group_subgroup_FPKM_sum(np_align_reads)

   
   def __load_Count(self):  # for ERCC count Stat
      l_breif_samp = [ self['sam_info']['samp_brief'][samp] for samp in self['sample'] ]
      
      genome_gtf = self['infile']['anno_file']
      Gtf_Info = m_gtf.GTFFeature( genome_gtf  )
      
      ERCC_info = m_cnt.CountInfo( self.HTS_k,l_breif_samp,"dexseq_ERCC_RGCPloyA",self.HTS )
      ERCC_info.load_mat()
      ERCC_info.cal_RPKM( Gtf_Info.gene,self.tophat )
      
      
      
            
   def __load_FPKM(self):
      Cuffnorm_file = "%s/genes.fpkm_table" % ( self.cuffnorm_ercc_k )
      f_Cuffnorm    = open( Cuffnorm_file,"r" )
      in_head        = f_Cuffnorm.readline()
      f_head         = in_head.split()
      l_samp_brief   = [ sam.split("_0")[0] for sam in f_head[1:] ]
      
      for brief_name in l_samp_brief:
         self.l_ERCC_FPKM[brief_name] = []
         self.l_RGCs_FPKM[brief_name] = []
         self.l_mRNA_FPKM[brief_name] = []
      
      for line in f_Cuffnorm:
         line = line.strip('\n')
         f    = line.split()
         gene = f[0]
         
         if gene[0:3] == "MIR":
            continue
         if gene[0:5] == "SNORD":
            continue
         
         if gene[0:5] == "ERCC-":
            self.l_ERCC_name.append( gene )
            for i,brief_name in enumerate( l_samp_brief ):
               self.l_ERCC_FPKM[brief_name].append( f[i+1] )
               
         elif gene[0:4] == "RGC-":
            self.l_RGCs_name.append( gene )
            for i,brief_name in enumerate( l_samp_brief ):
               self.l_RGCs_FPKM[brief_name].append( f[i+1] )
         
         else:
            self.l_mRNA_name.append( gene )
            for i,brief_name in enumerate( l_samp_brief ):
               self.l_mRNA_FPKM[brief_name].append( f[i+1] )
      
      f_Cuffnorm.close()
      
      self.__mRNA_refSeq_index()   
   
   def __mRNA_refSeq_index(self):
      self.mRNA_refSeq_index = []
      for i,gene in enumerate(self.l_mRNA_name):
         if gene[0:7] == "NONHSAG" or gene[0:5] == "XLOC_" or gene[0:9] == "novoXLOC_":
            continue
         self.mRNA_refSeq_index.append( i )
      self.mRNA_refSeq_index = np.array( self.mRNA_refSeq_index,dtype=int )
   
   def __load_MOLs(self):
      for samp in self['sample']:

         brief_name = self['sam_info']['samp_brief'][samp]
         if brief_name not in self.l_ERCC_MOLs:
            self.l_ERCC_MOLs[brief_name] = []
         if brief_name not in self.l_RGCs_MOLs:
            self.l_RGCs_MOLs[brief_name] = []
            
         for ERCC_id in self.l_ERCC_name:
            self.l_ERCC_MOLs[brief_name].append( self.ERCC_info['mol'][ERCC_id]*(6.02*10**23/10**18)*self['sam_info']['dilute'][samp] )
         
         self.l_RGCs_MOLs[ brief_name ].append( self['sam_info']['CRE_mols'][ samp ] )
         self.l_RGCs_MOLs[ brief_name ].append( self['sam_info']['GFP_mols'][ samp ] )
         self.l_RGCs_MOLs[ brief_name ].append( self['sam_info']['RFP_mols'][ samp ] )
         
   def __get_mRNA_MOLs(self):
      self.regression = {}
      for samp in self['sample']:
         brief_name = self['sam_info']['samp_brief'][samp]
         
         np_ERCC_FPKM = np.array( self.l_ERCC_FPKM[brief_name],dtype=float )
         np_ERCC_MOLs = np.array( self.l_ERCC_MOLs[brief_name],dtype=float )
         
         np_ERCC_FPKM = np.log2( np_ERCC_FPKM+0.1 )
         np_ERCC_MOLs = np.log2( np_ERCC_MOLs )

         mol_idx = ( np_ERCC_MOLs - np.log2(6.02*10**23/10**18) > -18 ) * (np_ERCC_FPKM>0)

         np_ERCC_FPKM = np_ERCC_FPKM[ mol_idx ]
         np_ERCC_MOLs = np_ERCC_MOLs[ mol_idx ]
         
         slope,intercept,r_value,p_value,slope_std_error =  scipy.stats.linregress( np_ERCC_MOLs,np_ERCC_FPKM )
         
         self.regression[brief_name] = {'slope'    :slope,           \
                                        'inter'    :intercept,       \
                                        'r_value'  :r_value,         \
                                        'p_value'  :p_value,         \
                                        'std_err'  :slope_std_error  \
         }
         
         np_mRNA_FPKM = np.array( self.l_mRNA_FPKM[brief_name],dtype=float )
         np_mRNA_FPKM = np.log2( np_mRNA_FPKM+0.1 )
         np_mRNA_MOLs = ( np_mRNA_FPKM - intercept ) / slope
         np_mRNA_MOLs = np.power( np_mRNA_MOLs, 2 )
         np_mRNA_MOLs[ np_mRNA_MOLs<0.0            ] = 0.0
         np_mRNA_MOLs[ np_mRNA_FPKM<=np.log2(0.1) ] = 0.0
         
         self.l_mRNA_MOLs[ brief_name ] = np_mRNA_MOLs   
   
   
   def plot_regression(self):
      pdfname = "All_samples.rpkm_cnt.pdf" 
      fig = plt.figure(figsize=(150,150))
      cnt = 0
      for samp in self['sample']:
         brief_name = self['sam_info']['samp_brief'][samp]
         
         slope          =  self.regression[brief_name]['slope']
         intercept      =  self.regression[brief_name]['inter']
         r_value        =  self.regression[brief_name]['r_value']
         p_value        =  self.regression[brief_name]['p_value']
         slope_std_error=  self.regression[brief_name]['std_err']
         
         np_ERCC_FPKM = np.array( self.l_ERCC_FPKM[brief_name],dtype=float )
         np_ERCC_MOLs = np.array( self.l_ERCC_MOLs[brief_name],dtype=float )/(6.02*10**23/10**18)
      
         np_ERCC_FPKM = np.log2( np_ERCC_FPKM+0.1 )
         np_ERCC_MOLs = np.log2( np_ERCC_MOLs )

         mol_idx = ( np_ERCC_MOLs > -18 ) * (np_ERCC_FPKM>0)

         np_ERCC_FPKM_use = np_ERCC_FPKM[ mol_idx ]
         np_ERCC_MOLs_use = np_ERCC_MOLs[ mol_idx ]
         
         MOLs_predict = np.linspace( -18,-5,100 )
         FPKM_predict = intercept + slope * (MOLs_predict + np.log2(6.02*10**23/10**18) )
      
         cnt += 1
         ax = plt.subplot( 20,20,cnt)
         ax.plot( np_ERCC_MOLs,     np_ERCC_FPKM      ,".r")
         ax.plot( np_ERCC_MOLs_use, np_ERCC_FPKM_use  ,"ro")
         ax.plot( MOLs_predict,     FPKM_predict           )
         title = brief_name + ' ERCC-mols vs RPKM in samples'
         ax.set_title(title,fontsize=12)
         ax.set_xlabel('Spike-in mols per reaction(attlmole,log2)')
         ax.set_ylabel('Mean read-coverage of spike-in molecules(RPKM,log2)')
         ax.text(-18,13.5,r'$ y = %f + %fx $' % (intercept,slope))
         ax.text(-18,13.0,r'$ Pearson R = %f $' % (r_value))
         ax.text(-18,12.5,r'$ p = %6.2e $' % (p_value))
         ax.text(-18,12.0,r'$ mRNA genes = %6.2e $' % (  np.sum( self.l_mRNA_MOLs[ brief_name ][ self.mRNA_refSeq_index ] )  ) )
         for tick in ax.xaxis.get_major_ticks():
            tick.tick1On = True
            tick.tick2On = False
         for tick in ax.yaxis.get_major_ticks():
            tick.tick1On = True
            tick.tick2On = False
         ax.set_xlim(-20,-5)
         ax.set_ylim(-4,15)
         
      plt.savefig(pdfname,format='pdf')
   
