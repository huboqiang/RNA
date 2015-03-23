from __future__ import division
import re,sys,os
import cPickle as pickle
import numpy as np
import matplotlib
import subprocess,time
matplotlib.use('Agg')
from matplotlib import pyplot as plt

import module_Matrix as m_mat

class CufflinksInfo(object):
   def __init__(self,root_dir,l_samp,out_dir):
      self.l_sample = l_samp
      self.l_file   = [ "%s/%s/isoforms.fpkm_tracking" % ( root_dir,samp ) for samp in l_samp ]
      self.out_dir  = out_dir
      self.gen_sam_FPKM = {}
      self.rep_group    = {}
   
   def load_samp(self):
      if not os.path.isdir( self.out_dir ):
         os.mkdir( self.out_dir )
      for i,samp in enumerate(self.l_sample):
         self.__read_cufflinks(i)
      self.__merge_FPKM_cufflinks()
   
   def element_group_subgroup_FPKM_sum(self):
      infile = "%s/merge.Repeat.FPKM.xls" % ( self.out_dir )
      
      m_matrix = m_mat.Matrix_info( infile, 6, "float" )
      m_matrix.load_mat()
      
      outfile_element  = "%s/merge.Repeat.SumFPKM.element.xls"  % ( self.out_dir )
      outfile_subgroup = "%s/merge.Repeat.SumFPKM.subgroup.xls" % ( self.out_dir )
      outfile_group    = "%s/merge.Repeat.SumFPKM.group.xls"    % ( self.out_dir )
      
      f_outfile_element  = open( outfile_element ,"w" )
      f_outfile_subgroup = open( outfile_subgroup,"w" )
      f_outfile_group    = open( outfile_group   ,"w" )
      
      print >>f_outfile_element ,"Element\t%s"    % ( "\t".join( m_matrix.rowname ) )
      print >>f_outfile_subgroup,"SubGroup\t%s"   % ( "\t".join( m_matrix.rowname ) )
      print >>f_outfile_group   ,"Group\t%s"      % ( "\t".join( m_matrix.rowname ) )
      
      l_repInfo  = m_matrix.colname
      np_element  = np.array( [ "%s" % ( inf.split('\t')[3] ) for inf in l_repInfo ],dtype="string" )
      np_subgroup = np.array( [ "%s" % ( inf.split('\t')[4] ) for inf in l_repInfo ],dtype="string" )
      np_group    = np.array( [ "%s" % ( inf.split('\t')[5] ) for inf in l_repInfo ],dtype="string" )
      
      for elem in sorted( set(np_element) ):
         idx = (np_element ==elem)
         sub_mat      = m_matrix.matrix[ idx,: ]
         sum_FPKM     = np.sum( sub_mat,axis=0 )
         str_sum_FPKM = np.array( sum_FPKM,dtype="string" )
         print >>f_outfile_element, "%s\t%s" % ( elem, "\t".join(str_sum_FPKM) )
      f_outfile_element.close()

      for subgroup in sorted( set(np_subgroup) ):
         idx = (np_subgroup ==subgroup)
         sub_mat      = m_matrix.matrix[ idx,: ]
         sum_FPKM     = np.sum( sub_mat,axis=0 )
         str_sum_FPKM = np.array( sum_FPKM,dtype="string" )
         print >>f_outfile_subgroup, "%s\t%s" % ( subgroup, "\t".join(str_sum_FPKM) )
      f_outfile_subgroup.close()


      for group in sorted( set(np_group) ):
         idx = (np_group ==group)
         sub_mat      = m_matrix.matrix[ idx,: ]
         sum_FPKM     = np.sum( sub_mat,axis=0 )
         str_sum_FPKM = np.array( sum_FPKM,dtype="string" )
         print >>f_outfile_group, "%s\t%s" % ( group, "\t".join(str_sum_FPKM) )
      f_outfile_group.close()
      
         
         
      
      
   def __read_cufflinks(self,i):
      infile = self.l_file[i]
      sam    = self.l_sample[i]
      print >>sys.stderr, "Reading %s" % (sam)
      f_infile = open( infile,"r" )
      f_infile.readline()
      for line in f_infile:
         line = line.strip('\n')
         f    = line.split()
         gene = f[0]
         fpkm = float(f[9])
         stat = f[12]
         group= f[3]
         
         if stat != "OK" or fpkm < 1:
            continue
         
         self.rep_group[ gene ] = group.split("_")[0]
         
         if gene not in self.gen_sam_FPKM:
            self.gen_sam_FPKM[ gene ] = {}
         self.gen_sam_FPKM[ gene ][ sam ] = fpkm
      f_infile.close()
      
   def __merge_FPKM_cufflinks(self):
      out_file_rep    = "%s/merge.Repeat.FPKM.xls" % ( self.out_dir )
      f_out_file_rep  = open( out_file_rep,"w" )
      out_info    = "Chrom\tBeg\tEnd\tElement\tSubGroup\tGroup\t%s" % ( "\t".join( self.l_sample ) )
      print >>f_out_file_rep, out_info
      
      pat = "(\w+)__(\w+)__(\w+):(\d+)-(\d+)"
      pattern = re.compile( pat )
      for repeat in sorted(self.gen_sam_FPKM):
         match   = pattern.search( repeat )
         if match:
            element  =      match.group(1)
            subgroup =      match.group(2)
            chrom    =      match.group(3)
            beg      = int( match.group(4) )
            end      = int( match.group(5) )
            group    = self.rep_group[ repeat ]
         
         out = "chr%s\t%d\t%d\t%s\t%s\t%s" % (chrom,beg,end,element, subgroup, group)
         for sam in self.l_sample:
            value = 0.0
            if sam in self.gen_sam_FPKM[repeat]:
               value = self.gen_sam_FPKM[repeat][sam]
            out += "\t%f" % ( value )
         print >>f_out_file_rep, out
      f_out_file_rep.close()
   

         
            
            
      
      