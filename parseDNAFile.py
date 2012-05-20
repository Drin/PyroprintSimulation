#!/usr/bin/python -tt
# Copyright 2010 Google Inc.
# Licensed under the Apache License, Version 2.0
# http://www.apache.org/licenses/LICENSE-2.0


"""
  So far we have done:
  Find all files with DNA info in DIR
  Parse all files and find all DNA Chunks
  in each chunk find the primer sequence (if it exists)
  store all the 104char segments after the primer
  find all unique sequences of 104char segments
  find all possible choose 7 combinations of unique sequences
  generate pyroprints (python)
  compare pyroprints
  CUDA pearson correlation
  Graph pyroprints
  """

import sys
import os
import itertools
import numpy as np
import math
import re 
from optparse import OptionParser
import ConfigParser
from extract_allele import extractAlleles

# detect gpu support
gpu_support = False
try:
    import biogpu.correlation
    gpu_support = True
except:
   pass


###PRIMER
#16-23 primer GGAACCTGCGGTTGGATCAC
#23-5 primer CGTGAGGCTTAACCTT

##Dispensation Sequence:
#23-5 AACACGCGA23(GATC)GAA
#16-23 CCTCTACTAGAGCG20(TCGA)TT

def main():
   #Detect gpu support
   if gpu_support:
      print "PyCUDA detected. GPU support enabled."
   else:
      print "PyCuda not detected. Exiting"
      sys.exit

   (dataDir, disp, primer, numIsolates) = handleArgs()

   alleles = extractAlleles(dataDir, disp, primer, numIsolates)
   for allele in alleles:
      print allele
   return

   print 'Storing alleles in constant memory'


   print 'Calculating Pearson Correlation'
   buckets = biogpu.correlation.pearson(ranges)
   print('Results:')
   for i in range(len(buckets)):
      print('\t[%d] (%.1f%%, %.1f%%) = %d' % (i, ranges[i][0] * 100.0, ranges[i][1] * 100.0, buckets[i]))
   print('\n')

'''
This method is parses command-line arguments. It returns a directory containing
DNA sequences (or other directories) a dispensation sequence to use on given
DNA sequences, a forward primer (which may need to be applied to the reverse
complement strand) and an upper bound on the number of isolates to generate in
silico.
'''
def handleArgs():
   '''
   parser = OptionParser(usage="Extracts alleles from input DNA sequences\n")
   parser.add_option("-p", "--path", dest="path", default="Genome Sequences/",
         help="Path containing DNA sequences")
   parser.add_option("--disp", dest="dispSeq", default="AACACGCGA23(GATC)GAA",
         help="Dispensation Sequence")
   parser.add_option("-n", dest="num", type="int", default=-1,
         help="Number of Isolates to generate")
   parser.add_option("-f", dest="file", default="config.cfg",
         help="File containing parameters")
   parser.add_option("--fp", dest="forPrimer", default="AACCTT",
         help="Forward Primer")
         '''

   parser = OptionParser(usage="Extracts alleles\n")
   parser.add_option("-p", "--path", dest="dir", default="Genome Sequences/rDNA plasmid sequences/23-5/", help="Path to Genome Sequence Folders")
   parser.add_option("-d", dest="DispSeq", default="AACACGCGA23(GATC)GAA", help="Dispensation order")
   parser.add_option("-m", "--max", dest="max", type="int", default=-1, help="Max number of combinations to generate")
   parser.add_option("-f", "--file", dest="file", help="File containing parameters", default="config.cfg")
   parser.add_option("--primer", dest="primer", default="AACCTT", help="Primer to use")

   (options, args) = parser.parse_args()

   if options.file:
      config = ConfigParser.RawConfigParser()
      config.read(options.file)
      dataPath = config.get("params", "path")

      if config.has_option("params", "max"):
         numIsolates = config.getint("params", "max")
      else:
         numIsolates = -1

      forwardPrimer = config.get("params", "primer")
      dispSeq = config.get("params", "DispSeq")

   else:
      #Use command line args
      dataPath = options.dir
      numIsolates = options.max
      dispSeq = options.DispSeq
      forwardPrimer = options.primer

   return (dataPath, dispSeq, forwardPrimer, numIsolates)


# This is the standard boilerplate that calls the main() function.
if __name__ == '__main__':
   main()
