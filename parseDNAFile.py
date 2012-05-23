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
import numpy
import math
import re 
from extract_allele import extractAlleles
from optparse import OptionParser
import ConfigParser
import pycuda.autoinit
import pycuda.compiler
import pycuda.gpuarray
import pycuda.driver

DEBUG = True
TESTING = True

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
   if DEBUG:
      print("Configuring simulation...\n")

   (dataDir, disp, primer, numIsolates, numRegions) = handleArgs()

   if DEBUG:
      print("Extracting allele information...\n")

   (alleles, num_alleles, num_pyro_peaks) = extractAlleles(dataDir, disp, primer, numIsolates)

   if DEBUG:
      print(("Running simulation for {0} alleles and {1} length " +
            "pyroprints...\n").format(num_alleles, num_pyro_peaks))

   if DEBUG:
      print("Preparing buckets for pearson correlation slices...\n")

   ranges = [(0.000, 1.000)]
   bucketSlices = generateRanges(0.9900, 1, 0.0001)
   ranges.extend([bucketSlice for bucketSlice in bucketSlices])

   if DEBUG:
      print('Preparing pyroprinted alleles for device constant memory...\n')

   alleles_c = numpy.zeros(shape=(num_alleles, num_pyro_peaks), dtype=numpy.uint8, order='C')
   for alleleNdx in range(num_alleles):
      numpy.put(alleles_c[alleleNdx], range(num_pyro_peaks), alleles[alleleNdx])

   if DEBUG:
      print("Loading CUDA kernel source code...\n")

   kernelFile = open(os.path.join(os.getcwd(), 'biogpu/pearson.cu'), 'r')
   kernel = pycuda.compiler.SourceModule(kernelFile.read())
   kernelFile.close()

   (const_ptr, size) = kernel.get_global("alleles")

   if DEBUG:
      print(("Transferring {0} bytes of allele pyroprints into device " +
            "constant memory...\n ").format(size))

   pycuda.driver.memcpy_htod(const_ptr, alleles_c)

   if DEBUG:
      print('Calculating pearson correlation for all pairwise combinations ' +
            'for {0} generated isolates...\n'.format(calcCombinations(num_alleles, numRegions)))

   if TESTING:
      buckets = biogpu.correlation.pearson(kernel, ranges, calcCombinations(10, numRegions), num_pyro_peaks)
   else:
      buckets = biogpu.correlation.pearson(kernel, ranges, calcCombinations(num_alleles, numRegions), num_pyro_peaks)

   print('Results:\n')
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
   parser = OptionParser(usage="Extracts alleles\n")
   parser.add_option("-p", "--path", dest="dir", default="Genome Sequences/", help="Path to Genome Sequence Folders")
   parser.add_option("-d", dest="DispSeq", default="AACACGCGA23(GATC)GAA", help="Dispensation order")
   parser.add_option("-m", "--max", dest="max", type="int", default=-1, help="Max number of combinations to generate")
   parser.add_option("-f", "--file", dest="file", help="File containing parameters", default="config.cfg")
   parser.add_option("--primer", dest="primer", default="TTGGATCAC", help="Primer to use")
   parser.add_option("--regionNum", dest="regionNum", type="int", default=7, help="Number of ITS Regions to simulate")

   (options, args) = parser.parse_args()

   if options.file:
      config = ConfigParser.RawConfigParser()
      config.read(options.file)
      dataPath = config.get("params", "path")

      if config.has_option("params", "max"):
         numIsolates = config.getint("params", "max")
      else:
         numIsolates = -1

      if config.has_option("params", "regionNum"):
         numRegions = config.getint("params", "regionNum")
      else:
         numRegions = 7

      forwardPrimer = config.get("params", "primer")
      dispSeq = config.get("params", "DispSeq")

   else:
      #Use command line args
      dataPath = options.dir
      numIsolates = options.max
      dispSeq = options.DispSeq
      forwardPrimer = options.primer
      numRegions = options.regionNum

   if TESTING:
      numIsolates = 3

   return (dataPath, dispSeq, forwardPrimer, numIsolates, numRegions)

"""
Generates bucket pearson correlation value slice ranges
"""
def generateRanges(start, stop, step):
   r = start
   s = start + step
   while r < stop:
      yield (r,s)
      r += step
      s += step

"""
Calculates the enumerations space for numAlleles number of alleles and
numRegions number of ITS Region copies for the bacteria to be simulated. More
formally, the enumeration space is represents all unique combinations (with
replacement) of numAlleles and numRegions.
"""
def calcCombinations(numAlleles, numRegions):
   return (math.factorial(numAlleles + numRegions - 1)/
          (math.factorial(numRegions) * math.factorial(numAlleles - 1)))

'''
If this file is called as a main file, then run the following
'''
if __name__ == '__main__':
   if gpu_support:
      print "PyCUDA detected. GPU support enabled."
      main()

   else:
      print "PyCuda not detected. Exiting"
