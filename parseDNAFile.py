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
import time
from extract_allele import extractAlleles
from optparse import OptionParser
import ConfigParser
import pycuda.autoinit
import pycuda.compiler
import pycuda.gpuarray
import pycuda.driver

DEBUG = True
TESTING = False
KERNELFILES = {'none': 'biogpu/pearson.cu',
               'globalMem': 'biogpu/ga_pearson.cu'} #,
               #'sharedMem': 'biogpu/sh_pearson.cu',
               #'texturedMem': 'biogput/text_pearson.cu'}

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
   startTime = time.time()

   if DEBUG:
      print("Configuring simulation...\n")

   (dataDir, disp, primer, maxAlleles, numRegions, memoryOpt) = handleArgs()

   if DEBUG:
      print("Extracting allele information...\n")

   (alleles, num_alleles, num_pyro_peaks) = extractAlleles(dataDir, disp, primer, maxAlleles)

   if DEBUG:
      print(("Running simulation for {0} alleles and {1} length " +
            "pyroprints...\n").format(num_alleles, num_pyro_peaks))

   if DEBUG:
      print("Preparing buckets for pearson correlation slices...\n")

   ranges = [(-1.000, 0.000), (0.000, 1.000)]
   bucketSlices = generateRanges(0.9900, 1.0000, 0.0001)
   ranges.extend([bucketSlice for bucketSlice in bucketSlices])

   if DEBUG:
      print('Preparing pyroprinted alleles for device constant memory...\n')

   alleles_c = numpy.zeros(shape=(num_alleles, num_pyro_peaks), dtype=numpy.uint8, order='C')
   for alleleNdx in range(num_alleles):
      numpy.put(alleles_c[alleleNdx], range(num_pyro_peaks), alleles[alleleNdx])

   if DEBUG:
      print("Loading CUDA kernel source code...\n")

   kernel_file = KERNELFILES[memoryOpt]

   kernelFile = open(os.path.join(os.getcwd(), kernel_file), 'r')
   kernel = pycuda.compiler.SourceModule(kernelFile.read())
   kernelFile.close()

   # default performance is to store alleles in constant memory
   if memoryOpt == "none":
      (const_ptr, size) = kernel.get_global("alleles")

      if DEBUG:
         print(("Transferring {0} bytes of allele pyroprints into device " +
            "constant memory...\n ").format(size))

      pycuda.driver.memcpy_htod(const_ptr, alleles_c)

   if DEBUG:
      print('Calculating pearson correlation for all pairwise combinations ' +
            'for {0} generated isolates...\n'.format(calcCombinations(num_alleles, numRegions)))

   if memoryOpt == "none":
      buckets = biogpu.correlation.pearson(kernel, ranges, memoryOpt,
                num_alleles, numRegions, calcCombinations(num_alleles,
                   numRegions), num_pyro_peaks, globalMem=None)
   elif memoryOpt == "globalMem":
      buckets = biogpu.correlation.pearson(kernel, ranges, memoryOpt,
                num_alleles, numRegions, calcCombinations(num_alleles, numRegions),
                num_pyro_peaks, globalMem=alleles_c)

   print("Elapsed Time: {0}\n".format(time.time() - startTime))
   print('Results:\n')
   for i in range(len(buckets)):
      print('\t[%d] (%.4f%%, %.4f%%) = %d' % (i, ranges[i][0] * 100.0, ranges[i][1] * 100.0, buckets[i]))
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
   parser.add_option("-p", "--path", dest="path", default="Genome Sequences/", help="Path to Genome Sequence Folders")
   parser.add_option("-d", dest="DispSeq", default="AACACGCGA23(GATC)GAA", help="Dispensation order")
   parser.add_option("-m", "--max", dest="maxAlleles", type="int", default=24, help="Max number of combinations to generate")
   parser.add_option("-f", "--file", dest="file", help="File containing parameters", default="config.cfg")
   parser.add_option("--primer", dest="primer", default="TTGGATCAC", help="Primer to use")
   parser.add_option("--numRegions", dest="numRegions", type="int", default=7, help="Number of ITS Regions to simulate")
   parser.add_option("--memory", dest="memory", default="none", help=("What " +
                     "memory to store alleles in. Possible options are " +
                     "'none', 'globalMem', 'sharedMem', 'texturedMem'"))

   (options, args) = parser.parse_args()

   if options.file:
      config = ConfigParser.RawConfigParser()
      config.read(options.file)
      dataPath = config.get("params", "path")

      if config.has_option("params", "maxAlleles"):
         maxAlleles = config.getint("params", "maxAlleles")
     
      if TESTING:
         maxAlleles = 4

      if config.has_option("params", "numRegions"):
         numRegions = config.getint("params", "numRegions")

      forwardPrimer = config.get("params", "primer")
      dispSeq = config.get("params", "DispSeq")
      memoryOpt = config.get("params", "memory")

   else:
      #Use command line args
      dataPath = options.path
      maxAlleles = options.maxAlleles
      dispSeq = options.DispSeq
      forwardPrimer = options.primer
      numRegions = options.numRegions
      memoryOpt = options.memory

   if memoryOpt not in KERNELFILES:
      print "parameter \"memory\" needs to be one of none, globalMem, or sharedMem"
      sys.exit()


   if DEBUG:
      print ("maxAlleles should be {0}\n".format(maxAlleles))
      print ("numRegions should be {0}\n".format(numRegions))

   return (dataPath, dispSeq, forwardPrimer, maxAlleles, numRegions, memoryOpt)

"""
Generates bucket pearson correlation value slice ranges
"""
def generateRanges(start, stop, step):
   r = start
   s = start + step
   while s < stop:
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
