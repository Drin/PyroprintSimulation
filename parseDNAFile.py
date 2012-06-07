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
from extract_allele import expandSequence
from optparse import OptionParser
import ConfigParser
import pycuda.compiler
import pycuda.driver

DEBUG = True
LEGACY = False
CORE_KERNEL_FILE = "biogpu/kernels.cu"
PEARSON_KERNEL_FILES = {'constantMem': 'biogpu/const_pearson.cu',
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
   if DEBUG:
      print("Configuring simulation...\n")

   (dataDir, disp, primer, maxAlleles, numRegions, memoryOpt, deviceId) = handleArgs()

   if DEBUG:
      print("Extracting allele information...\n")

   (alleles, num_alleles, num_pyro_peaks) = extractAlleles(dataDir, disp, primer, maxAlleles)

   if DEBUG:
      print(("Running simulation for {0} alleles and {1} length " +
            "pyroprints...\n").format(num_alleles, num_pyro_peaks))

   if DEBUG:
      print("Preparing buckets for pearson correlation slices...\n")

   ranges = [(-1.000, 1.000)]
   bucketSlices = generateRanges(0.9900, 1.0000, 0.0001)
   ranges.extend([bucketSlice for bucketSlice in bucketSlices])

   startTime = time.time()

   if LEGACY:
      print ("Running Legacy Code!\n")

      pycuda.driver.init()
      if (pycuda.driver.Device.count() >= deviceId):
         cudaDevice = pycuda.driver.Device(deviceId)
         cudaContext = cudaDevice.make_context()
       
         #find all combinations
         comboLimit = -1
         numCombos = 0
         allPyroPrints = []
         allCombinations = combinations_with_replacement(alleles, 7)

         if DEBUG:
            print ("Pyroprinting Sequences...\n")
        
         for oneCombo in allCombinations:
       
            if comboLimit > 0 and numCombos >= comboLimit:
               break 
       
            allPyroPrints.append(pyroprintData(oneCombo, expandSequence(disp)))
            numCombos += 1
           
            # Writing out every number is too much I/O, it slows down pyroprinting a lot.
            '''
            if numCombos % 10 == 0:
               sys.stdout.write("\rGenerating Combo %d" % numCombos)
               sys.stdout.flush()
            '''
       
         if DEBUG:
            print "\n" + str(numCombos) + " Pyroprints Generated"
         
         buckets = biogpu.correlation.legacyPearson(allPyroPrints, allPyroPrints, ranges)
      
         cudaContext.detach()

   else:
      if DEBUG:
         print('Preparing pyroprinted alleles for device constant memory...\n')

      alleleData_gpu = numpy.zeros(shape=(num_alleles, num_pyro_peaks), dtype=numpy.uint8, order='C')
      for alleleNdx in range(num_alleles):
         numpy.put(alleleData_gpu[alleleNdx], range(num_pyro_peaks), alleles[alleleNdx])

      if DEBUG:
         print("Loading CUDA kernel source code...\n")

      pearson_kernel = PEARSON_KERNEL_FILES[memoryOpt]

      if DEBUG:
         print("loading kernel file '{0}'\n".format(pearson_kernel))

      kernelFile = open(os.path.join(os.getcwd(), CORE_KERNEL_FILE), 'r')
      cudaSrc = kernelFile.read()
      kernelFile.close()

      kernelFile = open(os.path.join(os.getcwd(), pearson_kernel), 'r')
      cudaSrc = cudaSrc + kernelFile.read()
      kernelFile.close()

      #init the pycuda driver and now GPU interfacing begins
      pycuda.driver.init()

      if (pycuda.driver.Device.count() >= deviceId):
         if DEBUG:
            print("using single GPU '{0}'\n".format(deviceId))

         cudaDevice = pycuda.driver.Device(deviceId)
         cudaContext = cudaDevice.make_context()
         cudaModule = pycuda.compiler.SourceModule(cudaSrc)

         if DEBUG:
            print('Calculating pearson correlation for all pairwise combinations ' +
                  'for {0} generated isolates...\n'.format(calcCombinations(num_alleles, numRegions)))

         startTime = time.time()
         buckets = biogpu.correlation.pearson(cudaModule, ranges, memoryOpt,
                   num_alleles, numRegions, calcCombinations(num_alleles,
                   numRegions), num_pyro_peaks, alleleData=alleleData_gpu)

         cudaContext.detach()
      else:
         print("Error: invalid device Id '{0}'\n".format(deviceId))
         sys.exit()

   '''
   Not sure why this doesn't work just yet. Commenting out for testing though
   elif (pycuda.driver.Device.count() > 1):
      if DEBUG:
         print("Detected multiple GPUs!\n")

      gpuEnvs = [(None, None)] * pycuda.driver.Device.count()

      for envNdx in range(len(gpuEnvs)):
         cudaDevice = pycuda.driver.Device(envNdx)
         cudaContext = cudaDevice.make_context()
         cudaModule = pycuda.compiler.SourceModule(cudaSrc)

         gpuEnvs[envNdx] = (cudaContext, cudaModule)

         if DEBUG:
            print("Popping context from device {0}({1})\n".format(
                  pycuda.driver.Context.get_device().name(),
                  pycuda.driver.Context.get_device().pci_bus_id()))
         pycuda.driver.Context.pop()

      if DEBUG:
         print('Calculating pearson correlation for all pairwise combinations ' +
               'for {0} generated isolates...\n'.format(calcCombinations(num_alleles, numRegions)))

      buckets = biogpu.correlation.multi_pearson(gpuEnvs, ranges, memoryOpt,
                num_alleles, numRegions, calcCombinations(num_alleles,
                numRegions), num_pyro_peaks, alleleData=alleleData_gpu)
   '''

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
   parser.add_option("--device", dest="device", type="int", default=0, help="Selects the GPU to use for computation")
   parser.add_option("--memory", dest="memory", default="constantMem", help=("What " +
                     "memory to store alleles in. Possible options are " +
                     "'constantMem', 'globalMem', 'sharedMem', 'texturedMem'"))

   (options, args) = parser.parse_args()

   if options.file:
      config = ConfigParser.RawConfigParser()
      config.read(options.file)
      dataPath = config.get("params", "path")

      if config.has_option("params", "maxAlleles"):
         maxAlleles = config.getint("params", "maxAlleles")
      else:
         maxAlleles = options.maxAlleles
     
      if config.has_option("params", "numRegions"):
         numRegions = config.getint("params", "numRegions")

      if config.has_option("params", "device"):
         deviceId = config.getint("params", "device")
      else:
         deviceId = options.device


      if config.has_option("params", "memory"):
         memoryOpt = config.get("params", "memory")
      else:
         memoryOpt = options.memory

      forwardPrimer = config.get("params", "primer")
      dispSeq = config.get("params", "DispSeq")


   else:
      #Use command line args
      dataPath = options.path
      maxAlleles = options.maxAlleles
      dispSeq = options.DispSeq
      forwardPrimer = options.primer
      numRegions = options.numRegions
      memoryOpt = options.memory
      deviceId = options.device

   if memoryOpt not in PEARSON_KERNEL_FILES:
      print "parameter \"memory\" needs to be one of constantMem, globalMem, or sharedMem"
      sys.exit()


   if DEBUG:
      print ("maxAlleles should be {0}\n".format(maxAlleles))
      print ("numRegions should be {0}\n".format(numRegions))

   return (dataPath, dispSeq, forwardPrimer, maxAlleles, numRegions, memoryOpt, deviceId)

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
This is old legacy code for the purposes of testing current implementation against old implementation.
'''
def pyroprintData(oneCombo, dispSeq):
   sequence = oneCombo
   
   #Saved heights from all 7 sequences
   pyroData = [[],[],[],[],[],[],[]]
   #Final heights
   height = []
   #Current sequence position
   seqCount = 0
   #Current disposition
   dispCount = 0
   #Current height
   pyroCount = 0
   #Length of sequences
   length = [len(sequence[0]), len(sequence[1]), len(sequence[2]), len(sequence[3]), len(sequence[4]), len(sequence[5]), len(sequence[6])]
   #Sequence Counter
   t=0
   
   #Go through the 7 sequences and run through the disposition sequence getting the heights
   while t < 7:
      while seqCount < length[t]:
         if sequence[t][seqCount] == dispSeq[dispCount]:
            pyroCount += 1
            seqCount += 1
            if seqCount == length[t]:
               pyroData[t].append(pyroCount)
         elif (sequence[t][seqCount] != 'A') & (sequence[t][seqCount] != 'T') & (sequence[t][seqCount] != 'C') & (sequence[t][seqCount] != 'G'):
            seqCount += 1
            dispCount += 1
            pyroData[t].append(pyroCount)
         else:
            pyroData[t].append(pyroCount)
            pyroCount = 0
            dispCount += 1
      seqCount = 0
      dispCount = 0
      pyroCount = 0
      t += 1
   
   seqCount = 0

   #Get the max length of the heights (since they can be different - finish quicker/slower)
   maxVal = max(len(pyroData[0]),len(pyroData[1]),len(pyroData[2]),len(pyroData[3]),len(pyroData[4]),len(pyroData[5]),len(pyroData[6]))
   
   #Pad the heights that do not have 0's that need them for adding (to make all the lengths the same)
   x=0
   while x < 7:
      t = len(pyroData[x])
      while (len(dispSeq) - t) > 0:
         pyroData[x].append(0)
         t += 1
      x += 1
   
   #Get the final heights
   while seqCount < len(dispSeq):
      height.append( int(pyroData[0][seqCount]) + int(pyroData[1][seqCount]) +
                     int(pyroData[2][seqCount]) + int(pyroData[3][seqCount]) + 
                     int(pyroData[4][seqCount]) + int(pyroData[5][seqCount]) + 
                     int(pyroData[6][seqCount]))
      seqCount += 1
   
   return height

'''
More of Bob Legacy Code
'''
def combinations_with_replacement(iterable, r):
   # combinations_with_replacement('ABC', 2) --> AA AB AC BB BC CC
   pool = tuple(iterable)
   n = len(pool)
   if not n and r:
      return
   indices = [0] * r
   yield tuple(pool[i] for i in indices)
   while True:
      for i in reversed(range(r)):
         if indices[i] != n - 1:
            break
      else:
         return
      indices[i:] = [indices[i] + 1] * (r - i)
      yield tuple(pool[i] for i in indices)

'''
If this file is called as a main file, then run the following
'''
if __name__ == '__main__':
   if gpu_support:
      print "PyCUDA detected. GPU support enabled."
      main()

   else:
      print "PyCuda not detected. Exiting"
