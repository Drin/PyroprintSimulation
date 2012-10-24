#!/usr/bin/python -tt
# Copyright 2010 Google Inc.
# Licensed under the Apache License, Version 2.0
# http://www.apache.org/licenses/LICENSE-2.0

import os
import sys
import time
import math
import copy

import Queue
import numpy
import threading

import SimulationConfig
import SequenceParser
import Allele

from SimulationThread import SimulationThread

try:
   import pycuda.gpuarray
   import pycuda.driver

except:
   print("Could not load pyCUDA")
   sys.exit(1)

CORE_KERNEL_FILE = "biogpu/kernels.cu"

"""
Generates bucket pearson correlation value ranges in step incremements
"""
def generateRanges(start, stop, step):
   range_start = start
   range_end = start + step

   while (range_end < stop):
      yield (range_start, range_end)
      range_start = range_end
      range_end += step

"""
Calculates the enumerations space for numAlleles number of alleles and
numRegions number of ITS Region copies for the bacteria to be simulated. More
formally, the enumeration space is represents all unique combinations (with
replacement) of numAlleles and numRegions.
"""
def calcCombinations(numAlleles, numRegions):
   return (math.factorial(numAlleles + numRegions - 1)/
          (math.factorial(numRegions) * math.factorial(numAlleles - 1)))

class PyroprintSimulation(object):
   # begin - beginning of bucket range
   # end   - end of bucket range
   # step  - increment size between buckets
   def __init__(self, begin=0.9900, end=1.0000, step=0.0001):
      cli_desc = "Parameter configuration for a pyroprint simulation"
      self.config = SimulationConfig.parse_args(cli_desc)

      self.ranges = [(-1.000, 1.000)]
      self.ranges.extend([bucket for bucket in generateRanges(begin, end, step)])

   def prep_simulation(self):
      (allele_set, disp_len) = SequenceParser.extractAlleles(self.config)
      self.alleles = list(allele_set)
      self.num_alleles = len(self.alleles)
      self.pyro_len = max(disp_len, self.config.get('pyro_len'))
      self.num_isolates = calcCombinations(self.num_alleles,
                                           self.config.get('num_loci'))

      # Copy the alleles into a numpy array.
      self.alleles_cpu = numpy.zeros(shape=(self.num_alleles, self.pyro_len),
                                 dtype=numpy.uint8, order='C')

      for allele_ndx in range(self.num_alleles):
         numpy.put(self.alleles_cpu[allele_ndx], range(self.pyro_len),
                   self.alleles[allele_ndx].pyroprint)

   def get_state(self):
      state = ''

      for (paramKey, paramVal) in self.config.iteritems():
         state += "{0} => {1}\n".format(paramKey, paramVal)

      state += '\nAlleles:\n'

      for allele in self.alleles:
         state += "seq: {0}\n{1}\n".format(allele.sequence, allele.pyroprint)
         #state += "seq: {0}\n".format(allele.sequence)

      return state

   def load_cuda(self, *cuda_files):
      cuda_src = ''

      for src_file in cuda_files:
         if ('DEBUG' in os.environ):
            print("loading cuda from file '{0}'\n".format(src_file))

         open_file = open(os.path.join(os.getcwd(), src_file), 'r')
         cuda_src += open_file.read()
         open_file.close()

      self.cuda_src = cuda_src

   def has_loaded_cuda(self):
      return self.cuda_src is not None

   def is_ready(self):
      if (self.config is None): print("Simulation config not initialized.")
      if (self.cuda_src is None): print("CUDA code not loaded.")
      if (self.alleles is None): print("DNA not yet pyrosequenced.")

      return (self.config is not None and self.cuda_src is not None and
              self.alleles is not None)

   def run(self, num_threads=16, num_blocks=32):
      if (not self.is_ready()):
         print("Aborting attempt to run simulation before fully prepared.")
         return

      pycuda.driver.init()

      sim_threads = []
      sim_thread_num = (min(self.config.get('num_gpus'), pycuda.driver.Device.count()) or
                        pycuda.driver.Device.count())
      task_queues = [Queue.Queue() for thread_num in range(sim_thread_num)]
      task_results = [[0 for count in range(len(self.ranges))] for thread_num
                      in range(sim_thread_num)]
      progress_queue = Queue.Queue()

      if ('DEBUG' in os.environ):
         print("Num expected isolates: {0}".format(self.num_isolates))
         '''
         pyro1 = copy.deepcopy(self.alleles[0].pyroprint)
         pyro2 = copy.deepcopy(self.alleles[5].pyroprint)

         for x in range(1,7):
            for y in range(self.pyro_len):
               pyro1[y] += self.alleles[x].pyroprint[y]
               pyro2[y] += self.alleles[x+5].pyroprint[y]

         print("Pyro1: {0}\nPyro2: {1}".format(pyro1, pyro2))
         '''

      if ('NO_CUDA' in os.environ):
         sys.exit(0)
      startTime = time.time()

      tile_size = num_threads * num_blocks
      num_tiles = self.num_isolates / tile_size + 1
      num_task = 0
                            
      for tile_row in range(num_tiles):
         for tile_col in range(tile_row, num_tiles):
            shape = {'threads': num_threads, 'blocks': num_blocks}
            config = {'num': num_tiles, 'row': tile_row, 'col': tile_col}

            '''
            if ('DEBUG' in os.environ):
               print("adding task for thread %d" %
                     (num_task % sim_thread_num))
            '''
            task_queues[num_task % sim_thread_num].put({'config': config, 'shape': shape})
            num_task += 1

      buckets = [0 for count in range(len(self.ranges))]
      #task_result = {'lock': threading.RLock(), 'buckets': buckets}

      # start each thread
      for thread_ndx in range(sim_thread_num):
         sys.stdout.write('Starting thread {0}...\n'.format(thread_ndx))
         sys.stdout.flush()
         sim_thread = SimulationThread(thread_ndx, tile_size, num_threads,
                                       num_blocks, task_queues[thread_ndx],
                                       progress_queue,
                                       task_results[thread_ndx], self)
         sim_thread.start()
         sim_threads.append(sim_thread)

      # wait for each thread to finish
      for thread in sim_threads:
         thread.join()

      # aggregate results from each thread
      for result in task_results:
         for ndx in range(len(result)):
            buckets[ndx] += result[ndx]

      if ('TIMED' in os.environ):
         print("Elapsed Time: {0}\n".format(time.time() - startTime))

      print('Results:\n')
      for i in range(len(buckets)):
         print('\t[%d] (%.4f%%, %.4f%%) = %d' % (i, self.ranges[i][0] *
               100.0, self.ranges[i][1] * 100.0, buckets[i]))
      print('\n')

if __name__ == '__main__':
   simulation = PyroprintSimulation()
   simulation.prep_simulation()
   simulation.load_cuda(CORE_KERNEL_FILE, simulation.config.get('kernel_file'))

   if ('DEBUG' in os.environ):
      print("Simulation state: {0}\n".format(simulation.get_state()))
      print(("Simulation has {0} alleles and {1} pyroprint length...\n").format(
             simulation.num_alleles, simulation.pyro_len))

   simulation.run()
