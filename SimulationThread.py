#!/usr/bin/python -tt
# Copyright 2010 Google Inc.
# Licensed under the Apache License, Version 2.0
# http://www.apache.org/licenses/LICENSE-2.0

import os
import sys
import math
import time
import threading

import numpy
import pycuda.driver
import pycuda.compiler

class SimulationThread(threading.Thread):
   def __init__(self, device_id, tile_size, num_threads, num_blocks,
                task_queue, progress_queue, task_result, pyro_sim):
      super(SimulationThread, self).__init__()

      #given a device Id, create a new CUDA context. Since contexts must be
      #thread-local we initialize GPU Device contexts here
      self.pyro_sim = pyro_sim
      self.tile_size = tile_size
      self.num_threads = num_threads
      self.num_blocks = num_blocks

      self.cuda_device = pycuda.driver.Device(device_id)

      self.task_queue = task_queue
      self.progress_queue = progress_queue

      self.task_result = task_result

   def prep_alleles(self, cuda_module):
      alleles_gpu = None

      if (self.pyro_sim.config.get('memory') == "constant"):
         (const_ptr, size) = cuda_module.get_global("alleles")
         pycuda.driver.memcpy_htod(const_ptr, self.pyro_sim.alleles_cpu)

      elif (self.pyro_sim.config.get('memory') == "global" or
            self.pyro_sim.config.get('memory') == "textured"):
         alleles_gpu = pycuda.gpuarray.to_gpu(self.pyro_sim.alleles_cpu)

      return alleles_gpu

   def prep_gpu_buckets(self):
      # Copy the ranges into a numpy array.
      ranges_cpu = numpy.zeros(shape=(len(self.pyro_sim.ranges), 2),
                               dtype=numpy.float32, order='C')

      for bucketNdx in range(len(self.pyro_sim.ranges)):
          numpy.put(ranges_cpu[bucketNdx], range(2), self.pyro_sim.ranges[bucketNdx])

      # Create a zero-initialized array for the per-thread buckets.
      buckets_cpu = numpy.zeros(shape=(self.tile_size * self.tile_size *
                                       len(self.pyro_sim.ranges), 1),
                                dtype=numpy.uint64, order='C')
      buckets_gpu = pycuda.gpuarray.to_gpu(buckets_cpu)

      return (ranges_cpu, buckets_cpu, buckets_gpu)

   def run(self):
      cuda_context = self.cuda_device.make_context()
      cuda_module = pycuda.compiler.SourceModule(self.pyro_sim.cuda_src)

      (ranges_cpu, buckets_cpu, buckets_gpu) = self.prep_gpu_buckets()
      alleles_gpu = self.prep_alleles(cuda_module)

      # Compute correlations based on task_queue
      while (not self.task_queue.empty()):
         task_dict = self.task_queue.get()

         (shape, config) = (task_dict['shape'], task_dict['config'])

         pearson_kernel = cuda_module.get_function('pearson')
         pearson_kernel(buckets_gpu.gpudata,
                        pycuda.driver.In(ranges_cpu), 
                        numpy.uint32(len(self.pyro_sim.ranges)),
                        numpy.uint32(self.tile_size), 
                        numpy.uint32(config['row']), 
                        numpy.uint32(config['col']),
                        numpy.uint8(self.pyro_sim.num_alleles),
                        numpy.uint8(self.pyro_sim.config['num_loci']),
                        numpy.uint32(self.pyro_sim.num_isolates),
                        numpy.uint32(self.pyro_sim.pyro_len),
                        block = (shape['threads'], shape['threads'], 1),
                        grid  = (shape['blocks'],  shape['blocks']))

         '''
         progress = ((config['row'] * config['num'] + config['col']) * 100.0 /
                     (config['num'] * config['num']))
         sys.stdout.write('\rComputing correlations %.3f%%' % progress)
         sys.stdout.flush()
         '''

         self.task_queue.task_done()

      # Do a parallel reduction to sum all the buckets element-wise.
      reduction_kernel = cuda_module.get_function('reduction')
      reduction_kernel(buckets_gpu.gpudata,
                       numpy.uint32(len(self.pyro_sim.ranges)),
                       numpy.uint32(self.tile_size), numpy.uint32(self.num_blocks),
                       block=(self.num_threads, 1, 1),
                       grid=(self.tile_size, 1))

      # Copy buckets back from GPU.
      buckets_gpu.get(buckets_cpu)

      #reduce results from GPU buckets into cpu bucket totals
      #self.task_result['lock'].acquire()

      # Merge the results of the reduction from the first column of the matrix.
      # Merge results of reduction into the given array. Each thread has its
      # own result array to use
      for bucketNdx in range(len(self.pyro_sim.ranges)):
          for tileNdx in range(self.tile_size):
              bucket_index = ((self.tile_size * self.tile_size * bucketNdx) +
                              (tileNdx * self.tile_size))
              self.task_result[bucketNdx] += buckets_cpu[bucket_index]

      #self.task_result['lock'].release()

      #pop context so there are no problems of the GPU being used even after
      #this thread is finished
      pycuda.driver.Context.pop()

      #free context after it has been popped because otherwise this context
      #may be "leaked" in memory
      del cuda_context
