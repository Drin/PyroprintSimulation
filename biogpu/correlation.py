import os
import sys
import numpy
import pycuda.autoinit
import pycuda.compiler
import pycuda.gpuarray
import pycuda.driver
from math import factorial

   # CUDA parameters that seem to work well. The number of threads per tile
   # (the tile_size) should be a power of 2 for the parallel reduction to
   # work right!
def pearson(kernel, ranges, numIsolates, length_alleles, num_threads=16, num_blocks=64):
   pearson_cuda = kernel.get_function('pearson')
   reduction_cuda = kernel.get_function('reduction')

   num_buckets = len(ranges)
   tile_size = num_threads * num_blocks
   num_tiles = (numIsolates / tile_size + 1, numIsolates / tile_size + 1)

   # Copy the ranges into a numpy array.
   ranges_np = numpy.zeros(shape=(num_buckets, 2), dtype=numpy.float32, order='C')
   for bucketNdx in range(num_buckets):
       numpy.put(ranges_np[bucketNdx], range(2), ranges[bucketNdx])

   # Create a zero-initialized chunk of memory for the per-thread buckets and
   # copy it to the GPU.
   buckets = numpy.zeros(shape=(tile_size * tile_size * num_buckets, 1),
                      dtype=numpy.uint64, order='C')
   buckets_gpu = pycuda.gpuarray.to_gpu(buckets)

   # Do a kernel launch for each tile
   for s in range(num_tiles[0]):
       for t in range(num_tiles[1]):
           pearson_cuda(buckets_gpu.gpudata,
                        pycuda.driver.In(ranges_np), 
                        numpy.uint32(num_buckets),
                        numpy.uint32(tile_size), 
                        numpy.uint32(s), 
                        numpy.uint32(t),
                        numpy.uint32(numIsolates),
                        numpy.uint32(length_alleles),
                        block=(num_threads, num_threads, 1),
                        grid=(num_blocks, num_blocks))

           progress = (s * num_tiles[1] + t) * 100.0 / (num_tiles[0] * num_tiles[1])
           sys.stdout.write('\rComputing correlations %.3f%%' % progress)
           sys.stdout.flush()

   print('\rComputing correlations 100.000%')
   sys.stdout.write('Merging buckets...\n')
   sys.stdout.flush()

   # Do a parallel reduction to sum all the buckets element-wise.
   reduction_cuda(buckets_gpu.gpudata, numpy.uint32(num_buckets),
                  numpy.uint32(tile_size), numpy.uint32(num_blocks),
                  block=(num_threads, 1, 1),
                  grid=(tile_size, 1))

   # Copy buckets back from GPU.
   buckets_gpu.get(buckets)

   # Merge the results of the reduction from the first column of the matrix.
   merged = [0 for bucket in range(num_buckets)]

   for bucketNdx in range(num_buckets):
       for tileNdx in range(tile_size):
           bucket_index = ((tile_size * tile_size * bucketNdx) +
                          (tile_size * tileNdx) + 0)
           merged[bucketNdx] += buckets[bucket_index]

   print('finished merging buckets...\n')
   return merged
