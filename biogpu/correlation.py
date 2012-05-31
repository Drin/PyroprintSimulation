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
def pearson(kernel, ranges, testingOpts, num_alleles, alleles_per_isolate,
      num_isolates, length_alleles, num_threads=16, num_blocks=64,
      globalAlleles=None):
   pearson_cuda = kernel.get_function('pearson')
   reduction_cuda = kernel.get_function('reduction')

   num_buckets = len(ranges)
   tile_size = num_threads * num_blocks      # Threads per row/col of a tile
   num_tiles = num_isolates / tile_size + 1  # Tiles per row/col of matrix

   # Copy the ranges into a numpy array.
   ranges_np = numpy.zeros(shape=(num_buckets, 2), dtype=numpy.float32, order='C')
   for bucketNdx in range(num_buckets):
       numpy.put(ranges_np[bucketNdx], range(2), ranges[bucketNdx])

   # Create a zero-initialized chunk of memory for the per-thread buckets and
   # copy it to the GPU.
   buckets = numpy.zeros(shape=(tile_size * tile_size * num_buckets, 1), dtype=numpy.uint64, order='C')
   buckets_gpu = pycuda.gpuarray.to_gpu(buckets)

   # Do a kernel launch for each tile
   for tile_row in range(num_tiles):
      for tile_col in range(tile_row, num_tiles):
         if testingOpts == "none":
            pearson_cuda(buckets_gpu.gpudata,
                        pycuda.driver.In(ranges_np), 
                        numpy.uint32(num_buckets),
                        numpy.uint32(tile_size), 
                        numpy.uint32(tile_row), 
                        numpy.uint32(tile_col),
                        numpy.uint8(num_alleles),
                        numpy.uint8(alleles_per_isolate),
                        numpy.uint32(num_isolates),
                        numpy.uint32(length_alleles),
                        block=(num_threads, num_threads, 1),
                        grid=(num_blocks, num_blocks))
         else:
            pearson_cuda(buckets_gpu.gpudata,
                  pycuda.driver.In(ranges_np), 
                  pycuda.driver.In(globalAlleles),
                  numpy.uint32(num_buckets),
                  numpy.uint32(tile_size), 
                  numpy.uint32(tile_row), 
                  numpy.uint32(tile_col),
                  numpy.uint32(num_isolates),
                  numpy.uint32(length_alleles),
                  block=(num_threads, num_threads, 1),
                  grid=(num_blocks, num_blocks))



         progress = (tile_row * num_tiles + tile_col) * 100.0 / (num_tiles * num_tiles)
         sys.stdout.write('\rComputing correlations %.3f%%' % progress)
         sys.stdout.flush()

#print('\rComputing correlations 100.000%')
#sys.stdout.write('Merging buckets...\n')
#sys.stdout.flush()

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
