import os
import sys
import numpy as np
import pycuda.autoinit
import pycuda.compiler
import pycuda.gpuarray
import pycuda.driver
from math import factorial

def pearson(ranges, c, length_alleles):
   # Load the kernel and compile it.
   kernel_file = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                              'pearson.cu')
   f = open(kernel_file, 'r')
   kernel = pycuda.compiler.SourceModule(f.read())
   f.close()
   pearson_cuda = kernel.get_function('pearson')
   reduction_cuda = kernel.get_function('reduction')

   # CUDA parameters that seem to work well. The number of threads per tile
   # (the tile_size) should be a power of 2 for the parallel reduction to
   # work right!
   threads_per_block = 16
   blocks_per_tile = 64
   tile_size = threads_per_block * blocks_per_tile
   num_tiles = (c / tile_size + 1, c / tile_size + 1)

   # Copy the ranges into a numpy array.
   ranges_np = np.zeros(shape=(len(ranges), 2), dtype=np.float32, order='C')
   for i in range(len(ranges)):
       np.put(ranges_np[i], range(2), ranges[i])

   # Create a zero-initialized chunk of memory for the per-thread buckets and
   # copy it to the GPU.
   buckets = np.zeros(shape=(tile_size * tile_size * len(ranges), 1),
                      dtype=np.uint64, order='C')
   buckets_gpu = pycuda.gpuarray.to_gpu(buckets)

   # Do a kernel launch for each tile
   for s in range(num_tiles[0]):
       for t in range(num_tiles[1]):
           pearson_cuda(buckets_gpu.gpudata,
                        pycuda.driver.In(ranges_np), 
                        np.uint32(len(ranges)),
                        np.uint32(tile_size), 
                        np.uint32(s), 
                        np.uint32(t),
                        np.uint(c),
                        np.uint32(length_alleles),
                        block=(threads_per_block, threads_per_block, 1),
                        grid=(blocks_per_tile, blocks_per_tile))

           progress = (s * num_tiles[1] + t) * 100.0 / (num_tiles[0] * num_tiles[1])
           sys.stdout.write('\rComputing correlations %.3f%%' % progress)
           sys.stdout.flush()

   print('\rComputing correlations 100.000%')
   sys.stdout.write('Merging buckets... ')
   sys.stdout.flush()

   # Do a parallel reduction to sum all the buckets element-wise.
   reduction_cuda(buckets_gpu.gpudata, np.uint32(len(ranges)),
                  np.uint32(tile_size), np.uint32(blocks_per_tile),
                  block=(threads_per_block, 1, 1),
                  grid=(tile_size, 1))

   # Copy buckets back from GPU.
   buckets_gpu.get(buckets)

   # Merge the results of the reduction from the first column of the matrix.
   merged = [0 for k in range(len(ranges))]
   for k in range(len(ranges)):
       for i in range(tile_size):
           bucket_index = (tile_size * tile_size * k) + (tile_size * i) + 0
           merged[k] += buckets[bucket_index]
   print('done.')

   return merged

def cr(n, r):
   return factorial(n + r - 1)/(factorial(r) * factorial(n-1))
