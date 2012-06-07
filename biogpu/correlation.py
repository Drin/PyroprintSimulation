import os
import sys
import numpy
#import pycuda.autoinit
import pycuda.gpuarray
import pycuda.driver
import pycuda.tools
import pycuda.compiler

DEBUG = True

def legacyPearson(X, Y, ranges):
   # Check some preconditions to verify we're doing something sensical.
   # Doesn't cover all cases, but catches obvious mistakes.
   assert len(X[0]) == len(X[1]), 'Your sequences in X should all be the same length.'
   assert len(Y[0]) == len(Y[1]), 'Your sequences in Y should all be the same length.'
   assert len(X[0]) == len(Y[0]), 'Your sequences in X and Y should all be the same length.'

   n = len(X)
   m = len(Y)
   p = len(X[0])

   # Load the kernel and compile it.
   kernel_file = os.path.join(os.getcwd(), 'biogpu', 'bob_pearson.cu')

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
   num_tiles = (n / tile_size + 1, m / tile_size + 1)

   # Copy the ranges into a numpy array.
   ranges_np = numpy.zeros(shape=(len(ranges), 2), dtype=numpy.float32, order='C')
   for i in range(len(ranges)):
      numpy.put(ranges_np[i], range(2), ranges[i])

   # Create a zero-initialized chunk of memory for the per-thread buckets and
   # copy it to the GPU.
   buckets = numpy.zeros(shape=(tile_size * tile_size * len(ranges), 1),
                      dtype=numpy.uint64, order='C')
   buckets_gpu = pycuda.gpuarray.to_gpu(buckets)

   # Do a kernel launch for each tile, copying the appropriate chunks of the
   # input arrays into X and Y for each launch.
   for s in range(num_tiles[0]):
      for t in range(num_tiles[1]):
         num_A = tile_size
         remain_X = n - (s * tile_size)
         num_A = num_A if num_A < remain_X else remain_X

         A = numpy.zeros(shape=(num_A, p), dtype=numpy.float32, order='C')
         for i in range(num_A):
             numpy.put(A[i], range(p), X[(s * tile_size) + i])

         num_B = tile_size
         remain_Y = m - (t * tile_size)
         num_B = num_B if num_B < remain_Y else remain_Y

         B = numpy.zeros(shape=(num_B, p), dtype=numpy.float32, order='C')
         for j in range(num_B):
             numpy.put(B[j], range(p), Y[(t * tile_size) + j])

         pearson_cuda(buckets_gpu.gpudata,
                      pycuda.driver.In(ranges_np), numpy.uint32(len(ranges)),
                      pycuda.driver.In(A), pycuda.driver.In(B),
                      numpy.uint32(tile_size), numpy.uint32(s), numpy.uint32(t),
                      numpy.uint32(n), numpy.uint32(m), numpy.uint32(p),
                      block=(threads_per_block, threads_per_block, 1),
                      grid=(blocks_per_tile, blocks_per_tile))

         '''
         progress = (s * num_tiles[1] + t) * 100.0 / (num_tiles[0] * num_tiles[1])
         sys.stdout.write('\rComputing correlations %.3f%%' % progress)
         sys.stdout.flush()
         '''

   #print('\rComputing correlations 100.000%')
   sys.stdout.write('Merging buckets... ')
   sys.stdout.flush()

   # Do a parallel reduction to sum all the buckets element-wise.
   reduction_cuda(buckets_gpu.gpudata, numpy.uint32(len(ranges)),
                  numpy.uint32(tile_size), numpy.uint32(blocks_per_tile),
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

# CUDA parameters that seem to work well. The number of threads per tile
# (the tile_size) should be a power of 2 for the parallel reduction to
# work right!
def pearson(cudaModule, ranges, memoryOpt, num_alleles, alleles_per_isolate,
      num_isolates, length_alleles, alleleData=None, num_threads=16,
      num_blocks=64):
   pearson_cuda = cudaModule.get_function('pearson')
   reduction_cuda = cudaModule.get_function('reduction')

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

   alleles_gpu = None

   if memoryOpt == "constantMem":
      (const_ptr, size) = cudaModule.get_global("alleles")
      pycuda.driver.memcpy_htod(const_ptr, alleleData)

   elif memoryOpt == "globalMem" or memoryOpt == "texturedMem":
      alleles_gpu = pycuda.gpuarray.to_gpu(alleleData)

   buckets_gpu = pycuda.gpuarray.to_gpu(buckets)

   # Do a kernel launch for each tile
   for tile_row in range(num_tiles):
      for tile_col in range(tile_row, num_tiles):
         if memoryOpt == "constantMem":
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

         elif memoryOpt == "globalMem" or memoryOpt == "texturedMem":
            pearson_cuda(buckets_gpu.gpudata,
                  pycuda.driver.In(ranges_np), 
                  alleles_gpu.gpudata,
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

         '''
         progress = (tile_row * num_tiles + tile_col) * 100.0 / (num_tiles * num_tiles)
         sys.stdout.write('\rComputing correlations %.3f%%' % progress)
         sys.stdout.flush()
         '''

   #print('\rComputing correlations 100.000%')
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

def multi_pearson(gpuEnvs, ranges, memoryOpt, num_alleles, alleles_per_isolate,
      num_isolates, length_alleles, alleleData=None, num_threads=16,
      num_blocks=64):
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
   bucketTotals = numpy.zeros(shape=(tile_size * tile_size * num_buckets, 1), dtype=numpy.uint64, order='C')

   alleles_gpu_devs = [None] * len(gpuEnvs)
   buckets_gpu_devs = [None] * len(gpuEnvs)

   for envNdx in range(len(gpuEnvs)):
      (cudaContext, cudaModule) = gpuEnvs[envNdx]

      cudaContext.push()

      if DEBUG:
         print("Pushed context onto device {0}({1})\n".format(
               pycuda.driver.Context.get_device().name(),
               pycuda.driver.Context.get_device().pci_bus_id()))

      buckets_gpu_devs[envNdx] = pycuda.gpuarray.to_gpu(buckets)

      if memoryOpt == "constantMem":
         (const_ptr, size) = cudaModule.get_global("alleles")
         pycuda.driver.memcpy_htod(const_ptr, alleleData)

      elif memoryOpt == "globalMem" or memoryOpt == "texturedMem":
         alleles_gpu_devs[envNdx] = pycuda.gpuarray.to_gpu(alleleData)

      if DEBUG:
         print("Popping context from device {0}({1})\n".format(
               pycuda.driver.Context.get_device().name(),
               pycuda.driver.Context.get_device().pci_bus_id()))

      pycuda.driver.Context.pop()

   tileNdx = 0
   # Do a kernel launch for each tile
   for tile_row in range(num_tiles):
      for tile_col in range(tile_row, num_tiles):
         envNdx = tileNdx % len(gpuEnvs)
         tileNdx = tileNdx + 1

         (currContext, currModule) = gpuEnvs[envNdx]

         currContext.push()

         if DEBUG:
            print("\nPushed context onto device {0}({1})\n".format(
                  pycuda.driver.Context.get_device().name(),
                  pycuda.driver.Context.get_device().pci_bus_id()))

         pearson_cuda = cudaModule.get_function('pearson')

         if memoryOpt == "constantMem":
            pearson_cuda(buckets_gpu_devs[envNdx].gpudata,
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

         elif memoryOpt == "globalMem" or memoryOpt == "texturedMem":
            pearson_cuda(buckets_gpu_devs[envNdx].gpudata,
                  pycuda.driver.In(ranges_np), 
                  alleles_gpu_devs[envNdx].gpudata,
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

         if DEBUG:
            print("Popping context from device {0}({1})\n".format(
                  pycuda.driver.Context.get_device().name(),
                  pycuda.driver.Context.get_device().pci_bus_id()))

         pycuda.driver.Context.pop()

         '''
         progress = (tile_row * num_tiles + tile_col) * 100.0 / (num_tiles * num_tiles)
         sys.stdout.write('\rComputing correlations %.3f%%' % progress)
         sys.stdout.flush()
         '''

   #print('\rComputing correlations 100.000%')
   sys.stdout.write('Merging buckets...\n')
   sys.stdout.flush()

   for envNdx in range(len(gpuEnvs)):
      (cudaContext, cudaModule) = gpuEnvs[envNdx]
      cudaContext.push()

      if DEBUG:
         print("Pushed context onto device {0}({1})\n".format(
               pycuda.driver.Context.get_device().name(),
               pycuda.driver.Context.get_device().pci_bus_id()))
         print("performing reduction...\n")

      reduction_cuda = cudaModule.get_function('reduction')

      # Do a parallel reduction to sum all the buckets element-wise.
      reduction_cuda(buckets_gpu_devs[envNdx].gpudata, numpy.uint32(num_buckets),
                     numpy.uint32(tile_size), numpy.uint32(num_blocks),
                     block=(num_threads, 1, 1),
                     grid=(tile_size, 1))

      if DEBUG:
         print("finished reduction kernel...\n")

      # Copy buckets back from GPU.
      buckets_gpu_devs[envNdx].get(buckets)

      for bucketNdx in range(len(bucketTotals)):
         bucketTotals[bucketNdx] = bucketTotals[bucketNdx] + buckets[bucketNdx]

      if DEBUG:
         print("detaching context from device {0}({1})\n".format(
               pycuda.driver.Context.get_device().name(),
               pycuda.driver.Context.get_device().pci_bus_id()))

      pycuda.driver.Context.pop()
      cudaContext.detach()

   #clear caches from GPU
   pycuda.tools.clear_context_caches()

   # Merge the results of the reduction from the first column of the matrix.
   merged = [0 for bucket in range(num_buckets)]

   for bucketNdx in range(num_buckets):
      for tileNdx in range(tile_size):
         bucket_index = ((tile_size * tile_size * bucketNdx) +
                        (tile_size * tileNdx) + 0)
         merged[bucketNdx] += bucketTotals[bucket_index]

   print('finished merging buckets...\n')
   return merged
