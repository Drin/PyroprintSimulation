import os
import sys
import numpy
#import pycuda.autoinit
import pycuda.gpuarray
import pycuda.driver
import pycuda.tools

DEBUG = True

# CUDA parameters that seem to work well. The number of threads per tile
# (the tile_size) should be a power of 2 for the parallel reduction to
# work right!
def pearson(cudaModule, ranges, memoryOpt, num_alleles, alleles_per_isolate,
      num_isolates, length_alleles, alleleData=None, num_threads=16,
      num_blocks=32):
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


         progress = (tile_row * num_tiles + tile_col) * 100.0 / (num_tiles * num_tiles)
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

#TODO
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

   print('\rComputing correlations 100.000%')
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
