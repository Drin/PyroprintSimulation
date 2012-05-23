#include <stdint.h>
#include <string.h>

const int kNumAlleles = 24;
const int kAllelesPerIsolate = 7;

// n choose k
__device__ int comb(int n, int k) {
   if (n < k)
      return 0;

   uint64_t num = 1;
   uint64_t den = 1;
   
   for (int i = 1; i <= k; ++i)
      den *= i;
   
   for (int i = n-k+1; i <= n; ++i)
      num *= i;

   return (num/den);
}

// Gets isolate number |seq_num|, which is an index into a set of all isolates
// composed of |kAllelesPerIsolate| alleles, where each allele is one of 
// |kNumAlleles|. Order doesn't matter, and duplicates are allowed.
__device__ void get_isolate(int seq_num, uint8_t* alleles) {
   uint8_t cur_n = kNumAlleles;   // Shrinks every time we take a right turn
   uint8_t cur_row = kAllelesPerIsolate - 1;
   uint8_t cur_col = 0;

   // Init num to biggest "bucket"
   int num = comb(cur_n - cur_col - 1 + cur_row, cur_row);

   // Save in order to subtract from the "current" sequence number
   int old_num = 0;

   // Output |kAllelesPerIsolate| alleles
   int alleles_output = 0;
   while (1) {
      if (seq_num < num) {
         alleles[alleles_output] = cur_col + kNumAlleles - cur_n;
         if (++alleles_output >= kAllelesPerIsolate)
            break;

         seq_num -= old_num;
         old_num = 0;
         
         // New parameters for new location in tree
         cur_n -= cur_col;
         cur_row -= 1;

         num = comb(cur_n - 1 + cur_row, cur_row);

         cur_col = 0;
      } else {
         cur_col++;
         old_num = num;
         num += comb(cur_n - 1 - cur_col + cur_row, cur_row);
      }
   }
}

__constant__ uint8_t alleles[24 * 104];

__device__ void dump_bucket(uint64_t *buckets,
      uint32_t num_ranges, uint32_t tile_size,
      uint32_t src_i, uint32_t src_j,
      uint32_t dest_i, uint32_t dest_j) {
   // Element-wise sum for each in 0 -> num_ranges.
   for (uint32_t k = 0; k < num_ranges; k++) {
      uint32_t src_index = (tile_size * tile_size * k) +
         (tile_size * src_i) + src_j;
      uint32_t dest_index = (tile_size * tile_size * k) +
         (tile_size * dest_i) + dest_j;
      buckets[dest_index] += buckets[src_index];
   }
}

__global__ void reduction(uint64_t *buckets, uint32_t num_ranges,
      uint32_t tile_size, uint32_t chunk_size) {
   // Calculate <i, j> coords within the tile.
   uint32_t i = blockIdx.x; // row
   uint32_t j = threadIdx.x * chunk_size; // column

   // Each chunk leader reduces its chunk.
   for (uint32_t k = 1; k < chunk_size; k++) {
      dump_bucket(buckets, num_ranges, tile_size, i, j + k, i, j);
   }

   // Wait for all the threads in this row to finish.
   __syncthreads();

   // Reduce each chunk leader into the zeroth element of the row.
   if (j == 0) {
      for (uint32_t k = 1; k < blockDim.x; k++) {
         dump_bucket(buckets, num_ranges, tile_size, i, k * chunk_size, i, 0);
      }
   }
}

__global__ void pearson(uint64_t *buckets,
      float *ranges, uint32_t num_ranges,
      uint32_t tile_size, uint32_t s, uint32_t t,
      uint32_t c, uint32_t p) {
   // Calculate relative <i, j> coords within this tile.
   uint32_t i = blockIdx.y * blockDim.y + threadIdx.y; // row
   uint32_t j = blockIdx.x * blockDim.x + threadIdx.x; // column

   // Calculate the offsets based on the tile number.
   uint32_t i_offset = s * tile_size;
   uint32_t j_offset = t * tile_size;

   // Calculate the absolute <i, j> coords within the matrix.
   uint32_t i_abs = i_offset + i;
   uint32_t j_abs = j_offset + j;

   // Only compute values inside the bounds of the matrix.
   if (i_abs >= c || j_abs >= c)
      return;

   //TODO fix all these magic numbers
   // Generate isolate |i_abs| and |j_abs|
   uint8_t i_allele_indices[7]; 
   uint8_t j_allele_indices[7]; 
   get_isolate(i_abs, i_allele_indices);
   get_isolate(j_abs, j_allele_indices);

   // Allocate space for the two isolates
   uint16_t i_isolate[104]; 
   uint16_t j_isolate[104]; 
   memset(i_isolate, '\0', sizeof(uint16_t) * 104); 
   memset(j_isolate, '\0', sizeof(uint16_t) * 104); 

   // Add up each 7 alleles to generate the two isolates
   for (int num1 = 0; num1 < 7; ++num1) {
      for (int index = 0; index < 104; ++index) {
         i_isolate[index] += alleles[i_allele_indices[num1] * 104 + index];
         j_isolate[index] += alleles[j_allele_indices[num1] * 104 + index];
      }
   }
   
   // Initialize accumulators and the result.
   float sum_x, sum_y, sum_x2, sum_y2, sum_xy, coeff;
   sum_x = sum_y = sum_x2 = sum_y2 = sum_xy = coeff = 0.0f;

   // Compute the sums.
   for (int index = 0; index < 104; ++index) {
      float x = i_isolate[index];
      float y = j_isolate[index];

      sum_x += x;
      sum_y += y;
      sum_x2 += x * x;
      sum_y2 += y * y;
      sum_xy += x * y;
   }

   // Compute the Pearson coefficient using the "sometimes numerically
   // unstable" method because it's way more computationally efficient.
   coeff = (p * sum_xy - sum_x * sum_y) /
      sqrtf((p * sum_x2 - sum_x * sum_x) * (p * sum_y2 - sum_y * sum_y));

   // Dump it in the appropriate bucket. Buckets are allowed to overlap, so
   // we need to check all of them.
   for (uint32_t k = 0; k < num_ranges; k++) {
      float low = ranges[2 * k + 0];
      float high = ranges[2 * k + 1];
      if (coeff >= low && coeff < high) {
         uint32_t index = (tile_size * tile_size * k) +
            (tile_size * i) + j;
         buckets[index] += 1;
      }
   }
}
