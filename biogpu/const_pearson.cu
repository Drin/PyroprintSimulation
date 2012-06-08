#include <stdint.h>
#include <string.h>

const int kMaxAllelesPerIsolate = 7;
const int kMaxAlleles = 28;         // Maximum number of different alleles
const int kMaxAlleleLength = 104;   // Maximum length of each allele

__constant__ uint8_t alleles[kMaxAlleles * kMaxAlleleLength];

__global__ void pearson(uint64_t *buckets,
                        float *ranges, 
                        uint32_t num_ranges,
                        uint32_t tile_size, 
                        uint32_t tile_row, 
                        uint32_t tile_col, 
                        uint8_t num_alleles,         
                        uint8_t alleles_per_isolate, 
                        // Total number of isolates, dependent on num_alleles
                        // and alleles_per_isolate
                        uint32_t num_isolates,        
                        // Number of nucleotides per allele
                        uint32_t length_alleles) {
   // Calculate relative <i, j> coords within this tile.
   uint32_t i = blockIdx.x * blockDim.x + threadIdx.x; // column
   uint32_t j = blockIdx.y * blockDim.y + threadIdx.y; // row

   // Only compute values inside the bounds of the matrix.
   // And we don't want to compare isolates with themselves, or any comparisons
   // of a lower-numberes isolate to a higher-numbered one. Each pair of 
   // isolates (order doesn't matter) will only be compared once. This will
   // cause divergence only in the warps that lie along the main diagonal
   // of the comparison matrix.
   if (tile_col * tile_size + i >= num_isolates ||
       tile_row * tile_size + j >= num_isolates ||
       tile_col * tile_size + i <= tile_row * tile_size + j)
      return;


   // Generate isolate |i_abs| and |j_abs|
   uint8_t i_allele_indices[kMaxAllelesPerIsolate]; 
   uint8_t j_allele_indices[kMaxAllelesPerIsolate]; 
   get_isolate(tile_col * tile_size + i, i_allele_indices, num_alleles, alleles_per_isolate);
   get_isolate(tile_row * tile_size + j, j_allele_indices, num_alleles, alleles_per_isolate);

   // Initialize accumulators and the result.
   float sum_x = 0, sum_y = 0, sum_x2 = 0, sum_y2 = 0, sum_xy = 0;

   // Compute the sums.
   for (int index = 0; index < length_alleles; ++index) {
      uint16_t x = 0, y = 0;

      for (int alleleNdx = 0; alleleNdx < alleles_per_isolate; alleleNdx++) {
         x += alleles[i_allele_indices[alleleNdx] * length_alleles + index];
         y += alleles[j_allele_indices[alleleNdx] * length_alleles + index];
      }

      sum_x += x;
      sum_y += y;
      sum_x2 += x * x;
      sum_y2 += y * y;
      sum_xy += x * y;
   }

   // Compute the Pearson coefficient using the "sometimes numerically
   // unstable" method because it's way more computationally efficient.
   float coeff = (length_alleles * sum_xy - sum_x * sum_y) /
      sqrtf((length_alleles * sum_x2 - sum_x * sum_x) * 
            (length_alleles * sum_y2 - sum_y * sum_y));

   // Dump it in the appropriate bucket. 
   // Below is a commented-out comment that no longer applies. To re-implement
   // this feature, remove the break statement below.
   // //Buckets are allowed to overlap, so we need to check all of them.
   for (uint32_t k = 0; k < num_ranges; k++) {
      //float low = ranges[2 * k + 0];
      //float high = ranges[2 * k + 1];
      if (coeff >= ranges[2 * k] && coeff < ranges[2 * k + 1]) {
         uint32_t index = (tile_size * tile_size * k) +
            (tile_size * i) + j;
         buckets[index]++;
         //break;
      }
   }
}
