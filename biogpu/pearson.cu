#include <stdint.h>
#include <string.h>

const int kMaxAllelesPerIsolate = 7;
const int kMaxAlleles = 28;         // Maximum number of different alleles
const int kMaxAlleleLength = 104;   // Maximum length of each allele

__device__ int binomial_coefficient(int n, int k) {
   if (k < 0 || k > n)
      return 0;

   if (k > n - k)
      k = n - k;

   int c = 1;
   for (int i = 1; i <= k; ++i) {
      c *= n - k + i;
      c /= i;
   }

   return c;
}

// Gets isolate number |seq_num|, which is an index into a set of all isolates
// composed of |alleles_per_isolate| alleles, where each allele is one of 
// |num_alleles|. Order doesn't matter, and duplicates are allowed.
//
// Determining the isolate that corresponds to an index can be though of as
// traversing a tree similar to the following, where each node visited 
// represents an allele index. The branching factor is how many different 
// alleles there are, and the height is the number of alleles per isolate.
//
// The following drawing represents all possible 
// "order-doesn't-matter-with-replacement" isolates of length 3 (height = 3) 
// that can be generated from 4 isolates (branching factor = 4). The bottom
// line of numbers is the index of the isolate to be generated.
//
//
//                     _______________BEGIN________________________________
//                    /                        \                 \         \
//                   /                          \                 \         \
//      ____________0___________           ______1______          _2__       3
//     /         /      \       \         /       \     \        /    \      |
//    0          1        2      3       1        2      3      2      3     3
//  // \\       /|\      / \     |      /|\      / \     |     / \     |     |
// 0 1 2 3     1 2 3     2 3     3     1 2 3     2 3     3     2 3     3     3
//
// 0 1 2 3     4 5 6     7 8     9    ...   ...   ...   ...   ...   ...     19
// 
//
// Note that at each stage of tree traversal, some range of indices corresponds
// to each possible branch that can be traversed. At BEGIN, indices 0-9
// correspond to allele 0, indices 10-15 correspond to allele 1, and so on.
//
// The following table is from the perspective of BEGIN (before any traversal),
// and each number represents the number of isolates that correspond to that
// particular allele. 
//
// n=5   h=1   |   1     1     1     1     1
// n=5   h=2   |   5     4     3     2     1
// n=5   h=3   |   15    10    6     3     1
// n=5   h=4   |   35    20    10    4     1
// n=5   h=5   |   70    35    15    5     1
//
// Now, tilt your head 45 degrees to the right and look at the top right 
// number. Do you see it? Pascal's mother fucking Triangle. The number at row r, 
// column c of Pascal's Triangle is r choose c.
//
// So, given some input isolate index, we simply ask the question, "is this
// index <70?" If so, the first allele is 0. If not, add 35 to 70 (because 35
// is simply the number of indices that allele 1 maps to, not their values)
// and repeat. Once an allele is found, the tree is traversed and its height is,
// effectively reduced by one.
//
// Let's look at Pascal's Triangle and try to determine a mapping from this
// rectangle onto the triangle.
//
//                              /
//                           /  1  \                           0
//                        /  1     1  \                        1
//                     /  1     2     1  \                     2
//                  /  1     3     3     1  \                  3
//                \ 1     4     6     4     1 /                4
//               1   \ 5    10     10    5 /   1               5
//            1     6   \15     20    15/   6     1            6
//         1     7    21   \35     35/   21    7     1         7
//      1     8    28    56   \ 70/    56    28    8     1     8
//
// Let r = the row of the rectangle view, and let c = the column of the 
// rectangle view, and n = the number of different alleles (which equals the
// width of the rectangle).
//
// The Pascal's Triangle row that corresponds to the number at (r,c) is
// n - c - 1 + r. The Pascal's Triangle column is simply r.
__device__ void get_isolate(int seq_num, 
                            uint8_t* alleles, 
                            // Shrinks every time we take a right turn.
                            uint8_t num_alleles,
                            uint8_t alleles_per_isolate) {
   uint8_t orig_num_alleles = num_alleles;
   uint8_t cur_row = alleles_per_isolate - 1;
   uint8_t cur_col = 0;

   // Init num to biggest "bucket"
   int num = binomial_coefficient(num_alleles - cur_col - 1 + cur_row, cur_row);

   // Save in order to subtract from the "current" sequence number
   int old_num = 0;

   // Output |alleles_per_isolate| alleles
   int alleles_output = 0;
   while (1) {
      // Allele determined
      if (seq_num < num) {
         alleles[alleles_output] = cur_col + orig_num_alleles - num_alleles;
         if (++alleles_output >= alleles_per_isolate)
            break;

         seq_num -= old_num;
         old_num = 0;
         
         // New parameters for new location in tree
         num_alleles -= cur_col;
         cur_row--;

         num = binomial_coefficient(num_alleles - 1 + cur_row, cur_row);

         cur_col = 0;
      } 
      
      // Allele not in this "bucket" -- try the next one
      else {
         cur_col++;
         old_num = num;
         num += binomial_coefficient(num_alleles - 1 - cur_col + cur_row, cur_row);
      }
   }
}

__constant__ uint8_t alleles[kMaxAlleles * kMaxAlleleLength];

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
   uint32_t i = blockIdx.y * blockDim.y + threadIdx.y; // row
   uint32_t j = blockIdx.x * blockDim.x + threadIdx.x; // column

   // Calculate the absolute <i, j> coords within the matrix.
   uint32_t i_abs = tile_row * tile_size + i;
   uint32_t j_abs = tile_col * tile_size + j;

   // Only compute values inside the bounds of the matrix.
   if (i_abs >= num_isolates || j_abs >= num_isolates)
      return;

   // We don't want to compare isolates with themselves, or any comparisons
   // of a lower-numberes isolate to a higher-numbered one. Each pair of 
   // isolates (order doesn't matter) will only be compared once. This will
   // cause divergence only in the warps that lie along the main diagonal
   // of the comparison matrix.
   if (i_abs <= j_abs)
      return;

   // Generate isolate |i_abs| and |j_abs|
   uint8_t i_allele_indices[kMaxAllelesPerIsolate]; 
   uint8_t j_allele_indices[kMaxAllelesPerIsolate]; 
   get_isolate(i_abs, i_allele_indices, num_alleles, alleles_per_isolate);
   get_isolate(j_abs, j_allele_indices, num_alleles, alleles_per_isolate);

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
