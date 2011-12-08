#include <stdint.h>

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
                          uint32_t tile_size) {
    // Calculate <i, j> coords within the tile.
    uint32_t i = blockIdx.y * blockDim.y + threadIdx.y; // row
    uint32_t j = blockIdx.x * blockDim.x + threadIdx.x; // column

    // Start with 2x2 reductions.
    uint32_t size = 2;

    // Keep reducing until everything has been reduced into cell <0, 0>.
    while (size <= tile_size) {
        if (i % size == 0 && j % size == 0) {
            dump_bucket(buckets, num_ranges, tile_size,
                        i + size / 2, j, i, j);
            dump_bucket(buckets, num_ranges, tile_size,
                        i, j + size / 2, i, j);
            dump_bucket(buckets, num_ranges, tile_size,
                        i + size / 2, j + size / 2, i, j);
        }

        size *= 2;
    }
    
    // The results are in bucket <0, 0>!
}

__global__ void pearson(uint64_t *buckets,
                        float *ranges, uint32_t num_ranges,
                        float *A, float *B,
                        uint32_t tile_size, uint32_t s, uint32_t t,
                        uint32_t n, uint32_t m, uint32_t p) {
    // Calculate relative <i, j> coords within this tile.
    uint32_t i = blockIdx.y * blockDim.y + threadIdx.y; // row
    uint32_t j = blockIdx.x * blockDim.x + threadIdx.x; // column

    // Calculate the offsets based on the tile number.
    uint32_t i_offset = s * tile_size;
    uint32_t j_offset = t * tile_size;

    // Calculate the absolute <i, j> coords within the matrix.
    uint32_t i_abs = i_offset + i;
    uint32_t j_abs = j_offset + j;

    // Quick checks to bail out. Only compute values inside the bounds of the
    // matrix, and above the diagonal.
    if (i_abs >= n || j_abs >= n || i_abs >= j_abs) {
        return;
    }

    // Initialize accumulators and the result.
    float sum_x, sum_y, sum_x2, sum_y2, sum_xy, coeff;
    sum_x = sum_y = sum_x2 = sum_y2 = sum_xy = coeff = 0.0f;

    // Compute the sums.
    for (uint32_t k = 0; k < p; k++) {
        float x = A[i * m + k];
        float y = B[j * m + k];

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

    // Dump it in the appropriate bucket.
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
