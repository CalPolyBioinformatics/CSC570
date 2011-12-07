#include <stdint.h>

__global__ void pearson(int *buckets, int num_buckets,
                        float *A, int num_A, float *B, int num_B,
                        int s, int t, int n, int m) {
    // Calculate relative <i, j> coords within this tile.
    uint32_t i = blockIdx.y * blockDim.y + threadIdx.y; // row
    uint32_t j = blockIdx.x * blockDim.x + threadIdx.x; // column

    // Calculate the offsets based on the tile number.
    uint32_t i_offset = s * gridDim.y * blockDim.y;
    uint32_t j_offset = t * gridDim.x * blockDim.x;

    // Calculate the absolute <i, j> coords within the matrix.
    uint64_t i_abs = i_offset + i;
    uint64_t j_abs = j_offset + j;

    // Quick checks to bail out. Only compute values inside the bounds of the
    // matrix, and above the diagonal.
    if (i_abs >= n || j_abs >= n || i_abs >= j_abs) {
        return;
    }

    // Initialize accumulators and the result.
    float sum_x, sum_y, sum_x2, sum_y2, sum_xy, coeff;
    sum_x = sum_y = sum_x2 = sum_y2 = sum_xy = coeff = 0.0f;

    // Compute the sums.
    for (int k = 0; k < m; k++) {
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
    coeff = (m * sum_xy - sum_x * sum_y) /
            sqrtf((m * sum_x2 - sum_x * sum_x) * (m * sum_y2 - sum_y * sum_y));

    // Dump it in the appropriate bucket.
    int bucket = (int)(coeff * num_buckets);
    if (bucket >= num_buckets) {
        atomicAdd(&(buckets[num_buckets - 1]), 1);
    } else if (bucket >= 1) {
        atomicAdd(&(buckets[bucket - 1]), 1);
    }
}
