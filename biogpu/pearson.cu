__global__ void pearson(int *buckets, int num_buckets,
                        int *A, int num_A, int *B, int num_B,
                        int s, int t, int n, int m) {
    // calculate relative <i, j> coords within this tile
    int i = blockIdx.y * blockDim.y + threadIdx.y; // row
    int j = blockIdx.x * blockDim.x + threadIdx.x; // column

    // calculate the offsets based on the tile number
    int i_offset = s * gridDim.y * blockDim.y;
    int j_offset = t * gridDim.x * blockDim.x;

    // make sure this thread is inside the matrix
    if (i + i_offset >= n ||
        j + j_offset >= n) {
        return;
    }

    // initialize accumulators and result
    float sum_x, sum_y, sum_x2, sum_y2, sum_xy, coeff;
    sum_x = sum_y = sum_x2 = sum_y2 = sum_xy = coeff = 0.0f;

    // compute the sums
    for (int k = 0; k < m; k++) {
        int x = A[i * m + k];
        int y = B[j * m + k];

        sum_x += x;
        sum_y += y;
        sum_x2 += x * x;
        sum_y2 += y * y;
        sum_xy += x * y;
    }

    // compute the pearson coefficient using the "sometimes numerically
    // unstable" method because it's way more computationally efficient
    coeff = (m * sum_xy - sum_x * sum_y) /
            sqrtf((m * sum_x2 - sum_x * sum_x) * (m * sum_y2 - sum_y * sum_y));

    // dump it in the appropriate bucket
    int bucket = (int)(coeff * num_buckets);
    if (bucket >= num_buckets) {
        atomicAdd(&(buckets[num_buckets - 1]), 1);
    } else if (bucket >= 1) {
        atomicAdd(&(buckets[bucket - 1]), 1);
    }
}
