import pycuda.autoinit
import pycuda.driver as cuda
from pycuda.compiler import SourceModule
import numpy
from scipy.stats.stats import pearsonr

#mod = SourceModule("""
#__global__ void multiply_them(float *dest, float *a, float *b) {
#    const int i = threadIdx.x;
#    dest[i] = a[i] * b[i];
#}
#""")
#
#multiply_them = mod.get_function("multiply_them")
#
#a = numpy.random.randn(400).astype(numpy.float32)
#b = numpy.random.randn(400).astype(numpy.float32)
#
#dest = numpy.zeros_like(a)
#
#multiply_them(drv.Out(dest), drv.In(a), drv.In(b), block=(400, 1, 1), grid=(1, 1))
#
#print dest - a * b

##########################################

# Source data is an n x m matrix, where n = number of arrays and m = length of
# each array.

n = 6
m = 10

pearson_input = numpy.ndarray(shape=(n, m), dtype=numpy.int32, order='C')
for i in range(n):
    numpy.put(pearson_input[i], range(m),
              numpy.random.random_integers(0, 30, 10).astype(numpy.int32))

print('pearson_input:')
print(pearson_input)
print('\n');

pearson_matrix = numpy.zeros(shape=(n, n), dtype=numpy.float32, order='C')

print('pearson_matrix (empty):');
print(pearson_matrix);
print('\n');

for i in range(n):
    for j in range(n):
        coeff, _ = pearsonr(pearson_input[i], pearson_input[j])
        pearson_matrix[i][j] = coeff

print('pearson_matrix (computed):');
print(pearson_matrix);
print('\n');

kernel = SourceModule('''
/*
 * This kernel makes the following assumptions:
 *
 *  -> Input arrays for each sequence (of length SEQ_LEN) are a multiple of
 *     BLOCK_SIZE. Just pad the end with zeros, but still pass in the ACTUAL
 *     length in SEQ_LEN.
 */

#define BLOCK_SIZE 2

__global__ void pearson(float *corr_mat, int *global_A, int *global_B, int TILE_SIZE, int SEQ_LEN) {
    // keep the working set for each group of threads in shared memory
    __shared__ float shared_A[BLOCK_SIZE][BLOCK_SIZE];
    __shared__ float shared_B[BLOCK_SIZE][BLOCK_SIZE];

    int i, j, k;
    int offset;
    float a, b;
    float sum_a, sum_b, sum_a2, sum_b2, sum_ab, coeff;

    // calculate absolute <i, j> coords within the tile
    i = (blockIdx.y * blockDim.y) + threadIdx.y;
    j = (blockIdx.x * blockDim.x) + threadIdx.x;

    // initialize the accumulators
    sum_a = sum_b = sum_a2 = sum_b2 = sum_ab = 0.0f;

    // plow through the sequence in chunks of BLOCK_SIZE x BLOCK_SIZE threads
    for (offset = 0; offset < SEQ_LEN; offset += blockDim.x) {
        // copy the working set into shared memory
        shared_A[threadIdx.y][threadIdx.x] = global_A[(i * SEQ_LEN) + offset + threadIdx.x];
        shared_B[threadIdx.y][threadIdx.x] = global_B[(j * SEQ_LEN) + offset + threadIdx.x];
        
        // memory barrier
        __syncthreads();

        // accumulate the chunk
        for (k = 0; k < blockDim.x; k++) {
            a = shared_A[threadIdx.y][k];
            b = shared_B[threadIdx.x][k];
            sum_a += a;
            sum_b += b;
            sum_a2 += a * a;
            sum_b2 += b * b;
            sum_ab += a * b;
        }

        // memory barrier
        __syncthreads();
    }

    // compute the pearson correlation coefficient
    coeff = (float)((SEQ_LEN * sum_ab) - (sum_a * sum_b)) /
            (float)sqrtf(((SEQ_LEN * sum_a2) - (sum_a * sum_a)) * ((SEQ_LEN * sum_b2) - (sum_b * sum_b)));
    corr_mat[(i * TILE_SIZE) + j] = coeff;
}
''')

pearson_cuda_input = numpy.ndarray(shape=(n, m), dtype=numpy.int32, order='C')
for i in range(n):
    numpy.put(pearson_input[i], range(m),
              numpy.random.random_integers(0, 30, 10).astype(numpy.int32))

print('pearson_cuda_input:')
print(pearson_cuda_input)
print('\n');

pearson_cuda_matrix = numpy.zeros(shape=(n, n), dtype=numpy.float32, order='C')
print('pearson_cuda_matrix (empty):')
print(pearson_cuda_matrix)
print('\n');

pearson_cuda = kernel.get_function("pearson")
pearson_cuda(cuda.Out(pearson_cuda_matrix),
             cuda.In(pearson_cuda_input), cuda.In(pearson_cuda_input),
             numpy.int32(n), numpy.int32(m),
             block=(2, 2, 1), grid=(3, 3))

print('pearson_cuda_matrix (computed):');
print(pearson_cuda_matrix);
print('\n');
