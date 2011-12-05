import pycuda.autoinit
import pycuda.driver as cuda
import pycuda.compiler
import pycuda.gpuarray
import numpy
from scipy.stats.stats import pearsonr

def main():
    n = 6 # number of pyroprints
    m = 10 # pyroprint length

    pyroprints = numpy.zeros(shape=(n, m), dtype=numpy.int32, order='C')
    for i in range(n):
        numpy.put(pyroprints[i], range(m),
                  numpy.random.random_integers(0, 30, 10).astype(numpy.int32))

    print('fake pyroprints:')
    print(pyroprints)
    print('\n');

    print('=== computing with python ===')
    python_buckets = compute_python(pyroprints, 10000)
    print('buckets (abridged):')
    for i in range(10000):
        if python_buckets[i] > 0:
            print('\t[%d] = %d' % (i, python_buckets[i]))
    print('\n')

    print('=== computing with cuda ===')
    cuda_buckets = compute_cuda(pyroprints, 10000)
    print('buckets (abridged):')
    for i in range(10000):
        if cuda_buckets[i] > 0:
            print('\t[%d] = %d' % (i, cuda_buckets[i]))
    print('\n')

    print('done')

def compute_python(pyroprints, num_buckets):
    n = len(pyroprints)
    m = len(pyroprints[0])

    matrix = numpy.zeros(shape=(n, n), dtype=numpy.float32, order='C')
    buckets = numpy.zeros(shape=(num_buckets, 1), dtype=numpy.int32, order='C')

    for i in range(n):
        for j in range(n):
            coeff, _ = pearsonr(pyroprints[i], pyroprints[j])
            matrix[i][j] = coeff
            bucket = int(coeff * num_buckets)
            if bucket >= num_buckets:
                buckets[num_buckets - 1] += 1
            elif bucket >= 1:
                buckets[bucket - 1] += 1
    
    return buckets


def compute_cuda(pyroprints, num_buckets):
    kernel = pycuda.compiler.SourceModule('''
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
            // unstable" because it's waaaay more computationally efficient
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
    ''')
    pearson_kernel = kernel.get_function('pearson')

    n = len(pyroprints)
    m = len(pyroprints[0])
    
    block_size = 2
    tile_size = 2
    num_tiles = (n / (tile_size * block_size)) + 1

    buckets = numpy.zeros(shape=(num_buckets, 1), dtype=numpy.int32, order='C')
    buckets_gpu = pycuda.gpuarray.to_gpu(buckets)

    for s in range(num_tiles):
        for t in range(num_tiles):
            num_A = tile_size * block_size
            remain_A = n - (s * tile_size * block_size)
            num_A = num_A if num_A < remain_A else remain_A

            A = numpy.zeros(shape=(num_A, m), dtype=numpy.int32, order='C')
            for i in range(num_A):
                numpy.put(A[i], range(m), pyroprints[(s * tile_size * block_size) + i])

            num_B = tile_size * block_size
            remain_B = n - (t * tile_size * block_size)
            num_B = num_B if num_B < remain_B else remain_B

            B = numpy.zeros(shape=(num_B, m), dtype=numpy.int32, order='C')
            for i in range(num_B):
                numpy.put(B[i], range(m), pyroprints[(t * tile_size * block_size) + i])

            pearson_kernel(buckets_gpu.gpudata, numpy.int32(num_buckets),
                           cuda.In(A), numpy.int32(num_A),
                           cuda.In(B), numpy.int32(num_B),
                           numpy.int32(s), numpy.int32(t),
                           numpy.int32(n), numpy.int32(m),
                           block=(block_size, block_size, 1),
                           grid=(tile_size, tile_size))

    buckets_gpu.get(buckets)
    return buckets

if __name__ == '__main__':
    main()
