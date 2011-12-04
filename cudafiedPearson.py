import pycuda.autoinit
import pycuda.driver as cuda
from pycuda.compiler import SourceModule
import numpy
from scipy.stats.stats import pearsonr

# Source data is an n x m matrix, where n = number of arrays and m = length of
# each array.

n = 6
m = 10

pearson_input = numpy.zeros(shape=(n, m), dtype=numpy.int32, order='C')
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
    __global__ void pearson(float *corr_mat, int *input, int disp_len) {
        // calculate absolute <i, j> coords within the matrix
        int i = blockIdx.y * blockDim.y + threadIdx.y; // row
        int j = blockIdx.x * blockDim.x + threadIdx.x; // column

        // calculate the size of the matrix
        int n = gridDim.x * blockDim.x;

        // initialize accumulators and result
        float sum_x, sum_y, sum_x2, sum_y2, sum_xy, coeff;
        sum_x = sum_y = sum_x2 = sum_y2 = sum_xy = coeff = 0.0f;

        // compute the sums
        for (int k = 0; k < disp_len; k++) {
            int x = input[i * disp_len + k];
            int y = input[j * disp_len + k];

            sum_x += x;
            sum_y += y;
            sum_x2 += x * x;
            sum_y2 += y * y;
            sum_xy += x * y;
        }

        // compute the pearson coefficient using the "sometimes numerically
        // unstable" because it's waaaay more computationally efficient
        coeff = (disp_len * sum_xy - sum_x * sum_y) /
                sqrtf((disp_len * sum_x2 - sum_x * sum_x) * (disp_len * sum_y2 - sum_y * sum_y));
        corr_mat[i * n + j] = coeff;

    }
''')

pearson_cuda_input = numpy.zeros(shape=(n, m), dtype=numpy.int32, order='C')
for i in range(n):
    numpy.put(pearson_cuda_input[i], range(m), pearson_input[i])

print('pearson_cuda_input:')
print(pearson_cuda_input)
print('\n');

pearson_cuda_matrix = numpy.zeros(shape=(n, n), dtype=numpy.float32, order='C')
print('pearson_cuda_matrix (empty):')
print(pearson_cuda_matrix)
print('\n');

pearson_cuda = kernel.get_function("pearson")
pearson_cuda(cuda.Out(pearson_cuda_matrix), cuda.In(pearson_cuda_input),
             numpy.int32(m),
             block=(2, 2, 1), grid=(3, 3))

print('pearson_cuda_matrix (computed):');
print(pearson_cuda_matrix);
print('\n');

print('comparison: ');
print((pearson_matrix - pearson_cuda_matrix).round(6));
