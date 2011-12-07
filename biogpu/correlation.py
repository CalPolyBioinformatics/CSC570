import os
import sys
import numpy
import pycuda.autoinit
import pycuda.compiler
import pycuda.gpuarray
import pycuda.driver

# This expects the pyroprints argument to be a list of lists. It computes the
# Pearson correlation against itself. It also assumes each sub-list is the same
# length as every other sublist. Returns the list of buckets.
def pearson(pyroprints, num_buckets, show_progress = True):
    # We expect the kernel to be located in the package directory.
    kernel_file = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                               'pearson.cu')
    f = open(kernel_file, 'r')
    kernel = pycuda.compiler.SourceModule(f.read())
    f.close()
    pearson_kernel = kernel.get_function('pearson')

    n = len(pyroprints)
    m = len(pyroprints[0])
    
    block_size = 16
    tile_size = 64
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
                           pycuda.driver.In(A), numpy.int32(num_A),
                           pycuda.driver.In(B), numpy.int32(num_B),
                           numpy.int32(s), numpy.int32(t),
                           numpy.int32(n), numpy.int32(m),
                           block=(block_size, block_size, 1),
                           grid=(tile_size, tile_size))

            if show_progress:
                progress = (s * num_tiles + t) * 100.0 / (num_tiles * num_tiles)
                sys.stdout.write('\rComputing correlations %.3f%%' % progress)
                sys.stdout.flush()

    buckets_gpu.get(buckets)
    if show_progress:
        print('\rComputing correlations 100.000%')

    return buckets
