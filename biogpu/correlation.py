import os
import sys
import numpy as np
import pycuda.autoinit
import pycuda.compiler
import pycuda.gpuarray
import pycuda.driver

def pearson(X, Y, ranges):
    # Check some preconditions to verify we're doing something sensical.
    # Doesn't cover all cases, but catches obvious mistakes.
    assert len(X[0]) == len(X[1]), 'Your sequences in X should all be the same length.'
    assert len(Y[0]) == len(Y[1]), 'Your sequences in Y should all be the same length.'
    assert len(X[0]) == len(Y[0]), 'Your sequences in X and Y should all be the same length.'

    n = len(X)
    m = len(Y)
    p = len(X[0])

    # Load the kernel and compile it.
    kernel_file = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                               'pearson.cu')
    f = open(kernel_file, 'r')
    kernel = pycuda.compiler.SourceModule(f.read())
    f.close()
    pearson_cuda = kernel.get_function('pearson')
    reduction_cuda = kernel.get_function('reduction')

    # CUDA parameters that seem to work well. The number of threads per tile
    # (the tile_size) should be a power of 2 for the parallel reduction to
    # work right!
    threads_per_block = 16
    blocks_per_tile = 64
    tile_size = threads_per_block * blocks_per_tile
    num_tiles = (n / tile_size + 1, m / tile_size + 1)

    # Copy the ranges into a numpy array.
    ranges_np = np.zeros(shape=(len(ranges), 2), dtype=np.float32, order='C')
    for i in range(len(ranges)):
        np.put(ranges_np[i], range(2), ranges[i])

    # Create a zero-initialized chunk of memory for the per-thread buckets and
    # copy it to the GPU.
    buckets = np.zeros(shape=(tile_size * tile_size * len(ranges), 1),
                       dtype=np.uint64, order='C')
    buckets_gpu = pycuda.gpuarray.to_gpu(buckets)

    # Do a kernel launch for each tile, copying the appropriate chunks of the
    # input arrays into X and Y for each launch.
    for s in range(num_tiles[0]):
        for t in range(num_tiles[1]):
            num_A = tile_size
            remain_X = n - (s * tile_size)
            num_A = num_A if num_A < remain_X else remain_X

            A = np.zeros(shape=(num_A, p), dtype=np.float32, order='C')
            for i in range(num_A):
                np.put(A[i], range(p), X[(s * tile_size) + i])

            num_B = tile_size
            remain_Y = m - (t * tile_size)
            num_B = num_B if num_B < remain_Y else remain_Y

            B = np.zeros(shape=(num_B, p), dtype=np.float32, order='C')
            for j in range(num_B):
                np.put(B[j], range(p), Y[(t * tile_size) + j])

            pearson_cuda(buckets_gpu.gpudata,
                         pycuda.driver.In(ranges_np), np.uint32(len(ranges)),
                         pycuda.driver.In(A), pycuda.driver.In(B),
                         np.uint32(tile_size), np.uint32(s), np.uint32(t),
                         np.uint32(n), np.uint32(m), np.uint32(p),
                         block=(threads_per_block, threads_per_block, 1),
                         grid=(blocks_per_tile, blocks_per_tile))

            progress = (s * num_tiles[1] + t) * 100.0 / (num_tiles[0] * num_tiles[1])
            sys.stdout.write('\rComputing correlations %.3f%%' % progress)
            sys.stdout.flush()

    print('\rComputing correlations 100.000%')
    sys.stdout.write('Merging buckets... ')
    sys.stdout.flush()

    # Do a parallel reduction to sum all the buckets element-wise.
    reduction_cuda(buckets_gpu.gpudata, np.uint32(len(ranges)),
                   np.uint32(tile_size), np.uint32(blocks_per_tile),
                   block=(threads_per_block, 1, 1),
                   grid=(tile_size, 1))

    # Copy buckets back from GPU.
    buckets_gpu.get(buckets)

    # Merge the results of the reduction from the first column of the matrix.
    merged = [0 for k in range(len(ranges))]
    for k in range(len(ranges)):
        for i in range(tile_size):
            bucket_index = (tile_size * tile_size * k) + (tile_size * i) + 0
            merged[k] += buckets[bucket_index]
    print('done.')

    return merged
