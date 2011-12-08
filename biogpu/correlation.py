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
    num_tiles = (len(X) / tile_size + 1, len(Y) / tile_size + 1)

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
    for s in range(num_tiles[1]):
        for t in range(num_tiles[0]):
            num_A = tile_size
            remain_Y = len(Y) - (s * tile_size)
            num_A = num_A if num_A < remain_Y else remain_Y

            A = np.zeros(shape=(num_A, len(Y[0])), dtype=np.float32, order='C')
            for i in range(num_A):
                np.put(A[i], range(len(Y[0])), Y[(s * tile_size) + i])

            num_B = tile_size
            remain_X = len(X) - (t * tile_size)
            num_B = num_B if num_B < remain_X else remain_X

            B = np.zeros(shape=(num_B, len(X[0])), dtype=np.float32, order='C')
            for i in range(num_B):
                np.put(B[i], range(len(X[0])), X[(t * tile_size) + i])

            pearson_cuda(buckets_gpu.gpudata,
                         pycuda.driver.In(ranges_np), np.uint32(len(ranges)),
                         pycuda.driver.In(A), pycuda.driver.In(B),
                         np.uint32(tile_size), np.uint32(s), np.uint32(t),
                         np.uint32(len(X)), np.uint32(len(Y)), np.uint32(len(X[0])),
                         block=(threads_per_block, threads_per_block, 1),
                         grid=(blocks_per_tile, blocks_per_tile))

            progress = (s * num_tiles[1] + t) * 100.0 / (num_tiles[0] * num_tiles[1])
            sys.stdout.write('\rComputing correlations %.3f%%' % progress)
            sys.stdout.flush()

    print('\rComputing correlations 100.000%')
    sys.stdout.write('Merging buckets... ')
    sys.stdout.flush()
            
    reduction_cuda(buckets_gpu.gpudata, np.uint32(len(ranges)),
                   np.uint32(tile_size),
                   block=(threads_per_block, threads_per_block, 1),
                   grid=(blocks_per_tile, blocks_per_tile))

    # Copy buckets back from GPU.
    buckets_gpu.get(buckets)

    # Copy the results of the reduction out of cell <0, 0>.
    merged = []
    for i in range(len(ranges)):
        bucket_index = (tile_size * tile_size * i) + (tile_size * 0) + 0
        merged.append(buckets[bucket_index])
    print('done.')

    return merged
