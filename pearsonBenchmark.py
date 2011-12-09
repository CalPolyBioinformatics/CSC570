import sys
import numpy as np
from scipy.stats.stats import pearsonr
import biogpu.correlation
import time

def main():
    n = 2121 # length of X
    m = 1717 # length of Y
    p = 104 # pyroprint length
    ranges = [(0.00, 0.15),
              (0.15, 0.30),
              (0.30, 0.45),
              (0.45, 0.60),
              (0.60, 0.75),
              (0.75, 0.90),
              (0.90, 1.0)]

    A = np.zeros(shape=(n, p), dtype=np.float32, order='C')
    B = np.zeros(shape=(m, p), dtype=np.float32, order='C')
    for i in range(n):
        np.put(A[i], range(p), np.random.rand(p).astype(np.float32))

    for i in range(m):
        np.put(B[i], range(p), np.random.rand(p).astype(np.float32))

    print('X:')
    print(A)
    print('\nY:')
    print(B)
    print('\n')

    print('=== Computing with Python/SciPy ===')
    python_start = time.time()
    python_buckets = compute_python(A, B, ranges)
    python_end = time.time()

    print('Buckets:')
    for i in range(len(python_buckets)):
        print('\t[%d] (%.2f, %.2f) = %d' % (i, ranges[i][0], ranges[i][1], python_buckets[i]))
    print('\n')

    python_time = python_end - python_start
    print('Computed in %f seconds.\n' % python_time)

    print('=== Computing with CUDA ===')
    cuda_start = time.time()
    cuda_buckets = biogpu.correlation.pearson(A, B, ranges)
    cuda_end = time.time()

    print('Buckets:')
    for i in range(len(cuda_buckets)):
        print('\t[%d] (%.2f, %.2f) = %d' % (i, ranges[i][0], ranges[i][1], cuda_buckets[i]))
    print('\n')

    cuda_time = cuda_end - cuda_start
    print('Computed in %f seconds.\n' % cuda_time)

    speedup = python_time / cuda_time
    print('Speedup of %.2fx.\n' % speedup)

    print('Done.')

def compute_python(X, Y, ranges):
    # Check some preconditions to verify we're doing something sensical.
    # Doesn't cover all cases, but catches obvious mistakes.
    assert len(X[0]) == len(X[1]), 'Your sequences in X should all be the same length.'
    assert len(Y[0]) == len(Y[1]), 'Your sequences in Y should all be the same length.'
    assert len(X[0]) == len(Y[0]), 'Your sequences in X and Y should all be the same length.'

    n = len(X)
    m = len(Y)
    p = len(X[0])

    # Create the correlation matrix and buckets.
    matrix = np.zeros(shape=(n, m), dtype=np.float32, order='C')
    buckets = np.zeros(shape=(len(ranges), 1), dtype=np.uint64, order='C')

    for i in range(n):
        for j in range(m):
            # Compute the coefficient.
            coeff, _ = pearsonr(X[i], Y[j])
            matrix[i][j] = coeff

            # Drop it in appropriate bucket.
            for k in range(len(ranges)):
                low, high = ranges[k]
                if coeff >= low and coeff < high:
                    buckets[k] += 1

        progress = (i * m) * 100.0 / (n * m)
        sys.stdout.write('\rComputing correlations %.3f%%' % progress)
        sys.stdout.flush()

    print('\rComputing correlations 100.000%')
    return buckets

if __name__ == '__main__':
    main()
