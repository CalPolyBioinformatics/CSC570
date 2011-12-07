import sys
import numpy
from scipy.stats.stats import pearsonr
import biogpu.correlation
import time

def main():
    n = 512 # number of pyroprints
    m = 104 # pyroprint length
    #n = 10
    #m = 104

    pyroprints = numpy.zeros(shape=(n, m), dtype=numpy.float32, order='C')
    for i in range(n):
        numpy.put(pyroprints[i], range(m),
                  numpy.random.rand(m).astype(numpy.float32))

    print('Fake Pyroprints:')
    print(pyroprints)
    print('')

    print('=== Computing with Python/SciPy ===')
    python_start = time.time()
    python_buckets = compute_python(pyroprints, 10000)
    python_end = time.time()

    #print('Buckets (abridged):')
    #for i in range(10000):
    #    if python_buckets[i] > 0:
    #        print('\t[%d] = %d' % (i, python_buckets[i]))
    #print('\n')

    python_time = python_end - python_start
    print('Computed in %f seconds.\n' % python_time)

    print('=== Computing with CUDA ===')
    cuda_start = time.time()
    cuda_buckets = biogpu.correlation.pearson(pyroprints, 10000)
    cuda_end = time.time()

    #print('Buckets (abridged):')
    #for i in range(10000):
    #    if cuda_buckets[i] > 0:
    #        print('\t[%d] = %d' % (i, cuda_buckets[i]))
    #print('\n')

    cuda_time = cuda_end - cuda_start
    print('Computed in %f seconds.\n' % cuda_time)

    speedup = python_time / cuda_time
    print('Speedup of %.2fx.\n' % speedup)

    print('Done.')

def compute_python(pyroprints, num_buckets):
    n = len(pyroprints)
    m = len(pyroprints[0])

    matrix = numpy.zeros(shape=(n, n), dtype=numpy.float32, order='C')
    buckets = numpy.zeros(shape=(num_buckets, 1), dtype=numpy.int32, order='C')

    num_cells = n * n * 0.5 - n
    cell_count = 0
    for i in range(n):
        for j in range(n):
            if i >= j:
                continue
            coeff, _ = pearsonr(pyroprints[i], pyroprints[j])
            matrix[i][j] = coeff
            bucket = int(coeff * num_buckets)
            if bucket >= num_buckets:
                buckets[num_buckets - 1] += 1
            elif bucket >= 1:
                buckets[bucket - 1] += 1
            cell_count += 1

        progress = cell_count / num_cells * 100
        sys.stdout.write('\rComputing correlations %.3f%%' % progress)
        sys.stdout.flush()

    print('\rComputing correlations 100.000%')
    return buckets

if __name__ == '__main__':
    main()
