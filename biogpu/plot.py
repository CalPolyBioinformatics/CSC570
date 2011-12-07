import matplotlib.pyplot
import matplotlib.mlab

# Takes bucket output and writes a graph of the results
def correlation_ranges(buckets, output_file='correlation_ranges.png',
                       title='Correlation Ranges'):
    # TODO: redo all of this

    # merge into 6 bins as follows:
    #   0: 0.997 -> 1.000
    #   1: 0.990 -> 0.997
    #   2: 0.970 -> 0.990
    #   3: 0.900 -> 0.970
    #   4: 0.700 -> 0.900
    #   5: 0.000 -> 0.700

    bins = [0, 0, 0, 0, 0, 0]

    num_buckets = len(buckets)
    for i in range(num_buckets):
        current = i / float(num_buckets)
        if current > 0.997:
            bins[0] += buckets.item(i)
        elif current > 0.990:
            bins[1] += buckets.item(i)
        elif current > 0.970:
            bins[2] += buckets.item(i)
        elif current > 0.900:
            bins[3] += buckets.item(i)
        elif current > 0.700:
            bins[4] += buckets.item(i)
        else:
            bins[5] += buckets.item(i)
    
    bins.append(0) # visual spacer
    bins.reverse()

    fig = matplotlib.pyplot.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.bar(range(len(bins)), bins, align='center')
    ax.set_title(title)
    ax.set_ylabel('Number of Elements')
    ax.set_xticks(range(len(bins)))
    ax.set_xticklabels(['', '0% - 70%', '70% - 90%', '90% - 97%', '97% - 99%',
                        '99% - 99.7%', '99.7% - 100%'])
    ax.grid(True)
    fig.autofmt_xdate()
    print('Saving ' + output_file)
    fig.savefig(output_file)
