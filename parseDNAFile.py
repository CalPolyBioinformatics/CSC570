#!/usr/bin/python -tt
# Copyright 2010 Google Inc.
# Licensed under the Apache License, Version 2.0
# http://www.apache.org/licenses/LICENSE-2.0


"""
  So far we have done:
  Find all files with DNA info in DIR
  Parse all files and find all DNA Chunks
  in each chunk find the primer sequence (if it exists)
  store all the 104char segments after the primer
  find all unique sequences of 104char segments
  find all possible choose 7 combinations of unique sequences
  generate pyroprints (python)
  compare pyroprints
  Graph pyroprints
  TODO:
  generate pyroprints (CUDA) -- Bob: I'm not sure it makes a whole lot of sense
    to spend time writing a kernel for this. Since we're only talking about
    pyroprinting ~ 20k-50k strains it only takes like, a minute tops.
  calculate pearson correl (CUDA) -- Bob: This makes a lot of sense. The matrix
    is 50k x 50k and takes lightyears in Python right now.
  """

import sys
import os
import itertools
from scipy.stats.stats import pearsonr
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import math
import re 
from optparse import OptionParser
import ConfigParser
# detect gpu support
gpu_support = False
try:
    import biogpu.correlation
    gpu_support = True
except:
   pass

###PRIMER
#16-23 primer GGAACCTGCGGTTGGATCAC
#23-5 primer CGTGAGGCTTAACCTT

##Dispensation Sequence:
#23-5 AACACGCGA23(GATC)GAA
#16-23 CCTCTACTAGAGCG20(TCGA)TT

def main():
  
  #Parse Command line
  helpmsg = "Takes Genome Sequences from DIR and generates all posible pyroprints.\nThen produces a pearson correlation and displays the results"
  parser = OptionParser(usage=helpmsg)
  parser.add_option("-p", "--path", dest="dir", default="Genome Sequences/rDNA plasmid sequences/23-5/", help="Path to Genome Sequence Folders")
  parser.add_option("-d", dest="DispSeq", default="AACACGCGA23(GATC)GAA", help="Dispensation order")
  parser.add_option("-m", "--max", dest="max", type="int", default=-1, help="Max number of combinations to generate")
  parser.add_option("-f", "--file", dest="file", help="File containing parameters")
  parser.add_option("--primer", dest="primer", default="AACCTT", help="Primer to use")  

  (options, args) = parser.parse_args()
 
  if options.file:
    #Use file for args
    config = ConfigParser.RawConfigParser()
    config.read(options.file)
    dataPath = config.get("params", "path")
    
    if config.has_option("params", "max"):
       comboLimit = config.getint("params", "max")
    else:
       comboLimit = -1

    primerSequence = config.get("params", "primer")
    inDispSeq = config.get("params", "DispSeq")
  else:
    #Use command line args
     dataPath = options.dir
     comboLimit = options.max
     inDispSeq = options.DispSeq
     primerSequence = options.primer

  #Detect gpu support
  if gpu_support:
     print "PyCUDA detected. GPU support enabled."
  else:
     print('PyCUDA not detected. Calculating the Pearson correlation matrix for lots')
     print('of strains will take literally a gazillion years. Go grab a beer... or')
     print('two... or a billion.')

  #Build complete sequence string from input
  dispSeq = buildDispSeq(inDispSeq)
  
  #Find all files in Dir
  path = dataPath
  listing = os.listdir(path)
  allSequences =[]
  #TODO: parse both types of files. Cassettes and rDNA
  
  print "Fetching Files"
  
  for infile in listing:
    # Open File
    with open(dataPath + infile) as f:
        text = f.read()
        substring = ""

        if(text.find("ribosomal RNA")>0):
            for line in text:
                if ">" in line:
                    allSequences.append(substring)
                    substring = line
                else:
                    substring += line.replace("\n","")
        else:
            for line in text:
              substring += line.replace("\n","")

            allSequences.append(substring)
  
  seqList = []
  primer = primerSequence
  
  print "Generating Sequences"
  for sequence in allSequences:
    #find primer
    if primer in sequence:
      seqCount = 0
      dispCount = 0
      primerLoc = sequence.find(primer)
      #get next X dispensations(X = length of the dispensation sequence(def 104) - will make a pyroprint the length of the dispensation sequence)
      while dispCount < len(dispSeq):
        if sequence[primerLoc+len(primerSequence)+seqCount] == dispSeq[dispCount]:
          seqCount += 1
        elif (sequence[primerLoc+len(primerSequence)+seqCount] != 'A') & (sequence[primerLoc+len(primerSequence)+seqCount] != 'T') & (sequence[primerLoc+len(primerSequence)+seqCount] != 'C') &(sequence[primerLoc+len(primerSequence)+seqCount] != 'G'):
          seqCount += 1
          dispCount += 1
        else:
            dispCount += 1
      #add sequence to the list
      seqList.append(sequence[primerLoc+len(primerSequence):primerLoc+len(primerSequence)+seqCount])


  #find unique strings
  uniqueSequences = []
  for oneSeq in seqList:
    if oneSeq not in uniqueSequences:
      uniqueSequences.append(oneSeq)
  allCombinations = combinations_with_replacement(uniqueSequences,7)


  print "Pyroprinting Sequences"

  #find all combinations
  numCombos = 0
  allPyroPrints = []
  
  print "Total Combos " + str(len(uniqueSequences))
 
  for oneCombo in allCombinations:

    if comboLimit > 0 and numCombos >= comboLimit:
        break 

    allPyroPrints.append(pyroprintData(oneCombo, dispSeq))
    numCombos += 1

    # Writing out every number is too much I/O, it slows down pyroprinting a lot.
    if numCombos % 10 == 0:
        sys.stdout.write("\rGenerating Combo %d" % numCombos)
        sys.stdout.flush()

  print "\n" + str(numCombos) + " Pyroprints Generated"
  
  allPCorrs = [] 
  smallestPCor = 1000 
  largestPCor = 0 

  print 'Calculating Pearson Correlation'

  if gpu_support:
    buckets = biogpu.correlation.pearson(allPyroPrints, 10000)
    #cuda_plot(buckets)
  else:
    for i in range(0, len(allPyroPrints)-1):
      for j in range(i+1,len(allPyroPrints)-1):
        
        currPearsonR = pearsonr(allPyroPrints[i],allPyroPrints[j])[0]
        if(math.isnan(currPearsonR)!=True):
           if(currPearsonR < smallestPCor):
              smallestPCor = currPearsonR

           if(currPearsonR > largestPCor):
              largestPCor = currPearsonR

           allPCorrs.append(currPearsonR)

        mu, sigma = 100, 15
        #x = mu + sigma * np.random.randn(10000)
        x = allPCorrs
        fig = plt.figure()
        ax = fig.add_subplot(111)

        # the histogram of the data
        #n, bins, patches = ax.hist(x, 50, normed=1, facecolor='green', alpha=0.75)
        n, bins, patches = ax.hist(x)
        # hist uses np.histogram under the hood to create 'n' and 'bins'.
        # np.histogram returns the bin edges, so there will be 50 probability
        # density values in n, 51 bin edges in bins and 50 patches.  To get
        # everything lined up, we'll compute the bin centers
        bincenters = 0.5*(bins[1:]+bins[:-1])
        # add a 'best fit' line for the normal PDF
        y = mlab.normpdf( bincenters, mu, sigma)
        l = ax.plot(bincenters, y, 'r--', linewidth=1)

        ax.set_xlabel('Correlation Range')
        ax.set_ylabel('Number of Correlation')
        ax.set_title('Pearson Correlation of Data')
        rangePearsonCor = largestPCor - smallestPCor
        largestN = 0 ;
        for i in range(0, len(n)):
          print n[i]
          if n[i] > largestN:
            largestN = n[i] 
        ax.set_xlim(smallestPCor - (.1*rangePearsonCor), largestPCor + (.1*rangePearsonCor))
        ax.set_ylim(0, largestN*1.1)
        ax.grid(True)

        #plt.show()
        fname = 'pyroprintHisto.png'
        print 'Saving frame', fname
        fig.savefig(fname)

# Builds the whole dispensation order string from the string seqExp
# seqExp should be in the format of [StartSeq](NumRepeat)[RepeatedSeq]
def buildDispSeq(seqExp):
        seq = re.findall('[a-zA-Z]+|\d+\([a-zA-Z]+\)',seqExp)
        
        complete = ''
        for item in seq:
                if re.match('\d',item):
                        loopinfo = re.split('\(|\)',item)
                        count = int(loopinfo[0])
                        chars = loopinfo[1]
                        i = 0
                        while i < count:
                                complete += chars
                                i += 1
                else:
                        complete += item

        return complete


#Generates a pyroprint.
#Input: oneCombo - A single combination of 7 sequences, used to generate the pyroprint.
#Input: dispSeq  - The current dispensation sequence, can be changed from the command line.
#This is done by iterating through the entire dispensation sequence and getting the "heights" at that dispensation.
#The heights are a count of how many of the next characters in the sequence are the current dispensation.
#All 7 heights are added up to get the final simulated pyroprint.
def pyroprintData(oneCombo, dispSeq):
  sequence = oneCombo
  
  #Saved heights from all 7 sequences
  pyroData = [[],[],[],[],[],[],[]]
  #Final heights
  height = []
  #Current sequence position
  seqCount = 0
  #Current disposition
  dispCount = 0
  #Current height
  pyroCount = 0
  #Length of sequences
  length = [len(sequence[0]), len(sequence[1]), len(sequence[2]), len(sequence[3]), len(sequence[4]), len(sequence[5]), len(sequence[6])]
  #Sequence Counter
  t=0
  
  #Go through the 7 sequences and run through the disposition sequence getting the heights
  while t < 7:
    while seqCount < length[t]:
      if sequence[t][seqCount] == dispSeq[dispCount]:
        pyroCount += 1
        seqCount += 1
        if seqCount == length[t]:
          pyroData[t].append(pyroCount)
      elif (sequence[t][seqCount] != 'A') & (sequence[t][seqCount] != 'T') & (sequence[t][seqCount] != 'C') & (sequence[t][seqCount] != 'G'):
        seqCount += 1
        dispCount += 1
        pyroData[t].append(pyroCount)
      else:
        pyroData[t].append(pyroCount)
        pyroCount = 0
        dispCount += 1
    seqCount = 0
    dispCount = 0
    pyroCount = 0
    t += 1
  
  seqCount = 0

  #Get the max length of the heights (since they can be different - finish quicker/slower)
  maxVal = max(len(pyroData[0]),len(pyroData[1]),len(pyroData[2]),len(pyroData[3]),len(pyroData[4]),len(pyroData[5]),len(pyroData[6]))
  
  #Pad the heights that do not have 0's that need them for adding (to make all the lengths the same)
  x=0
  while x < 7:
    t = len(pyroData[x])
    while (len(dispSeq) - t) > 0:
      pyroData[x].append(0)
      t += 1
    x += 1
  
  #Get the final heights
  while seqCount < len(dispSeq):
    height.append( int(pyroData[0][seqCount]) + int(pyroData[1][seqCount]) + int(pyroData[2][seqCount]) + int(pyroData[3][seqCount]) + int(pyroData[4][seqCount]) + int(pyroData[5][seqCount]) + int(pyroData[6][seqCount]))
    seqCount += 1
  
  return height

#Get the length of the iterator
def getIterLength(iterator):
  temp = list(iterator)
  result = len(temp)
  iterator = iter(temp)
  return result

#def expandDispSeq(dispSeq):
#  #cctctactagagcg20(tcga)tt
#  #aacacgcga23(gatc)gaa
#  #number_regex = re.compile('[0-9]*')
#  multiplier = ''
#  multStart = 0 
#  multEnd = 0
#  multSeq = ''
#  result = ''
#  p = re.compile('[0-9]+')
#  print p.findall(dispSeq)
#  print 'result'
#  print result
#  return result
#
#
#  for i in range(0, len(dispSeq)):
#    if dispSeq[i:i+1].isdigit():
#      multiplierStart = i 
#
#    
#  for char in dispSeq:
#    if char.isalpha():
#      result.append(char)
#    if char.isdigit():
#      multiplier.appent(char)
#    if char.expandDispSeq
#      print "Seq"
#  #TODO finish this code
  
# Produces all the Chose(total,r) combonations
# Used to produce all posible pryoprints
def combinations_with_replacement(iterable, r):
  # combinations_with_replacement('ABC', 2) --> AA AB AC BB BC CC
  pool = tuple(iterable)
  n = len(pool)
  if not n and r:
    return
  indices = [0] * r
  yield tuple(pool[i] for i in indices)
  while True:
    for i in reversed(range(r)):
      if indices[i] != n - 1:
        break
    else:
      return
    indices[i:] = [indices[i] + 1] * (r - i)
    yield tuple(pool[i] for i in indices)

# Takes output from pearson_kernel and displays a graph of the results
def cuda_plot(buckets):
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

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.bar(range(len(bins)), bins, align='center')
    ax.set_title('Pearson Correlation of Data')
    ax.set_ylabel('Number of Pyroprints')
    ax.set_xticks(range(len(bins)))
    ax.set_xticklabels(['', '0% - 70%', '70% - 90%', '90% - 97%', '97% - 99%',
                        '99% - 99.7%', '99.7% - 100%'])
    ax.grid(True)
    fig.autofmt_xdate()
    print('Saving pyroprintHisto.png')
    fig.savefig('pyroprintHisto.png')


# This is the standard boilerplate that calls the main() function.
if __name__ == '__main__':
  main()
