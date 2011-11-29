#!/usr/bin/python -tt
# Copyright 2010 Google Inc.
# Licensed under the Apache License, Version 2.0
# http://www.apache.org/licenses/LICENSE-2.0


"""
  So far I have done:
  Find all files with DNA info in DIR
  Parse all files and find all DNA Chunks
  in each chunk find the primer sequence (if it exists)
  store all the 104char segments after the primer
  find all unique sequences of 104char segments
  find all possible choose 7 combinations of unique sequences
  TODO:
  generate pyroprints (python)
  generate pyroprints (CUDA)
  compare pyroprints
  """

import sys
import os
import itertools
from scipy.stats.stats import pearsonr
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab


# Define a main() function that prints a little greeting.
def main():
  #Find all files in Dir
  path = 'Genome Sequences/ecoli36cassettes/'
  listing = os.listdir(path)
  allSequences =[]
  for infile in listing:
    #print "current file is: " + 'Genome Sequences/ecoli36cassettes/' + infile
    # Open File
    f = open('Genome Sequences/ecoli36cassettes/' + infile)
    #print("currentFile: "+infile)
    substring = ""
    for line in f:
      if ">" in line:
        #print("NewSection:")
        allSequences.append(substring)
        #print(entireFile)
        substring = line
      else:
        substring += line.replace("\r\r\n","")
    #print(substring)
    f.close()
  #TODO: Parse Primer from file Primer.txt
  #16-23 primer GGAACCTGCGGTTGGATCAC
  #23-5 primer CGTGAGGCTTAACCTT
  
  seqList = []
  primer = "GGAACCTGCGGTTGGATCAC"
  
  for sequence in allSequences:
    #print("sequence" + sequence)
    #find primer
    if primer in sequence:
      primerLoc = sequence.find(primer)
      #print('Found Primer At:' + str(primerLoc))
      #print("Found Primer in Sequence:" + sequence[primerLoc:primerLoc+20])
      #print("Found Sequence after Primer:" + sequence[primerLoc+20:primerLoc+20+140])
      #get next 104 chars
      seqList.append(sequence[primerLoc+20:primerLoc+20+104])
  uniqueSequences = []
  numUniqueSeq = 0
  #find unique strings
  for oneSeq in seqList:
    #print("OneSeq::")
    #print(oneSeq)
    if oneSeq not in uniqueSequences:
      #print("found non unique seq:")
      #print(oneSeq)
      numUniqueSeq = numUniqueSeq + 1
      uniqueSequences.append(oneSeq)
  
  # else:
  #uniqueSequences.append(oneSeq)
  print 'num unique seq:'
  print numUniqueSeq
  #print uniqueSequences
  allCombinations = combinations_with_replacement(uniqueSequences,7)
  #allCombinations = itertools.combinations_with_replacement(uniqueSequences,7)
  #print ("Number of Unique Combinations: " + len(allCombinations))
  #print(len(allCombinations))
  numCombos = 0
  allPyroPrints = []
  
  #print(getIterLength(allCombinations))

  for oneCombo in allCombinations:
    allPyroPrints.append(pyroprintData(oneCombo))
    numCombos = numCombos +1
    #print('numCombos')
    moduloResults = numCombos%1000    
    if moduloResults ==0:
      print numCombos
    if numCombos > 10000:
      break
      
  allPCorrs = [] 
  smallestPCor = 1000 
  largestPCor = 0 
  
  for i in range(0, len(allPyroPrints)-1):
    # for onePyroPrints in allPyroPrints:
    #if numCombos%1000 = 0 
    #print scipy.stats.pearsonr(onePyroPrints,lastPyroprint)
    #print "New Pyroprints"
    #print len(allPyroPrints[i])
    #print len(allPyroPrints[i+1])
    if(len(allPyroPrints[i]) == len(allPyroPrints[i+1])):
      currPearsonR = pearsonr(allPyroPrints[i],allPyroPrints[i+1])[0]
      if(currPearsonR < smallestPCor):
          smallestPCor = currPearsonR
      if(currPearsonR > largestPCor):
          largestPCor = currPearsonR
    allPCorrs.append(currPearsonR)  
  print allPCorrs




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
  #print "onePyroPrint + lastPyroPrint"
    #lastPyroprint = onePyroPrints
    #print onePyroPrints  
  #print numCombos%1000
  #print(numCombos)

  #print(numCombos)
  # for i in range(0, 100):
  # print(allCombinations[0])


#perform Pearson Corralation on pyroprints

#save all sequences

#???

#Profit!

def pyroprintData(oneCombo):
  sequence = oneCombo
  #Current disposition sequence (more can be added/changed)
  dispSeq = "ATCG"
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
      elif (sequence[t][seqCount] != 'A') & (sequence[t][seqCount] != 'T') & (sequence[t][seqCount] != 'C') & (sequence[t][seqCount] != 'G'):
        seqCount += 1
      else:
        pyroData[t].append(pyroCount)
        pyroCount = 0
        if dispCount == 3:
          dispCount = 0
        else:
          dispCount += 1
    seqCount = 0
    dispCount = 0
    pyroCount = 0
    t += 1
  
  seqCount = 0
  #Get the max length of the heights (since they can be different - finish quicker/slower)
  maxVal = max(len(pyroData[0]),len(pyroData[1]),len(pyroData[2]),len(pyroData[3]),len(pyroData[4]),len(pyroData[5]),len(pyroData[6]))
  
  #Pad the heights that do not have 0's that need them for adding
  x=0
  while x < 7:
    t = len(pyroData[x])
    while (maxVal - t) > 0:
      pyroData[x].append(0)
      t += 1
    x += 1
  
  #Get the final heights
  while seqCount < maxVal:
    height.append( int(pyroData[0][seqCount]) + int(pyroData[1][seqCount]) + int(pyroData[2][seqCount]) + int(pyroData[3][seqCount]) + int(pyroData[4][seqCount]) + int(pyroData[5][seqCount]) + int(pyroData[6][seqCount]))
    seqCount += 1
  
  #Print out the heights
  '''x=0
    while x<len(height):
    print x, dispSeq[dispCount], ':', height[x]
    x += 1
    if dispCount == 3:
    dispCount = 0
    else:
    dispCount += 1
    '''
  #print height
  #print '-->height'
  return height

def getIterLength(iterator):
  temp = list(iterator)
  result = len(temp)
  iterator = iter(temp)
  return result


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




# This is the standard boilerplate that calls the main() function.
if __name__ == '__main__':
  main()
