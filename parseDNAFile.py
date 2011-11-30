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
import math




def main():
  #TODO: Get Dispensation Seq from Command Line
  #Current disposition sequences (more can be added/changed)
  dispSeq = "ATCG"
  dispSeq1 = "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"
  dispSeq2 = "aacacgcgagatcgatcgatcgatcgatcgatcgatcgatcgatcgatcgatcgatcgatcgatcgatcgatcgatcgatcgatcgatcgatcgatcgatcgaa"
  dispSeq3 = "cctctactagagcgtcgatcgatcgatcgatcgatcgatcgatcgatcgatcgatcgatcgatcgatcgatcgatcgatcgatcgatcgatcgatt"
  dispSeq2 = dispSeq2.upper()
  dispSeq3 = dispSeq3.upper()
  
  
  #Find all files in Dir
  #TODO: get Dir from command line arg
  path = 'Genome Sequences/ecoli36cassettes/'
  listing = os.listdir(path)
  allSequences =[]
  #TODO: parse both types of files. Cassettes and rDNA
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
  #TODO: get Primer from Command Line Arg
  #16-23 primer GGAACCTGCGGTTGGATCAC
  #23-5 primer CGTGAGGCTTAACCTT
  
  seqList = []
  primer = "GGAACCTGCGGTTGGATCAC"
  
  for sequence in allSequences:
    #print("sequence" + sequence)
    #find primer
    if primer in sequence:
      #Current sequence position
      seqCount = 0
      #Current disposition
      dispCount = 0
      primerLoc = sequence.find(primer)
      #print('Found Primer At:' + str(primerLoc))
      #print("Found Primer in Sequence:" + sequence[primerLoc:primerLoc+20])
      #print("Found Sequence after Primer:" + sequence[primerLoc+20:primerLoc+20+140])
      #get next 104 dispensations
      while dispCount < len(dispSeq1):
        #print sequence[primerLoc+20+seqCount], dispSeq1[dispCount], seqCount, dispCount
        if sequence[primerLoc+20+seqCount] == dispSeq1[dispCount]:
          seqCount += 1
        elif (sequence[primerLoc+20+seqCount] != 'A') & (sequence[primerLoc+20+seqCount] != 'T') & (sequence[primerLoc+20+seqCount] != 'C') &(sequence[primerLoc+20+seqCount] != 'G'):
          seqCount += 1
          dispCount += 1
        else:
            dispCount += 1
        seqList.append(sequence[primerLoc+20:primerLoc+20+seqCount])
        #seqList.append(sequence[primerLoc+20:primerLoc+20+104])
  
  #find unique strings
  uniqueSequences = []
  for oneSeq in seqList:
    if oneSeq not in uniqueSequences:
      uniqueSequences.append(oneSeq)
  allCombinations = combinations_with_replacement(uniqueSequences,7)
  
  #find all combinations
  numCombos = 0
  allPyroPrints = []
  for oneCombo in allCombinations:
    allPyroPrints.append(pyroprintData(oneCombo))
    numCombos = numCombos +1
    #print('numCombos')
    moduloResults = numCombos%100    
    if moduloResults ==0:
      print numCombos
    if numCombos > 10:
      break
      
  allPCorrs = [] 
  smallestPCor = 1000 
  largestPCor = 0 
  print 'TRY TO FIND PEARSON CORS'
  for i in range(0, len(allPyroPrints)-1):
    for j in range(i+1,len(allPyroPrints)-1):
      print 'allPyroPrints[i]'
      print allPyroPrints[i]
      print 'allPyroPrints[j]'
      print allPyroPrints[j]
      print 'i'
      print i 
      print 'j'
      print j 
      currPearsonR = pearsonr(allPyroPrints[i],allPyroPrints[j])[0]
      if(math.isnan(currPearsonR)!=True):
         if(currPearsonR < smallestPCor):
            smallestPCor = currPearsonR
         if(currPearsonR > largestPCor):
            largestPCor = currPearsonR
         allPCorrs.append(currPearsonR)
      # for onePyroPrints in allPyroPrints:
      #if numCombos%1000 = 0 
      #print scipy.stats.pearsonr(onePyroPrints,lastPyroprint)
      #print "New Pyroprints"
      #print len(allPyroPrints[i])
      #print len(allPyroPrints[i+1])
      '''if(len(allPyroPrints[i]) == len(allPyroPrints[i+1])):
        currPearsonR = pearsonr(allPyroPrints[i],allPyroPrints[i+1])[0]
        if(currPearsonR < smallestPCor):
          smallestPCor = currPearsonR
        if(currPearsonR > largestPCor):
          largestPCor = currPearsonR
      allPCorrs.append(currPearsonR)  '''
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

def pyroprintData(oneCombo):
  sequence = oneCombo

  #Current disposition sequences (more can be added/changed)
  dispSeq = "ATCG"
  dispSeq1 = "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"
  dispSeq2 = "aacacgcgagatcgatcgatcgatcgatcgatcgatcgatcgatcgatcgatcgatcgatcgatcgatcgatcgatcgatcgatcgatcgatcgatcgatcgaa"
  dispSeq3 = "cctctactagagcgtcgatcgatcgatcgatcgatcgatcgatcgatcgatcgatcgatcgatcgatcgatcgatcgatcgatcgatcgatcgatt"
  dispSeq2 = dispSeq2.upper()
  dispSeq3 = dispSeq3.upper() 
  
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
  #print len(sequence[0]), len(sequence[1]), len(sequence[2]), len(sequence[3]), len(sequence[4]), len(sequence[5]), len(sequence[6])
  
  #Sequence Counter
  t=0
  
  #Go through the 7 sequences and run through the disposition sequence getting the heights
  while t < 7:
    while seqCount < length[t]:
      #print sequence[t][seqCount], dispSeq1[dispCount]
      if sequence[t][seqCount] == dispSeq1[dispCount]:
        pyroCount += 1
        seqCount += 1
        if seqCount == length[t]:
          pyroData[t].append(pyroCount)
      elif (sequence[t][seqCount] != 'A') & (sequence[t][seqCount] != 'T') & (sequence[t][seqCount] != 'C') & (sequence[t][seqCount] != 'G'):
        seqCount += 1
        dispCount += 1
        if seqCount == length[t]:
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
  #print len(pyroData[0]),len(pyroData[1]),len(pyroData[2]),len(pyroData[3]),len(pyroData[4]),len(pyroData[5]),len(pyroData[6])
  
  
  #Pad the heights that do not have 0's that need them for adding
  x=0
  while x < 7:
    t = len(pyroData[x])
    while (len(dispSeq1) - t) > 0:
      pyroData[x].append(0)
      t += 1
    x += 1
  
  
  #Get the final heights
  while seqCount < len(dispSeq1):
    height.append( int(pyroData[0][seqCount]) + int(pyroData[1][seqCount]) + int(pyroData[2][seqCount]) + int(pyroData[3][seqCount]) + int(pyroData[4][seqCount]) + int(pyroData[5][seqCount]) + int(pyroData[6][seqCount]))
    seqCount += 1
  
  
  #print len(height)
  
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
