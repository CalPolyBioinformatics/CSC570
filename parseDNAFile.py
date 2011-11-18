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
#import itertools
#import pearsonr

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
  #find unique strings
  for oneSeq in seqList:
    #print("OneSeq::")
    #print(oneSeq)
    if oneSeq not in uniqueSequences:
      #print("found non unique seq:")
      #print(oneSeq)
      uniqueSequences.append(oneSeq)
  # else:
  #uniqueSequences.append(oneSeq)
  
  #print uniqueSequences
  allCombinations = combinations_with_replacement(uniqueSequences,7)
  #print ("Number of Unique Combinations: " + len(allCombinations))
  #print(len(allCombinations))
  numCombos = 0
  allPyroPrints = []
  
  #print(getIterLength(allCombinations))
  
  for oneCombo in allCombinations:
    allPyroPrints.append(pyroprintData(oneCombo))
    numCombos = numCombos +1
    #print('numCombos')
    moduloResults = numCombos%100
    if moduloResults ==0:
      print numCombos
    if numCombos > 1000:
      break
  
  lastPyroprint = 0 ;
  for onePyroPrints in allPyroPrints:
    #if numCombos%1000 = 0 
    #print pearsonr(onePyroPrints,lastPyroprint)
    lastPyroprint = onePyroPrints
    print onePyroPrints
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
  #for oneEle in oneCombo:
  #print(oneEle)
  #Open a file to make a pyroprint of the data
  #fd = open('Genome Sequences/ecoli36cassettes/E.-coli-536-16s.txt')
  #7 sequences will be stored here
  #### sequence = ["", "", "", "", "", "", ""]
  sequence = oneCombo
  #print("Sequence: ")
  
  #grab all the sequences from the file and close the file
  #t=-1
  # for line in fd:
  #if ">" in line:
  # t += 1
  # sequence = sequence
  #else:
  # sequence[t] += line.replace("\r\r\n","")
  # fd.close()
  
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
