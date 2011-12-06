#!/usr/bin/python -tt
# Copyright 2010 Google Inc.
# Licensed under the Apache License, Version 2.0
# http://www.apache.org/licenses/LICENSE-2.0

# Google's Python Class
# http://code.google.com/edu/languages/google-python-class/

# Projecct Part 2 - Pyroprinting
# Generates pyroprints from 7 sequences that are located in a file
# Depending on how/what we want to pyroprint, slight changes should be made
# Expects a file with 7 sequences starting with the '>' line (which is ignored)
# Sequences can be different lengths and of any length

import sys
import os
import itertools

def main():
  #Open a file to make a pyroprint of the data
  fd = open('Genome Sequences/ecoli36cassettes/E.-coli-536-16s.txt')
  #7 sequences will be stored here
  sequence = ["", "", "", "", "", "", ""]
  #grab all the sequences from the file and close the file
  t=-1
  for line in fd:
    if ">" in line:
       t += 1
       sequence = sequence
    else:
      sequence[t] += line.replace("\r\r\n","")
  fd.close()  

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
        seqCount  += 1
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
  x=0
  while x<len(height):
    print x, dispSeq[dispCount], ':', height[x]
    x += 1
    if dispCount == 3:
      dispCount = 0
    else:
      dispCount += 1    

# This is the standard boilerplate that calls the main() function.
if __name__ == '__main__':
  main()
