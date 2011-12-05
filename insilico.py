#!/usr/bin/python

import pymysql
import string
import random
import sys
import math
import pyroPrint
import cPickle

DEFAULT_SLOTS = 4000 # Default "length" of pyroprints to keep track of
DEFAULT_SAMPLE_SIZE = 200
SEEN = 1 # constant for the index of the "seen" count (used in average)
HEIGHT = 0 # constant for the index of the average height of the pyroprint
DISP = 0 # constant for the index of the dispensation seq in the pyroprint struct
HVAL = 1 # constant for the index of the heights generated from the pyroprint

# Database specific configs
global con
ghost = 'abra.csc.calpoly.edu'
gport = 3306
guser = 'csc570'
gpasswd = 'ilovedata570'
gdb_name = 'cplop'
# For pulling samples from the database to build the table
primer_len = 14
short_disp = 'CCTCTACTAGAGCG20(TCGA)TT'
pyro_ds_name = "ModTCGA-2c"
global pyro_dis_seq
pyro_dis_seq = []
number_samples = 20
global num_opts
opteron_size = 104

# Dispensation Sequence Derived from File:
# '12-15-10b rep plasmid ratios-ModGATC-2controls.xls'
graph_disp_seq = "AACACGCGAGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGAA"

def main():
   global con
   
   opts = parseData(sys.argv[1])
   
   con = pymysql.connect(host = ghost, port = gport, user = guser, passwd = gpasswd, db = gdb_name)

   getPyroDispensationSeq()
   sdata = getSampleData()

   (phcomps, stdDevs) = getPeakHeightCompensations()
   con.close()

   tbl = averageHeights(pyro_dis_seq, sdata, phcomps, stdDevs)
   printTable(tbl)

   graph(opts, tbl, pyro_dis_seq)

def parseData(filename):
    global num_opts
    num_opts = 0
    # Open file, assign to a list
    fd = open(filename, "r")
    fileInput = fd.readlines()
    fd.close()

    # Define a new list to just hold the opterons
    optList = []

    # FOR each string terminated by a newline
    for opt in fileInput :
        opt = opt.strip()
        # IF the string is not an empty string
        if len(opt) > 0 :
            # IF this is an opteron
            if  opt[0] != '#' :
                # Add to the list
                optList.append(opt)
                num_opts += 1

    return optList

def getPeakHeightCompensations():
   nsamps = number_samples * 20
   # peak height dictionary for calculating standard dev and averages
   # 
   phc_dict = {"A" : [], "T" : [], "C" : [], "G" : []}
   samplePyroIds = []

   sql = "SELECT DISTINCT H.pyroId FROM Histograms as H"
   sql += " JOIN Pyroprints as P ON (H.pyroId = P.pyroId) "
   sql += "WHERE dsName = '" + short_disp + "'"

   cur = con.cursor()

   cur.execute(sql)

   pyroids = []
   for (r,) in cur:
      pyroids.append(int(r))

   plen = len(pyroids) - 1

   for i in range(nsamps):
      pid = pyroids[random.randint(0, plen)]

      while pid in samplePyroIds:
         pid = pyroids[random.randint(0, plen)]

      samplePyroIds.append(pid)

      # get the first primer_len number of nucleotides peak heights pairs with pyroId
      sql = "SELECT nucleotide, pHeight FROM Histograms WHERE pyroId = " + str(pid)
      sql += " LIMIT " + str(primer_len)

      cur.execute(sql)

      # FOR each nucleotide-peak height pair
      for (n,h) in cur:
         # IF the height was less than one, junk data and discard
         if h < 1:
            continue
         else:
            try:
               phc_dict[n].append(h / (num_opts * (int(math.ceil(h)) / num_opts)))
            except ZeroDivisionError:
               # if the peak height is smaller than the number of opterons
               # it could be the case that the unit value for this nucleotide
               #  is less than 1, in this case we do a straight division
               phc_dict[n].append(h / num_opts)

   cur.close()

   # the average unit value for each nucleotide
   avgNucs = {'A' : 0, 'T' : 0, 'C' : 0, 'G' : 0}
   # the standard deviation for the unit value for each nucleotide
   stdDev_Nucs = {'A' : 0, 'T' : 0, 'C' : 0, 'G' : 0}

   # Get the Average Compensated Peak Height Value for each Nucleotide
   for key in phc_dict:
      for h in phc_dict[key]:
         avgNucs[key] += h
      avgNucs[key] /= nsamps

   # Get the Standard Deviation for Each Complensated Peak Height for Each Nucleotide
   for ph_nuc in phc_dict:
      sum = 0

      # get the average for the current nucleotide
      for h in phc_dict[key]:
         sum += math.pow((h - avgNucs[ph_nuc]), 2)
      avg = sum / nsamps

      # calculate and store the standard dev for this nucleotide
      stdDev_Nucs[ph_nuc] = math.sqrt(avg)

      # No longer need the list for the dictionary this nucleotide
      # Change the value at index ph_nuc to just the avg
      phc_dict[ph_nuc] = avg

   return (phc_dict, stdDev_Nucs)

def getSampleData():
   samples = []
   samplePyroIds = []
   cur = con.cursor()
   sql = "SELECT DISTINCT H.pyroId FROM Histograms as H"
   sql += " JOIN Pyroprints as P ON (H.pyroId = P.pyroId) "
   sql += "WHERE dsName = '" + short_disp + "' "
   # sql += "ORDER BY RAND() LIMIT " + number_samples

   cur.execute(sql)

   rws = cur.fetchall()
   pyroids = []
   for (r,) in rws:
      pyroids.append(int(r))

   plen = len(pyroids) - 1

   for i in range(number_samples):
      pid = pyroids[random.randint(0, plen)]

      while pid in samplePyroIds:
         pid = pyroids[random.randint(0, plen)]

      samplePyroIds.append(pid)

      sql = "SELECT pHeight FROM Histograms WHERE pyroId = " + str(pid)

      cur.execute(sql)

      sample = []

      for (r,) in cur:
         sample.append(float(r))

      samples.append(sample)

   cur.close()

   return samples

def getPyroDispensationSeq():
   global pyro_dis_seq
   sql = "SELECT dispSeq FROM Dispensation WHERE dsName = '" + pyro_ds_name + "'"

   cur = con.cursor()
   cur.execute(sql)

   tmp, = cur.fetchone()

   for c in tmp:
      pyro_dis_seq.append(c)

   cur.close()

def printTable(table):
   print("Disp Seq Table\n")

   for dispSeq in table:
      print("| "+dispSeq+" : ")
      for i in range(len(table[dispSeq])):
         print(" | "+str(table[dispSeq][i])+" | ")
      print("\n")

def blankSlot():
   heightEntry = [0,0]
   return {"A": heightEntry, "G": heightEntry, "C": heightEntry, "T": heightEntry}

def createTableEntry():
   entry = []
   for i in range(DEFAULT_SAMPLE_SIZE):
      entry.append(blankSlot())
   return entry

def avgHeights(table, dispSeq, smpl, unitHeights, stdDev):
   dSeq = "".join(dispSeq)

   for idx in range(len(smpl)):
      letter = dispSeq[idx%len(dispSeq)]
      if dSeq in table:
         currAvg = table[dSeq][idx][letter][HEIGHT] # current average
         totalSeen = table[dSeq][idx][letter][SEEN] # total number of values counted
      else:
         currAvg = 0
         totalSeen = 0
         table[dSeq] = createTableEntry()

      totalVal = currAvg * totalSeen # undo averaging
      numNucleotides = round(smpl[idx]/unitHeights[letter])
      calcUnitAvg = smpl[idx]/numNucleotides
      if (abs(calcUnitAvg - unitHeights[letter]) > stdDev[letter]):
         numNucleotides = numNucleotides - 1
      currAvg = (totalVal + smpl[idx]/numNucleotides) / (totalSeen + numNucleotides) # add a height and average
      # Update values in list
      table[dispSeq][idx][letter][HEIGHT] = currAvg
      table[dispSeq][idx][letter][SEEN] = totalSeen + numNucleotides

# Takes in a generated pyroprint and keeps a running average of the heights in
# a given slot for a given dispensation sequence
# @param pyroprint - a tuple containing lists of the dispensation seq and the height
def averageHeights(dispSeq, sampleList, unitHeights, stdDev):
   dispIdx = 0

   dataFile = open("pyroprintTrends.pkl", "wb+") # load the persistent averages
   try:
      table = cPickle.load(dataFile)
   except (EOFError, IOError): # file doesn't exist or is in wrong format
      table = {}

   for sample in sampleList:
      avgHeights(table, dispSeq, sample, unitHeights, stdDev)

   cPickle.dump(table, dataFile) # write out the changes
   print("Average heights recomputed based on new data and written to file!\n")

   return table

def graph(opterons, tbl, pyro_dis_seq):
    dispSeq = list(graph_disp_seq)
    outputList = []
    fd = open("pyroprint.out", "w")

    # FOR each position
    for i in range(0, opteron_size) :

        # Declare a nucleotide counter
        nucleotidesCounted = 0

        # Get the next nucleotide from the dispensation sequence
        dispSeqNuc = dispSeq[i % 4]

        # FOR each opteron
        for opt in opterons :
            # Set a temporary index to the current position (i.e. i)
            j = i

            # WHILE each nucleotide matches the dispSeqNucleotide
            while opt[j] == dispSeqNuc :
                # Increment the number of nucleotides counted
                nucleotidesCounted += 1
                # Move the temporary index ahead one position
                j += 1

        # Multiply the number of nucelotides counted by the comp. peak height
        ph = tbl[pyro_dis_seq][i][dispSeqNuc][HEIGHT] * nucleotidesCounted

        # Add a tuple to the output list in format: (position, dispensation nucleotide, peak height value)
        outputList.append( str(i) + ", " +  dispSeqNuc + ", "  + str(ph) + "\n")

    # Dump contents into an output file
    cPickle.dump(outputList, fd)

if __name__ == "__main__":
   main()
