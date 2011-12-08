#!/usr/bin/python

import pymysql
import string
import random
import sys
import math
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
db_num_opts = 7      # The number of opterons in a pyrorun from the database

# For Parsing Opteron Files and Building the Graph 
global num_opts
opteron_size = 104

# Dispensation Sequence Derived from File:
# '12-15-10b rep plasmid ratios-ModGATC-2controls.xls'
graph_disp_seq = "AACACGCGAGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGAA"

def main():
   global con

   try:
      opts = parseData(sys.argv[1])
   except IndexError:
      print "Usage: ./insilico.py <opterons.opt>"
      return

   con = pymysql.connect(host = ghost, port = gport, user = guser, passwd = gpasswd, db = gdb_name)

   getPyroDispensationSeq()
   sdata = getSampleData()

   # Get the Unit Heights for each nucleotide (phcomps) and the standard dev for
   #     each nucleotide's unit value
   (phcomps, stdDevs) = getPeakHeightCompensations()
   con.close()

   tbl = averageHeights(pyro_dis_seq, sdata, phcomps, stdDevs)
   # printTable(tbl)

   graph(opts, tbl, pyro_dis_seq)

# Parses the file of opterons in filename.
# @param filename - A newline-separated list of comments and opterons
# RETURNS: A list of opertons
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
                optList.append(opt[0:104])
                num_opts += 1

    return optList

# Gets the list of unit value compensated peak height values for each
#     of the nucleotides and the standard deviation of unit values for
#     each nucleotide
# RETURNS: a tuple of the average unit peak height values and the
#     standard deviations
def getPeakHeightCompensations():
   nsamps = number_samples * 20  # Increase the number of samples to ensure
                                 # and adequate sample
   # peak height dictionary for calculating standard dev and averages
   # phc_dict is a LIST of unit compensaited peak height values
   phc_dict = {"A" : [], "T" : [], "C" : [], "G" : []}
   samplePyroIds = []      # A list of pyroprint samples used

   # GET a list of all unique pyroIds from the Histograms
   sql = "SELECT DISTINCT H.pyroId FROM Histograms as H"
   sql += " JOIN Pyroprints as P ON (H.pyroId = P.pyroId) "
   sql += "WHERE dsName = '" + short_disp + "'"

   cur = con.cursor()

   cur.execute(sql)

   # FOR each returned result
   pyroids = []
   for (r,) in cur:
      # APPEND the pyroId to the list of possible pyroids
      pyroids.append(int(r))

   plen = len(pyroids) - 1

   # FOR the number of samples to get
   for i in range(nsamps):
      # SELECT a random pyroId from the list
      pid = pyroids[random.randint(0, plen)]

      # WHILE the candidate pyroId is in the samples list
      while pid in samplePyroIds:
         # SELECT a random pyroId from the list
         pid = pyroids[random.randint(0, plen)]

      # APPEND the qualified pyroId to the list of sample pyroids
      samplePyroIds.append(pid)

      # get the first primer_len number of nucleotides peak heights pairs with pyroId
      sql = "SELECT nucleotide, pHeight FROM Histograms WHERE pyroId = " + str(pid)
      sql += " LIMIT " + str(primer_len)

      cur.execute(sql)

      j = 1
      (basen,baseh) = cur.fetchone()
      baseh = float(baseh)
      phc_dict[basen].append(baseh / db_num_opts)
      # FOR each nucleotide-peak height pair
      for (n,h) in cur:
         # IF the height was less than one, junk data and discard
         if h < 1:
            continue
         else:
            try:
               # Compute the unit value for this peak height given the number of opterons
               # NOTE: you must convert the peak height 'h' to a float!
               multi = 1
               while (baseh * multi) / float(h) < 1:
#                  print str(multi + 1) + " " + str((baseh * multi) / (float(h)))
                  multi += 1

#               print (str(j) + " " + n + "(" + str(float(h)) + ") / 7 * " + str(multi) + ": " 
#                     + str(float(h) / (db_num_opts * multi)))
#               print "  " + basen + " baseh(" + str(baseh) + ") / h(" + str(float(h)) + "): " + str(baseh / float(h))
               phc_dict[n].append(float(h) / (db_num_opts * multi));
               # phc_dict[n].append(float(h) / (db_num_opts * (int(math.ceil(h)) / db_num_opts)))
            except ZeroDivisionError:
               # if the peak height is smaller than the number of opterons
               # it could be the case that the unit value for this nucleotide
               #  is less than 1, in this case we do a straight division
               phc_dict[n].append(float(h) / db_num_opts)

         j+= 1

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

   # Return the average for each nucleotide and the standard deviations
   return (avgNucs, stdDev_Nucs)

# Connects to the database and retreives the samples to be used when
#     building the table
# RETURNS: a list of samples where a single sample is a list of peak heights
#     for a single pyro run
def getSampleData():
   samples = []         # A list of pyroprint samples from the DB
   samplePyroIds = []   # A list of pyroprints added to the samples
   cur = con.cursor()
   sql = "SELECT DISTINCT H.pyroId FROM Histograms as H"
   sql += " JOIN Pyroprints as P ON (H.pyroId = P.pyroId) "
   sql += "WHERE dsName = '" + short_disp + "' "
   # sql += "ORDER BY RAND() LIMIT " + number_samples

   cur.execute(sql)

   # GET a list of distinct pyroIds to randomly select from
   rws = cur.fetchall()
   pyroids = []
   for (r,) in rws:
      pyroids.append(int(r))

   plen = len(pyroids) - 1

   # FOR the number of samples to be retrieved
   for i in range(number_samples):
      # GET a random pyroId
      pid = pyroids[random.randint(0, plen)]

      # WHILE this pyroId has already been selected
      while pid in samplePyroIds:
         # GET another random pyroId
         pid = pyroids[random.randint(0, plen)]

      # APPEND a valid pyroId (one not already selected) to the list of samples
      samplePyroIds.append(pid)

      # GET the peak heights from the Histograms with pyroId pid
      sql = "SELECT pHeight FROM Histograms WHERE pyroId = " + str(pid)

      cur.execute(sql)

      # Create a new empty sample list
      sample = []

      # FOR each peak height returned
      for (r,) in cur:
         # APPEND the peak height to the sample list
         # NOTE: you must convert the returned value to a float since
         #     PyMySQL returns a Decimal type
         sample.append(float(r))

      # APPEND the sample to the list of samples
      samples.append(sample)

   cur.close()

   return samples

# Connects to the database and gets the dispSeq from the Dispensation Table
#     where the dsName matches pyro_ds_name
# SETS: the global pyro_dis_seq to the database result
# Required because there are no relations in the DB to do this conversion
#     automatically
def getPyroDispensationSeq():
   global pyro_dis_seq
   sql = "SELECT dispSeq FROM Dispensation WHERE dsName = '" + pyro_ds_name + "'"

   cur = con.cursor()
   cur.execute(sql)

   # The extra comma (,) is required here WITHOUT the parens because pyMySQL
   #     returns a tuple no matter the number of columns selected
   tmp, = cur.fetchone()

   # FOR each character returned, build the pyro_dis_seq list
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

# Creates a blank slot to put in a sample entry 
# Each entry contains the current unit height (idx 0) and the total number 
# of nucleotides seen (idx 1)
# RETURNS a new slot
def blankSlot():
   heightEntry = [0,0]
   return {"A": heightEntry, "G": heightEntry, "C": heightEntry, "T": heightEntry}

# Creates a list to use as a "sample". The size of the empty sample created is 
# DEFAULT_SAMPLE_SIZE which is defined at the top of this file 
# RETURNs a list of dictionaries of lists (a new sample)
def createTableEntry():
   entry = []
   for i in range(DEFAULT_SAMPLE_SIZE):
      entry.append(blankSlot())
   return entry

# Handles the unit height calculation for a single sample. Calculates the 
# average height for each position in a sample and keeps count of the total 
# number of nucleotides of a particular letter at the position
# When the algorithm completes, the entry in the table for the sample is updated 
# with the calculated values  
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
      if numNucleotides < 1:
            numNucleotides = 1

      calcUnitAvg = smpl[idx]/numNucleotides

      if (abs(calcUnitAvg - unitHeights[letter]) > stdDev[letter]):
         numNucleotides = numNucleotides - 1

      currAvg = (totalVal + smpl[idx]/numNucleotides) / (totalSeen + numNucleotides) # add a height and average

      # Update values in list
      table[dSeq][idx][letter][HEIGHT] = currAvg
      table[dSeq][idx][letter][SEEN] = totalSeen + numNucleotides

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

# Parses the list of opterons and writes to pyroprint.csv with a list
# of peak height values starting with position 0 and going until 
# position 103. 
# @param opterons - A list of opertons
# @param tbl - A 4D table used in calculating the peak height values
#              for the list of opterons
# @param pyro_dis_seq - Used to index into the correct dispensation
#                       sequence for the 4D table
def graph(opterons, tbl, pyro_dis_seq):
    dispSeq = list(graph_disp_seq)
    dSeq = "".join(pyro_dis_seq)
    fd = open("pyroprint.csv", "w")

    # FOR each position
    for i in range(0, opteron_size):

        # Get the next nucleotide from the dispensation sequence
        nuc = dispSeq[i]

        # Declare a nucleotide counter initialized to 0
        nucleotidesCounted = 0

        # FOR each opteron
        for j in range(0, num_opts):

            # WHILE each nucleotide matches the dispseq nucleotide
            while opterons[j][0] == nuc:
                # Increment the number of nucleotides counted
                nucleotidesCounted += 1
                opterons[j] = opterons[j][1:len(opterons[j])]
                      
        # Multiply the number of nucelotides counted by the comp. peak height
        ph = tbl[dSeq][i][nuc][HEIGHT] * nucleotidesCounted
      
        # Add a tuple to the output list in format: (position, dispensation nucleotide, peak height value)
        # Uncomment and use to include position and nucleotide character
            # fd.write(str(i) + ", " + str(nuc) + ", ")
        fd.write(str(ph) + ",")

if __name__ == "__main__":
   main()
