#!/usr/bin/python

import pymysql
import string
import random
import sys
import math
import pyroPrint
import cPickle

DEFAULT_SLOTS = 4000 # Default "length" of pyroprints to keep track of
SEEN = 1 # constant for the index of the "seen" count (used in average)
HEIGHT = 0 # constant for the index of the average height of the pyroprint
DISP = 0 # constant for the index of the dispensation seq in the pyroprint struct
HVAL = 1 # constant for the index of the heights generated from the pyroprint

global con
ghost = 'abra.csc.calpoly.edu'
gport = 3306
guser = 'csc570'
gpasswd = 'ilovedata570'
gdb_name = 'cplop'
primer_len = 14
hist_ds_name = 'CCTCTACTAGAGCG20(TCGA)TT'
pyro_ds_name = "ModTCGA-2c"
pyro_dis_seq = []
number_samples = 20
num_opts = 7

def main():
   global con
   con = pymysql.connect(host = ghost, port = gport, user = guser, passwd = gpasswd, db = gdb_name)
   
   opts = parseData(sys.argv[1])
   getPyroDispensationSeq()
   sdata = getSampleData()
   
   con.close()

   tbl = averageHeights("".join(pyro_dis_seq), sdata)
   
   (phcomps, stdDevs) = getPeakHeightCompensations()
   
#   graph(opts, tbl, phcomps)

def parseData(filename):
    global num_opts
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
   phc_dict = {"A" : [], "T" : [], "C" : [], "G" : []}
   samplePyroIds = []

   sql = "SELECT DISTINCT H.pyroId FROM Histograms as H"
   sql += " JOIN Pyroprints as P ON (H.pyroId = P.pyroId) "
   sql += "WHERE dsName = '" + hist_ds_name + "'"

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

      sql = "SELECT nucleotide, pHeight FROM Histograms WHERE pyroId = " + str(pid)

      cur.execute(sql)

      j = 0
      for (n,h) in cur:
         if h < 1:
            continue
         elif j < primer_len:
            phc_dict[n].append(h/(num_opts * (int(math.ceil(h)) / num_opts)))
            j += 1
         else:
            break

   cur.close()

   avgNucs = {'A' : 0, 'T' : 0, 'C' : 0, 'G' : 0}
   stdDev_Nucs = {'A' : 0, 'T' : 0, 'C' : 0, 'G' : 0}
   # Get the Average Compensated Peak Height Value for each Nucleotide
   for key in phc_dict:
      for h in phc_dict[key]:
         avgNucs[key] += h
      avgNucs[key] /= nsamps

   # Get the Standard Deviation for Each Complensated Peak Height for Each Nucleotide
   for ph_nuc in phc_dict:
      sum = 0
      for h in phc_dict[key]:
         sum += math.pow((h - avgNucs[ph_nuc]), 2)
      sum /= nsamps
      
      stdDev_Nucs[ph_nuc] = math.sqrt(sum)

   print stdDev_Nucs
      
   return (phc_dict, stdDev_Nucs)

def getSampleData():
   samples = []
   samplePyroIds = []
   cur = con.cursor()
   sql = "SELECT DISTINCT H.pyroId FROM Histograms as H"
   sql += " JOIN Pyroprints as P ON (H.pyroId = P.pyroId) " 
   sql += "WHERE dsName = '" + hist_ds_name + "'"
 
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
         sample.append(r)

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

def avgHeights(table, dispSeq, smpl):
   for idx in range(0, len(smpl)):
      letter = dispSeq[idx]
      try:
         currAvg = table[dispSeq][idx][letter][HEIGHT] # current average
         totalSeen = table[dispSeq][idx][letter][SEEN] # total number of values counted      
      except (KeyError):
         currAvg = 0
         totalSeen = 0
  
      totalVal = currAvg * totalSeen # undo averaging
      currAvg = (totalVal + smpl[idx]) / (totalSeen + 1) # add a height and average
      # Update values in list  
      try:
         table[dispSeq][idx][letter][HEIGHT] = currAvg 
         table[dispSeq][idx][letter][SEEN] = totalSeen + 1
      except (KeyError):
         print table[dispSeq]

# Takes in a generated pyroprint and keeps a running average of the heights in
# a given slot for a given dispensation sequence 
# @param pyroprint - a tuple containing lists of the dispensation seq and the height
def averageHeights(dispSeq="GATCGTAC", sampleList=[]):
   dispIdx = 0

   dataFile = open("pyroprintTrends.pkl", "wb+") # load the persistent averages
   try:
      table = cPickle.load(dataFile)
   except (EOFError, IOError): # file doesn't exist or is in wrong format
      table = {}   

   for sample in sampleList:
      avgHeights(table, dispSeq, sample)
         
   cPickle.dump(table, dataFile) # write out the changes
   print("Average heights recomputed based on new data and written to file!\n")

   return table


# def graph(opterons, phcomps):
   


if __name__ == "__main__":
   main()
