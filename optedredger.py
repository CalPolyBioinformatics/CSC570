#!/usr/bin/python
import pymysql
import string
import sys

primer = "CGTGAGGCTTAACCTT"
revcomprimer = "AAGGTTAAGCCTCACG"
seqFile1 = ""
ratioSeq1 = ""
seqFile2 = ""
ratioSeq2 = ""
fileDir = "./Genome Sequences/rDNA plasmid sequences/23-5"
outfile = ""

# A quick tool to dredge the opterons from .seq files
# Takes in three params:
#  - the ratio of sequences
#  - the two sequence files
def main():
   global seqFile1
   global ratioSeq1
   global seqFile2
   global ratioSeq2
   global fileDir
   global outfile

   if len(sys.argv) < 5 or sys.argv[1] == "help":
      printUsage()
      return

   # Parse input params
   i = 1
   while i < len(sys.argv):
      if sys.argv[i] == "-r":
         try:
            (r1, split, r2) = sys.argv[i + 1].partition(":")
            ratioSeq1 = int(r1)
            ratioSeq2 = int(r2)
            i += 2
         except ValueError:
            print "exception"
            printUsage()
            return
      elif sys.argv[i] == "-s":
         (seqFile1, split, seqFile2) = sys.argv[i + 1].partition(":")
         outfile = seqFile1 + "-" + seqFile2 + ".opts"
         i += 2
      elif sys.argv[i] == "-d":
         fileDir = argv[i + 1]
         i += 2
      elif sys.argv[i] == "-o":
         outfile = argv[i + 1]
         i += 2
      else:
         printUsage()
         return

   seq1 = findSeq(seqFile1)
   seq2 = findSeq(seqFile2)

   # Print the found sequences to an output file using the ratios
   fd = open(outfile, "w")
   for i in range(ratioSeq1):
      fd.write("# " + str(i + 1) + " - " + seqFile1 + "\n")
      fd.write(seq1 + "\n")

   for i in range(ratioSeq2):
      fd.write("# " + str(i + 1) + " - " + seqFile2 + "\n")
      fd.write(seq2 + "\n")

   fd.close()

   return

def findSeq(seqFile):
   try:
      fd = open(fileDir + "/23-5 " + seqFile + ".seq", "r")
   except:
      print ("File Not Found: " +
            fileDir + "/23-5 " + seqFile + ".seq")
      sys.exit()

   seq = ""

   for line in fd:
      seq += line.strip()

   # FIRST try to find the sequence after the original primer
   (pre, p, end) = seq.partition(primer)

   # IF didn't find using the original primer try the reverse compliment
   if len(p) == 0 or len(end) == 0:
      (pre, p, end) = seq.partition(revcomprimer)

      if len(p) == 0 or len(end) == 0:
         print ("Runtime Error: Could not find the primer or the reverse " +
                " compliment of the primer")
         sys.exit()

      # REVERSE the string pre
      pre = pre[::-1]

      # Compliment the string
      pre = compliment(pre)

      return pre
   else:
      return end

def compliment(seg):
   # Compliment the string
   slen = len(seg)
   for i in range(slen):
      if seg[i] == "A":
         seg = seg[0:i] + "T" + seg[i + 1:slen]
      elif seg[i] == "T":
         seg = seg[0:i] + "A" + seg[i + 1:slen]
      elif seg[i] == "C":
         seg = seg[0:i] + "G" + seg[i + 1:slen]
      elif seg[i] == "G":
         seg = seg[0:i] + "C" + seg[i + 1:slen]

   return seg

def printUsage():
   print ("Usage: " + sys.argv[0] + " -r <X:Y> -s " +
         "<sequence 1 name>:<sequence 2 name> [-d <sequence file directory>] [-o <outfile>]")
   print ("   -r : The ratio of sequence 1 to sequence 2 in the final result")
   print ("   -s : A \"ratio\" of sequence names to be used.\n" +
         "         Example: Dg03-5:Hu01-3\n" +
         "         NOTE: \'23-5 \' and \'.seq\' are automatically added")
   print ("   -d : An optional parameter to change the default location to search " +
         "for the sequence files.\n" +
         "         The default location is: " + fileDir)
   print ("   -o : Optional parameter to change the default output file.\n" +
          "        The default output file name is \"<seq 1 name>-<seq 2 name>.opts\"")

if __name__ == "__main__":
   main()
