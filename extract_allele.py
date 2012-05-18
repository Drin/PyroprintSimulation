import sys
import os
import itertools
import math
import re 
from optparse import OptionParser
import ConfigParser


DEBUG = True

'''
A Method which extracts alleles from a given set of DNA sequences. An allele is
a DNA sequence that is unique in content given particular pyroprinting
parameters. These alleles can be represented as DNA sequences, or alternatively
as pyroprints (histograms). Biologically speaking, an allele suggests a genetic
indication of a different strain (I think?).
'''
def extractAlleles():
   #Parse arguments and retrieve appropriate variables
   (dataDir, disp, primer, numIsolates) = handleArgs()

   #Expand an abbreviated "dispensation sequence"
   #e.g. 2(ATCG) expands to ATCGATCG
   dispSeq = expandSequence(disp)
  
   #Find files that contain actual DNA sequence data
   validSequenceFiles = findSequenceFiles(dataDir)
   #From DNA sequence files, extract the actual DNA sequences
   allSequences = extractFileSequences(validSequenceFiles)

   #pyroprint the DNA sequences to discover alleles amongst all of the given
   #DNA. seqList is a tuple containing:
   #(<fileName>, <allele DNA sequence>, <allele histogram>)
   seqList = pyroprintSequences(allSequences, dispSeq, primer)

   alleles = []
   alleleFiles = []
   allelePeaks = []

   for (seqFile, seqStr, seqPeaks) in seqList:
      if seqStr not in alleles:
         alleles.append(seqStr)
         alleleFiles.append(seqFile)
         allelePeaks.append(seqPeaks)

   for alleleFile, allele, peak in map(None, alleleFiles, alleles, allelePeaks):
      print "allele '{0}' from file '{1}'\n\thas pyroprint '{2}'\n".format(allele,
            alleleFile, peak)

   return allelePeaks

'''
This method is parses command-line arguments. It returns a directory containing
DNA sequences (or other directories) a dispensation sequence to use on given
DNA sequences, a forward primer (which may need to be applied to the reverse
complement strand) and an upper bound on the number of isolates to generate in
silico.
'''
def handleArgs():
   parser = OptionParser(usage="Extracts alleles from input DNA sequences\n")
   parser.add_option("--path", dest="path", default="Genome Sequences/",
                     help="Path containing DNA sequences")
   parser.add_option("--disp", dest="dispSeq", default="AACACGCGA23(GATC)GAA",
                     help="Dispensation Sequence")
   parser.add_option("-n", dest="num", type="int", default=-1,
                     help="Number of Isolates to generate")
   parser.add_option("-f", dest="file", default="parameters.config",
                     help="File containing parameters")
   parser.add_option("--fp", dest="forPrimer", default="AACCTT",
                     help="Forward Primer")

   (options, args) = parser.parse_args()
 
   if options.file:
      config = ConfigParser.RawConfigParser()
      config.read(options.file)
      dataPath = config.get("params", "path")
    
      if config.has_option("params", "num"):
         numIsolates = config.getint("params", "num")
      else:
         numIsolates = -1

      forwardPrimer = config.get("params", "forPrimer")
      dispSeq = config.get("params", "dispSeq")

   else:
      #Use command line args
      dataPath = options.path
      numIsolates = options.num
      dispSeq = options.dispSeq
      forwardPrimer = options.primer

   return (dataPath, dispSeq, forwardPrimer, numIsolates)


def findSequenceFiles(dataDir):
   validSequenceFiles = []
   allSequenceFiles = os.listdir(dataDir)

   for entry in allSequenceFiles:
      if os.path.isdir(dataDir + "/" + entry):
         for subEntry in os.listdir(dataDir + "/" + entry):
            allSequenceFiles.append(entry + "/" + subEntry)

      elif os.path.isfile(dataDir + "/" + entry):
         if entry.find(".seq") > 0 or entry.find(".txt") > 0:
            validSequenceFiles.append(dataDir + entry)

   return validSequenceFiles 

def pyroprintSequences(allSequences, dispSeq, primer):
   seqList = []

   for (seqFile, seq) in allSequences:
      peakVals = [0] * len(dispSeq)
      peakNdx = 0
      seqCount = 0
      dispCount = 0

      if primer not in seq and primer not in reverseComplSeq(seq):
         continue

      primerLoc = seq.find(primer)
      if (primerLoc < 0):
         seq = reverseComplSeq(seq)
         primerLoc = seq.find(primer)

      while (dispCount < len(dispSeq)):
         if (seq[primerLoc+len(primer)+seqCount] == dispSeq[dispCount]):
            seqCount += 1
            peakVals[peakNdx] += 1

         elif ((seq[primerLoc+len(primer)+seqCount] != 'A') and
               (seq[primerLoc+len(primer)+seqCount] != 'T') and
               (seq[primerLoc+len(primer)+seqCount] != 'C') and
               (seq[primerLoc+len(primer)+seqCount] != 'G')):
            seqCount += 1
            dispCount += 1
            peakNdx += 1

         else:
            dispCount += 1
            peakNdx += 1

      seqList.append((seqFile, seq[primerLoc+len(primer):primerLoc+len(primer)+seqCount], peakVals))
   return seqList

def extractFileSequences(sequenceFiles):
   allSequences = []
  
   for sequenceFile in sequenceFiles:
      with open(sequenceFile) as f:
         text = f.read()
         substring = ""

         if (text.find("ribosomal RNA") > 0):
            for line in text:
               if ">" in line:
                  if (substring != ""):
                     allSequences.append(substring)
                  substring = ""
               else:
                  substring += line.replace("\n","")
         else:
            for line in text:
               substring += line.replace("\n","")

         allSequences.append((sequenceFile, substring))

   return allSequences

def complement(char):
   if (char == 'A'):
      return 'T'
   elif (char == 'T'):
      return 'A'
   elif (char == 'C'):
      return 'G'
   elif (char == 'G'):
      return 'C'
   return char 

def reverseComplSeq(seq):
   reverseCompl = ""

   for char in seq:
      reverseCompl = complement(char) + reverseCompl

   return reverseCompl

def reverseSeq(seq):
   reverseSeq = ""
   for char in seq:
      reverseSeq = char + reverseSeq

   return reverseSeq

# Builds the whole dispensation order string from the string seqExp
# seqExp should be in the format of [StartSeq](NumRepeat)[RepeatedSeq]
def expandSequence(seqExp):
   seq = re.findall('[a-zA-Z]+|\d+\([a-zA-Z]+\)', seqExp)
        
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

if __name__ == '__main__':
   extractAlleles()
