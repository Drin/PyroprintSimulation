import sys
import os
import itertools
import math
import re 
from optparse import OptionParser
import ConfigParser

def extractAlleles():
   DEBUG = True

   if DEBUG:
      print "parsing args..\n"

   (dataDir, disp, primer, numIsolates) = handleArgs()

   if DEBUG:
      print "expanding dispensation sequence..\n"

   dispSeq = expandSequence(disp)
  
   if DEBUG:
      print "retrieving sequence files..\n"

   validSequenceFiles = findSequenceFiles(dataDir)
   allSequences = extractFileSequences(validSequenceFiles)

   if DEBUG:
      print "retrieved all DNA sequences.\n"

   seqList = pyroprintSequences(allSequences, dispSeq, primer)

   uniqueSequences = []
   uniqueSeqFiles = []
   for (seqFile, seqStr) in seqList:
      if seqStr not in uniqueSequences:
         uniqueSequences.append(seqStr)
         uniqueSeqFiles.append(seqFile)

   for alleleFile, allele in map(None, uniqueSeqFiles, uniqueSequences):
      print "allele '{0}' from file '{1}'\n".format(allele, alleleFile)

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

         elif ((seq[primerLoc+len(primer)+seqCount] != 'A') and
               (seq[primerLoc+len(primer)+seqCount] != 'T') and
               (seq[primerLoc+len(primer)+seqCount] != 'C') and
               (seq[primerLoc+len(primer)+seqCount] != 'G')):
            seqCount += 1
            dispCount += 1

         else:
            dispCount += 1

      seqList.append((seqFile, seq[primerLoc+len(primer):primerLoc+len(primer)+seqCount]))
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
