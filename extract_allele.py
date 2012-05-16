import sys
import os
import itertools
import math
import re 
from optparse import OptionParser
import ConfigParser

def extractAlleles():
   (dataDir, disp, primer, numIsolates) = handleArgs()

   dispSeq = expandSequence(disp)
  
   validSequenceFiles = findSequenceFiles(dataDir)
   allSequences = extractFileSequences(validSequenceFiles)

   for seq in allSequences:
      print "{0} is a valid file\n".format(seq)

   sys.exit(0)
  
   seqList = []
   primer = primerSequence
  
   print "Generating Sequences"
   for sequence in allSequences:
      #find primer
      if primer in sequence:
         seqCount = 0
         dispCount = 0
         primerLoc = sequence.find(primer)
         #get next X dispensations(X = length of the dispensation sequence(def 104) -
         #will make a pyroprint the length of the dispensation sequence)
         while dispCount < len(dispSeq):
            if sequence[primerLoc+len(primerSequence)+seqCount] == dispSeq[dispCount]:
               seqCount += 1

            elif (sequence[primerLoc+len(primerSequence)+seqCount] != 'A') and (sequence[primerLoc+len(primerSequence)+seqCount] != 'T') and (sequence[primerLoc+len(primerSequence)+seqCount] != 'C') and (sequence[primerLoc+len(primerSequence)+seqCount] != 'G'):
               seqCount += 1
               dispCount += 1

            else:
               dispCount += 1

         #add sequence to the list
         seqList.append(sequence[primerLoc+len(primerSequence):primerLoc+len(primerSequence)+seqCount])

   #find unique strings
   uniqueSequences = []
   for oneSeq in seqList:
      if oneSeq not in uniqueSequences:
         uniqueSequences.append(oneSeq)

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
         #print "{0} is a directory\n".format(dataDir + "/" + entry)
         for subEntry in os.listdir(dataDir + "/" + entry):
            #print "{0} is a sub entry\n".format(subEntry)
            allSequenceFiles.append(entry + "/" + subEntry)

      elif os.path.isfile(dataDir + "/" + entry):
         #print "{0} is not a directory\n".format(dataDir + "/" + entry)
         if entry.find(".seq") > 0 or entry.find(".txt") > 0:
            #print "{0} is a valid sequence file\n".format(entry)
            validSequenceFiles.append(dataDir + "/" + entry)

   return validSequenceFiles 

def extractFileSequences(sequenceFiles):
   allSequences = []
  
   for pathEntry in os.listdir(sequenceFiles):
      with open(pathEntry) as f:
         text = f.read()
         substring = ""

         if (text.find("ribosomal RNA") > 0):
            for line in text:
               if ">" in line:
                  allSequences.append(substring)
                  substring = line
               else:
                  substring += line.replace("\n","")
         else:
            for line in text:
               substring += line.replace("\n","")
               allSequences.append(substring)

   return allSequences

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
