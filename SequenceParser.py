#!/usr/bin/python -tt

import sys
import os
import itertools
import math
import re 


NUCL_COMPL = {'A' : 'T',
              'T' : 'A',
              'C' : 'G',
              'G' : 'C',
             }

'''
A Method which extracts alleles from a given set of DNA sequences. An allele is
a DNA sequence that is unique in content given particular pyroprinting
parameters. These alleles can be represented as DNA sequences, or alternatively
as pyroprints (histograms). Biologically speaking, an allele suggests a genetic
indication of a different strain (I think?).
'''
def extractAlleles(dataDir, disp, primer, numIsolates):
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

   print ("numAlleles: {0}\n".format(numIsolates))

   # find unique strings
   for (seqFile, seqStr, seqPeaks) in seqList:
      if (numIsolates == -1 or len(alleles) < numIsolates):
         if (seqStr not in alleles):
            alleles.append(seqStr)
            alleleFiles.append(seqFile)
            allelePeaks.append(seqPeaks)

   if ('DEBUG' in os.environ):
      for alleleFile, allele, peak in map(None, alleleFiles, alleles, allelePeaks):
         print "allele '{0}' from file '{1}'\n\thas pyroprint '{2}'\n".format(allele, alleleFile, peak)

   return (allelePeaks, len(allelePeaks), len(dispSeq))

def findSequenceFiles(data_path):
   validSequenceFiles = []
   allSequenceFiles = os.listdir(data_path)

   for entry in allSequenceFiles:
      if os.path.isdir(data_path + "/" + entry):
         for subEntry in os.listdir(data_path + "/" + entry):
            allSequenceFiles.append(entry + "/" + subEntry)

      elif os.path.isfile(data_path + "/" + entry):
         if entry.find(".seq") > 0 or entry.find(".txt") > 0:
            validSequenceFiles.append(data_path + entry)

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
   return NUCL_COMPL.get(char.upper())

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

################################################################################
#
# This method expands the abbreviated DNA sequence short_seq into a
# verbose/explicit DNA sequence. An abbreviated DNA sequence is a sequence
# where contiguous repeats (e.g. ATGATGATG) are shortened to
# '<numRepeats>(<repeatSeq>)' (e.g. 3(ATG)).
#
################################################################################
def expandSequence(short_seq):
   expanded_seq = ''
   seq = re.findall('[a-zA-Z]+|\d+\([a-zA-Z]+\)', short_seq)

   for item in seq:
      if (re.match('\d', item)):
         loopinfo = re.split('\(|\)', item)

         repeat_count = int(loopinfo[0])
         repeat_seq = loopinfo[1]

         for repeat in range(repeat_count):
            expanded_seq += repeat_seq

      else:
         expanded_seq += item

   return expanded_seq

if __name__ == '__main__':
   extractAlleles()
