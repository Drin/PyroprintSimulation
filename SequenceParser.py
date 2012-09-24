#!/usr/bin/python -tt

import sys
import os
import itertools
import math
import re 

from Allele import Allele

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
def extractAlleles(configuration):
   #Expand an abbreviated "dispensation sequence"
   #e.g. 2(ATCG) expands to ATCGATCG
   expanded_seq = expandSequence(configuration.get('disp'))

   #Find files that contain actual DNA sequence data
   validSequenceFiles = findSequenceFiles(configuration.get('data_path'))
   #From DNA sequence files, extract the actual DNA sequences
   allSequences = extractFileSequences(validSequenceFiles)

   #pyroprint the DNA sequences to discover alleles amongst all of the given
   #DNA. seqList is a tuple containing:
   #(<fileName>, <allele DNA sequence>, <allele histogram>)
   seqList = pyroprintSequences(allSequences, expanded_seq, configuration)

   if ('DEBUG' in os.environ):
      print ("numAlleles: {0}\n".format(configuration.get('alleles')))

   alleles = set()

   # find unique strings
   for (seqFile, seqStr, seqPeaks) in seqList:
      if ((configuration.get('alleles') == -1 or
          len(alleles) < configuration.get('alleles')) and
          seqStr not in alleles):
         alleles.add(Allele(seqStr, seqFile, seqPeaks))

   if ('DEBUG' in os.environ):
      for allele in alleles:
         print ("allele '{0}' from file '{1}'\n\thas pyroprint '{2}'\n".format(
                allele.sequence, allele.src_file, allele.pyroprint))

   return (alleles, len(expanded_seq))

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


def pyroprintSequences(allSequences, dispSeq, config):
   seqList = []

   pyro_len = config.get('pyro_len') if config.get('pyro_len') > 0 else len(dispSeq)

   for (seqFile, seq) in allSequences:
      peakVals = [0] * pyro_len
      peakNdx = 0
      seqCount = 0
      dispCount = 0

      if (config.get('primer') not in seq and
          config.get('primer') not in reverseComplSeq(seq)):
         continue

      primerLoc = seq.find(config.get('primer'))
      if (primerLoc < 0):
         seq = reverseComplSeq(seq)
         primerLoc = seq.find(config.get('primer'))

      while (dispCount < pyro_len):
         if (seq[primerLoc + len(config.get('primer')) + seqCount] ==
             dispSeq[dispCount]):
            seqCount += 1
            peakVals[peakNdx] += 1

         elif ((seq[primerLoc+len(config.get('primer'))+seqCount] != 'A') and
               (seq[primerLoc+len(config.get('primer'))+seqCount] != 'T') and
               (seq[primerLoc+len(config.get('primer'))+seqCount] != 'C') and
               (seq[primerLoc+len(config.get('primer'))+seqCount] != 'G')):
            seqCount += 1
            dispCount += 1
            peakNdx += 1

         else:
            dispCount += 1
            peakNdx += 1

      seqList.append((seqFile, seq[primerLoc + len(config.get('primer')) :
                                   primerLoc + len(config.get('primer')) + seqCount],
                                   peakVals))
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
