#!/usr/bin/python -tt

import sys
import os
import itertools
import math
import re 

from Sequence import Sequence

NUCL_COMPL = {'A' : 'T',
              'T' : 'A',
              'C' : 'G',
              'G' : 'C',
              'N' : 'N',
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
   if ('DEBUG' in os.environ):
      print("{0} sequence files".format(len(validSequenceFiles)))
   #From DNA sequence files, extract the actual DNA sequences
   allSequences = extractFileSequences(validSequenceFiles)

   #pyroprint the DNA sequences to discover alleles amongst all of the given
   #DNA. seqList is a tuple containing:
   #(<fileName>, <allele DNA sequence>, <allele histogram>)
   seqList = pyroprintSequences(allSequences, expanded_seq, configuration)

   if ('DEBUG' in os.environ):
      print ("numAlleles: {0}\n".format(configuration.get('alleles')))

   uniqueSeqs = set()

   # find unique strings
   for seq_obj in seqList:
      (seqFile, seqStr, seqPeaks) = (seq_obj.src_file, seq_obj.get_allele(),
                                     seq_obj.pyro)

      if ((configuration.get('alleles') == -1 or
          len(uniqueSeqs) < configuration.get('alleles')) and
          seqStr not in uniqueSeqs):
         uniqueSeqs.add(seq_obj)

   if ('DEBUG' in os.environ):
      for seq_obj in uniqueSeqs:
         print ("allele '{0}' from file '{1}'\n\thas pyroprint '{2}'\n".format(
                seq_obj.sequence, seq_obj.src_file, seq_obj.pyro))

   return (uniqueSeqs, len(expanded_seq))

def findSequenceFiles(data_path):
   validSequenceFiles = []
   allSequenceFiles = os.listdir(data_path)

   for entry in allSequenceFiles:
      if os.path.isdir(data_path + "/" + entry):
         for subEntry in os.listdir(data_path + "/" + entry):
            allSequenceFiles.append(entry + "/" + subEntry)

      elif os.path.isfile(data_path + "/" + entry):
         if entry.find(".seq") > 0 or entry.find(".txt") > 0:
            validSequenceFiles.append(os.path.join(data_path, entry))

   return validSequenceFiles 


'''
   Returns the same list passed in (allSequences) but with each sequence object
   having an additional field set -- pyroprint.
'''
def pyroprintSequences(allSequences, dispSeq, config):
   pyro_len = config.get('pyro_len') if config.get('pyro_len') > 0 else len(dispSeq)

   for seq_obj in allSequences:
      (seqFile, seq) = (seq_obj.src_file, seq_obj.sequence)

      peakVals = [0] * pyro_len
      (peakNdx, seqCount, dispCount) = (0, 0, 0)

      if (config.get('primer') not in seq and
          config.get('primer') not in reverseComplSeq(seq)):
         continue

      primerLoc = seq.find(config.get('primer'))
      if (primerLoc < 0):
         seq = reverseComplSeq(seq)
         primerLoc = seq.find(config.get('primer'))

      primer_end = primerLoc + len(config.get('primer'))

      while (dispCount < pyro_len):
         if (seq[primer_end + seqCount] == dispSeq[dispCount]):
            seqCount += 1
            peakVals[peakNdx] += 1

         elif ((seq[primer_end + seqCount] != 'A') and
               (seq[primer_end + seqCount] != 'T') and
               (seq[primer_end + seqCount] != 'C') and
               (seq[primer_end + seqCount] != 'G')):
            seqCount += 1
            dispCount += 1
            peakNdx += 1

         else:
            dispCount += 1
            peakNdx += 1

      seq_obj.set_pyroprint(peakVals)
      seq_obj.set_allele(seq[primer_end : primer_end + seqCount])
   return allSequences

'''
   Returns a list of Sequence objects (see Sequence.py)
'''
def extractFileSequences(sequenceFiles):
   allSequences = []

   for sequenceFile in sequenceFiles:
      with open(sequenceFile) as f:
         text = f.read()
         substring = ""
         matches = re.search("^>(\d+) .* (.*):", text)

         if (matches is not None):
            for line in text:
               if ">" in line:
                  if (substring != ""):
                     allSequences.append(Sequence(sequenceFile, substring,
                                                  matches.group(2),
                                                  matches.group(1)))
                  substring = ""
               else:
                  substring += line.replace("\n","")
         else:
            for line in text:
               substring += line.replace("\n","")

         if (matches is not None):
            allSequences.append(Sequence(sequenceFile, substring,
                                         matches.group(2), matches.group(1)))
         else:
            allSequences.append(Sequence(sequenceFile, substring))

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
