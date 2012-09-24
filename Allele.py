#!/usr/bin/python -tt

class Allele(object):
   def __init__(self, dna_seq, src_file, pyroprint):
      self.sequence = dna_seq
      self.src_file = src_file
      self.pyroprint = pyroprint

   def __hash__(self):
      return self.sequence.__hash__()

   def __eq__(self, otherAllele):
      if (isinstance(otherAllele, Allele)):
         return self.sequence.__eq__(otherAllele.sequence)

   def __str__(self):
      return self.sequence

   def __repr__(self):
      return self.sequence
