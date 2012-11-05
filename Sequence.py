#!/usr/bin/python -tt

import os
import sys

class Sequence(object):
   def __init__(self, src_file, seq, strain_id, loci_id):
      self.src_file = src_file
      self.sequence = seq
      self.strain_id = strain_id
      self.loci_id = loci_id

   def set_pyroprint(self, pyroprint):
      self.pyro = pyroprint

   def set_allele(self, allele):
      self.allele = allele

   def __hash__(self):
      return self.sequence.__hash__()

   def __eq__(self, otherSeq):
      if (isinstance(otherSeq, Sequence)):
         return self.sequence.__eq__(otherSeq.sequence)

   def __str__(self):
      return self.sequence

   def __repr__(self):
      return self.sequence
