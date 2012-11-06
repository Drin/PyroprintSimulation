#!/usr/bin/python -tt

import os
import sys

class Sequence(object):
   def __init__(self, src_file, seq, strain_id=None, loci_id=None):
      self.src_file = src_file
      self.sequence = seq
      self.strain_id = strain_id
      self.loci_id = loci_id
      self.allele = None
      self.pyro = None

   def set_pyroprint(self, pyroprint):
      self.pyro = pyroprint

   def get_pyroprint(self):
      return self.pyro

   def set_allele(self, allele):
      self.allele = allele

   def get_allele(self):
      return self.allele

   def __hash__(self):
      return self.sequence.__hash__()

   def __eq__(self, otherSeq):
      if (isinstance(otherSeq, Sequence)):
         return self.sequence.__eq__(otherSeq.sequence)

   def __str__(self):
      return self.sequence

   def __repr__(self):
      return self.sequence
