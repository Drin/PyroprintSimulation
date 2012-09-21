#!/usr/bin/python -tt
# Copyright 2010 Google Inc.
# Licensed under the Apache License, Version 2.0
# http://www.apache.org/licenses/LICENSE-2.0

import os
import sys
import math
import time

import numpy
import pycuda.compiler
import pycuda.driver

import SimulationConfig
import extractAlleles

try:
   import biogpu.correlation

except:
   print("Could not load pyCUDA")
   sys.exit(1)

CORE_KERNEL_FILE = "biogpu/kernels.cu"

class PyroprintSimulation():
   def __init__(self):
      cli_desc = "Parameter configuration for a pyroprint simulation"
      self.configuration = SimulationConfig.parse_args(cli_desc)

   """
   Generates bucket pearson correlation value ranges in step incremements
   """
   def generateRanges(start, stop, step):
      range_start = start
      range_end = start + step

      while (range_end < stop):
         yield (range_start, range_end)
         range_start = range_end
         range_end += step

   """
   Calculates the enumerations space for numAlleles number of alleles and
   numRegions number of ITS Region copies for the bacteria to be simulated. More
   formally, the enumeration space is represents all unique combinations (with
   replacement) of numAlleles and numRegions.
   """
   def calcCombinations(numAlleles, numRegions):
      return (math.factorial(numAlleles + numRegions - 1)/
             (math.factorial(numRegions) * math.factorial(numAlleles - 1)))

if __name__ == '__main__':
   simulation = PyroprintSimulation()
   print("Simulation state: {0}\n".format(simulation.get_state()))
