#!/usr/bin/python -tt
# Copyright 2010 Google Inc.
# Licensed under the Apache License, Version 2.0
# http://www.apache.org/licenses/LICENSE-2.0

import os
import sys
import math
import time
import threading

import numpy
import pycuda.compiler
import pycuda.driver

from extract_allele import extractAlleles

try:
   import biogpu.correlation
except:
   print >>sys.stderr, ("Could not load pycuda")

CORE_KERNEL_FILE = "biogpu/kernels.cu"
PEARSON_KERNEL_FILES = {'constantMem': 'biogpu/const_pearson.cu',
                        'globalMem': 'biogpu/ga_pearson.cu'} #,
                       #'sharedMem': 'biogpu/sh_pearson.cu',
                       #'texturedMem': 'biogput/text_pearson.cu'}
