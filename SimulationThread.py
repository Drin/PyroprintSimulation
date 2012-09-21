#!/usr/bin/python -tt
# Copyright 2010 Google Inc.
# Licensed under the Apache License, Version 2.0
# http://www.apache.org/licenses/LICENSE-2.0

import os
import sys
import math
import time
import threading

class SimulationThread(threading.Thread):
   def __init__(self, device_id, cudaModule, workload):
      self.device_id = device_id
      self.module = cudaModule
      self.workload = workload

   def run(self):
      #given a device Id, create a new CUDA context. Since contexts must be
      #thread-local we have to do this instead of passing in contexts to be
      #used
      cudaDevice = pycuda.driver.Device(self.device_id)
      cudaContext = cudaDevice.make_context()

      self.workload(self.module)

      #pop context so there are no problems of the GPU being used even after
      #this thread is finished
      pycuda.driver.Context.pop()

      #free context after it has been popped because otherwise this context
      #may be "leaked" in memory
      del cudaContext
