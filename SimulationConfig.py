#!/usr/bin/python -tt

import argparse
import ConfigParser

################################################################################
#
#  Constants
#
################################################################################
PEARSON_KERNEL_FILES = {'constant': 'biogpu/const_pearson.cu',
                        'global'  : 'biogpu/ga_pearson.cu',
                        'shared'  : 'biogpu/sh_pearson.cu',
                        'texture' : 'biogpu/text_pearson.cu',
                       }

# Lengthy default values for input parameters
DEFAULT_DISP = 'AACACGCGA23(GATC)GAA'
DEFAULT_PATH = 'Genome Sequences/'

#  Help descriptions for each input parameter
FILE_HELP   = 'Configuration file specifying configuration values'
PATH_HELP   = 'Path containing genomic sequence data for use in this simulation'
LOCI_HELP   = 'The number of ITS Regions to simulate for each isolate'
DEV_HELP    = 'The GPU card to be used for computation in this simulation'

MEM_HELP    = ('What memory to store alleles in. Possible options are ' +
               "'constantMem', 'globalMem', 'sharedMem', 'texturedMem'")
ALLELE_HELP = ('The max number of alleles to use for each isolate in this ' +
               'simulation. Defaults to the number of alleles observed.')
DISP_HELP   = ('The order in which reagents will be dispensed during the ' +
               'pyrosequencing simulation')
PRIMER_HELP = ('The forward primer to use during the PCR phase before the ' +
               'pyrosequencing phase.')

'''
This method parses command-line arguments. It returns a dictionary containing
the following keys (each key has a description of what it's value holds):
 * data_path - Directory containing DNA sequences for use in this simulation
 * disp      - Sequence/ordering for reagent dispensation when pyrosequencing
 * primer    - A forward primer (the reverse complement of the reverse primer)
 * alleles   - A limit to the number of isolates to generate in this simulation
 * num_loci  - The number of loci (replicate ITS regions) to simulate for each
               isolate
 * device    - A value indicating a specific GPU device to use for computation
 * memory    - Type of GPU memory to store alleles in during the simulation
'''
def parse_args(programDesc):
   #parser = argparse.ArgumentParser(description=programDesc, usage="Extracts alleles")
   parser = argparse.ArgumentParser(description=programDesc)

   parser.add_argument("-f",         dest="file",      default="config.cfg", help=FILE_HELP)
   parser.add_argument("-p",         dest="data_path", default=DEFAULT_PATH, help=PATH_HELP)
   parser.add_argument("-d",         dest="disp",      default=DEFAULT_DISP, help=DISP_HELP)
   parser.add_argument("-n",         dest="alleles",   default=-1, type=int, help=ALLELE_HELP)
   parser.add_argument("--primer",   dest="primer",    default="TTGGATCAC",  help=PRIMER_HELP)
   parser.add_argument("--num_loci", dest="num_loci",  default=7,  type=int, help=LOCI_HELP)
   parser.add_argument("--device",   dest="device",    default=-1, type=int, help=DEV_HELP)
   parser.add_argument("--mem_type", dest="memory",    default="constant",   help=MEM_HELP)

   args = vars(parser.parse_args())

   if (args.get('file') is not None):
      config = ConfigParser.RawConfigParser()
      config.read(args.get('file'))

      args['data_path'] = config.get("params", "data_path")
      args['disp']      = config.get("params", "disp")
      args['primer']    = config.get("params", "primer")

      if config.has_option("params", "alleles"):
         args['alleles'] = config.getint("params", "alleles")

      if config.has_option("params", "num_loci"):
         args['num_loci'] = config.getint("params", "num_loci")

      if config.has_option("params", "device"):
         args['device'] = config.getint("params", "device")

      if config.has_option("params", "memory"):
         args['memory'] = config.get("params", "memory")

   if (args.get('memory') not in PEARSON_KERNEL_FILES):
      print "Invalid memory option provided"
      sys.exit()
   else:
      args['kernel_file'] = PEARSON_KERNEL_FILES.get(args.get('memory'))

   return args
