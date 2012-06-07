#!/bin/bash
#export CUDA_PROFILE=1
#export CUDA_PROFILE_CSV=1
#export CUDA_PROFILE_CONFIG=cuda_profile.cfg

PROFILING_SETTINGS="CUDA_PROFILE=1 CUDA_PROFILE_CSV=1 CUDA_PROFILE_CONFIG=cuda_profile.cfg"

TEST1_PARAMS="--memory=globalMem --device=0 -m -1"
TEST2_PARAMS="--memory=globalMem --device=1 -m -1"
TEST3_PARAMS="--memory=constantMem --device=2 -m -1"
TEST4_PARAMS="--memory=constantMem --device=3 -m -1"

TEST1_OUT="0_global_full.out"
TEST2_OUT="1_global_full.out"
TEST3_OUT="0_constant_full.out"
TEST4_OUT="1_constant_full.out"

if [ ! -f $TEST1_OUT ] && [ ! -f $TEST2_OUT ] && [ ! -f $TEST3_OUT ] && [ ! -f $TEST4_OUT ]; then
   #python parseDNAFile.py $TEST1_PARAMS > $TEST1_OUT &
   #python parseDNAFile.py $TEST2_PARAMS > $TEST2_OUT &
   #python parseDNAFile.py $TEST3_PARAMS > $TEST3_OUT 
   python parseDNAFile.py $TEST4_PARAMS > $TEST4_OUT

elif [ ! -f $TEST4_OUT ]; then
   python parseDNAFile.py $TEST4_PARAMS > $TEST4_OUT

else
   echo "Error: Not running, don't want to clobber old results..."

fi
