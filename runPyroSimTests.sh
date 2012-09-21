#!/bin/bash
#export CUDA_PROFILE=1
#export CUDA_PROFILE_CSV=1
#export CUDA_PROFILE_CONFIG=cuda_profile.cfg

PROFILING_SETTINGS="CUDA_PROFILE=1 CUDA_PROFILE_CSV=1 CUDA_PROFILE_CONFIG=cuda_profile.cfg"
MEM_TYPE="constant"
NUM_ALLELES="14"
DEVICE_NUM=0

TEST0_PARAMS="--memory=${MEM_TYPE}Mem --device=${DEVICE_NUM} -m $NUM_ALLELES"
TEST1_PARAMS="--memory=${MEM_TYPE}Mem --device=${DEVICE_NUM} -m $NUM_ALLELES"
TEST2_PARAMS="--memory=${MEM_TYPE}Mem --device=${DEVICE_NUM} -m $NUM_ALLELES"
TEST3_PARAMS="--memory=${MEM_TYPE}Mem --device=${DEVICE_NUM} -m $NUM_ALLELES"
TEST4_PARAMS="--memory=${MEM_TYPE}Mem --device=${DEVICE_NUM} -m $NUM_ALLELES"

OUT_0="${MEM_TYPE}_${NUM_ALLELES}_0.out"
LOG_0="${MEM_TYPE}_${NUM_ALLELES}_0.log"

OUT_1="${MEM_TYPE}_${NUM_ALLELES}_1.out"
LOG_1="${MEM_TYPE}_${NUM_ALLELES}_1.log"

OUT_2="${MEM_TYPE}_${NUM_ALLELES}_2.out"
LOG_2="${MEM_TYPE}_${NUM_ALLELES}_2.log"

OUT_3="${MEM_TYPE}_${NUM_ALLELES}_3.out"
LOG_3="${MEM_TYPE}_${NUM_ALLELES}_3.log"

OUT_4="${MEM_TYPE}_${NUM_ALLELES}_4.out"
LOG_4="${MEM_TYPE}_${NUM_ALLELES}_4.log"

if [ ! -f $OUT_0 ] && [ ! -f $OUT_1 ] && [ ! -f $OUT_2 ] && [ ! -f $OUT_3 ] && [ ! -f $OUT_4 ]; then
   python2.7 parseDNAFile.py $TEST0_PARAMS > $OUT_0
   mv cuda_profile_0.log $LOG_0
   python2.7 parseDNAFile.py $TEST1_PARAMS > $OUT_1
   mv cuda_profile_0.log $LOG_1
   python2.7 parseDNAFile.py $TEST2_PARAMS > $OUT_2
   mv cuda_profile_0.log $LOG_2
   python2.7 parseDNAFile.py $TEST3_PARAMS > $OUT_3
   mv cuda_profile_0.log $LOG_3
   python2.7 parseDNAFile.py $TEST4_PARAMS > $OUT_4
   mv cuda_profile_0.log $LOG_4
else
   echo "Error: Not running, don't want to clobber old results..."
fi

MEM_TYPE="global"
NUM_ALLELES="14"

TEST0_PARAMS="--memory=${MEM_TYPE}Mem --device=${DEVICE_NUM} -m $NUM_ALLELES"
TEST1_PARAMS="--memory=${MEM_TYPE}Mem --device=${DEVICE_NUM} -m $NUM_ALLELES"
TEST2_PARAMS="--memory=${MEM_TYPE}Mem --device=${DEVICE_NUM} -m $NUM_ALLELES"
TEST3_PARAMS="--memory=${MEM_TYPE}Mem --device=${DEVICE_NUM} -m $NUM_ALLELES"
TEST4_PARAMS="--memory=${MEM_TYPE}Mem --device=${DEVICE_NUM} -m $NUM_ALLELES"

OUT_0="${MEM_TYPE}_${NUM_ALLELES}_0.out"
LOG_0="${MEM_TYPE}_${NUM_ALLELES}_0.log"

OUT_1="${MEM_TYPE}_${NUM_ALLELES}_1.out"
LOG_1="${MEM_TYPE}_${NUM_ALLELES}_1.log"

OUT_2="${MEM_TYPE}_${NUM_ALLELES}_2.out"
LOG_2="${MEM_TYPE}_${NUM_ALLELES}_2.log"

OUT_3="${MEM_TYPE}_${NUM_ALLELES}_3.out"
LOG_3="${MEM_TYPE}_${NUM_ALLELES}_3.log"

OUT_4="${MEM_TYPE}_${NUM_ALLELES}_4.out"
LOG_4="${MEM_TYPE}_${NUM_ALLELES}_4.log"

if [ ! -f $OUT_0 ] && [ ! -f $OUT_1 ] && [ ! -f $OUT_2 ] && [ ! -f $OUT_3 ] && [ ! -f $OUT_4 ]; then
   python2.7 parseDNAFile.py $TEST0_PARAMS > $OUT_0
   mv cuda_profile_0.log $LOG_0
   python2.7 parseDNAFile.py $TEST1_PARAMS > $OUT_1
   mv cuda_profile_0.log $LOG_1
   python2.7 parseDNAFile.py $TEST2_PARAMS > $OUT_2
   mv cuda_profile_0.log $LOG_2
   python2.7 parseDNAFile.py $TEST3_PARAMS > $OUT_3
   mv cuda_profile_0.log $LOG_3
   python2.7 parseDNAFile.py $TEST4_PARAMS > $OUT_4
   mv cuda_profile_0.log $LOG_4
else
   echo "Error: Not running, don't want to clobber old results..."
fi
