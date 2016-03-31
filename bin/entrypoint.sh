#!bin/bash

export LIBRARY_PATH=/opt/AMDAPPSDK-3.0/lib/x86_64
export LD_LIBRARY_PATH=/opt/AMDAPPSDK-3.0/lib/x86_64
export OPENCL_VENDOR_PATH=/etc/OpenCL/vendors
export DATA_PATH=/data

cd /sibernetic

git checkout development

#if the $BRANCH_SIBERNETIC environment variable is set, check out the latest in that branch and use it
if [ -n "$BRANCH_SIBERNETIC" ]; then
git pull
git checkout $BRANCH_SIBERNETIC
fi

#If the REBUILD_SIBERNETIC environment variable is set to TRUE, rebuild sibernetic
if [ -n "$REBUILD_SIBERNETIC" ]; then 
make clean
make all
fi

#If the SIBERNETIC_NAME environment variable is set, use it to define the name of the data directory
if [ -n "$SIBERNETIC_NAME" ]; then 
export DATA_PATH=/data-$SIBERNETIC_NAME
fi

export PYTHONPATH=./src

if [ -n "$OCL_FILE" ] && [-n "$CONF" ] && [-n "$TIMELIMIT" ]; then
./Release/Sibernetic -f $CONF timelimit=$TIMELIMIT -no_g oclsourcepath=src/$OCL_FILE -l_to lpath=$DATA_PATH >>$DATA_PATH/log.out 2>&1
if [ -n "$OCL_FILE" ] && [-n "$TIMELIMIT" ]; then
./Release/Sibernetic -f worm timelimit=$TIMELIMIT -no_g oclsourcepath=src/$OCL_FILE -l_to lpath=$DATA_PATH >>$DATA_PATH/log.out 2>&1
if [ -n "$OCL_FILE" ] && [-n "$CONF" ]; then
./Release/Sibernetic -f $CONF -no_g oclsourcepath=src/$OCL_FILE -l_to lpath=$DATA_PATH >>$DATA_PATH/log.out 2>&1
if [ -n "$TIMELIMIT" ] && [-n "$CONF" ]; then
./Release/Sibernetic -f $CONF -no_g timelimit=$TIMELIMIT -l_to lpath=$DATA_PATH >>$DATA_PATH/log.out 2>&1
elif [-n "$OCL_FILE" ]; then
./Release/Sibernetic -f worm -no_g oclsourcepath=src/$OCL_FILE -l_to lpath=$DATA_PATH >>$DATA_PATH/log.out 2>&1
elif [-n "$CONF" ]; then
./Release/Sibernetic -f $CONF -no_g -l_to lpath=$DATA_PATH >>$DATA_PATH/log.out 2>&1
elif [-n "$TIMELIMIT" ]; then
./Release/Sibernetic -f worm -no_g -l_to timelimit=$TIMELIMIT lpath=$DATA_PATH >>$DATA_PATH/log.out 2>&1
else
./Release/Sibernetic -f worm -no_g -l_to lpath=$DATA_PATH >>$DATA_PATH/log.out 2>&1
