#!/bin/bash

script_name=$BASH_SOURCE
currentDir=$PWD

##setup gcc and python
/bin/grep ' release [2345]\.' /etc/redhat-release >/dev/null 2>&1 && \
  echo "WARNING: This version of the HSG7 ROOT build will not work on an SLC5 machine. Please use a different machine (eg. lxplus.cern.ch)" >&2

# ROOT and GCC
# export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
# source $ATLAS_LOCAL_ROOT_BASE/user/atlasLocalSetup.sh
# asetup AnalysisBase,21.2.68,here

lsetup 'root 6.08.06-HiggsComb-x86_64-slc6-gcc49-opt'

# fix numpy problem.
#lsetup "sft --cmtconfig=x86_64-slc6-gcc49-opt releases/LCG_88/pyanalysis/2.0,releases/LCG_88/lapack/3.5.0,releases/LCG_88/blas/20110419"

# # additional lib.
# export PYTHONPATH=/afs/cern.ch/user/x/xju/public/python2.7.15/lib/python2.7/site-packages:$PYTHONPATH
# export PATH=$PATH:/afs/cern.ch/user/x/xju/public/python2.7.15/bin

# increase stack size - needed for large workspaces
ulimit -S -s unlimited
