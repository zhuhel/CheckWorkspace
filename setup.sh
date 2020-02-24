#!/bin/bash

script_name=$BASH_SOURCE
currentDir=$PWD

##setup gcc and python
/bin/grep ' release [2345]\.' /etc/redhat-release >/dev/null 2>&1 && \
  echo "WARNING: This version of the HSG7 ROOT build will not work on an SLC5 machine. Please use a different machine (eg. lxplus.cern.ch)" >&2

# first setup HZZWorkspace

lsetup "sft --cmtconfig=x86_64-slc6-gcc49-opt releases/LCG_88/pyanalysis/2.0,releases/LCG_88/lapack/3.5.0,releases/LCG_88/blas/20110419"

# increase stack size - needed for large workspaces
ulimit -S -s unlimited
