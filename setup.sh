#!/bin/bash
script_name=$BASH_SOURCE
currentDir=$PWD

##setup gcc and python
/bin/grep ' release [2345]\.' /etc/redhat-release >/dev/null 2>&1 && \
  echo "WARNING: This version of the HSG7 ROOT build will not work on an SLC5 machine. Please use a different machine (eg. lxplus.cern.ch)" >&2

# ROOT and GCC
export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
source $ATLAS_LOCAL_ROOT_BASE/user/atlasLocalSetup.sh
#lsetup "root 6.08.06-HiggsComb-x86_64-slc6-gcc49-opt"
lsetup "root 6.18.00-x86_64-centos7-gcc8-opt"

# PYTHON
PYTHONDIR=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/x86_64/python/2.7.13-x86_64-slc6-gcc49/2.7.13-597a5/x86_64-slc6-gcc49-opt
export PATH=${PYTHONDIR}/bin:$PATH

# pip
export PATH=$PATH:/afs/cern.ch/user/x/xju/public/python/bin

# additional lib.
#export PYTHONPATH=$PYTHONPATH:$PYTHONDIR/lib:/afs/cern.ch/user/x/xju/public/python/lib/python2.7/site-packages/

# increase stack size - needed for large workspaces
ulimit -S -s unlimited
