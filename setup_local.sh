#!/bin/bash

currentDir=$PWD
##setup gcc and python
/bin/grep ' release [2345]\.' /etc/redhat-release >/dev/null 2>&1 && \
  echo "WARNING: This version of the HSG7 ROOT build will not work on an SLC5 machine. Please use a different machine (eg. lxplus.cern.ch)" >&2

# first setup HZZWorkspace
cd /afs/cern.ch/user/x/xju/work/code/HZZWorkspace/build
export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
source $ATLAS_LOCAL_ROOT_BASE/user/atlasLocalSetup.sh
asetup AnalysisBase,21.2.89,here && source x86_64-centos7-gcc8-opt/setup.sh
cd $currentDir

export PYTHONPATH=$PYTHONPATH:/afs/cern.ch/work/x/xju/Install/python2p7/lib/python2.7/site-packages
export PATH=/afs/cern.ch/work/x/xju/Install/python2p7/bin:$PATH

# increase stack size - needed for large workspaces
ulimit -S -s unlimited