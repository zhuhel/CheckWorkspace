# CheckWorkspace

A collection of scripts to check workspaces.
```
pip install -e . --prefix=/afs/cern.ch/work/x/xju/Install/python2p7 --process-dependency-links
```

```bash
source setup.sh
```

* *draw_pulls.py* is to plot pull distributions of nuisance parameters
* *print_np.py* is to print nuisance paramters and parameter of interest in the workspace
* *check_ws.py* is to check the workspaces... A typical example to run the combined workspace of 4l+llvv is:


### Examples
* prepare to run, assuming workspace is *workspace.root*, inside which signal POIs are *XS_ggF, XS_VBF*, 
mass parameter *mH*.
RooWorkspace name is *combined*, data is *obsData*.
```bash
mkdir run
cd run
```

* Check expected number of signal and background events in the workspace at 500 GeV.
```
python ../src/check_ws.py workspace.root --out_dir numbers -w combined -d obsData --poi_name XS_ggF --fixVar "XS_VBF=0, mH=500" --nBins 60
```
* background-only fit, with correlation matrix
```bash
python ../src/check_ws.py workspace.root --out_dir bkg_only_fit -w combined -d obsData --poi_name XS_ggF --fixVar "XS_ggF=0,XS_VBF=0,mH=500" --afterFit --matrix --conditionalFit --nBins 60
```
* background-only fit, with qqZZ free. assuming qqZZ scaling factor is *mu_qqZZ*.
```bash
python ../src/check_ws.py workspace.roo --out_dir bkg_only_fit_qqZZFree -w combined -d obsData --poi_name XS_ggF --fixVar "XS_ggF=0,XS_VBF=0,mH=500" --floatVar "mu_qqZZ" --afterFit --matrix --conditionalFit --nBins 60
```
* signal+background fit
```bash
python ../src/check_ws.py workspace.roo --out_dir sig_plus_bkg_fit -w combined -d obsData --poi_name XS_ggF --fixVar "XS_VBF=0,mH=500" --afterFit --matrix --nBins 60
```

To obtain full help:
```
cd src
python check_ws.py --help

Usage: check_ws.py [options] file_name out_name

check yields and shape for WS

Options:
  --version             show program's version number and exit
  -h, --help            show this help message and exit
  -w WS_NAME, --wsname=WS_NAME
  -m MC_NAME, --mcname=MC_NAME
  -d DATA_NAME, --dataname=DATA_NAME
                        name of observed data
  --poi_name=POI_NAME   name of POI
  --fixVar=FIXVAR       set variables as constant:  mu=1,lumi=1
  --floatVar=FLOATVAR   set variables float:  mu_ggF,mu_VBF
  --noPlot              don't make plots
  --nBins=N_BINS        setup binning of the observable
  --xMax=X_MAX          max value of the observable
  --xMin=X_MIN          min value of the observable
  --logY                if log scale for y-axis
  --signalScale=SIG_SCALE
                        scale factor applied to data
  --afterFit            make plots and yields after fit
  --conditionalFit      perform conditional fit
  --lumi=LUMI           which luminosity used
  --matrix              plot covariance matrix
  --debug               in debug mode
  --out_dir=OUT_DIR     name of output directory
```  
 
