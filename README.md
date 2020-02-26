# CheckWorkspace

## Instructions

### Prerequisites

Install the HZZWorkspace using the `AnalysisBase,21.2.89` to avoid the issue of not able to use `numpy`.
Assuming the HZZWorkspace is build in the directory: `HZZWORKSPACE-BUILD-DIR`

### Installation

* run the following to properly setup `python` and `root`:

```bash
source setup.sh HZZWORKSPACE-BUILD-DIR
```

* select a directory you want to install this package, assuming you chose `YOUR-OWN-DIR`,
which could be something like `/afs/cern.ch/work/x/xju/Install/python2p7`. Then run the following:

```bash
export PYTHONPATH=$PYTHONPATH:YOUR-OWN-DIR/lib/python2.7/site-packages

pip install -e . --prefix=YOUR-OWN-DIR --process-dependency-links
```

This installs all the exectuables specified in the `setup.py` in an __editable__ model (thanks to the `-e` option), in which case all changes you would make during developments will take effects immediately (i.e. no need of installing again).

* then run the following so that your bash can find these executables:

```bash
export PATH=YOUR-OWN-DIR/bin:$PATH
```

Note that it is convenient to put the two `export` command lines in your `setup.sh`, otherwise **they have to be executed every time after the python and root is setup**.

There are 3 executables in the `YOUR-OWN-DIR/bin`, which can be run anywhere. Type `which check_ws` to see if your bash can find the script.

Thanks to the magic header comments `#!/usr/bin/env python`, bash can automatically identify the script as a `python` script and use `python` to interpet it. No additional typing of `python` is needed!

* *draw_pulls* is to plot pull distributions of nuisance parameters
* *print_np* is to print nuisance paramters and parameter of interest in the workspace
* *check_ws* is to check the workspaces... A typical example to run the combined workspace of 4l+llvv is:

### Examples

* prepare to run, assuming workspace is *workspace.root*, inside which the signal POIs are *XS_ggF, XS_VBF* and the mass parameter is *mH*. In addition, RooWorkspace name is *combined*, data is *obsData*.

```bash
mkdir run
cd run
```

* Check expected number of signal and background events in the workspace at 500 GeV (pre-fit).

```bash
check_ws workspace.root --out_dir numbers -w combined -d obsData --poi_name XS_ggF --fixVar "XS_VBF=0, mH=500" --nBins 60
```

* Check expected number of signal and background events in the workspace at 500 GeV (pre-fit) **with background uncertainties**.

```bash
check_ws workspace.root --out_dir numbers -w combined -d obsData --poi_name XS_ggF --fixVar "XS_VBF=0, mH=500" --nBins 60 --add-bkg-sys --signalScale 0.
```

* background-only fit, with correlation matrix

```bash
check_ws workspace.root --out_dir bkg_only_fit -w combined -d obsData --poi_name XS_ggF --fixVar "XS_ggF=0,XS_VBF=0,mH=500" --afterFit --matrix --conditionalFit --nBins 60 --add-bkg-sys
```

* background-only fit, with qqZZ free. assuming qqZZ scaling factor is *mu_qqZZ*.

```bash
check_ws workspace.root --out_dir bkg_only_fit_qqZZFree -w combined -d obsData --poi_name XS_ggF --fixVar "XS_ggF=0,XS_VBF=0,mH=500" --floatVar "mu_qqZZ" --afterFit --matrix --conditionalFit --nBins 60 --add-bkg-sys
```

* signal+background fit

```bash
check_ws workspace.root --out_dir sig_plus_bkg_fit -w combined -d obsData --poi_name XS_ggF --fixVar "XS_VBF=0,mH=500" --afterFit --matrix --nBins 60 --add-bkg-sys
```

To obtain full help:

```text
check_ws

Usage: check_ws.py [options] file_name out_name

check yields and shape for WS

Options:
  -o OUT_DIR, --out_dir OUT_DIR
                        name of output directory
  --postfix POSTFIX     postfix used in naming histograms
  -v, --verbose         print info for debuging

  -w WS_NAME, --wsname=WS_NAME
  -m MC_NAME, --mcname=MC_NAME
  -d DATA_NAME, --dataname=DATA_NAME
                        name of observed data
  --poi_name POI_NAME   name of POI

  --fixVar FIXVAR       set variables as constant: mu=1,lumi=1
  --floatVar FLOATVAR   set variables float: mu_ggF,mu_VBF
  --signalScale SIGNALSCALE
                        scale factor applied to data
  --add-bkg-sys         add systematic uncertainties on background

  --afterFit            make plots and yields after fit
  --save-fitted-to SAVE_FITTED_TO
                        save fitted results to a ROOT file
  --save-snapshot-to SAVE_SNAPSHOT_TO
                        snapshot name of fitted results
  --conditionalFit      in conditional fit, POI is set to constant
  --matrix              plot covariance matrix

  --lumi LUMI           which luminosity used
  --plot-type {Internal,Preliminary, }
                        Internal or preliminary
  --noPlot              do not make plots
  --nBins NBINS         setup binning of the observable
  --xMax XMAX           max value of the observable
  --xMin XMIN           min value of the observable
  --logY                if log scale for y-axis
  --signal-mc SIGNAL_MC
                        signal MC names
```
