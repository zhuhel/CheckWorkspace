# CheckWorkspace

A collection of scripts to check workspaces.

```bash
source setup.sh
```

* *draw_pulls.py* is to plot pull distributions of nuisance parameters
* *print_np.py* is to print nuisance paramters and parameter of interest in the workspace
* *check_ws.py* is to check the workspaces... A typical example to run the combined workspace of 4l+llvv is:

```bash
cd src
python check_ws.py workspace.root bkg_only_fit -w combWS -d combData --poi_name XS_ggF --fixVar "XS_ggF=0,mH=500" --afterFit --matrix --conditionalFit --nBins 60
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
  --sigPdfName=SIGPDF   signal pdf name
  --matrix              plot covariance matrix
  --debug               in debug mode
```  
 
