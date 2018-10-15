# CheckWorkspace

A collection of scripts to check workspaces.

```bash
source setup.sh
```

*print_np.py* is to print nuisance paramters and parameter of interest in the workspace
*check_ws.py* is to check the workspaces... A typical example to run the combined workspace of 4l+llvv is:
```
python check_ws.py workspace.root bkg_only_fit -w combWS -d combData --poi_name XS_ggF --fixVar "XS_ggF=0,mH=500" --afterFit --matrix --conditionalFit --nBins 60
```
