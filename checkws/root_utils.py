# -*- coding: utf-8 -*-
import ROOT
import numpy as np

def TH2_to_numpy(h2d):
    array = []
    x_range = range(1, h2d.GetXaxis().GetNbins()+1)
    y_range = range(h2d.GetYaxis().GetNbins(), 0, -1)

    for iy in y_range:
        x_array = [h2d.GetBinContent(ix, iy) for ix in x_range]
        array.append(x_array)

    return np.array(array)


def get_diagnal(array):
    if array.shape[0] == array.shape[1]:
        return [ array[i][i] for i in range(array.shape[0]) ]
    else:
        return None


def redefine_range(obs, nbins, max_var, min_var):
    obs_nbins = nbins if obs.getBins() > nbins else obs.getBins()
    obs_max = max_var if obs.getMax() < max_var else obs.getMax()
    obs_min = min_var if obs.getMin() > min_var else obs.getMin()
    
    obs.setMax(obs_max)
    obs.setMin(obs_min)
    obs.setBins(obs_nbins)

def integral_and_error(hist):
    error = ROOT.double(0.0)
    evts = hist.IntegralAndError(1, hist.GetNbinsX(), error)
    return (evts, error)

def print_hist(hist, add_sys):
    if add_sys:
        evts, error = integral_and_error(hist)
        out_str = "{:.2f} \pm {:.2f}".format(evts, error)
    else:
        out_str = "{:.2f}".format(hist.Integral())
    return out_str