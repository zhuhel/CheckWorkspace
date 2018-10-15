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
