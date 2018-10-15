#!/usr/bin/env python

import ROOT
from array import array

class maker:
    """
    make ROOT objects from python objets
    """
    def __init__(self):
        pass

    @staticmethod
    def graph(name, x, y):
        gr = ROOT.TGraph(len(x), array('f', x), array('f', y))
        gr.SetName(name)
        return gr

    @staticmethod
    def graph_error(name, x, xe, y, ye):
        """
        Symmetric errors
        """
        gr = ROOT.TGraphErrors(
            len(x), array('f', x), array('f', y),
            array('f', xe), array('f', ye)
        )
        gr.SetName(name)
        return gr

    @staticmethod
    def unequal_bin_hist(hist_name, bin_list):
        """
        create TH1F using a list as x-axis
        """
        nbins = len(bin_list) -1
        h1 = ROOT.TH1F(hist_name, hist_name, nbins, array('f', bin_list))
        return h1
