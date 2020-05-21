#!/usr/bin/env python
__author__ = "Heling Zhu"
__version__ = "0.2"
import ROOT
import AtlasStyle
import sys

import os
import errno

#####################
## Import Module  ###
#####################
from array import array
import glob
from math import sqrt,fabs,sin,log
from ROOT import TFile,TTree,TChain,TBranch,TH1,TH1F,TList
from ROOT import TLorentzVector,TGraphAsymmErrors,TMath
from ROOT import THStack,TCanvas,TLegend,TColor,TPaveText,TPad
from ROOT import gStyle,gDirectory,gROOT


def plot_matrix(corr_hist, out_name, labels, **kw):
    m_can = ROOT.TCanvas("canvas", "canvas", 1200, 800)
    m_can.SetRightMargin(0.1)
    m_can.SetLeftMargin(0.18)
    gROOT.SetBatch()

    if len(labels)>0: h_new = remove_bins(corr_hist, labels)
    else: h_new = remove_weak_correlation_bins(corr_hist, 0.05) 

    #ROOT.gStyle.SetPalette(55)
    # palette = array('d', [15, 20, 23, 30, 32])
    # ROOT.gStyle.SetPalette(5, palette)
    set_palette()
    h_new.Draw("colz text")
    h_new.SetMarkerSize(1.2)
    ROOT.gStyle.SetPaintTextFormat("3.2f")
    h_new.GetXaxis().SetLabelSize(0.035)
    h_new.GetXaxis().LabelsOption('u')
    h_new.GetYaxis().SetLabelSize(0.035)
    h_new.GetZaxis().SetLabelSize(0.035)
    h_new.GetZaxis().SetRangeUser(-1, 1)
    h_new.GetXaxis().SetTickSize(0.)
    h_new.GetYaxis().SetTickSize(0.)

    m_can.SaveAs(out_name+".eps")
    m_can.SaveAs(out_name+".pdf")

def set_palette(**kw):
    red_list = [1., 1., 0.]
    white_list = [0., 1.0, 0.0]
    blue_list = [0., 1., 1.]
    length_list = [0.0, 0.5, 1.0]
    nb = 50
    ROOT.TColor.CreateGradientColorTable(
        3,
        array('d', length_list),
        array('d', red_list),
        array('d', white_list),
        array('d', blue_list),
        nb)


def remove_bins(h2d, labels, **kw):
    empty_xbins = []
    for xbin in range(h2d.GetNbinsX()):
        is_empty = True
        xlabel = h2d.GetXaxis().GetBinLabel(xbin+1)
        if xlabel in labels:
            is_empty = False
        if is_empty:
            empty_xbins.append(xbin+1)

    final_bins = h2d.GetNbinsX() - len(empty_xbins)
    org_bins = h2d.GetNbinsX()
    print "original bins", org_bins
    print "final bins: ", final_bins
    #print empty_xbins
    if True:
        h2d_new = ROOT.TH2D("reduced_correlation", "reduced_correlation",
                            final_bins, 0.5, final_bins+0.5,
                            final_bins, 0.5, final_bins+0.5)
        new_x = 0
        for xbin in range(h2d.GetNbinsX()):
            if (xbin+1) in empty_xbins:
                continue
            new_x += 1
            new_y = 0
            for ybin in range(h2d.GetNbinsY()):
                if (org_bins - ybin) in empty_xbins:
                    continue
                new_y += 1
                value = h2d.GetBinContent(xbin+1, ybin+1)
                # print xbin,ybin,new_x,new_y, value

                h2d_new.SetBinContent(new_x, new_y, value)
                h2d_new.GetXaxis().SetBinLabel(new_x, h2d.GetXaxis().GetBinLabel(xbin+1))
                h2d_new.GetYaxis().SetBinLabel(new_y, h2d.GetYaxis().GetBinLabel(ybin+1))

        return h2d_new

def remove_weak_correlation_bins(h2d, thre=0.05, **kw):
    empty_xbins = []
    for xbin in range(h2d.GetNbinsX()):
        is_empty = True
        for ybin in range(h2d.GetNbinsY()):
            value = h2d.GetBinContent(xbin+1, ybin+1)
            if abs(value) > thre and abs(value - 1) > 1E-3 :
                is_empty = False
                break
        if is_empty:
            empty_xbins.append(xbin+1)

    final_bins = h2d.GetNbinsX() - len(empty_xbins)
    org_bins = h2d.GetNbinsX()
    print "original bins", org_bins
    print "final bins: ", final_bins
    #print empty_xbins
    if True:
        h2d_new = ROOT.TH2D("reduced_correlation", "reduced_correlation",
                            final_bins, 0.5, final_bins+0.5,
                            final_bins, 0.5, final_bins+0.5)
        new_x = 0
        for xbin in range(h2d.GetNbinsX()):
            if (xbin+1) in empty_xbins:
                continue
            new_x += 1
            new_y = 0
            for ybin in range(h2d.GetNbinsY()):
                if (org_bins - ybin) in empty_xbins:
                    continue
                new_y += 1
                value = h2d.GetBinContent(xbin+1, ybin+1)
                # print xbin,ybin,new_x,new_y, value

                h2d_new.SetBinContent(new_x, new_y, value)
                h2d_new.GetXaxis().SetBinLabel(new_x, h2d.GetXaxis().GetBinLabel(xbin+1))
                h2d_new.GetYaxis().SetBinLabel(new_y, h2d.GetYaxis().GetBinLabel(ybin+1))

        return h2d_new

        return None


###################
## Main Function ##
###################
if __name__ == "__main__":

    if len(sys.argv) == 2:
        inputdir=sys.argv[1]
    else:
        raise RuntimeError('One and only one arg needed: root file path')

    inputf=inputdir+"/correlation.root"
    #labels=["mu_llll_qqZZ", "mu_llll_ggZZ", "mu_llvv_qqZZ", "mu_llvv_ggZZ"]
    labels=[]
    #out_name=inputdir+"/out_corrMatrix_scaleF"
    out_name=inputdir+"/correlation_matrix"

    fin=TFile(inputf);
    #corr_hist=gDirectory.Get("correlation_matrix")
    corr_hist=gDirectory.Get("corr")

    plot_matrix(corr_hist, out_name, labels)

