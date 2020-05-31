#!/usr/bin/env python
import ROOT
from array import array
import math
from checkws.root_utils import sumw2

def add_text(x, y, color, text, size=0.05, font=42):
    l = ROOT.TLatex()
    l.SetTextSize(size)
    l.SetNDC()
    l.SetTextColor(color)
    l.SetTextFont(font)
    l.DrawLatex(x, y, text)


def add_line(hist, y_val, color=1, style=2, option="x"):
    x_low = hist.GetBinLowEdge(hist.GetXaxis().GetFirst())
    x_hi = hist.GetBinLowEdge(hist.GetXaxis().GetLast()+1)
    y_low = hist.GetBinLowEdge(hist.GetYaxis().GetFirst())
    y_hi = hist.GetBinLowEdge(hist.GetYaxis().GetLast()+1)
    line = ROOT.TLine()
    line.SetLineColor(color)
    line.SetLineStyle(style)
    line.SetLineWidth(2)
    if option.lower() == "x":
        line.DrawLine(x_low, y_val, x_hi, y_val)
    else:
        line.DrawLine(y_val, y_low, y_val, y_hi)


def make_legend(x1, y1, x2, y2):
    legend = ROOT.TLegend(x1, y1, x2, y2)
    legend.SetBorderSize(0)
    legend.SetFillColor(0)
    legend.SetTextFont(42)
    legend.SetTextSize(0.04)
    return legend
    

def make_error_band(hist, center, width, add_stats=True, scale=1.):
    x = array('d')
    y = array('d')
    up = array('d')
    down = array('d')
    weight = hist.Integral()/hist.GetEntries()
    for i in range(hist.GetXaxis().GetNbins()):
        ibin = i+1
        content = hist.GetBinContent(ibin)/weight
        if add_stats and content < 1E-10:
            new_width = width
        elif add_stats:
            tot_variance = width**2 + scale/content
            new_width = math.sqrt(tot_variance)
        else:
            new_width = width

        x.append(hist.GetXaxis().GetBinCenter(ibin))
        y.append(center)
        up.append(center + new_width)
        down.append(max(center - new_width, 0))

    n = len(x)
    grband = ROOT.TGraph(2*n)
    for i in range(n):
        grband.SetPoint(i, x[i], up[i])
        grband.SetPoint(n+i, x[n-i-1], down[n-i-1])

    grband.SetFillStyle(3013)
    grband.SetFillColor(16)
    #grband.Draw("F SAME")
    return grband

def make_self_ratio_band(hist):
    return make_error_band(hist, 1., 0, scale=2.)

def make_error_band2(hist, center, width, scale=1.):
    x = array('d')
    y = array('d')
    up = array('d')
    down = array('d')
    sumw2(hist)
    sum_w2 = hist.GetSumw2()
    for i in range(hist.GetXaxis().GetNbins()):
        ibin = i+1
        content = hist.GetBinContent(ibin)
        if content < 1E-10:
            new_width = width
        else:
            tot_variance = width**2 + scale*sum_w2[ibin]
            new_width = math.sqrt(tot_variance)

        x.append(hist.GetXaxis().GetBinCenter(ibin))
        y.append(center)
        up.append(center + new_width)
        down.append(max(center - new_width, 0))

    n = len(x)
    grband = ROOT.TGraph(2*n)
    for i in range(n):
        grband.SetPoint(i, x[i], up[i])
        grband.SetPoint(n+i, x[n-i-1], down[n-i-1])

    grband.SetFillStyle(3013)
    grband.SetFillColor(16)
    return grband

def make_error_band3(hist, center, width, scale=1.):
    x = array('d')
    y = array('d')
    up = array('d')
    down = array('d')
    name = hist.GetName()
    name = name + "_relative_syst"
    ## hist
    if isinstance(hist, ROOT.TH1):
      sumw2(hist)
      h_err = ROOT.TGraphAsymmErrors(hist)
      h_err.SetName(name)
      h_err.SetTitle(name)
      for i in range(hist.GetXaxis().GetNbins()):
        ibin = i+1
        x = hist.GetBinCenter(ibin)
        content = hist.GetBinContent(ibin)
        e_up = hist.GetBinErrorUp(ibin)
        e_down = hist.GetBinErrorLow(ibin)
        r_up = 0
        r_down = 0
        if content!=0:
          r_up = e_up/content
          r_down = e_down/content

        h_err.SetPoint(i, x, 1.)
        ex_up = hist.GetBinWidth(ibin)*0.5
        ex_dn = ex_up
        h_err.SetPointEXhigh(i, ex_up)
        h_err.SetPointEXlow(i, ex_dn)
        h_err.SetPointEYhigh(i, r_up)
        h_err.SetPointEYlow(i, r_down)

        
    ## Tgraph
    if isinstance(hist, ROOT.TGraphAsymmErrors):
      h_err = hist.Clone()
      h_err.SetName(name)
      h_err.SetTitle(name)
      for i in range(hist.GetN()):
        ibin = i
        x, y = ROOT.Double(), ROOT.Double()
        hist.GetPoint(ibin, x, y)
        ex_up, ex_dn, ey_up, ey_dn = hist.GetErrorXhigh(ibin), hist.GetErrorXlow(ibin), hist.GetErrorYhigh(ibin), hist.GetErrorYlow(ibin)
        r_up = 0
        r_dn = 0
        content=y
        if content!=0:
          r_up = ey_up/content
          r_dn = ey_dn/content
        h_err.SetPoint(i, x, 1.)
        h_err.SetPointEXhigh(ibin, ex_up)
        h_err.SetPointEXlow(ibin, ex_dn)
        h_err.SetPointEYhigh(ibin, r_up)
        h_err.SetPointEYlow(ibin, r_dn)

    h_err.SetLineColor(1)
    h_err.SetMarkerStyle(0)
    h_err.SetFillStyle(3004)
    h_err.SetFillColor(1)
    return h_err
