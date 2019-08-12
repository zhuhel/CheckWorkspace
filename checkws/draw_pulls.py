#!/usr/bin/env python

from ploter import Ploter
from adder import adder
import ROOT
ROOT.gROOT.SetBatch()
import AtlasStyle
from optparse import OptionParser

usage = "%prog [options]"
version = "%prog 1.0"
parser = OptionParser(usage=usage, description="check yields and shape for WS", version=version)
parser.add_option("-m", '--message', dest='message', default='vvqq bkg')
parser.add_option("-d", '--out_dir', dest='out_dir', default='.')
(options,args) = parser.parse_args()

corfile = options.out_dir + "/correlation.root"
f1 = ROOT.TFile.Open(corfile)
res = f1.Get("nll_res")

fitted_list = res.floatParsFinal()
iterator = fitted_list.createIterator()
obj = iterator()

pulls = []
n_skipped = 0
while obj:
    # skip gamma terms, not of interest
    np_name = obj.GetName()

    # reduced np name
    np_name = np_name.replace("alpha_","")
    np_name = np_name.replace("Sys","")
    np_name = np_name.replace("JET_SR1_","")
    if "gamma_stat" in np_name:
        pass
    else:
        nom = obj.getVal()
        error_hi = obj.getErrorHi()
        error_lo = obj.getErrorLo()
        if abs(nom)  < 1E-2 and abs(error_hi-1) < 1E2:
            n_skipped += 1
        else:

            np_info = (np_name, nom, error_hi, error_lo)
            pulls.append(np_info)
        #print np_name, obj.getVal(), obj.getError()
    obj = iterator()

print "skipped NPs:", n_skipped
sorted_list = sorted(pulls, key=lambda x: (abs(x[2])+abs(x[3]))/2)
n_nps = len(sorted_list)


can = ROOT.TCanvas("canvas", "canvas", 800, 600)
can.SetBottomMargin(0.35);
can.SetGridy()
ps = Ploter()
print "nuisances:",n_nps

dummy = ROOT.TH1F("dummy", "dummy", n_nps, 0, n_nps)
dummy.GetYaxis().SetRangeUser(-2.8, 2.8)
dummy.GetXaxis().SetLabelSize(0.020)
x_axis = dummy.GetXaxis()
for ibin in range(n_nps):
    x_axis.SetBinLabel(ibin+1, sorted_list[ibin][0])

dummy.GetXaxis().LabelsOption('v')
dummy.Draw()
gr, gr_one, gr_two = ps.create_graph_pulls(sorted_list)
gr_one.SetFillStyle(1001)
gr_one.SetFillColor(ROOT.kGreen)
gr_two.SetFillStyle(1001)
gr_two.SetFillColor(ROOT.kYellow)
gr_two.Draw("3")
gr_one.Draw("3")
gr.Draw("P")
dummy.Draw("same")
dummy.Draw("Axis Same")

adder.add_text(0.2, 0.88, 1, options.message)

can.SaveAs("pulls.pdf")

f1.Close()
