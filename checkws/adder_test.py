#!/usr/bin/env python

import ROOT
import random
random.seed(100)

def test_ratio():
    """
    TH1.Divide, assuming the errors are uncorrelated.
    In case of the ratio of the histgoram itself,
    the error in ratio = sqrt(2/bin-content).

    SetBinError changes Sumw2 vector
    """
    nbins = 10
    h1 = ROOT.TH1F("h1", "h1", nbins, -5, 5)
    nentries = 1000
    mu = 0.0
    sigma = 1.0
    for ientry in range(nentries):
        h1.Fill(random.gauss(mu, sigma))

    h1.Sumw2()
    # h1_sum_w2 = h1.GetSumw2()
    # for ibin in range(nbins):
    #     print(ibin+1, h1_sum_w2[ibin+1])

    # for ibin in range(nbins):
    #     h1.SetBinError(ibin+1, 1.)

    # for ibin in range(nbins):
    #     print(ibin+1, h1_sum_w2[ibin+1])

    ratio = h1.Clone("ratio")
    ratio.Divide(h1)
    sum_w2 = ratio.GetSumw2()
    for ibin in range(nbins):
        print("{} {:.2f} {:.2f} {:.2f} {:.3f}".format(
            ibin+1,
            h1.GetBinContent(ibin+1),
            sum_w2[ibin+1],
            ratio.GetBinContent(ibin+1),
            ratio.GetBinError(ibin+1)))

if __name__ == "__main__":
    test_ratio()