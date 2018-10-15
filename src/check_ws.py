#!/usr/bin/env python
"""
plot the shape from workspace, written in py2
"""

import ROOT
from optparse import OptionParser
from ROOT import RooFit
import math

import numpy as np
import os
import sys

from ploter import Ploter
from adder import adder

#sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)), '..')))

import root_utils as rtu

class WSReader:
    def __init__(self, file_name, out_name, options):
        self.fin = ROOT.TFile.Open(file_name)
        if not self.fin:
            print file_name,"does not exit"
            exit(1)
        self.out_name = out_name
        self.options = options
        self.colors = [2, 4, 64, 206, 95, 28, 29, 209, 5, 8, 6, 432, 433, 434, 435, 436]
        self.plotted_id = 0
        self.conditional_snapshot_name = "unconditionalFit"

        self.ps = Ploter("Internal", options.lumi)

        # for plotting the histograms
        self.hist_list = []
        self.tag_list = []
        if os.path.exists('correlation.root'):
            f1 = ROOT.TFile.Open("correlation.root")
            self.corr = f1.Get("correlation_matrix")
            self.corr.SetDirectory(0)
            f1.Close()
            print "Correlation matrix is ready"
            print self.corr.GetXaxis().GetNbins()

    def get_fit_name(self):
        if self.options.cond_fit:
            fit_name = "Cond"
        else:
            fit_name= "UnCond"
        return fit_name

    def get_fitted_file_name(self):
        file_name = os.path.basename(self.fin.GetName())
        file_name = options.poi_name+"_"+file_name
        fixed_names = self.get_fix_var_name()
        if fixed_names:
            file_name = fixed_names +file_name

        fitted_file_name = self.get_fit_name()+"_"+file_name
        return fitted_file_name

    def get_fix_var_name(self):
        if self.options.fixVar is not None:
            out = ""
            for poi_str in self.options.fixVar.split(','):
                # name,value = poi_str.split('=')
                # out += name+value+"_"
                out += "_".join(poi_str.split('='))
            return out
        else:
            return None

    def open_ws(self):
        ##  decide which workspace to use
        # if already fitted, reuse the existing workspace
        if self.options.after_fit and os.path.isfile(self.get_fitted_file_name()):
                fin_fitted = ROOT.TFile.Open(self.get_fitted_file_name())
                if fin_fitted:
                    print "find fitted workspace, use this one"
                    self.ws = fin_fitted.Get(self.options.ws_name)
                    self.ws.loadSnapshot(self.conditional_snapshot_name)
        else:
            self.ws = self.fin.Get(self.options.ws_name)

        ws = self.ws
        self.mc = ws.obj(self.options.mc_name)
        #self.simPdf = ws.obj(self.options.simpdfname)
        self.simPdf = self.mc.GetPdf()
        self.data = ws.obj(self.options.data_name)
        if not self.data:
            print self.options.data_name,"is not there"
            exit(1)
        else:
            self.data.Print()
        self.poi = ws.var(self.options.poi_name)
        if not self.poi:
            print "POI:", self.options.poi_name,"does not exit"
            exit(1)
        else:
            self.poi.Print()

        self.poi.setVal(self.options.sig_scale)
        self.observables = self.mc.GetObservables()


        self.categories = self.simPdf.indexCat()
        if self.data is not None:
            self.data_lists = self.data.split(self.categories, True)

    def fix_var(self):
        if self.options.fixVar:
            for poi_str in self.options.fixVar.split(','):
                name,value = poi_str.split('=')
                obj = self.ws.var(name)
                if obj:
                    print "set",name,"to",float(value)
                    obj.setVal(float(value))
                    obj.setConstant(True)
                else:
                    print name,"cannot be found"

    def get_hist(self, pdf, cat_name, obs, events, tag):
        #hist = pdf.createHistogram( "hist_"+tag, obs,
        #                          RooFit.Binning(obs.getBins(), obs.getMin(), obs.getMax()) )
        hist = pdf.createHistogram( "hist"+cat_name+tag, obs, RooFit.IntrinsicBinning(),
                                   RooFit.Extended(True))
        if hist.Integral() > 1E-6:
            hist.Scale(events/hist.Integral())
        else:
            print hist.GetName(),"MISSING"

        self.hist_list.append(hist)
        self.tag_list.append(tag)

    def plot_on(self, pdf, line_style, events, tag, leg_opt="L", do_norm=True):
        if do_norm:
            opt = RooFit.Normalization(events, ROOT.RooAbsReal.NumEvent)
        else:
            opt = ROOT.RooCmdArg.none()

        pdf.plotOn(self.frame, RooFit.LineStyle(line_style),
                   RooFit.LineColor(self.colors[self.plotted_id]),
                   opt
                  )
        self.legend.AddEntry(self.frame.getObject(self.plotted_id),
                             tag+"={:.2f}".format(events), leg_opt)
        self.plotted_id += 1

    def fit(self):
        """
        do the fit...and save a snapshot
        """
        if not hasattr(self, "simPdf") or\
           not hasattr(self, "data"):
            print "pdf nor data is missing, cannot fit"
            return None

        if self.ws.loadSnapshot(self.conditional_snapshot_name):
            return True

        if not self.options.cond_fit:
            self.poi.setConstant(False)

        if self.ws.var('muZZ'):
            self.ws.var('muZZ').setConstant(False)

        nuisance = self.mc.GetNuisanceParameters()
        nuisance.Print("v")
        nll = self.simPdf.createNLL(
            self.data,
            RooFit.Constrain(nuisance),
            RooFit.GlobalObservables(self.mc.GetGlobalObservables())
        )
        nll.enableOffsetting(True)
        minim = ROOT.RooMinimizer(nll)
        minim.optimizeConst(2)
        minim.setStrategy(2)
        #minim.setProfile()
        status = minim.minimize(
            "Minuit2",
            ROOT.Math.MinimizerOptions.DefaultMinimizerAlgo()
        )
        if status != 0:
            print "status is not zero!"

        self.ws.saveSnapshot(self.conditional_snapshot_name, self.ws.allVars())
        self.ws.writeToFile(self.get_fitted_file_name())

        # after performed the Fit, plot the coefficient matrix
        if self.options.matrix:
            fit_res = minim.save()
            corr_hist = fit_res.correlationHist()
            self.corr = corr_hist.Clone("corr")
            self.corr.SetDirectory(0)

            self.ps.plot_correlation(corr_hist, self.out_name+"_correlation_matrix", 0.05)
            fout = ROOT.TFile.Open("correlation.root", 'recreate')
            corr_hist.Write()
            fit_res.SetName("nll_res")
            fit_res.Write()
            fout.Close()


    def get_sys_for_NPs(self):
        """
        assuming all nuisance parameters are in best-fitted state...
        Fix one of them each time, then evaluate the total bakcground yields
        """
        ##self.sys_dict = {}
        self.sys_list = []
        itr = self.mc.GetNuisanceParameters()
        var = itr()
        while var:
            name = var.GetName()
            error = var.getError()

    def redefine_range(self, obs_first):
        #obs_first = obs.first()
        obs_nbins = obs_first.getBins()
        obs_max = obs_first.getMax()
        obs_min = obs_first.getMin()

        # re-define the range and binning
        if obs_nbins > self.options.n_bins:
            obs_nbins = self.options.n_bins

        if obs_max < self.options.x_max:
            obs_max = self.options.x_max

        if obs_min > self.options.x_min:
            obs_min = self.options.x_min

        obs_first.setMax(obs_max)
        obs_first.setMin(obs_min)
        obs_first.setBins(obs_nbins)

    def get_outplot_name(self, out_name):
        name_ = self.out_name+"_"+out_name
        if self.options.is_logY:
            name_ += "_Log"

        if self.options.after_fit:
            name_ += "_"+self.get_fit_name()

        return name_

    def get_output_histogram_name(self):
        return self.out_name+"_histograms.root"

    def loop_categories(self):
        m_debug = self.options.debug
        #f_out = ROOT.TFile.Open(self.out_name+"_histograms.root", "recreate")
        f_out = ROOT.TFile.Open(self.get_output_histogram_name(), "recreate")

        if self.options.after_fit:
            self.fit()

        iter_category = ROOT.TIter(self.categories.typeIterator())
        obj = iter_category()

        no_plot = options.no_plot
        yield_out_str = "\n"+self.fin.GetName()+"\n"

        nuisances = self.mc.GetNuisanceParameters()

        while obj:
            del self.hist_list[:]
            del self.tag_list[:]

            cat_name = obj.GetName()
            yield_out_str += "In category: {}\n".format(cat_name)
            self.plotted_id = 0 # for plots

            pdf = self.simPdf.getPdf(cat_name)
            obs = pdf.getObservables( self.observables )
            obs_var = obs.first()
            #redefine the range
            self.redefine_range(obs_var)

            ## get data first...
            hist_data = None
            if self.data is not None:
                data_ch = self.data_lists.At(obj.getVal())
                num_data =  data_ch.sumEntries()
                hist_data = ROOT.TH1F("data_"+cat_name, "data", obs_var.getBins(), obs_var.getMin(), obs_var.getMax())
                hist_data.SetXTitle(obs_var.GetName())
                data_ch.fillHistogram(hist_data, ROOT.RooArgList(obs_var))

            # get signal + background
            hist_splusb = ROOT.TH1F("splusb_"+cat_name, "data", obs_var.getBins(), obs_var.getMin(), obs_var.getMax())
            pdf.Print()
            hist_splusb = pdf.createHistogram("splusb_"+cat_name, obs_var, RooFit.IntrinsicBinning(), RooFit.Extended(True))
            spb_evts = pdf.expectedEvents(obs)
            if hist_splusb.Integral() > 1E-6:
                hist_splusb.Scale(spb_evts/hist_splusb.Integral())

            hist_splusb.SetLineColor(206)
            print "nbins for s-plus-b:", hist_splusb.GetNbinsX()

            # get background-only events
            old_poi_val = self.poi.getVal()
            print "OLD POI: " , old_poi_val
            self.poi.setVal(0.0)
            hist_bonly = ROOT.TH1F("bonly_"+cat_name, "data", obs_var.getBins(), obs_var.getMin(), obs_var.getMax())
            hist_bonly = pdf.createHistogram("bonly_"+cat_name, obs_var, RooFit.IntrinsicBinning(), RooFit.Extended(True))
            bonly_evts = pdf.expectedEvents(obs)
            if hist_bonly.Integral() > 1E-6:
                hist_bonly.Scale(bonly_evts/hist_bonly.Integral())


            # get background uncertainties
            bkg_sys = 0
            if self.corr:
                sys_list = []
                for ix in range(1, self.corr.GetXaxis().GetNbins()+1):
                    name = self.corr.GetXaxis().GetBinLabel(ix)
                    nom = self.ws.var(name).getVal()
                    error = self.ws.var(name).getError()
                    self.ws.var(name).setVal(nom + error)
                    up_bkg = pdf.expectedEvents(obs)
                    self.ws.var(name).setVal(nom - error)
                    down_bkg = pdf.expectedEvents(obs)
                    bkg_error = (up_bkg - down_bkg)/2.
                    sys_list.append(bkg_error)
                    self.ws.var(name).setVal(nom)

                print "Total Variations: ", len(sys_list)
                print sys_list
                sys_np = np.array(sys_list)
                correlation_matrix = rtu.TH2_to_numpy(self.corr)
                total_bkg_sys = math.sqrt(np.matmul(sys_np,
                                                    np.matmul(correlation_matrix, sys_np.transpose())
                                                   )
                                         )
                bkg_sys = total_bkg_sys/bonly_evts


            # get signal only spectrum
            if self.options.sigPDF:
                pdf_name = self.options.sigPDF+"_"+cat_name+"_cbga"
                sig_pdf = self.ws.obj(pdf_name.replace('Cat',''))
                hist_sonly = ROOT.TH1F("sonly_"+cat_name, "data", obs_var.getBins(), obs_var.getMin(), obs_var.getMax())
                hist_sonly = sig_pdf.createHistogram("sonly_"+cat_name, obs_var, RooFit.IntrinsicBinning(), RooFit.Extended(False))
                if hist_sonly.Integral() > 1E-6:
                    hist_sonly.Scale((spb_evts-bonly_evts)/hist_sonly.Integral())
                hist_sonly.SetLineColor(206)
            else:
                hist_sonly = hist_splusb.Clone("signalOnly_"+cat_name)
                hist_sonly.Add(hist_bonly, -1)

            #self.poi.setVal(old_poi_val)
            if "RooProdPdf" in pdf.ClassName():
                # break down each component for the PDF
                pdf_list = pdf.pdfList()
                this_pdf = None
                for pdf_index in range(pdf_list.getSize()):
                    ipdf = pdf_list[pdf_index]
                    if "RooGaussian" in ipdf.ClassName():
                        continue
                this_pdf = ipdf

                if this_pdf:
                    this_pdf_class_name = this_pdf.ClassName()
                    if "RooRealSumPdf" in this_pdf_class_name or\
                       "RooAddPdf" in this_pdf_class_name:
                        is_sum_pdf = "RooRealSumPdf" in this_pdf_class_name
                        if is_sum_pdf:
                            func_list = this_pdf.funcList()
                        else:
                            func_list = this_pdf.pdfList()
                        coeff_list = this_pdf.coefList()
                        total = 0

                        # define the function of calculate the number of events for each component
                        if is_sum_pdf:
                            nevts_func = lambda k:func_list[k].createIntegral(obs).getVal()*coeff_list[k].getVal()
                        else:
                            nevts_func = lambda k:coeff_list[k].getVal()

                        for func_index in sorted(
                            range(func_list.getSize()), key=nevts_func,
                            reverse=True
                        ):
                            func = func_list[func_index]
                            sum_ch = nevts_func(func_index)
                            total += sum_ch
                            if not no_plot and sum_ch > 1E-5:
                                simple_name = func.GetName().split('_')[2]
                                self.get_hist(func, cat_name, obs_var, sum_ch, simple_name) ## the normalization is included in tags
                            if m_debug:
                                print "{} {:.2f}".format(func.GetName(), sum_ch)
                            yield_out_str += "{} {:.2f}\n".format(func.GetName(), sum_ch)
                        yield_out_str += "total yields {:.2f}\n".format(total)
                    else:
                        print "no baseline pdf avaiable!"
                        continue

            else:
                print pdf.ClassName(),"should be RooProdPdf"

            self.poi.setVal(old_poi_val)

            if not no_plot:
                # make the plots
                self.ps.color(self.hist_list)
                sum_bkg, hs = self.ps.stack(self.hist_list)
                sum_bkg.SetLineColor(8)

                self.ps.prepare_2pad_canvas("canvas", 600, 600)
                self.ps.pad2.cd()
                self.ps.add_ratio_panel([hist_data, sum_bkg], "Data/MC", 0.50, 2.0)
                # add uncertainty band for background...
                if bkg_sys > 0:
                    print "background systematic: ", bkg_sys
                    gr = adder.add_band(sum_bkg, 1., bkg_sys, False)
                    #gr.Print()
                    gr.Draw("CF SAME")


                self.ps.pad1.cd()
                self.ps.get_offset(hist_splusb)

                self.ps.set_y_range([hist_data, hist_splusb], options.is_logY)
                #this_hist = self.ps.set_y_range(hist_data, hist_splusb, options.is_logY)
                hist_data.SetXTitle(obs_var.GetName())
                hist_data.SetYTitle("Events / 20 GeV")

                if hist_data:
                    legend = self.ps.get_legend(len(self.hist_list) + 3)
                    legend.AddEntry(hist_data, "Data {:.0f}".format(hist_data.Integral()), "LP")
                    hist_data.SetMarkerStyle(20)
                    hist_data.SetMarkerSize(1.2)
                    hist_data.Draw("EP")
                    #hist_splusb.Draw("HIST same")
                    hs.Draw("HIST same")
                    hist_data.Draw("AXISsame")
                    hist_data.Draw("EPsame")
                    hist_sonly.Draw("HIST same")
                else:
                    legend = self.ps.get_legend(len(self.hist_list) + 2)
                    #hist_splusb.Draw("HIST")
                    hs.Draw("HIST")
                    hist_sonly.Draw("HIST same")

                legend.AddEntry(hist_splusb, "S({:.1f})+B {:.1f}" .format(self.poi.getVal(), hist_splusb.Integral()), "L")
                sum_bkg.SetLineColor(0)
                sum_bkg.SetFillColor(0)
                legend.AddEntry(sum_bkg, "Total Bkg {:.0f}".format(sum_bkg.Integral()), "F")

                for hist, tag in zip(self.hist_list, self.tag_list):
                    legend.AddEntry(hist, tag+" {:.1f}".format(hist.Integral()), "F")

                legend.Draw("same")
                self.ps.add_atlas()
                if "2017" in cat_name:
                    self.ps.add_lumi(43.7)
                else:
                    self.ps.add_lumi(36.1)
                out_plot_name = self.get_outplot_name(cat_name)
                self.ps.can.SaveAs(out_plot_name+".pdf")

            # start next category
            obj = iter_category()
            # save the histograms..
            f_out.cd()
            hist_splusb.Write()
            if hist_data:
                hist_data.Write()
            for hist in self.hist_list:
                hist.Write()

        print yield_out_str
        with open("yield.log", 'a') as f:
            f.write(yield_out_str)

        f_out.Close()


if __name__ == "__main__":
    usage = "%prog [options] file_name out_name"
    version="%prog 1.0"
    parser = OptionParser(usage=usage, description="check yields and shape for WS", version=version)
    parser.add_option("-w", "--wsname", dest='ws_name', default='combined')
    parser.add_option("-m", '--mcname', dest='mc_name', default='ModelConfig')
    parser.add_option("-d", '--dataname', dest='data_name', default='obsData', help="name of observed data")
    parser.add_option("--poi_name", dest='poi_name', help="name of POI", default="SigXsecOverSM")

    parser.add_option("--fixVar", help="set variables as constant:  mu=1,lumi=1", default=None)
    parser.add_option("--noPlot", dest='no_plot', help="don't make plots", default=False, action="store_true")
    parser.add_option("--nBins", dest='n_bins', help="setup binning of the observable", default=10000, type="int")
    parser.add_option("--xMax", dest='x_max', help="max value of the observable", default=0, type="float")
    parser.add_option("--xMin", dest='x_min', help="min value of the observable", default=100000, type="float")
    parser.add_option("--logY", dest='is_logY', help="if log scale for y-axis", default=False, action="store_true")
    parser.add_option("--signalScale", dest='sig_scale', help="scale factor applied to data", default=10., type="float")
    parser.add_option("--afterFit", dest='after_fit', help="make plots and yields after fit", default=False, action="store_true")
    parser.add_option("--conditionalFit", dest='cond_fit', help="perform conditional fit", default=False, action="store_true")

    parser.add_option("--lumi", dest='lumi', help="which luminosity used",  default=36.1, type='float')
    parser.add_option("--sigPdfName", dest='sigPDF', help="signal pdf name",  default=None)

    parser.add_option("--matrix", dest='matrix', help="plot covariance matrix",  default=False, action='store_true')
    parser.add_option("--debug", dest='debug', help="in debug mode", action="store_true", default=False)

    (options,args) = parser.parse_args()

    if len(args) < 1:
        print parser.print_help()
        exit(1)

    if len(args) > 1:
        out_name = args[1]
    else:
        out_name = "test"

    ROOT.gROOT.SetBatch()
    import AtlasStyle

    ws_reader = WSReader(args[0], out_name, options)
    ws_reader.open_ws()
    ws_reader.fix_var()
    ws_reader.loop_categories()
