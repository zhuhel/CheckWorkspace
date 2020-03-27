"""
Read Workspaces and make histograms/yields
"""
from __future__ import print_function
import ROOT
from ROOT import RooFit

ROOT.gErrorIgnoreLevel = ROOT.kInfo
ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.ERROR)

import os

class WSReader:
    def __init__(self, 
                file_name,  # input workspace
                out_dir,    # output directory
                ws_name='combined',
                mc_name='ModelConfig',
                data_name='obsData',
                poi_name='XS_ggF',
                ):
        if not os.path.exists(file_name):
            print(file_name, "not there")

        self.fin = ROOT.TFile.Open(file_name)
        self.ws = self.fin.Get(ws_name)
        self.mc = self.ws.obj(mc_name)
        self.simPdf = self.mc.GetPdf()
        self.categories = self.simPdf.indexCat()

        self.data = self.ws.obj(data_name)
        if not self.data:
            print("{} not there".format(data_name))
            self.data_lists = None
        else:
            self.data.Print()
            self.data_lists = self.data.split(self.categories, True)
        
        self.poi = self.ws.var(poi_name)
        if not self.poi:
            print("POI: {} not there".format(poi_name))
        else:
            self.poi.Print()

        self.observables = self.mc.GetObservables()
        # hardcoded outputs
        self.out_dir = out_dir
        self.corr_root_path = os.path.join(self.out_dir, 'correlation.root')
        self.mc.GetParametersOfInterest().Print()

        # figure out the category anmes and mc samples
        self.category_names = []
        self.mc_samples = []
        iter_category = ROOT.TIter(self.categories.typeIterator())
        obj = iter_category()
        while obj:
            cat_name = obj.GetName()
            self.category_names.append(cat_name)
            pdf = self.simPdf.getPdf(cat_name)
            # break down the PDF into small components
            if "RooProdPdf" in pdf.ClassName():
                pdf_list = pdf.pdfList()
                this_pdf = None
                # skip all Gaussian constraints
                for ipdf in pdf_list:
                    if "RooGaussian" in ipdf.ClassName():
                        continue
                    else:
                        this_pdf = ipdf
                        break
                if this_pdf:
                    this_pdf_class_name = this_pdf.ClassName()
                    if "RooRealSumPdf" in this_pdf_class_name or\
                       "RooAddPdf" in this_pdf_class_name:
                        is_sum_pdf = "RooRealSumPdf" in this_pdf_class_name
                        if is_sum_pdf:
                            func_list = this_pdf.funcList()
                        else:
                            func_list = this_pdf.pdfList()
                        for func_index in range(func_list.getSize()):
                            func = func_list[func_index]
                            simple_name = func.GetName().split('_')[2]
                            self.mc_samples.append(simple_name)
                    else:
                        continue
            elif "RooAddPdf" in pdf.ClassName():
                func_list = pdf.pdfList()
                for func_index in range(func_list.getSize()):
                    func = func_list[func_index]
                    simple_name = func.GetName().split('_')[2]
                    self.mc_samples.append(simple_name)
            else:
                print(pdf.ClassName()," should be either RooProdPdf or RooAddPdf")
            obj = iter_category()
        self.mc_samples = list(set(self.mc_samples))

    def fix_var(self, var_name, var_value):
        obj = self.ws.var(var_name)
        if obj:
            print("set {} to {:.4f}".format(var_name, var_value))
            obj.setVal(var_value)
            obj.setConstant(True)
        else:
            print(var_name, "cannot be found")

    def fix_var_str(self, var_str):
        for poi_str in var_str.split(','):
            name, value = poi_str.strip().split('=')
            self.fix_var(name, float(value))

    def float_var(self, var_name):
        obj = self.ws.var(var_name.strip())
        if obj:
            print("set {} floating".format(var_name))
            obj.setConstant(False)
        else:
            print(var_name, "cannot be found")

    def float_var_str(self, var_str):
        map(self.float_var, var_str.split(','))

    @staticmethod
    def create_hist_from_pdf(pdf, hist_name, obs):
        hist = pdf.createHistogram(hist_name, obs, RooFit.IntrinsicBinning(), RooFit.Extended(True))
        events = pdf.expectedEvents(ROOT.RooArgSet(obs))
        if hist.Integral() > 1E-6:
            if hist.GetSumw2 is None:
                hist.Sumw2(True)
            hist.Scale(events/hist.Integral())
        else:
            print(hist.GetName(), "MISSING")
        return hist


    def get_hist(self, pdf, obs, events, hist_name):
        hist = pdf.createHistogram(hist_name, obs, RooFit.IntrinsicBinning(),
                                   RooFit.Extended(True))
        if hist.Integral() > 1E-6:
            if hist.GetSumw2 is None:
                hist.Sumw2(True)
            hist.Scale(events/hist.Integral())
        else:
            print(hist.GetName(), "MISSING")
        return hist

    def fit(self, snapshot_name, do_matrix=True, save_to_file=None):
        """
        do the fit...and save a snapshot
        """
        if not self.data:
            print("data is missing")
            return

        if self.ws.loadSnapshot(snapshot_name):
            print("Reuse fitted results saved at snapshot: {}".format(snapshot_name))
            return True

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
            print("status is not zero!")

        self.ws.saveSnapshot(snapshot_name, self.ws.allVars())
        if save_to_file:
            print("debug", self.out_dir, save_to_file)
            #self.ws.writeToFile(os.path.join(self.out_dir, save_to_file))

        # after performed the Fit, plot the coefficient matrix
        if do_matrix:
            fit_res = minim.save()
            corr_hist = fit_res.correlationHist()
            self.corr = corr_hist.Clone("corr")
            self.corr.SetDirectory(0)
            fout = ROOT.TFile.Open(self.corr_root_path, 'recreate')
            self.corr.Write()
            fit_res.SetName("nll_res")
            fit_res.Write()
            fout.Close()


    def loop_categories(self, postfix='nominal', force_update=False):
        """
        make histograms for each component, save them into file [optional]
        print and save the yields
        """
        outname = os.path.join(self.out_dir, "histograms_{}.root".format(postfix))
        if os.path.exists(outname):
            print("{} is there, ".format(outname))
            if not force_update:
                return outname
            else:
                print("forcing updates. May the observable range be changed?")

        f_out = ROOT.TFile.Open(outname, "recreate")

        iter_category = ROOT.TIter(self.categories.typeIterator())
        obj = iter_category()

        yield_out_str = "\n"+self.fin.GetName()+"\n"
        self.hist_file_name = outname
        while obj:
            all_histograms = []
            cat_name = obj.GetName()
            yield_out_str += "In category: {}\n".format(cat_name)

            pdf = self.simPdf.getPdf(cat_name)
            obs = pdf.getObservables( self.observables )
            obs_var = obs.first()
            obs_bins = (obs_var.getBins(), obs_var.getMin(), obs_var.getMax())

            ## get data first...
            hist_data = None
            if self.data:
                data_ch = self.data_lists.At(obj.getVal())
                hist_data = ROOT.TH1F("data-"+cat_name+"-"+postfix, "data", *obs_bins)
                hist_data.SetXTitle(obs_var.GetName())
                data_ch.fillHistogram(hist_data, ROOT.RooArgList(obs_var))
                all_histograms.append(hist_data)

            # get signal + background
            print("get signal + background")
            hist_splusb = ROOT.TH1F("splusb-"+cat_name+'-'+postfix,
                                    "signal + background", *obs_bins)
            hist_splusb = self.create_hist_from_pdf(pdf, hist_splusb.GetName(), obs_var)
            hist_splusb.SetLineColor(206)
            all_histograms.append(hist_splusb)
            hist_splusb.Print()

            # get background-only events
            print("get background-only events")
            old_poi_val = self.poi.getVal()
            self.poi.setVal(1e-10)
            hist_bonly = ROOT.TH1F("bonly-"+cat_name+"-"+postfix, "background only ", *obs_bins)
            hist_bonly = self.create_hist_from_pdf(pdf, hist_bonly.GetName(), obs_var)
            bonly_evts = hist_bonly.Integral()
            all_histograms.append(hist_bonly)
            hist_bonly.Print()

            # get signal only spectrum
            print("get signal only spectrum")
            hist_sonly = hist_splusb.Clone("signalOnly-"+cat_name+"-"+postfix)
            hist_sonly.Add(hist_bonly, -1)
            hist_sonly.SetLineColor(205)
            all_histograms.append(hist_sonly)
            hist_sonlypdf = ROOT.TH1F("sonlyPdf-"+cat_name+"-"+postfix, "signal only from PDF", *obs_bins)
            self.poi.setVal(old_poi_val)

            # break down the PDF into small components
            print ("break down the PDF into small components")
            if "RooProdPdf" in pdf.ClassName():
                pdf_list = pdf.pdfList()
                this_pdf = None
                # skip all Gaussian constraints
                for ipdf in pdf_list:
                    if "RooGaussian" in ipdf.ClassName():
                        continue
                    else:
                        this_pdf = ipdf
                        break

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
                            simple_name = func.GetName().split('_')[2]
                            if sum_ch > 1E-5 and "signal" not in func.GetName().lower():
                                hist_name = simple_name+"-"+cat_name+"-"+postfix
                                histogram = self.get_hist(func, obs_var, sum_ch, hist_name)
                                all_histograms.append(histogram)

                            ## signal pdf
                            if sum_ch > 1E-5 and "signal" in func.GetName().lower():
                                hist_name = "Signal-"+cat_name+"-"+postfix
                                hist_sonlypdf=self.get_hist(func, obs_var, sum_ch, hist_name)
                                all_histograms.append(hist_sonlypdf)
                            yield_out_str += "{} {:.2f}\n".format(func.GetName(), sum_ch)
                        yield_out_str += "total yields {:.2f}\n".format(total)
                    else:
                        print("no baseline pdf avaiable!")
                        continue

            elif "RooAddPdf" in pdf.ClassName():
                func_list = pdf.pdfList()
                coeff_list = pdf.coefList()
                total = 0
                nevts_func = lambda k:coeff_list[k].getVal()
                for func_index in sorted(
                    range(func_list.getSize()), key=nevts_func,
                    reverse=True
                ):
                    func = func_list[func_index]
                    sum_ch = nevts_func(func_index)
                    total += sum_ch
                    simple_name = func.GetName().split('_')[2]
                    if sum_ch > 1E-5 and "signal" not in func.GetName().lower():
                        hist_name = simple_name+"-"+cat_name+"-"+postfix
                        histogram = self.get_hist(func, obs_var, sum_ch, hist_name)
                        all_histograms.append(histogram)
                    yield_out_str += "{} {:.2f}\n".format(func.GetName(), sum_ch)
                yield_out_str += "total yields {:.2f}\n".format(total)
            else:
                print(pdf.ClassName()," should be either RooProdPdf or RooAddPdf")
            self.poi.setVal(old_poi_val)
            # start next category
            obj = iter_category()
            # save the histograms..
            f_out.cd()
            for histograms in all_histograms:
                histograms.Write()

        # print(yield_out_str)
        with open(os.path.join(self.out_dir, "yields_{}.txt".format(postfix)), 'w') as f:
            f.write(yield_out_str)

        f_out.Close()
        return outname
