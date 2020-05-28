"""
Read Workspaces and make histograms/yields
"""
from __future__ import print_function
import ROOT
from ROOT import RooFit

ROOT.gErrorIgnoreLevel = ROOT.kInfo
ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.ERROR)

import os
import array
from math import fabs, pow, sqrt
from root_utils import integral_and_error

thispath=os.path.abspath(__file__)
thispath=os.path.dirname(thispath)
ROOT.gInterpreter.Declare('#include "'+"{}/RooExpandedFitResult.h".format(thispath)+'"')
ROOT.gROOT.LoadMacro("{}/RooExpandedFitResult.C".format(thispath))

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

        self.fit_res=None # fit result
        # event yield
        self.dict_yield={}

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

    def fix_shapeNPs(self, ):
        """ set all shape sys NPs as constants """
        nuisances = self.mc.GetNuisanceParameters()
        itr = ROOT.TIter(nuisances.createIterator())
        var = itr()
        while var:
            name = var.GetName()
            if "alpha" in name and "Shape" in name:
              var.setConstant(True)
              print("WARNING! set shape NP as constant", name)
              var.Print()
            var = itr()

    def fix_allNPs(self, ):
        """ set all sys NPs as constants """
        nuisances = self.mc.GetNuisanceParameters()
        itr = ROOT.TIter(nuisances.createIterator())
        var = itr()
        while var:
            name = var.GetName()
            if "alpha" in name:
              var.setConstant(True)
              print("WARNING! set NP as constant", name)
              var.Print()
            var = itr()

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

    def translate_binning(self, hist, cate='ggF'):
        """
        for 2l2v mainly, where variable binning used
        """
        #TODO: hardcoded binning
        dbin_ggF_llvv = [0., 50., 100., 150., 200., 250., 300., 350., 400., 450., 500., 550., 600., 650., 700., 750., 800., 850., 900., 950., 1000., 1100., 1200., 1300., 1400., 1600., 1800., 2000.]
        dbin_VBF_llvv = [0., 100., 420., 540., 820., 1700.]
       
        if 'ggF' in cate:
            dbin = dbin_ggF_llvv
        elif 'VBF' in cate:
            dbin = dbin_VBF_llvv

        if not hist: return hist
        name = hist.GetName()
        name = name + "_newbinning"
         
        if isinstance(hist, ROOT.TH1): 
          if not hist or hist.GetNbinsX() > 28: return hist
          new_hist = ROOT.TH1F(name, name, ( len(dbin)-1), array.array('d',dbin))
          for i in range(len(dbin)-1):
              #print("translate_binning:", hist.GetName(), i, hist.GetBinContent(i+1))
              new_hist.SetBinContent(i+1, hist.GetBinContent(i+1))
              new_hist.SetBinError(i+1, hist.GetBinError(i+1))
              new_hist.SetXTitle("mT")
         
          return new_hist

        ## Tgraph
        if isinstance(hist, ROOT.TGraphAsymmErrors):
          new_hist = hist.Clone()
          new_hist.SetName(name)
          new_hist.SetTitle(name)
          for i in range(hist.GetN()):
            ibin = i
            x, y = ROOT.Double_t(), ROOT.Double_t()
            hist.GetPoint(ibin, x, y)
            ex_up, ex_dn, ey_up, dy_dn = hist.GetErrorXhigh(ibin), hist.GetErrorXlow(ibin), hist.GetErrorYhigh(ibin), hist.GetErrorYlow(ibin)
            # do the translation
            n_x, n_ex_dn, n_ex_up, = x, ex_dn, ex_up
            if x >0 and x < len(dbin)-1:
              n_x = (dbin[i-1] + dbin[i])*0.5
              n_ex_dn = n_x - dbin[i-1]
              n_ex_up = dbin[i] - n_x
            else:
              if x<0: # underflow
                n_x = dbin[0] - 1
                n_ex_dn, n_ex_up=1, 1
              else: # overlfow
                n_x = dbin[-1] + 1
                n_ex_dn, n_ex_up=1, 1
   
            r_up = 0
            r_down = 0
            content=y
            if content!=0:
              r_up = e_up/content
              r_down = e_down/content
            new_hist.SetPoint(i, n_x, 1.)
            new_hist.SetPointEXhigh(ibin, n_ex_up)
            new_hist.SetPointEXlow(ibin, n_ex_dn)
            new_hist.SetPointEYhigh(ibin, ey_up)
            new_hist.SetPointEYlow(ibin, ey_dn)

        return new_hist

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

    @staticmethod
    def create_err_from_pdf(pdf, fitres, hist, hist_name, obs, doIntegralOnly=1, nBins=0):
        """
        get the error band. the uncertainties are propagted via a RooFitResult
        """
        print("get error for {}".format(hist_name))
        if not hist:
            hist = self.create_hist_from_pdf(pdf, hist_name, obs)
 
        [h_errors, sampleYield, sampleError]=[None, 0, 0]

        obs.Print()
        obs_binning=obs.getBinning()
        #ba=obs_binning.array()
        #for i in range(obs.getBins()): print(ba[i])
        hist.Print("All")

        #get the error band
        if not doIntegralOnly:
          frame=obs.frame()
          pdf.plotOn(frame, ROOT.RooFit.VisualizeError(fitres,1), ROOT.RooFit.FillColor(ROOT.kBlack), ROOT.RooFit.LineColor(ROOT.kBlack), ROOT.RooFit.Normalization(hist.Integral(), ROOT.RooAbsReal.NumEvent))
          h_errors_tmp=frame.getCurve()
          if not h_errors_tmp:
            print("Error=> could not get the postfit error for {}, {}\n".format(hist_name, pdf.GetName()))
          else:
            #find all the high point in tgraph
            yuptmp,ydowntmp=[], []
            upYield=0
            downYield=0
            ne=h_errors_tmp.GetN()
            nescan=ne/2
            if ne>=4*(obs.getBins()): ## TODO: why is it like this??
              step=2 
            else:  step=1
            print("check ne={}, step={}".format(ne, step))
            
            #for ib in range(1, h_errors_tmp.GetN()/2, 2):
            for ib in range(1, nescan, step):
              x1=ROOT.Double(0)
              y1=ROOT.Double(0)
              h_errors_tmp.GetPoint(ib,x1,y1)
              yuptmp.append(y1)
              upYield+=y1
        
            #find all low point in tgraph
            #for ib in range(h_errors_tmp.GetN()-2, h_errors_tmp.GetN()/2, -2):
            for ib in range(h_errors_tmp.GetN()-2, nescan, -step):
              x1=ROOT.Double(0)
              y1=ROOT.Double(0)
              h_errors_tmp.GetPoint(ib,x1,y1)
              ydowntmp.append(y1)
              downYield+=y1
        
            #finally make arrays storing required values
            xs,ys,yup,ydown,xerr=[], [], [], [], []
            for ib in range(1, hist.GetNbinsX()):
              x=hist.GetBinCenter(ib)
              y=hist.GetBinContent(ib)
        
              xs.append(x)
              ys.append(y)
              yup.append(yuptmp[ib]-y)
              ydown.append(y-ydowntmp[ib])
              xerr.append(hist.GetBinCenter(ib)-hist.GetBinLowEdge(ib))
            
            #Now make the graph
            axs=array.array('f', xs)
            ays=array.array('f', ys)
            ayup=array.array('f', yup)
            aydown=array.array('f', ydown)
            axerr=array.array('f', xerr)
        
            h_errors=ROOT.TGraphAsymmErrors(len(xs), axs, ays, axerr, axerr, aydown, ayup)
            h_errors.SetName("{}_err".format(hist_name))
            h_errors.SetTitle("{}_err".format(hist_name))
            h_errors.SetLineColor(ROOT.kBlack)
            h_errors.SetFillColor(ROOT.kBlack)
            h_errors.SetMarkerStyle(0)
            h_errors.SetLineStyle(1)
            h_errors.SetFillStyle(3004)
        
          postFitYield = pdf.expectedEvents(ROOT.RooArgSet(obs))
          print("\texpectedEvents={:.4f}\n".format(postFitYield))
          if postFitYield==0:
            postFitYield=hist.Integral()
          print("\tPostfit Integral={:.4f} + {:.4f} - {:.4f}\n".format(postFitYield,upYield-postFitYield,postFitYield-downYield))

          sampleYield=postFitYield
          sampleError=(fabs( upYield-postFitYield ) + fabs( postFitYield-downYield) )*0.5

        else:
          h_errors=None

          #Find the error on the yield
          if nBins>0:
            sampleIntegral= pdf.createIntegral(ROOT.RooArgSet(obs))
            sampleYield=sampleIntegral.getVal()*nBins
            sampleError=sampleIntegral.getPropagatedError(fitres)*nBins

            ##do scaling
            nevents=hist.Integral()
            sf=sampleYield/nevents
            sampleYield=nevents
            sampleError*=sf

            print("\tPostfit Integral={:.4f} +/- {:.4f} (nBins= {:.4f}, sf= {:.4f})\n".format(sampleYield, sampleError, nBins, sf))

        return [h_errors, sampleYield, sampleError]


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

        if os.path.isfile(self.corr_root_path):
            print("Reuse fitted results saved at root file: {}".format(self.corr_root_path))
            fout = ROOT.TFile.Open(self.corr_root_path)
            self.corr = fout.Get("corr")
            self.corr.SetDirectory(0)
            self.fit_res = fout.Get("nll_res")
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
        #minim.setStrategy(2)
        minim.setStrategy(1)
        #minim.setProfile()
        status = minim.minimize(
            "Minuit2",
            ROOT.Math.MinimizerOptions.DefaultMinimizerAlgo()
        )
        if status != 0:
            minim.setStrategy(1)
            status = minim.minimize(
                "Minuit2",
                ROOT.Math.MinimizerOptions.DefaultMinimizerAlgo()
            )
        assert status == 0, "Fit did not converge..change to a different strategy"

        self.ws.saveSnapshot(snapshot_name, self.ws.allVars())
        if save_to_file:
            print("debug", self.out_dir, save_to_file)
            #self.ws.writeToFile(os.path.join(self.out_dir, save_to_file))

        # after performed the Fit, plot the coefficient matrix
        if do_matrix:
            self.fit_res = minim.save()
            corr_hist = self.fit_res.correlationHist()
            self.corr = corr_hist.Clone("corr")
            self.corr.SetDirectory(0)
            fout = ROOT.TFile.Open(self.corr_root_path, 'recreate')
            self.corr.Write()
            self.fit_res.SetName("nll_res")
            self.fit_res.Write()
            fout.Close()

    @staticmethod
    def resetError(m_w, mc, parList, vetoList=None):
      """
      taken from CommonStatTools/HistFactoryInspector.C
      RooWorkspace* m_w, RooStats::ModelConfig* mc, const RooArgList &parList, const RooArgList &vetoList= RooArgList()
    
      For the given workspace,
      find the input systematic with
      the given name and shift that
      systematic by 1-sigma
      """
    
      new_parList=ROOT.RooArgList(parList)
      it = ROOT.TIter(new_parList.createIterator())
    
      globs = ROOT.RooArgSet(mc.GetGlobalObservables())
      tmpvar=ROOT.RooRealVar()
    
      it = it.Begin()
      arg = it()
      while arg:
        #arg.Print()
        UncertaintyName=""
        # Ignore constant parameters
        if (arg.InheritsFrom("RooRealVar") and (not arg.isConstant())):
          UncertaintyName = arg.GetName()
        else:
          arg = it()
          continue
    
        # Ignore parameters without constraint (Floattig normalization factors, for eg)
        hasConstraint=False
        itr = ROOT.TIter(globs.createIterator())
        itr = itr.Begin()
        tmpvar = itr()
        while tmpvar:
          varName=tmpvar.GetName()
          #tmpvar.Print()
          if UncertaintyName.replace("alpha_", "") in varName:
            if "_In" in varName or "nom_" in varName:
              hasConstraint=True
              break
          tmpvar = itr()
    
        if (not hasConstraint): 
          arg = it()
          print("Warning=> not constraint, ignore", UncertaintyName)
          continue
    
        if vetoList and (vetoList.FindObject(UncertaintyName) != 0):
          arg = it()
          continue
    
        var = m_w.var(UncertaintyName)
    
        # Initialize
        val_hi  = 5
        val_low = -5
        sigma   = 0.
        resetRange=False
    
        # If it is a stat uncertainty (gamma)
        if "gamma" in UncertaintyName:
    
          # Get the constraint and check its type:
          constraint     = m_w.obj((UncertaintyName + "_constraint"))
          ConstraintType = "RooPoisson"
          if constraint:
            ConstraintType = constraint.IsA().GetName()
    
          if (ConstraintType == "RooGaussian"):
            sigmaVar = m_w.obj((UncertaintyName + "_sigma"))
            sigma    = sigmaVar.getVal()
    
            # Symmetrize shifts
            val_hi     = 1 + sigma
            val_low    = 1 - sigma
            resetRange = True
          elif (ConstraintType == "RooPoisson"):
            nom_gamma     = m_w.obj(("nom_" + UncertaintyName))
            nom_gamma_val = nom_gamma.getVal()
    
            sigma      = 1 / ROOT.TMath.Sqrt(nom_gamma_val)
            val_hi     = 1 + sigma
            val_low    = 1 - sigma
            resetRange = True
    
        # Assume it's standard (gaussian) uncertainty by default
        else:
          # Assume the values are +1, -1
          val_hi     = 1.0
          val_low    = -1.0
          sigma      = 1.0
    
        var.setError(fabs(sigma))
        #var.Print()
        if (resetRange):
          minrange = var.getMin()
          maxrange = var.getMax()
          newmin   = var.getVal() - 6. * sigma
          newmax   = var.getVal() + 6. * sigma
          if (minrange < newmin): var.setMin(newmin)
          if (newmax < maxrange): var.setMax(newmax)
    
        arg = it()
    
      new_parList.Print()
      return new_parList  
    
    def create_prefitres(self,):

      nuis = self.mc.GetNuisanceParameters()
      #nuis.add(self.mc.GetParametersOfInterest()) ## doesn't matter for prefit anyway
    
      np=ROOT.RooArgList(nuis)
    
      np=self.resetError(self.ws, self.mc, np)
      rfrexp = ROOT.RooExpandedFitResult(np)
      rfrexp.Print()
      self.fit_res=rfrexp

    def loop_categories(self, postfix='nominal', force_update=False, use_fitres=True):
        """
        make histograms for each component, save them into file [optional]
        print and save the yields
        """
        outname = os.path.join(self.out_dir, "histograms_{}.root".format(postfix))
        if os.path.exists(outname):
            f_out = ROOT.TFile.Open(outname, "READ")
            lkeys=f_out.GetListOfKeys()
            if len(lkeys)>0 and (not force_update):
                print("{} is there, ".format(outname))
                f_out.Close()
                return outname
            else:
                print("forcing updates. May the observable range be changed?")

        f_out = ROOT.TFile.Open(outname, "recreate")

        iter_category = ROOT.TIter(self.categories.typeIterator())
        obj = iter_category()

        if use_fitres:
          if not self.fit_res:
            if os.path.isfile(self.corr_root_path):
                print("Reuse fitted results saved at root file: {}".format(self.corr_root_path))
                fout_rfr = ROOT.TFile.Open(self.corr_root_path)
                self.fit_res = fout_rfr.Get("nll_res")
          if self.fit_res:
            print("Load fit results")
            self.fit_res.Print()
            if not self.ws.loadSnapshot("vars_final"):
              np = self.mc.GetNuisanceParameters()
              np.add(self.mc.GetParametersOfInterest())
              snp_init="nuisance_norminal"
              self.ws.loadSnapshot(snp_init)
              print("Prefit Values:\n")
              np.Print("v")
              #Get Paramters from fit
              fpf = self.fit_res.floatParsFinal().Clone()  #floating param from fit result
              cp = self.fit_res.constPars() #constant param from fit result
              fpf.add(cp) # add all const parameters of the RooFitResult to the floating ones
          
              #assign np to fitted values
              np.assignValueOnly(fpf) # only sets the central value. Should be ok as VisualizeError uses errors from the RooFitResult, according to ROOT doc
              print("Postfit Values:\n")
              np.Print("v")
              self.ws.saveSnapshot("vars_final", np)
            self.ws.loadSnapshot("vars_final")

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
                if cat_name not in self.dict_yield: self.dict_yield[cat_name]={}
                self.dict_yield[cat_name]["data"]=[hist_data.Integral(), 0]

            # get signal + background
            print("get signal + background")
            hist_splusb = ROOT.TH1F("splusb-"+cat_name+'-'+postfix,
                                    "signal + background", *obs_bins)
            hist_splusb = self.create_hist_from_pdf(pdf, hist_splusb.GetName(), obs_var)
            hist_splusb.SetLineColor(206)
            all_histograms.append(hist_splusb)
            hist_splusb.Print()
            if cat_name not in self.dict_yield: self.dict_yield[cat_name]={}
            ## get the error
            if self.fit_res:
              [hist_splusb_err, y, err] = self.create_err_from_pdf(pdf, self.fit_res, hist_splusb, hist_splusb.GetName(), obs_var, 0, 0)
              if hist_splusb_err: all_histograms.append(hist_splusb_err)
              self.dict_yield[cat_name]["splusb"]=[y, err]
            else:
              (y, err) = integral_and_error(hist_splusb)
              self.dict_yield[cat_name]["splusb"]=[y, err]

            # get background-only events
            print("get background-only events")
            old_poi_val = self.poi.getVal()
            self.poi.setVal(1e-10)
            hist_bonly = ROOT.TH1F("bonly-"+cat_name+"-"+postfix, "background only ", *obs_bins)
            hist_bonly = self.create_hist_from_pdf(pdf, hist_bonly.GetName(), obs_var)
            bonly_evts = hist_bonly.Integral()
            all_histograms.append(hist_bonly)
            hist_bonly.Print()
            ## get the error
            if self.fit_res:
              [hist_bonly_err, y, err] = self.create_err_from_pdf(pdf, self.fit_res, hist_bonly, hist_bonly.GetName(), obs_var, 0, 0)
              if hist_bonly_err: all_histograms.append(hist_bonly_err)
              self.dict_yield[cat_name]["bonly"]=[y, err]
            else:
              (y, err) = integral_and_error(hist_bonly)
              self.dict_yield[cat_name]["splusb"]=[y, err]

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
                print("Look at RooProdPdf pdf", pdf.GetName())
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
                        total_err = 0

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
                            y, err = 0, 0
                            total += sum_ch
                            simple_name = func.GetName().split('_')[2]
                            print("  Look at component {}, yield {:.3f}".format(func.GetName(), sum_ch))
                            if sum_ch > 1E-5 and "signal" not in func.GetName().lower():
                                hist_name = simple_name+"-"+cat_name+"-"+postfix
                                histogram = self.get_hist(func, obs_var, sum_ch, hist_name)
                                all_histograms.append(histogram)
                                ## get the error
                                if self.fit_res:
                                  if is_sum_pdf:
                                    [hist_err, y, err] = self.create_err_from_pdf(func, self.fit_res, histogram, histogram.GetName(), obs_var, 1, coeff_list[0].getVal())
                                  else:
                                    ### need to consider the syst eff for both the normalization and shape
                                    pname = "RooAddpdf_"+simple_name+"_"+cat_name
                                    #samp_pdf = ROOT.RooAddPdf(pname, pname, ROOT.RooArgList(func), ROOT.RooArgList(coeff_list[func_index])) ## much slower
                                    samp_pdf = func
                                    [hist_err, y, err] = self.create_err_from_pdf(samp_pdf, self.fit_res, histogram, histogram.GetName(), obs_var, 1, 1)
                     
                                  if hist_err: all_histograms.append(hist_err)
                                  self.dict_yield[cat_name][simple_name]=[y, err]
                                else:
                                  (y, err) = integral_and_error(histogram)
                                  self.dict_yield[cat_name][simple_name]=[y, err]


                            ## signal pdf
                            if sum_ch > 1E-5 and "signal" in func.GetName().lower():
                                hist_name = "Signal-"+cat_name+"-"+postfix
                                hist_sonlypdf=self.get_hist(func, obs_var, sum_ch, hist_name)
                                all_histograms.append(hist_sonlypdf)
                                ## get the error
                                if self.fit_res:
                                  [hist_err, y, err] = self.create_err_from_pdf(func, self.fit_res, hist_sonlypdf, hist_sonlypdf.GetName(), obs_var, 1, 1)
                                  if hist_err: all_hist_sonlypdfs.append(hist_err)
                                  if cat_name not in self.dict_yield: self.dict_yield[cat_name]={}
                                  self.dict_yield[cat_name]["signal"]=[y, err]
                                else:
                                  (y, err) = integral_and_error(hist_sonlypdf)
                                  self.dict_yield[cat_name]["signal"]=[y, err]

                            yield_out_str += "{} {:.2f} {:.2f} +/- {:.3f}\n".format(func.GetName(), sum_ch, y, err)
                            total_err = sqrt( pow(total_err, 2) + pow(err, 2) )
                        yield_out_str += "total yields {:.2f} +/- {:.3f}\n".format(total, total_err)
                    else:
                        print("no baseline pdf avaiable!")
                        continue

            elif "RooAddPdf" in pdf.ClassName():
                print("Look at RooAddPdf pdf", pdf.GetName())
                func_list = pdf.pdfList()
                coeff_list = pdf.coefList()
                total = 0
                total_err = 0
                nevts_func = lambda k:coeff_list[k].getVal()
                for func_index in sorted(
                    range(func_list.getSize()), key=nevts_func,
                    reverse=True
                ):
                    func = func_list[func_index]
                    sum_ch = nevts_func(func_index)
                    total += sum_ch
                    simple_name = func.GetName().split('_')[2]
                    print("  Look at component", func.GetName())
                    if sum_ch > 1E-5 and "signal" not in func.GetName().lower():
                        hist_name = simple_name+"-"+cat_name+"-"+postfix
                        histogram = self.get_hist(func, obs_var, sum_ch, hist_name)
                        all_histograms.append(histogram)
                        ## get the error
                        if self.fit_res:
                          [hist_err, y, err] = self.create_err_from_pdf(func, self.fit_res, histogram, histogram.GetName(), obs_var, 1, 1)
                          if hist_err: all_histograms.append(hist_err)
                          self.dict_yield[cat_name][simple_name]=[y, err]
                        else:
                          (y, err) = integral_and_error(histogram)
                          self.dict_yield[cat_name][simple_name]=[y, err]

                    yield_out_str += "{} {:.2f} {:.2f} +/- {:.3f}\n".format(func.GetName(), sum_ch, y, err)
                    total_err = sqrt( pow(total_err, 2) + pow(err, 2) )
                yield_out_str += "total yields {:.2f} +/- {:.3f}\n".format(total, total_err)
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
