import ROOT
import os

script_dir = os.path.dirname(os.path.abspath(__file__))
ROOT.gROOT.LoadMacro(script_dir+"/AtlasStyle.C")
ROOT.SetAtlasStyle()
