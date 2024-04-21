
import ROOT
import numpy as np
import argparse
import uproot
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("filename", help="filename (.root)")
args = parser.parse_args()

inc_file = ROOT.TFile("reso_eta_inclusive_2018.root")
ROOT.gROOT.cd()
inc = inc_file.Get("h_inc").Clone()
inc_file.Close()

lowbrem_file = ROOT.TFile("reso_eta_lowbrem_2018.root")
ROOT.gROOT.cd()
lowbrem = lowbrem_file.Get("h_lowbrem").Clone()
lowbrem_file.Close()

unconv_frac_file = ROOT.TFile("unconv_frac.root")
ROOT.gROOT.cd()
unconv_frac = unconv_frac_file.Get("golden1").Clone()
unconv_frac_file.Close()

f = ROOT.TFile(args.filename)
tree = f.Get("Events")

h = ROOT.TH1F("sumw", "sumw", 3, 0, 3)

tree.Draw("1>>sumw", "genWeight", "goff")

if "_intqg_" in args.filename:
  gamma = ".".join(args.filename.split("intqg_")[1].split("_")[:2])
  df = pd.read_csv("total_XS_intqg_plain100k.csv")
  xsec = df[df.gamma==float(gamma)]["xsec"].iloc[0]
  open("norm_intqg.csv", "a").write(f"{gamma},{xsec*1e3*1470/h.Integral()}\n")
if "_int_" in args.filename:
  gamma = ".".join(args.filename.split("int_")[1].split("_")[:2])
  df = pd.read_csv("total_XS_int_plain100k.csv")
  xsec = df[df.gamma==float(gamma)]["xsec"].iloc[0]
  open("norm_int.csv", "a").write(f"{gamma},{xsec*1e3*1470/h.Integral()}\n")

vartype_list = ["pt", "eta", "phi", "en_gen", "en", "has_conv"]
part_list = ["lead", "sublead"]
var_list = [f"{part}_{vartype}" for vartype in vartype_list for part in part_list]

v = {var: [] for var in var_list}
v.update({"mass": [], "weight": [], "higgs_pt": []})

j=0
for entry in tree:
  if j%1000==0: print(j)
  j+=1
  if entry.nGenIsolatedPhoton != 2: continue
  if np.abs(entry.GenIsolatedPhoton_eta[0]) > 2.5 or np.abs(entry.GenIsolatedPhoton_eta[1]) > 2.5:
    continue
  if np.abs(entry.GenIsolatedPhoton_eta[0]) > 1.479 and np.abs(entry.GenIsolatedPhoton_eta[0]) < 1.57:
    continue
  if np.abs(entry.GenIsolatedPhoton_eta[1]) > 1.479 and np.abs(entry.GenIsolatedPhoton_eta[1]) < 1.57:
    continue
  vects = ROOT.TLorentzVector(), ROOT.TLorentzVector()
  for i, part in enumerate(part_list):
    v[f"{part}_eta"].append(entry.GenIsolatedPhoton_eta[i])
    v[f"{part}_phi"].append(entry.GenIsolatedPhoton_phi[i])
    pt, has_conv = entry.GenIsolatedPhoton_pt[i], 0
    v[f"{part}_pt"].append(pt)
    v[f"{part}_en"].append(pt/np.arctan(2*np.exp(-entry.GenIsolatedPhoton_eta[i])))
    v[f"{part}_en_gen"].append(entry.GenIsolatedPhoton_pt[i]/np.arctan(2*np.exp(-entry.GenIsolatedPhoton_eta[i])))
    v[f"{part}_has_conv"].append(has_conv)
    vects[i].SetPtEtaPhiM(pt, entry.GenIsolatedPhoton_eta[i], entry.GenIsolatedPhoton_phi[i], 0)
  v[f"higgs_pt"].append((vects[0]+vects[1]).Pt())
  v["mass"].append((vects[0]+vects[1]).M() + ROOT.gRandom.Gaus(0, 1.7))
  v["weight"].append(entry.genWeight)

f.Close()

outf = uproot.recreate(f"{args.filename.split('.root')[0]}_smear.root")
outf["tree"] = v
outf.close()
