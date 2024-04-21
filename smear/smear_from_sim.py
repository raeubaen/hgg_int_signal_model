import ROOT
import numpy as np
import argparse
import uproot
import pandas as pd

def smear_pt(pt, eta, inc, lowbrem, unconv_frac):
  theta = 2*np.arctan(np.exp(-eta))
  en = pt/np.sin(theta)
  inc_relerr = inc.GetBinContent(inc.FindBin(np.abs(eta)))
  lowbrem_relerr = lowbrem.GetBinContent(lowbrem.FindBin(np.abs(eta)))
  unconv_probability = unconv_frac.GetBinContent(unconv_frac.FindBin(eta))
  has_conv = (np.random.uniform() > unconv_probability)
  if has_conv: en += np.random.normal(0, en*inc_relerr)
  else: en += np.random.normal(0, en*lowbrem_relerr)
  return (124.7/125*en*np.sin(theta), has_conv)

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

vartype_list = ["fullpt", "pt", "pt_gen", "eta", "phi", "en_gen", "en", "has_conv"]
part_list = ["lead", "sublead"]
var_list = [f"{part}_{vartype}" for vartype in vartype_list for part in part_list]

v = {var: [] for var in var_list}
v.update({"mass": [], "weight": [], "higgs_pt": []})

j=0
for entry in tree:
  if j%1000==0: print(j)
  j+=1
  if entry.nPhoton != 2:
    continue
  if entry.Photon_genPartIdx[0]==-1 or entry.Photon_genPartIdx[1]==-1:
    continue
  if np.abs(entry.GenPart_eta[entry.Photon_genPartIdx[0]]) > 2.5 or np.abs(entry.GenPart_eta[entry.Photon_genPartIdx[1]]) > 2.5:
    continue
  if np.abs(entry.GenPart_eta[entry.Photon_genPartIdx[0]]) > 1.479 and np.abs(entry.GenPart_eta[entry.Photon_genPartIdx[0]]) < 1.57:
    continue
  if np.abs(entry.GenPart_eta[entry.Photon_genPartIdx[1]]) > 1.479 and np.abs(entry.GenPart_eta[entry.Photon_genPartIdx[1]]) < 1.57:
    continue
  vects = ROOT.TLorentzVector(), ROOT.TLorentzVector()
  for i, part in enumerate(part_list):
    ptgen = entry.GenPart_pt[entry.Photon_genPartIdx[i]]
    etagen = entry.GenPart_eta[entry.Photon_genPartIdx[i]]
    v[f"{part}_pt_gen"].append(ptgen)
    v[f"{part}_fullpt"].append(entry.Photon_pt[i])
    v[f"{part}_eta"].append(etagen)
    v[f"{part}_phi"].append(entry.GenPart_phi[entry.Photon_genPartIdx[i]])
    pt, has_conv = smear_pt(ptgen, etagen, inc, lowbrem, unconv_frac)
    v[f"{part}_pt"].append(pt)
    v[f"{part}_en"].append(pt/np.arctan(2*np.exp(-etagen)))
    v[f"{part}_en_gen"].append(ptgen/np.arctan(2*np.exp(-etagen)))
    v[f"{part}_has_conv"].append(has_conv)
    vects[i].SetPtEtaPhiM(pt, etagen, entry.GenPart_phi[entry.Photon_genPartIdx[i]], 0)
  v[f"higgs_pt"].append((vects[0]+vects[1]).Pt())
  v["mass"].append((vects[0]+vects[1]).M())
  v["weight"].append(entry.genWeight)

f.Close()

outf = uproot.recreate(f"{args.filename.split('.root')[0]}_smear.root")
outf["tree"] = v
outf.close()
