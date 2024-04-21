import glob
import ROOT
import os
import sys

files = glob.glob("../../eos/sherpa_int/*_100k*.root")

to_discard = []

for td in to_discard:
  files.remove(f"../../eos/sherpa_int/sherpa_int_{td}_NANOAODGEN_100k.root")

ROOT.gROOT.LoadMacro("root_logon.C")

for file in files:
  print(file)
  try:
    g_rel = float("_".join(file.split("/")[-1].split("_")[2:4]).replace("_", "."))
  except IndexError:
    print("skipping - not OK")
    continue
  print(g_rel)
  cnv = ROOT.TCanvas(f"c_{g_rel:.1f}")
  c = ROOT.TChain("Events")
  c.Add(file)
  h = ROOT.TH1F(f"h_{g_rel:.1f}", f"#Gamma/#Gamma_{{SM}} = {g_rel:.1f}", 200, 124, 126)
  c.Draw(f"sqrt(2*GenIsolatedPhoton_pt[0]*GenIsolatedPhoton_pt[1]*(cosh(GenIsolatedPhoton_eta[0]-GenIsolatedPhoton_eta[1]) - cos(GenIsolatedPhoton_phi[0]-GenIsolatedPhoton_phi[1])))>>h_{g_rel:.1f}", "(nGenIsolatedPhoton==2) * genWeight")
  h.GetXaxis().SetTitle("Diphoton mass [GeV]")
  h.GetYaxis().SetTitle("Entries / 10 MeV")
  cnv.SaveAs(f"cnv_{g_rel:.1f}.root")
  cnv.SaveAs(f"cnv_{g_rel:.1f}.pdf")
  del c

os.system("pdfunite $(ls -1 cnv_*.pdf | sort --version-sort) all_gammas.pdf")
