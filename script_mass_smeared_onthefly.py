import glob
import ROOT
import os
import sys

files = glob.glob("int/*_100k*_smear.root")

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
  c = ROOT.TChain("tree")
  c.Add(file)
  h = ROOT.TH1F(f"h_{g_rel:.1f}", f"#Gamma/#Gamma_{{SM}} = {g_rel:.1f}", 70, 105, 140)
  c.Draw(f"mass>>h_{g_rel:.1f}", "weight", "E")
  h.GetXaxis().SetTitle("Diphoton mass [GeV]")
  h.GetYaxis().SetTitle("Entries / 0.5 GeV")
  cnv.SaveAs(f"cnv_smear_{g_rel:.1f}.root")
  cnv.SaveAs(f"cnv_smear_{g_rel:.1f}.pdf")
  del c

os.system("pdfunite $(ls -1 cnv_smear*.pdf | sort --version-sort) all_gammas_smear.pdf")
