import ROOT
import glob
import os
import pandas as pd
import flashgg_fit
import attridict
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-hist", action="store_true")
parser.add_argument("-fitmass", action="store_true")
parser.add_argument("-fitdeltam", action="store_true")

args = parser.parse_args()

nointfile = "gluglu_2018_hgg_noint_fullsim_smear.root"

gammas = ("0_1", "1_0", "1_5", "2_0", "2_5", "3_0", "3_5", "4_0", "4_5", "5_0", "5_5", "6_0", "6_5", "7_0", "7_5", "8_0", "8_5", "9_0", "9_5", "10_0", "20_0", "30_0", "40_0", "50_0")

dfbins = pd.read_csv("higgsptbins.csv")
nbins = len(dfbins)-1

if args.hist: os.system(f"./norm_and_dump_histos_couplingscale.sh {nointfile} noint 1 4.379")

if args.fitmass: os.system("echo gamma,ptbin,mass_w,mass_h,mass_h_err,chi2_over_ndf > mean_nointqg.csv")

for g in gammas:

  g = g.replace("_", ".")

  dfint = pd.read_csv("norm_int.csv")
  normint = dfint[dfint.gamma==float(g)]["norm"].iloc[0]

  #dfintqg = pd.read_csv("norm_intqg.csv")
  #normintqg = dfint[dfint.gamma==float(g.replace("_", "."))]["norm"].iloc[0]

  if args.hist:
    intfile = glob.glob(f"int/*{g.replace('.', '_')}*smear*")[0]
    #intqgfile = glob.glob(f"intqg/*{g.replace('.', '_')}*smear*")[0]

    os.system(f"./norm_and_dump_histos_couplingscale.sh {intfile} int {g} {normint}")
    #os.system(f"./norm_and_dump_histos_couplingscale.sh {intqgfile} intqg {g} {normintqg}")

  else: g = g.replace("_", ".")

  for i in range(nbins):
    if args.hist:
      os.system(f"hadd -f histos/mass_couplingscale_nointqg_{g}_hptbin{i}.root histos/mass_couplingscale_int_{g}_hptbin{i}.root histos/mass_couplingscale_noint_1_hptbin{i}.root")
      #os.system(f"hadd -f histos/mass_couplingscale_nointqg_{g}_hptbin{i}.roothistos/mass_couplingscale_intqg_{g}_hptbin{i}.root  histos/mass_couplingscale_int_{g}_hptbin{i}.root histos/mass_couplingscale_noint_1_hptbin{i}.root")

    if args.fitmass:
      options = attridict({
        "inputFile": f"histos/mass_couplingscale_nointqg_{g}_hptbin{i}.root",
        "proc": f"ggH_plus_intgg_gammamod_{g}", "cat": f"hptbin{i}",
        "year": '2018', "nGauss": 3, "outputDir": "histos"
      })

      mass_w, mass_h, mass_h_err, chi2_over_ndf = flashgg_fit.run(options)

      os.system(f"echo {g},{i},{mass_w:.3f},{mass_h:.3f},{mass_h_err:.4f},{chi2_over_ndf:.3f} >> mean_nointqg.csv")

if args.fitdeltam:
  os.system("echo ptbin,p0,p0err,p1,p1err > sqrt_fit_dm_gamma.csv")

  for i in range(nbins):
    os.system(f"./plot_means_hmeanstat_vs_hmeananal.sh mean_nointqg.csv {i} dmfit/media_fittata_ptbin_{i} >> sqrt_fit_dm_gamma.csv")
    #os.system(f"./plot_means_noerr.sh mean_nointqg.csv {i} media_fittata_incl_{i} >> sqrt_fit_dm_gamma.csv")
