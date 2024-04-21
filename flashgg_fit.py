import ROOT
import glob
from flashgg_fit_tool import *
from flashgg_plottingTools import *
from collections import OrderedDict as od
from optparse import OptionParser

def run(opt):
  ROOT.gROOT.SetBatch(True)
  ROOT.gStyle.SetOptStat(0)


  # Open ROOT file storing datasets
  f = opt.inputFile
  fin = ROOT.TFile(f)

  pm = "mass"

  MHLow = '120'
  MHHigh = '130'
  massPoints = '125'
  nBins = 20 #nBins for fit
  MHPolyOrder = 0 # dependence of fit params on MH, set to 0 if using one mass point
  minimizerMethod = 'TNC'
  minimizerTolerance = 1e-8
  nGauss = opt.nGauss
  useDCB = False

  xvar = ROOT.RooRealVar("mass","mass", int(MHLow), int(MHHigh))

  # MH var
  MH = ROOT.RooRealVar("MH","m_{H}", int(MHLow), int(MHHigh))
  MH.setError(1e-4)

  # Create dict to store datasets: key=mass point, value = ROOT.TH1()
  datasets = od()
  datasets['125'] = fin.Get("%s_%s"%(pm, opt.cat))

  # Build ssf object + pdfs
  ssf = SimultaneousFit("name",opt.proc,opt.cat,datasets,xvar.Clone(),MH,MHLow,MHHigh,massPoints,nBins,MHPolyOrder,minimizerMethod,minimizerTolerance)
  ssf.buildNGaussians(nGauss)

  # Run fits and build mean + sigma splines
  chi2_over_ndf = ssf.runFit()

  print("FIT COMPLETED")

  #ssf.buildSplines()

  fracs = [0, 0, 0]
  means = [0, 0, 0]
  sigmas = [0, 0, 0]

  pdfs = od()

  '''
  ssf.Pdfs['final'].getComponents().Print()

  ssf.Pdfs['final'].getComponents().find("dm_g0").Print()
  print(ssf.Pdfs['final'].getComponents().find("dm_g0").getAsymErrorHi())
  print(ssf.Pdfs['final'].getComponents().find("mean_g0").getAsymErrorLo())

  ssf.Pdfs['final'].getComponents().find("frac_g0_constrained").Print()
  print(ssf.Pdfs['final'].getComponents().find("frac_g0_constrained").getAsymErrorHi())
  print(ssf.Pdfs['final'].getComponents().find("frac_g0_constrained").getAsymErrorLo())

  ssf.Pdfs['final'].getComponents().find("ggH_plus_intgg_gamma20_hptbin0_recursive_fraction_gaus_g2_3").Print()
  print(ssf.Pdfs['final'].getComponents().find("ggH_plus_intgg_gamma20_hptbin0_recursive_fraction_gaus_g2_3").getAsymErrorHi())
  print(ssf.Pdfs['final'].getComponents().find("ggH_plus_intgg_gamma20_hptbin0_recursive_fraction_gaus_g2_3").getAsymErrorLo())

  '''


  # Individual Gaussian histograms
  for k,v in ssf.Pdfs.items():
    if k == 'final': continue
    pdfs[k] = v
  if len(pdfs.keys())!=1:
    pdfItr = 0
    for k,v in pdfs.items():
      means[pdfItr] = ssf.Pdfs['final'].getComponents().getRealValue(f"mean_g{pdfItr}")
      sigmas[pdfItr] = ssf.Pdfs['final'].getComponents().getRealValue(f"sigma_g{pdfItr}")
      if pdfItr == 0:
        if "gaus" in k: fracs[pdfItr] = ssf.Pdfs['final'].getComponents().getRealValue("frac_g0_constrained")
        else: fracs[pdfItr] = ssf.Pdfs['final'].getComponents().getRealValue("frac_constrained")
      else:
        fracs[pdfItr] = ssf.Pdfs['final'].getComponents().getRealValue("%s_%s_recursive_fraction_%s_%i"%(ssf.proc,ssf.cat,k,pdfItr+1))
      pdfItr += 1

  #print(fracs, means, sum([f*m for f,m in zip(fracs, means)]))

  exp = ROOT.TMath.Exp
  sqrt = ROOT.TMath.Sqrt
  erf = ROOT.Math.erf

  unw_nums = lambda s, m: s/sqrt(2*ROOT.TMath.Pi()) * (exp(-(120-m)**2/(2*s**2)) - exp(-(130-m)**2/(2*s**2)) ) + m * (0.5*erf((130-m)/(sqrt(2)*s)) - 0.5*erf((120-m)/(sqrt(2)*s)))
  unw_dens = lambda s, m: 0.5*erf((130-m)/(sqrt(2)*s)) - 0.5*erf((120-m)/(sqrt(2)*s))

  mean_w = sum([f*unw_nums(s, m) for f, s, m in zip(fracs, sigmas, means)])/sum([f*unw_dens(s, m) for f, s, m in zip(fracs, sigmas, means)])

  h = ssf.Pdfs['final'].generateBinned(xvar, 1e6).createHistogram("mass")

  #print("\n\n\n\n\n\n\n\n\n", sum([f*m for f,m in zip(fracs, means)]), mean_w, h.GetMean(), "\n\n\n\n\n\n\n\n")

  plotPdfComponents(ssf,_outdir=opt.outputDir,_extension='total_',_proc=ssf.proc,_cat=ssf.cat)

  return (mean_w, h.GetMean(), h.GetMeanError(), chi2_over_ndf)

if __name__ == "__main__":
  parser = OptionParser()
  parser.add_option('--inputFile', dest='inputFile', default = "", help='Input flashgg workspace to fit')
  parser.add_option('--proc', dest='proc', default='ggH_plus_intgg', help="Name of signal process")
  parser.add_option('--cat', dest='cat', default='hptbin0', help="Name of analysis category")
  parser.add_option('--year', dest='year', default='2018', help="Year")
  parser.add_option('--nGauss', dest='nGauss', default=3, type='int', help="Number of gaussians")
  parser.add_option('--outputDir', dest='outputDir', default='.', help="Plot output directory")
  opt, args = parser.parse_args()
  print(run(opt))
