import json
from math import exp, sqrt
import numpy as np
import ROOT as r

r.gStyle.SetOptStat(0)
r.gROOT.ProcessLine(".L crystalBall.C")

data_points = []
signal_points = []
with open("unblinded_data_rho.json", "r") as fin:
    data_points = json.load(fin)
with open("signal_rho.json", "r") as fin:
    signal_points = json.load(fin)
N_sig = len(signal_points)*0.08330784*5e-3

# Background Exponential fit
bkg_fnc = r.TF1("bkg_fnc", r.Exp, 100, 180, 2)
bkg_fnc.SetParameter(0, len(data_points)*8/7)
bkg_fnc.SetParameter(1, 30)
# bkg_fnc = r.TF1("bkg_fnc", "[0]*(1-[1])/(pow(180, 1-[1])-pow(100, 1-[1]))*pow(x, -[1])", 100, 180)
# bkg_fnc.SetParameter(1, 0.1)
bkg_fnc.SetNpx(300)

# Signal Single Crystal Ball fit
sig_fnc = r.TF1("sig_fnc", r.CrystalBall, 100, 180, 5)
sig_fnc.SetParameter(0, 1)
sig_fnc.SetParameter(1, 125)
sig_fnc.SetParameter(2, 2)
sig_fnc.SetParameter(3, 2)
sig_fnc.SetParameter(4, 1)
fid = r.TFile("/nfs-7/userdata/bemarsh/rare-higgs/btree/v1/2018/wh.root")
t = fid.Get("Events")
h_sig =r.TH1D("h_sig","",160,100,180)
t.Draw("mcand_photon_mass[best_rhoCand_idx]>>h_sig", "fabs(mcand_mass_kaon[best_rhoCand_idx]-1.02)<0.02 && lepton_pt>35", "goff")
h_sig.Scale(1.0/ h_sig.GetMaximum())
h_sig.Fit("sig_fnc", "N", "goff", 100, 180)
sig_fnc.SetParameter(0, N_sig * sig_fnc.GetParameter(0) / sig_fnc.Integral(100,180))
sig_mean = sig_fnc.GetParameter(1)
sig_width = sig_fnc.GetParameter(2)
sig_alpha = sig_fnc.GetParameter(3)
sig_n = sig_fnc.GetParameter(4)
sig_fnc.SetNpx(300)

# Fill plot
h_blinded =r.TH1D("h_blinded","",80,100,180)
for d in data_points:
    h_blinded.Fill(d)

# Blind plot
# for i in range(21,31):
#     h_blinded.SetBinContent(i, 0)
#     h_blinded.SetBinError(i, 10000)

# Fit
r.REJECT=False
h_blinded.Fit(bkg_fnc, "NML", "goff", 100, 180)
r.REJECT=False

# Plot fit
gr_fit = r.TGraphErrors()
gr_err1 = r.TGraphErrors()
gr_err2 = r.TGraphErrors()
for i,x in enumerate(np.linspace(100, 180, 200)):
    gr_err1.SetPoint(i, x, bkg_fnc.Eval(x))
    gr_err2.SetPoint(i, x, bkg_fnc.Eval(x))
r.TVirtualFitter.GetFitter().GetConfidenceIntervals(gr_err1, 0.68)
r.TVirtualFitter.GetFitter().GetConfidenceIntervals(gr_err2, 0.95)

# Reset bin error
# for i in range(21,31):
#     h_blinded.SetBinError(i, 0)

# Canvas
c2 = r.TCanvas("c2","c2", 700, 600)

# Style
h_blinded.SetMarkerStyle(20)
h_blinded.SetMarkerColor(r.kBlack)
h_blinded.SetLineColor(r.kBlack)

bkg_fnc.SetLineWidth(2)
bkg_fnc.SetLineStyle(1)
bkg_fnc.SetLineColor(r.kRed)
gr_err1.SetFillColor(r.kGreen+1)
gr_err2.SetFillColor(r.kOrange-2)

bkg_fnc.SetLineColor(r.kBlack)
bkg_fnc.SetLineStyle(7)

h_blinded.Draw("PE")
gr_err2.Draw("3 SAME")
gr_err1.Draw("3 SAME")
sig_fnc.Draw("L SAME")
bkg_fnc.Draw("L SAME")
bkg_fnc.Draw("L SAME")
h_blinded.Draw("PE SAME")
c2.SaveAs("fit_rho.pdf")

raw_input()

# Roofit - expect seg fault
ws = r.RooWorkspace('w', "unblinded_workspace_rho")

M = r.RooRealVar("M","M",100,180)
data = r.RooDataSet("data_obs","data_obs",r.RooArgSet(M))

for i in data_points:
    M.setVal(i)
    data.add(r.RooArgSet(M))
data.Print()

x_scale = r.RooRealVar("x_scale","x_scale",-1,-0.001)
x_scale.setVal(-1.0 / bkg_fnc.GetParameter(1))
bkg_exp = r.RooExponential("bkg_exp", "bkg_exp", M, x_scale)

bkg_norm = bkg_fnc.Integral(100,180)
bkg_exp_norm = r.RooRealVar("bkg_exp_norm","bkg_exp_norm", 0, 2*bkg_norm)
bkg_exp_norm.setVal(bkg_norm)

sig_alpha = r.RooRealVar("sig_alpha","sig_alpha", 1, 3)
sig_n = r.RooRealVar("sig_n","sig_n", 1, 3)
sig_width = r.RooRealVar("sig_width","sig_width", 0.1, 10)
sig_mean = r.RooRealVar("sig_mean","sig_mean", 120, 130)
# sig_pdf = r.RooGaussian("sig_pdf","sig_pdf", M, sig_mean, sig_width)
sig_pdf = r.RooCBShape("sig_pdf", "sig_pdf", M, sig_mean, sig_width, sig_alpha, sig_n)
sig_mean.setVal(sig_fnc.GetParameter(1))
sig_width.setVal(sig_fnc.GetParameter(2))
sig_alpha.setVal(sig_fnc.GetParameter(3))
sig_n.setVal(sig_fnc.GetParameter(4))
sig_norm = sig_fnc.Integral(100,180)
sig_pdf_norm = r.RooRealVar("sig_pdf_norm", "sig_pdf_norm", 0, 2*sig_norm)
sig_pdf_norm.setVal(sig_norm)
sig_pdf_norm.setConstant(True)
sig_alpha.setConstant(True)
sig_n.setConstant(True)

getattr(ws, 'import')(data)
getattr(ws, 'import')(bkg_exp)
getattr(ws, 'import')(bkg_exp_norm)
getattr(ws, 'import')(sig_pdf)
getattr(ws, 'import')(sig_pdf_norm)
getattr(ws, 'import')(sig_alpha)
getattr(ws, 'import')(sig_n)

ws.writeToFile("unblinded_workspace_rho.root")
