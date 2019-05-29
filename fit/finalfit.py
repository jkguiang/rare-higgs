import ROOT as r
import json
from ppmUtils import ConvertToPoissonGraph

sig_plus_bkg = True
analysis = "phi"
hads = "#pi" if analysis == "rho" else "K"

r.gROOT.SetBatch(1)
r.gStyle.SetOptStat(0)
r.gROOT.ProcessLine(".L crystalBall.C")

data_points = []
signal_points = []
with open("unblinded_data{}.json".format("_"+analysis), "r") as fin:
    data_points = json.load(fin)
with open("signal{}.json".format("_"+analysis), "r") as fin:
    signal_points = json.load(fin)
N_sig = len(signal_points)*0.08330784*5e-3

h_data = r.TH1D("h_data", "", 80,100,180)
for dp in data_points:
    h_data.Fill(dp)
# Data
g_data = r.TGraphAsymmErrors()
ConvertToPoissonGraph(h_data, g_data, drawZeros=False)

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
t.Draw("mcand_photon_mass[best_"+analysis+"Cand_idx]>>h_sig", "fabs(mcand_mass_kaon[best_"+analysis+"Cand_idx]-1.02)<0.02 && lepton_pt>35", "goff")
h_sig.Scale(1.0/ h_sig.GetMaximum())
h_sig.Fit("sig_fnc", "N", "goff", 100, 180)
sig_fnc.SetParameter(0, N_sig * sig_fnc.GetParameter(0) / sig_fnc.Integral(100,180))
sig_mean = sig_fnc.GetParameter(1)
sig_width = sig_fnc.GetParameter(2)
sig_alpha = sig_fnc.GetParameter(3)
sig_n = sig_fnc.GetParameter(4)
sig_fnc.SetNpx(1000)

# Fit function
tfile = r.TFile("fitDiagnostics{}.root".format("_"+analysis))
fit_sigbkg = tfile.Get("shapes_fit_s/total_overall")
if sig_plus_bkg:
    fit_bkg = tfile.Get("shapes_fit_s/total_background")
else:
    fit_bkg = tfile.Get("shapes_fit_b/total_background")
fit_bkg.Scale(fit_bkg.GetNbinsX()/h_data.GetNbinsX())

# Plot fit
gr_fit = r.TGraphErrors()
gr_sigbkg = r.TGraphErrors()
gr_err1 = r.TGraphErrors()
gr_err2 = r.TGraphErrors()
gr_fit.SetPoint(0, 100, fit_bkg.GetBinContent(1))
gr_sigbkg.SetPoint(0, 100, fit_sigbkg.GetBinContent(1))
gr_err1.SetPoint(0, 100, fit_bkg.GetBinContent(1))
gr_err2.SetPoint(0, 100, fit_bkg.GetBinContent(1))
gr_err1.SetPointError(0, 0, fit_bkg.GetBinError(1))
gr_err2.SetPointError(0, 0, 2*fit_bkg.GetBinError(1))
nbins = fit_bkg.GetNbinsX()
for i in range(1, nbins+1):
    gr_fit.SetPoint(i, (i-0.5)/100*80.+100, fit_bkg.GetBinContent(i))
    gr_sigbkg.SetPoint(i, (i-0.5)/100*80.+100, fit_sigbkg.GetBinContent(i))
    gr_err1.SetPoint(i, (i-0.5)/100*80.+100, fit_bkg.GetBinContent(i))
    gr_err2.SetPoint(i, (i-0.5)/100*80.+100, fit_bkg.GetBinContent(i))
    gr_err1.SetPointError(i, 0, fit_bkg.GetBinError(i))
    gr_err2.SetPointError(i, 0, 2*fit_bkg.GetBinError(i))
gr_fit.SetPoint(gr_fit.GetN(), 180, fit_bkg.GetBinContent(nbins))
gr_sigbkg.SetPoint(gr_sigbkg.GetN(), 180, fit_sigbkg.GetBinContent(nbins))
gr_err1.SetPoint(gr_err1.GetN(), 180, fit_bkg.GetBinContent(nbins))
gr_err2.SetPoint(gr_err2.GetN(), 180, fit_bkg.GetBinContent(nbins))
gr_err1.SetPointError(gr_err1.GetN()-1, 0, fit_bkg.GetBinError(nbins))
gr_err2.SetPointError(gr_err2.GetN()-1, 0, 2*fit_bkg.GetBinError(nbins))

# Use h_data as formatting device
h_data.Reset()
g_dummy = r.TGraph()
g_dummy.SetLineColor(r.kWhite)

# Canvas
c2 = r.TCanvas("c2","c2", 700, 600)
c2.SetTicky(1)

# --> Style <-- #
# x-axis
h_data.GetXaxis().SetTitle("m_{"+hads+"^{+}"+hads+"^{#minus}#gamma} (GeV)")
h_data.GetXaxis().SetTitleSize(0.04)
h_data.GetXaxis().SetLabelOffset(0.01)
h_data.GetXaxis().SetTitleOffset(1.1)
# y-axis
h_data.GetYaxis().SetTitle("Events / 1 GeV")
h_data.GetYaxis().SetTitleSize(0.04)
h_data.GetYaxis().SetRangeUser(0,10)
# Marker and Line styles
gr_fit.SetLineWidth(2)
gr_fit.SetLineColor(r.kRed)
if sig_plus_bkg:
    gr_fit.SetLineStyle(2)
    gr_sigbkg.SetLineWidth(2)
    gr_sigbkg.SetLineColor(r.kRed)
sig_fnc.SetLineColor(r.kBlue)
g_data.SetMarkerStyle(20)
g_data.SetMarkerColor(r.kBlack)
g_data.SetLineColor(r.kBlack)
# Fill
gr_err1.SetFillColor(r.kGreen+1)
gr_err2.SetFillColor(r.kOrange-2)
sig_fnc.SetFillColorAlpha(r.kBlue, 0.35)
sig_fnc.SetFillStyle(1001)
# --> Plot <-- #
h_data.Draw()
gr_err2.Draw("3 SAME")
gr_err1.Draw("3 SAME")
gr_fit.Draw("L SAME")
if sig_plus_bkg: gr_sigbkg.Draw("L SAME")
else: sig_fnc.Draw("LF SAME")
g_data.Draw("SAME PZ")
h_data.Draw("SAME AXIS")
# --> Legend <-- #
# Options
l = r.TLegend(0.12, 0.7, 0.6, 0.88)
l.SetBorderSize(0)
l.SetNColumns(2)
l.SetFillStyle(0)
# Draw statements
l.AddEntry(g_data, "Data", "ep")
if sig_plus_bkg:
    l.AddEntry(gr_fit, "B component", "l")
    l.AddEntry(gr_sigbkg, "S+B fit", "l")
    l.AddEntry(gr_err1, "#pm1#sigma", "f")
    l.AddEntry(g_dummy, "", "l")
    l.AddEntry(gr_err2, "#pm2#sigma", "f")
else:
    l.AddEntry(gr_fit, "B-Only Fit", "l")
    l.AddEntry(sig_fnc, "H#rightarrow#"+analysis+"#gamma", "lf")
    l.AddEntry(gr_err1, "#pm1#sigma", "f")
    l.AddEntry(g_dummy, "#(){BR = 5#times10^{-3}}  ", "l")
    l.AddEntry(gr_err2, "#pm2#sigma", "f")
l.Draw()
# --> Text <-- #
# Options
tl = r.TLatex()
tl.SetNDC(1)
tl.SetTextSize(.045)
tl.SetTextFont(42)
# Draw statements
tl.SetTextAlign(11)
tl.DrawLatex(0.11,0.91,"#bf{CMS} #it{Preliminary}")
tl.SetTextAlign(31)
tl.DrawLatex(0.89,0.91,"137 fb^{-1} (13 TeV)")
if sig_plus_bkg:
    fit_result = tfile.Get("fit_s")
    mu = fit_result.floatParsFinal().at(1).getValV()
    error_up = fit_result.floatParsFinal().at(1).getErrorHi()
    error_down = fit_result.floatParsFinal().at(1).getErrorLo()
    tl.SetTextColor(r.kBlue+2)
    tl.SetTextAlign(13)
    tl.SetTextSize(0.04)
    tl.DrawLatex(0.6,0.82, "#hat{{ #scale[0.8]{{BR}} }} = {0:.1f}^{{+{1:.1f} }}_{{#minus{2:.1f} }} #times 10^{{-3}}".format(mu*5, abs(error_up*5), abs(error_down*5)))
    tl.DrawLatex(0.605,0.87, "H#rightarrow#"+analysis+"#gamma")

# Save
c2.SaveAs("finalfit{0}{1}.pdf".format("_"+analysis, "_bgOnlyFit" if not sig_plus_bkg else ""))
