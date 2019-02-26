{
  gSystem->Load("/home/users/jguiang/projects/mt2/MT2Analysis/CORE/CMS3_CORE.so");

  gROOT->ProcessLine(".L SanityCheck.C+");

  // MONTE CARLO
  TChain *priv_mc = new TChain("Events"); 
  priv_mc->Add("/hadoop/cms/store/user/bemarsh/ProjectMetis/WH_HtoRhoGammaPhiGamma_privateMC_102x_MINIAOD_v1/merged_ntuple_1.root");
  SanityCheck(priv_mc, "priv_mc.root");
  // END MONTE CARLO

}
