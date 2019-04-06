{
  gSystem->Load("/home/users/jguiang/projects/mt2/MT2Analysis/CORE/CMS3_CORE.so");

  gROOT->ProcessLine(".L mcTree.C+");
  gROOT->ProcessLine(".L BuildMCTree.C+");

  // MONTE CARLO
  TChain *priv_mc = new TChain("Events"); 
  priv_mc->Add("/hadoop/cms/store/group/snt/run2_mc2018_private/WH_HtoRhoGammaPhiGamma_privateMC_102x_MINIAOD_v1/merged_ntuple_1.root");
  BuildMCTree(priv_mc, "mcTree.root");
  // END MONTE CARLO

}
