{
  TString sampleName = "/hadoop/cms/store/group/snt/run2_mc2018/WGToLNuG_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1_MINIAODSIM_CMS4_V10-02-04/merged_ntuple_1.root";
  gSystem->Load("./CMS3_CORE.so");

  gROOT->ProcessLine(".L mcTree.C++");
  gROOT->ProcessLine(".L BuildMCTree.C++");

  // MONTE CARLO
  // TChain *priv_mc = new TChain("Events"); 
  // priv_mc->Add(sampleName);
  // BuildMCTree(priv_mc, "WGamma_mcTree.root", sampleName); // ScanChain(TChain, outName, sampleName, ...)
  // END MONTE CARLO

}
