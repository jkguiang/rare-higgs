// -*- C++ -*-

// C++
#include <iostream>
#include <vector>

// ROOT
#include "TBenchmark.h"
#include "TChain.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TROOT.h"
#include "TTreeCache.h"
#include "TString.h"

// CORE
#include "CORE/CMS3.h"
#include "CORE/ElectronSelections.h"
#include "CORE/MuonSelections.h"
#include "CORE/PhotonSelections.h"
#include "CORE/IsolationTools.h"
#include "CORE/TriggerSelections.h"
#include "CORE/MetSelections.h"

// Tools
#include "CORE/Tools/JetCorrector.h"
#include "CORE/Tools/goodrun.h"
#include "CORE/Tools/jetcorr/FactorizedJetCorrector.h"
#include "CORE/Tools/datasetinfo/getDatasetInfo.h"

// Custom
#include "babyTree.h"
#include "magicAngles.h"

// Namespaces
using namespace std;
using namespace tas;

/**
 * Initialize TTree
 * @params: none
 * @return: none
 */
BabyTree::BabyTree() {
    /* --> TTree Setup <-- */
    t = new TTree("tree", "tree");
    /* --> Meta Branches Setup <-- */
    // Event info
    b_run = t->Branch("run", &run, "run/I");
    b_lumi = t->Branch("lumi", &lumi, "lumi/I");
    b_event = t->Branch("event", &event, "event/I");
    b_scale1fb = t->Branch("scale1fb", &scale1fb, "scale1fb/F");
    b_isGold = t->Branch("isGold", &isGold, "isGold/I");
    b_isHEM = t->Branch("isHEM", &isHEM, "isHEM/I");
    // MET
    b_met_pt = t->Branch("met_pt", &met_pt, "met_pt/F"); 
    b_met_phi = t->Branch("met_phi", &met_phi, "met_phi/F"); 
    b_rawMet_pt = t->Branch("rawMet_pt", &rawMet_pt, "rawMet_pt/F"); 
    b_rawMet_phi = t->Branch("rawMet_phi", &rawMet_phi, "rawMet_phi/F");
    // Triggers
    b_HLT_singleMu = t->Branch("HLT_singleMu", &HLT_singleMu, "HLT_singleMu/I"); 
    b_HLT_singleEl = t->Branch("HLT_singleEl", &HLT_singleEl, "HLT_singleEl/I"); 
    // Filters
    b_passFilters = t->Branch("passFilters", &passFilters, "passFilters/I");
    /* --> Special Branches Setup <-- */
    // dR(gen, reco)
    b_genRecoGamma_dR = t->Branch("genRecoGamma_dR", &genRecoGamma_dR, "genRecoGamma_dR/F");
    b_genRecoPhi_dR = t->Branch("genRecoPhi_dR", &genRecoPhi_dR, "genRecoPhi_dR/F");
    b_genRecoRho_dR = t->Branch("genRecoRho_dR", &genRecoRho_dR, "genRecoRho_dR/F");
    // Magic Angles
    b_recoMagAng_cosThetaStar = t->Branch("recoMagAng_cosThetaStar", &recoMagAng_cosThetaStar, "recoMagAng_cosThetaStar/F");
    b_recoMagAng_cosTheta1 = t->Branch("recoMagAng_cosTheta1", &recoMagAng_cosTheta1, "recoMagAng_cosTheta1/F");
    b_recoMagAng_cosTheta2 = t->Branch("recoMagAng_cosTheta2", &recoMagAng_cosTheta2, "recoMagAng_cosTheta2/F");
    b_recoMagAng_Phi = t->Branch("recoMagAng_Phi", &recoMagAng_Phi, "recoMagAng_Phi/F");
    b_recoMagAng_Phi1 = t->Branch("recoMagAng_Phi1", &recoMagAng_Phi1, "recoMagAng_Phi1/F");
    b_recoMagAng_m1 = t->Branch("recoMagAng_m1", &recoMagAng_m1, "recoMagAng_m1/F");
    b_recoMagAng_m2 = t->Branch("recoMagAng_m2", &recoMagAng_m2, "recoMagAng_m2/F");
    b_genMagAng_cosThetaStar = t->Branch("genMagAng_cosThetaStar", &genMagAng_cosThetaStar, "genMagAng_cosThetaStar/F");
    b_genMagAng_cosTheta1 = t->Branch("genMagAng_cosTheta1", &genMagAng_cosTheta1, "genMagAng_cosTheta1/F");
    b_genMagAng_cosTheta2 = t->Branch("genMagAng_cosTheta2", &genMagAng_cosTheta2, "genMagAng_cosTheta2/F");
    b_genMagAng_Phi = t->Branch("genMagAng_Phi", &genMagAng_Phi, "genMagAng_Phi/F");
    b_genMagAng_Phi1 = t->Branch("genMagAng_Phi1", &genMagAng_Phi1, "genMagAng_Phi1/F");
    b_genMagAng_m1 = t->Branch("genMagAng_m1", &genMagAng_m1, "genMagAng_m1/F");
    b_genMagAng_m2 = t->Branch("genMagAng_m2", &genMagAng_m2, "genMagAng_m2/F");
    /* --> Gen Branches Setup <-- */
    // Gen W
    b_genW_pt = t->Branch("genW_pt", &genW_pt, "genW_pt/F");
    b_genW_eta = t->Branch("genW_eta", &genW_eta, "genW_eta/F");
    b_genW_phi = t->Branch("genW_phi", &genW_phi, "genW_phi/F");
    b_genW_mass = t->Branch("genW_mass", &genW_mass, "genW_mass/F");
    // Gen Leptons
    b_genWLepton_id = t->Branch("genWLepton_id", &genWLepton_id, "genWLepton_id/I");
    b_genWLepton_pt = t->Branch("genWLepton_pt", &genWLepton_pt, "genWLepton_pt/F");
    b_genWLepton_eta = t->Branch("genWLepton_eta", &genWLepton_eta, "genWLepton_eta/F");
    b_genWLepton_phi = t->Branch("genWLepton_phi", &genWLepton_phi, "genWLepton_phi/F");
    // Gen Higgs
    b_genHiggs_pt = t->Branch("genHiggs_pt", &genHiggs_pt, "genHiggs_pt/F");
    b_genHiggs_eta = t->Branch("genHiggs_eta", &genHiggs_eta, "genHiggs_eta/F");
    b_genHiggs_phi = t->Branch("genHiggs_phi", &genHiggs_phi, "genHiggs_phi/F");
    b_genHiggs_mass = t->Branch("genHiggs_mass", &genHiggs_mass, "genHiggs_mass/F");
    // Gen Mesons from Higgs
    b_genHiggsMeson_id = t->Branch("genHiggsMeson_id", &genHiggsMeson_id, "genHiggsMeson_id/I");
    b_genHiggsMeson_pt = t->Branch("genHiggsMeson_pt", &genHiggsMeson_pt, "genHiggsMeson_pt/F");
    b_genHiggsMeson_eta = t->Branch("genHiggsMeson_eta", &genHiggsMeson_eta, "genHiggsMeson_eta/F");
    b_genHiggsMeson_phi = t->Branch("genHiggsMeson_phi", &genHiggsMeson_phi, "genHiggsMeson_phi/F");
    b_genHiggsMeson_mass = t->Branch("genHiggsMeson_mass", &genHiggsMeson_mass, "genHiggsMeson_mass/F");
    b_genHiggsMesonGamma_dR = t->Branch("genHiggsMesonGamma_dR", &genHiggsMesonGamma_dR, "genHiggsMesonGamma_dR/F");
    // Gen K+, K-
    b_genKm_pt = t->Branch("genKm_pt", &genKm_pt, "genKm_pt/F");
    b_genKm_phi = t->Branch("genKm_phi", &genKm_phi, "genKm_phi/F");
    b_genKm_eta = t->Branch("genKm_eta", &genKm_eta, "genKm_eta/F");
    b_genKp_pt = t->Branch("genKp_pt", &genKp_pt, "genKp_pt/F");
    b_genKp_phi = t->Branch("genKp_phi", &genKp_phi, "genKp_phi/F");
    b_genKp_eta = t->Branch("genKp_eta", &genKp_eta, "genKp_eta/F");
    b_genKpKm_dR = t->Branch("genKpKm_dR", &genKpKm_dR, "genKpKm_dR/F");
    // Gen Photons
    b_genGamma_pt = t->Branch("genGamma_pt", &genGamma_pt, "genGamma_pt/F");
    b_genGamma_phi = t->Branch("genGamma_phi", &genGamma_phi, "genGamma_phi/F");
    b_genGamma_eta = t->Branch("genGamma_eta", &genGamma_eta, "genGamma_eta/F");
    /* --> Reco Branches Setup <-- */
    // Reco Higgs
    b_recoHiggs_mass = t->Branch("recoHiggs_mass", &recoHiggs_mass, "recoHiggs_mass/F");
    // Reco Mesons
    b_recoMeson_nCands = t->Branch("recoMeson_nCands", &recoMeson_nCands, "recoMeson_nCands/I");
    // Reco Phi
    b_recoPhi_mass = t->Branch("recoPhi_mass", &recoPhi_mass, "recoPhi_mass/F");
    b_recoPhi_pt = t->Branch("recoPhi_pt", &recoPhi_pt, "recoPhi_pt/F");
    b_recoPhi_eta = t->Branch("recoPhi_eta", &recoPhi_eta, "recoPhi_eta/F");
    b_recoPhi_phi = t->Branch("recoPhi_phi", &recoPhi_phi, "recoPhi_phi/F");
    b_recoPhi_relIso = t->Branch("recoPhi_relIso", &recoPhi_relIso, "recoPhi_relIso/F");
    b_recoPhiGamma_dR = t->Branch("recoPhiGamma_dR", &recoPhiGamma_dR, "recoPhiGamma_dR/F");
    // Reco K+, K-
    b_recoKm_pt = t->Branch("recoKm_pt", &recoKm_pt, "recoKm_pt/F");
    b_recoKm_eta = t->Branch("recoKm_eta", &recoKm_eta, "recoKm_eta/F");
    b_recoKm_phi = t->Branch("recoKm_phi", &recoKm_phi, "recoKm_phi/F");
    b_recoKm_relIso = t->Branch("recoKm_relIso", &recoKm_relIso, "recoKm_relIso/F");
    b_recoKp_pt = t->Branch("recoKp_pt", &recoKp_pt, "recoKp_pt/F");
    b_recoKp_eta = t->Branch("recoKp_eta", &recoKp_eta, "recoKp_eta/F");
    b_recoKp_phi = t->Branch("recoKp_phi", &recoKp_phi, "recoKp_phi/F");
    b_recoKp_relIso = t->Branch("recoKp_relIso", &recoKp_relIso, "recoKp_relIso/F");
    b_recoKpKm_dR = t->Branch("recoKpKm_dR", &recoKpKm_dR, "recoKpKm_dR/F");
    // Reco Rho
    b_recoRho_mass = t->Branch("recoRho_mass", &recoRho_mass, "recoRho_mass/F");
    b_recoRho_pt = t->Branch("recoRho_pt", &recoRho_pt, "recoRho_pt/F");
    b_recoRho_eta = t->Branch("recoRho_eta", &recoRho_eta, "recoRho_eta/F");
    b_recoRho_phi = t->Branch("recoRho_phi", &recoRho_phi, "recoRho_phi/F");
    b_recoRho_relIso = t->Branch("recoRho_relIso", &recoRho_relIso, "recoRho_relIso/F");
    b_recoRhoGamma_dR = t->Branch("recoRhoGamma_dR", &recoRhoGamma_dR, "recoRhoGamma_dR/F");
    // Reco Pi+, Pi-
    b_recoPim_pt = t->Branch("recoPim_pt", &recoPim_pt, "recoPim_pt/F");
    b_recoPim_eta = t->Branch("recoPim_eta", &recoPim_eta, "recoPim_eta/F");
    b_recoPim_phi = t->Branch("recoPim_phi", &recoPim_phi, "recoPim_phi/F");
    b_recoPim_relIso = t->Branch("recoPim_relIso", &recoPim_relIso, "recoPim_relIso/F");
    b_recoPip_pt = t->Branch("recoPip_pt", &recoPip_pt, "recoPip_pt/F");
    b_recoPip_eta = t->Branch("recoPip_eta", &recoPip_eta, "recoPip_eta/F");
    b_recoPip_phi = t->Branch("recoPip_phi", &recoPip_phi, "recoPip_phi/F");
    b_recoPip_relIso = t->Branch("recoPip_relIso", &recoPip_relIso, "recoPip_relIso/F");
    b_recoPipPim_dR = t->Branch("recoPipPim_dR", &recoPipPim_dR, "recoPipPim_dR/F");
    // Reco Photons
    b_recoGamma_pt = t->Branch("recoGamma_pt", &recoGamma_pt, "recoGamma_pt/F");
    b_recoGamma_phi = t->Branch("recoGamma_phi", &recoGamma_phi, "recoGamma_phi/F");
    b_recoGamma_eta = t->Branch("recoGamma_eta", &recoGamma_eta, "recoGamma_eta/F");
    b_recoGamma_relIso = t->Branch("recoGamma_relIso", &recoGamma_relIso, "recoGamma_relIso/F");
    b_genRecoGamma_isMatch = t->Branch("genRecoGamma_isMatch", &genRecoGamma_isMatch, "genRecoGamma_isMatch/I");
    b_minGammaParton_dR = t->Branch("minGammaParton_dR", &minGammaParton_dR, "minGammaParton_dR/F");
    // Reco Leptons
    b_recoWLepton_id = t->Branch("recoWLepton_id", &recoWLepton_id, "recoWLepton_id/I");
    b_recoWLepton_pt = t->Branch("recoWLepton_pt", &recoWLepton_pt, "recoWLepton_pt/F");
    b_recoWLepton_eta = t->Branch("recoWLepton_eta", &recoWLepton_eta, "recoWLepton_eta/F");
    b_recoWLepton_phi = t->Branch("recoWLepton_phi", &recoWLepton_phi, "recoWLepton_phi/F");
    b_recoWLepton_nLep = t->Branch("recoWLepton_nLep", &recoWLepton_nLep, "recoWLepton_nLep/I");
    b_recoGammaWLepton_dR = t->Branch("recoGammaWLepton_dR", &recoGammaWLepton_dR, "recoGammaWLepton_dR/F");
}

/**
 * Reset TTree branch values
 * @params: none
 * @return: none
 */
void BabyTree::Reset() {
    // Meta
    run = -999;
    lumi = -999;
    event = -999;
    scale1fb = -999;
    isGold = -999;
    isHEM = -999;
    met_pt = -999; 
    met_phi = -999;
    rawMet_pt = -999;
    rawMet_phi = -999;
    HLT_singleMu = -999; 
    HLT_singleEl = -999; 
    passFilters = -999;
    // Special
    genRecoGamma_dR = -999;
    genRecoPhi_dR = -999;
    genRecoRho_dR = -999;
    recoMagAng_cosThetaStar = -999;
    recoMagAng_cosTheta1 = -999;
    recoMagAng_cosTheta2 = -999;
    recoMagAng_Phi = -999;
    recoMagAng_Phi1 = -999;
    recoMagAng_m1 = -999;
    recoMagAng_m2 = -999;
    genMagAng_cosThetaStar = -999;
    genMagAng_cosTheta1 = -999;
    genMagAng_cosTheta2 = -999;
    genMagAng_Phi = -999;
    genMagAng_Phi1 = -999;
    genMagAng_m1 = -999;
    genMagAng_m2 = -999;
    // Gen
    genW_pt = -999;
    genW_eta = -999;
    genW_phi = -999;
    genW_mass = -999;
    genWLepton_id = -999;
    genWLepton_pt = -999;
    genWLepton_eta = -999;
    genWLepton_phi = -999;
    genHiggs_pt = -999;
    genHiggs_eta = -999;
    genHiggs_phi = -999;
    genHiggs_mass = -999;
    genHiggsMeson_id = -999;
    genHiggsMeson_pt = -999;
    genHiggsMeson_eta = -999;
    genHiggsMeson_phi = -999;
    genHiggsMeson_mass = -999;
    genHiggsMesonGamma_dR = -999;
    genGamma_pt = -999;
    genGamma_phi = -999;
    genGamma_eta = -999;
    genKm_pt = -999;
    genKm_phi = -999;
    genKm_eta = -999;
    genKp_pt = -999;
    genKp_phi = -999;
    genKp_eta = -999;
    genKpKm_dR = -999;
    // Reco
    recoHiggs_mass = -999;
    recoMeson_nCands = -999;
    recoPhi_mass = -999;
    recoPhi_pt = -999;
    recoPhi_eta = -999;
    recoPhi_phi = -999;
    recoPhi_relIso = -999;
    recoPhiGamma_dR = -999;
    recoKm_pt = -999;
    recoKm_eta = -999;
    recoKm_phi = -999;
    recoKm_relIso = -999;
    recoKp_pt = -999;
    recoKp_eta = -999;
    recoKp_phi = -999;
    recoKp_relIso = -999;
    recoKpKm_dR = -999;
    recoRho_mass = -999;
    recoRho_pt = -999;
    recoRho_eta = -999;
    recoRho_phi = -999;
    recoRho_relIso = -999;
    recoRhoGamma_dR = -999;
    recoPim_pt = -999; 
    recoPim_eta = -999;
    recoPim_phi = -999;
    recoPim_relIso = -999;
    recoPip_pt = -999;
    recoPip_eta = -999;
    recoPip_phi = -999;
    recoPip_relIso = -999;
    recoPipPim_dR = -999;
    recoGamma_pt = -999;
    recoGamma_phi = -999;
    recoGamma_eta = -999;
    recoGamma_relIso = -999;
    genRecoGamma_isMatch = -999;
    minGammaParton_dR = -999;
    recoWLepton_id = -999;
    recoWLepton_pt = -999;
    recoWLepton_eta = -999;
    recoWLepton_phi = -999;
    recoWLepton_nLep = -999;

    return;
}

/**
 * Fill config struct from sample name
 * @params: sample name
 * @return: none
 */
void BabyTree::MakeConfig(TString sample) {
    // Get year, determine if data or MC
    config.isData = (sample.Contains("Run20"));
    if (config.isData) {
        if (sample.Contains("Run2016")) config.year = 2016;
        else if (sample.Contains("Run2017")) config.year = 2017;
        else if (sample.Contains("Run2018")) config.year = 2018;
    }
    else if (sample.Contains("WH_HtoRhoGammaPhiGamma")) config.year = 2018;
    else if (sample.Contains("RunIISummer16") && sample.Contains("94X")) config.year = 2016;
    else if (sample.Contains("RunIIFall17") && sample.Contains("94X")) config.year = 2017;
    else if (sample.Contains("RunIIAutumn18") && sample.Contains("102X")) config.year = 2018;
    // Set CORE global year variable
    gconf.year = config.year;
    gconf.ea_version = (config.year == 2016) ? 1 : 3;
    // Get intergrated lumi
    config.lumi = (config.year == 2016) ? 35.92 : (config.year == 2017) ? 41.53 : 59.74;
    // Get golden JSON
    if (config.year == 2016) config.json = "jsons/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON_snt.txt";
    if (config.year == 2017) config.json = "jsons/Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON_v1_snt.txt";
    if (config.year == 2018) config.json = "jsons/Cert_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18_JSON_snt.txt";
    // Get JEC tag for dynamic file collection
    string jecTag, jecEra;
    if (config.year == 2018) {
        if (config.isData) {
            if (sample.Contains("Run2018A")) jecEra = "RunA";
            else if (sample.Contains("Run2018B")) jecEra = "RunB";
            else if (sample.Contains("Run2018C")) jecEra = "RunC";
            else if (sample.Contains("Run2018D")) jecEra = "RunD";
            jecTag = "Autumn18_"+jecEra+"_V8_DATA";
        }
        else jecTag = "Autumn18_V8_MC";
    }
    if (config.year == 2017) {
        if (config.isData) {
            if (sample.Contains("2017B")) jecEra = "17Nov2017B";
            else if (sample.Contains("2017C")) jecEra = "17Nov2017C";
            else if (sample.Contains("2017D")) jecEra = "17Nov2017DE";
            else if (sample.Contains("2017E")) jecEra = "17Nov2017DE";
            else if (sample.Contains("2017F")) jecEra = "17Nov2017F";
            jecTag = "Fall17_"+jecEra+"_V32_DATA";
        }
        else jecTag = "Fall17_17Nov2017_V32_MC";
    }
    if (config.year == 2016) {
        if (config.isData) {
            if (sample.Contains("Run2016B")) jecEra = "07Aug2017BCD";
            else if (sample.Contains("Run2016C")) jecEra = "07Aug2017BCD";
            else if (sample.Contains("Run2016D")) jecEra = "07Aug2017BCD";
            else if (sample.Contains("Run2016E")) jecEra = "07Aug2017EF";
            else if (sample.Contains("Run2016F")) jecEra = "07Aug2017EF";
            else if (sample.Contains("Run2016G")) jecEra = "07Aug2017GH";
            else if (sample.Contains("Run2016H")) jecEra = "07Aug2017GH";
            jecTag = "Summer16_"+jecEra+"_V11_DATA";
        }
        else jecTag = "Summer16_07Aug2017_V11_MC";
    }
    cout << "JEC Tag: " << jecTag << endl;
    // Collect jet correction files
    vector<string> jetCorr_files;
    jetCorr_files.push_back("jetCorrections/"+jecTag+"_L1FastJet_AK4PFchs.txt");
    jetCorr_files.push_back("jetCorrections/"+jecTag+"_L2Relative_AK4PFchs.txt");
    jetCorr_files.push_back("jetCorrections/"+jecTag+"_L3Absolute_AK4PFchs.txt");
    jetCorr_files.push_back("jetCorrections/"+jecTag+"_L2L3Residual_AK4PFchs.txt");
    // Make jet corrector
    config.jetCorr = makeJetCorrector(jetCorr_files);
    return;
}

/**
 * Get dR between two four-vectors
 * @params: phi, eta from two four-momenta
 * @return: dR value
 */
float BabyTree::dR(float phi1, float phi2, float eta1, float eta2) {
    float dphi = abs(phi2 - phi1);
    if (dphi > M_PI){
        dphi = 2*M_PI - dphi;
    }
    return sqrt(pow((dphi), 2) + pow((eta2 - eta1), 2));
}

/**
 * Fill relevant branches with gen-level information
 * @params: none
 * @return: none
 */
void BabyTree::FillGenBranches() {

    /* --> Find and save relevant decay modes <-- */

    // Decay mode tracking information
    enum modes { H_to_PhiGamma = 0, H_to_RhoGamma = 1, H_to_KKGamma = 2, Phi_to_KK = 3, W_to_ElNu = 4, W_to_MuNu = 5 };
    int mothers[6] = { 0,0,0,0,0,0 };
    vector<int> daughters[6];

    // START Loop over gen-level data ------------
    for (unsigned int i = 0; i < genps_id().size(); i++) {
        // Get PDG ID of daugher, mother
        int thisID = genps_id().at(i);
        int thisMotherID = genps_id_mother().at(i);
        int thisMotherIdx = genps_idx_mother().at(i);
        // Veto duplicates
        if (!genps_isLastCopy().at(i)) continue;
        if (!genps_isLastCopy()[thisMotherIdx]) continue;
        // Veto weird photons
        if (thisID == 22 && !genps_isHardProcess().at(i)) continue;
        // Fill branches for mother particles
        if (thisID == 25) {
            // Fill Higgs boson branches
            genHiggs_pt = genps_p4().at(i).pt();
            genHiggs_eta = genps_p4().at(i).eta();
            genHiggs_phi = genps_p4().at(i).phi(); 
        }
        if (abs(thisID) == 24) {
            // Fill W boson branches
            genW_pt = genps_p4().at(i).pt();
            genW_eta = genps_p4().at(i).eta();
            genW_phi = genps_p4().at(i).phi(); 
        }
        // Add particle id, idx, p4 to relevant systems
        if (abs(thisMotherID) == 24) {
            if (abs(thisID) == 11) {
                mothers[W_to_ElNu] = thisMotherIdx;
                daughters[W_to_ElNu].push_back(i);
            }
            if (abs(thisID) == 13) {
                mothers[W_to_MuNu] = thisMotherIdx;
                daughters[W_to_MuNu].push_back(i);
            }
            if (abs(thisID) == 12) {
                daughters[W_to_ElNu].push_back(i);
            }
            if (abs(thisID) == 14) {
                daughters[W_to_MuNu].push_back(i);
            }
        }
        if (thisMotherID == 25) {
            if (thisID == 22 || thisID == 333) {
                mothers[H_to_PhiGamma] = thisMotherIdx;
                daughters[H_to_PhiGamma].push_back(i);
            }
            if (thisID == 22 || thisID == 113) {
                mothers[H_to_RhoGamma] = thisMotherIdx;
                daughters[H_to_RhoGamma].push_back(i);
            }
            if (thisID == 22) {
                mothers[H_to_KKGamma] = thisMotherIdx;
                daughters[H_to_KKGamma].push_back(i);
            }
        }
        if (thisMotherID == 333) {
            if (abs(thisID) == 321 && genps_id_mother()[thisMotherIdx] == 25) {
                // Fill Phi->KK Mother/Daughters
                mothers[Phi_to_KK] = thisMotherIdx;
                daughters[Phi_to_KK].push_back(i);
                // Fill H->KK,Gamma Daughters
                daughters[H_to_KKGamma].push_back(i);
            }
        }
    } // END Loop over gen-level data ------------

    // General variables used for filling tree and code readability
    LorentzVector p4_sum;
    int decay;

    /* --> Fill W --> el/mu, nu <-- */
    LorentzVector lepton_p4, nu_p4;
    if (daughters[W_to_MuNu].size() == 2 || daughters[W_to_ElNu].size() == 2) {
        decay = (daughters[W_to_MuNu].size() == 1) ? W_to_MuNu : W_to_ElNu;
        p4_sum -= p4_sum; // Reset p4 sum
        for (unsigned int j = 0; j < daughters[decay].size(); j++) {
            // Get daughter
            int idx = daughters[decay].at(j);
            int id = genps_id().at(idx);
            LorentzVector p4 = genps_p4().at(idx);
            // Fill leptons from W branches
            if (abs(id) == 11 || abs(id) == 13) {
                lepton_p4 = p4;
                genWLepton_id = id;
                genWLepton_pt = p4.pt();
                genWLepton_eta = p4.eta();
                genWLepton_phi = p4.phi();
                p4_sum += p4;
            }
            else {
                // Save neutrino kinematics for magic angles calculation
                nu_p4 = p4;
            }
        }
        genW_mass = p4_sum.M();
    }

    /* --> Fill H --> phi/rho, gamma <-- */
    LorentzVector gamma_p4;
    if (daughters[H_to_PhiGamma].size() == 2 || daughters[H_to_RhoGamma].size() == 2) {
        decay = (daughters[H_to_PhiGamma].size() == 2) ? H_to_PhiGamma : H_to_RhoGamma;
        p4_sum -= p4_sum; // Reset p4 sum
        for (unsigned int j = 0; j < daughters[decay].size(); j++) {
            // Get daughter
            int idx = daughters[decay].at(j);
            int id = genps_id().at(idx);
            LorentzVector p4 = genps_p4().at(idx);
            // Fill photon branches
            if (id == 22) {
                gamma_p4 = p4;
                genGamma_pt = p4.pt();
                genGamma_eta = p4.eta();
                genGamma_phi = p4.phi();
                p4_sum += p4;
            }
            else {
                // Fill mesons from Higgs branches
                genHiggsMeson_id = id;
                genHiggsMeson_pt = p4.pt();
                genHiggsMeson_eta = p4.eta();
                genHiggsMeson_phi = p4.phi();
                p4_sum += p4;
            }
        }
        genHiggsMesonGamma_dR = dR(genHiggsMeson_phi, genGamma_phi, genHiggsMeson_eta, genGamma_eta);
        genHiggs_mass = p4_sum.M();
    }

    /* --> Fill phi --> K+, K- <-- */
    LorentzVector phi_p4;
    if (daughters[Phi_to_KK].size() == 2) {
        decay = Phi_to_KK;
        p4_sum -= p4_sum; // Reset p4 sum
        for (unsigned int j = 0; j < daughters[decay].size(); j++) {
            // Get daughter
            int idx = daughters[decay].at(j);
            int id = genps_id().at(idx);
            LorentzVector p4 = genps_p4().at(idx);
            // Fill Kaon from Phi branches
            if (id == 321) {
                genKp_pt = p4.pt();
                genKp_phi = p4.phi();
                genKp_eta = p4.eta();
                p4_sum += p4;
            }
            else if (id == -321) {
                genKm_pt = p4.pt();
                genKm_phi = p4.phi();
                genKm_eta = p4.eta();
                p4_sum += p4;
            }
        }
        genKpKm_dR = dR(genKp_phi, genKm_phi, genKp_eta, genKm_eta);
        phi_p4 = p4_sum;
        genHiggsMeson_mass = phi_p4.M();
    }

    /* --> Fill Magic Angles <-- */
    if ((daughters[W_to_MuNu].size() == 2 || daughters[W_to_ElNu].size() == 2) && daughters[H_to_PhiGamma].size() == 2) {
        // Output from MELA
        MagicAngles mAngles;
        // Variables for MELA
        TLorentzVector lep_p4(lepton_p4.Px(),lepton_p4.Py(),lepton_p4.Pz(),lepton_p4.E());
        TLorentzVector mes_p4(phi_p4.Px(),phi_p4.Py(),phi_p4.Pz(),phi_p4.E());
        TLorentzVector gam_p4(gamma_p4.Px(),gamma_p4.Py(),gamma_p4.Pz(),gamma_p4.E());
        // Get magic angles
        mAngles = getAngles(nu_p4.pt(), nu_p4.phi(), genWLepton_id, lep_p4, mes_p4, gam_p4, nu_p4.Pz());
        // Fill branches
        genMagAng_cosThetaStar = mAngles.angles[0];
        genMagAng_cosTheta1 = mAngles.angles[1];
        genMagAng_cosTheta2 = mAngles.angles[2];
        genMagAng_Phi = mAngles.angles[3];
        genMagAng_Phi1 = mAngles.angles[4];
        genMagAng_m1 = mAngles.angles[5];
        genMagAng_m2 = mAngles.angles[6];
    }

    return;
}

/**
 * Fill relevant branches with reco information
 * @params: none
 * @return: none
 */
void BabyTree::FillRecoBranches() {

    /* --> Fill MET <-- */
    // Uncorrected
    rawMet_pt = evt_pfmet();
    rawMet_phi = evt_pfmetPhi();
    // Corrected
    pair<float, float> t1met = getT1CHSMET_fromMINIAOD(config.jetCorr, 0, 0, false, 0);
    met_pt = t1met.first;
    met_phi = t1met.second;

    /* --> Mesons, Hadrons <-- */
    // Store K+, K- indexes
    vector<int> posHadronIdxs;
    vector<int> negHadronIdxs;
    // START Loop over pfcands -----------------------------
    for (unsigned int i = 0; i < pfcands_p4().size(); i++) {
        // Find Kaons (and Pions)
        int thisID = pfcands_particleId().at(i);
        if (thisID == 211) posHadronIdxs.push_back(i);
        else if (thisID == -211) negHadronIdxs.push_back(i);
    } // END Loop over pfcands -----------------------------

    // Find Phi candidates
    vector<float> mesonCandsFromK_mass;
    vector<float> mesonCandsFromPi_mass;
    vector<int> mesonCands_posHadronIdx;
    vector<int> mesonCands_negHadronIdx;
    // START Loop over neg hadrons -------------------------
    for (unsigned int i = 0; i < posHadronIdxs.size(); i++) {
        int posHadron_i = posHadronIdxs.at(i);
        float posHadron_pt = pfcands_p4()[posHadron_i].pt();
        // Only consider high-pt pos hadrons
        if (posHadron_pt < 10) continue;
        // START Loop over neg hadrons ---------------------
        for (unsigned int j = 0; j < negHadronIdxs.size(); j++) {
            int negHadron_i = negHadronIdxs.at(j);
            float negHadron_pt = pfcands_p4()[negHadron_i].pt();
            // Only consider high-pt neg hadrons
            if (negHadron_pt < 10) continue;
            // Get Kaon pair p4
            LorentzVector posHadron_p4 = pfcands_p4()[posHadron_i];
            LorentzVector negHadron_p4 = pfcands_p4()[negHadron_i];
            // Get Kaon pair eta, phi
            float posHadron_phi = posHadron_p4.phi();
            float negHadron_phi = negHadron_p4.phi();
            float posHadron_eta = posHadron_p4.eta();
            float negHadron_eta = negHadron_p4.eta();
            // Ignore large eta
            if (abs(posHadron_eta) > 2.4 || abs(negHadron_eta) > 2.4) continue;
            // Get hadron isolation
            float posHadron_iso = pfcands_trackIso().at(posHadron_i);
            float negHadron_iso = pfcands_trackIso().at(negHadron_i);
            bool posHadron_dzCheck = (pfcands_dz().at(posHadron_i) >= 0 && pfcands_dz().at(posHadron_i) < 0.1);
            bool negHadron_dzCheck = (pfcands_dz().at(negHadron_i) >= 0 && pfcands_dz().at(negHadron_i) < 0.1);
            posHadron_iso = (negHadron_dzCheck || pfcands_fromPV().at(negHadron_i) > 1) ? posHadron_iso - negHadron_p4.pt() : posHadron_iso;
            negHadron_iso = (posHadron_dzCheck || pfcands_fromPV().at(posHadron_i) > 1) ? negHadron_iso - posHadron_p4.pt() : negHadron_iso;
            // Get meson candidate isolation
            float meson_iso = max(posHadron_iso, negHadron_iso);
            float meson_pt = (posHadron_p4+negHadron_p4).pt();
            // Relative Isloation cut
            if (meson_iso/meson_pt > 0.06) continue;
            // Get dR between hadrons
            float posNegHadron_dR = dR(posHadron_phi, negHadron_phi, posHadron_eta, negHadron_eta);
            if (posNegHadron_dR < 0.1) {
                // Pi hypothesis
                float mesonFromPi_mass = (posHadron_p4 + negHadron_p4).M();    
                // K hypothesis
                posHadron_p4.SetE(pow(posHadron_p4.P2()+pow(.493677,2), 0.5)); // Assumed by tracker to be a
                negHadron_p4.SetE(pow(negHadron_p4.P2()+pow(.493677,2), 0.5)); // Pi, so must manually set E 
                float mesonFromK_mass = (posHadron_p4 + negHadron_p4).M();    
                // Fill meson candidate vectors
                mesonCandsFromK_mass.push_back(mesonFromK_mass);
                mesonCandsFromPi_mass.push_back(mesonFromPi_mass);
                mesonCands_posHadronIdx.push_back(posHadron_i);
                mesonCands_negHadronIdx.push_back(negHadron_i);
            }
        } // END Loop over neg hadrons --------------------
    } // END Loop over neg hadrons ------------------------

    // Find best Phi, Rho candidates
    int bestPhiCand = 0;
    int bestRhoCand = 0;
    float phiMass = 1.019;
    float rhoMass = 0.775;
    // Check if any meson candidates found
    if (mesonCands_posHadronIdx.size() > 0) {
        // START Loop over meson cands from K --------------
        for (unsigned int i = 0; i < mesonCandsFromK_mass.size(); i++) {
            float thisDiff = abs(phiMass - mesonCandsFromK_mass.at(i));
            float bestDiff = abs(phiMass - mesonCandsFromK_mass.at(bestPhiCand));
            if (thisDiff < bestDiff) {
                bestPhiCand = i;
            }
        } // END Loop over meson cands from K --------------
        // START Loop over meson cands from Pi -------------
        for (unsigned int i = 0; i < mesonCandsFromPi_mass.size(); i++) {
            float thisDiff = abs(rhoMass - mesonCandsFromPi_mass.at(i));
            float bestDiff = abs(rhoMass - mesonCandsFromPi_mass.at(bestRhoCand));
            if (thisDiff < bestDiff) {
                bestRhoCand = i;
            }
        } // END Loop over meson cands from Pi -------------
    }

    /* --> Photons <-- */
    vector<int> goodPhotons;
    // START Loop over photons -----------------------------
    for (unsigned int i = 0; i < photons_p4().size(); i++) {
        LorentzVector photon_p4 = photons_p4().at(i);
        float photon_pt = photon_p4.pt();
        // Define cuts for readability
        bool photonPtCut = (photon_pt > 20);
        bool photonEtaCut = (abs(photon_p4.eta()) < 2.5);
        bool photonTightIDCut = isMediumPhotonPOG_Fall17V2(i);
        bool photonIsoCut = (photons_recoChargedHadronIso().at(i)/photon_pt < 0.06);
        // Store 'good' photons
        if (photonPtCut && photonEtaCut && photonTightIDCut && photonIsoCut) {
            // Check for overlap
            bool isOverlap = false;
            for (unsigned int j = 0; j < els_p4().size(); j++) {
                if (els_p4().at(j).pt() > 10) {
                    LorentzVector el_p4 = els_p4().at(j);
                    isOverlap = dR(photon_p4.phi(), el_p4.phi(), photon_p4.eta(), el_p4.eta()) < 0.2;
                    if (isOverlap) break;
                }
            }
            if (!isOverlap) goodPhotons.push_back(i);
        }
    } // END Loop over photons -----------------------------

    // Find best (highest Pt) photon
    int bestPhoton = 0;
    for (unsigned int i = 0; i < goodPhotons.size(); i++) {
        float thisPhoton_pt = photons_p4().at(goodPhotons.at(i)).pt();
        float bestPhoton_pt = photons_p4().at(goodPhotons.at(bestPhoton)).pt();
        if (thisPhoton_pt > bestPhoton_pt) {
            bestPhoton = i;
        }
    }

    // Find best gen photon match
    int bestMatch = -1;
    float minPhotonParton_dR = 999.; 
    if (!evt_isRealData() && goodPhotons.size() > 0) {
        float bestMatch_dR = 0.1;
        float bestMatch_eta = 999;
        float bestMatch_phi = 999;
        for (unsigned int i = 0; i < genps_p4().size(); i++) {
            LorentzVector thisPhoton_p4 = genps_p4().at(i);
            float thisPhoton_pt = thisPhoton_p4.pt();
            float thisPhoton_eta = thisPhoton_p4.eta();
            float thisPhoton_phi = thisPhoton_p4.phi();
            LorentzVector recoPhoton_p4 = photons_p4().at(goodPhotons.at(bestPhoton));
            float recoPhoton_pt = recoPhoton_p4.pt();
            float recoPhoton_eta = recoPhoton_p4.eta();
            float recoPhoton_phi = recoPhoton_p4.phi();
            // Only consider final-state photons
            if (genps_id().at(i) != 22) continue;
            if (genps_status().at(i) != 1) continue; 
            if (fabs(genps_id_simplemother().at(i)) > 22  && genps_id_simplemother().at(i) != 2212) continue; // pions etc // but keep photons from the leading proton
            // Pre-dR cut (saves on computation time)
            if (fabs(recoPhoton_eta - thisPhoton_eta) > 0.1) continue;
            // Pt cut
            if (recoPhoton_pt > 2*thisPhoton_pt || recoPhoton_pt < 0.5*thisPhoton_pt) continue;
            // Get/compare dR
            float thisMatch_dR = dR(thisPhoton_eta, recoPhoton_eta, thisPhoton_phi, recoPhoton_phi);
            if (thisMatch_dR < bestMatch_dR) {
                bestMatch = i;
                bestMatch_dR = thisMatch_dR;
                bestMatch_eta = thisPhoton_eta;
                bestMatch_phi = thisPhoton_phi;
            }
            // Find closest parton
            if (bestMatch != -1) {
                for(unsigned int i = 0; i < genps_p4().size(); i++){
                    if (genps_status().at(i) != 22 && genps_status().at(i) != 23) continue;
                    if (fabs(genps_id().at(i)) > 21) continue;
                    LorentzVector thisParton_p4 = genps_p4().at(i);
                    float thisParton_eta = thisParton_p4.eta();
                    float thisParton_phi = thisParton_p4.phi();
                    float thisPartonPhoton_dR = dR(thisParton_eta, bestMatch_eta, thisParton_phi, bestMatch_phi);
                    if (thisPartonPhoton_dR < minPhotonParton_dR) {
                        minPhotonParton_dR = thisPartonPhoton_dR;
                    }
                }
            }
        }
    }

    /* --> Leptons <-- */
    vector<int> goodLeptonIdxs;
    vector<int> goodLeptonIDs;
    // START Loop over electrons ---------------------------
    for (unsigned int i = 0; i < els_p4().size(); i++) {
        // Define cuts for readability
        bool elsPtCut = (els_p4().at(i).pt() > 20);
        bool elsEtaCut = (abs(els_p4().at(i).eta()) < 2.4);
        bool elsIDCut = (electronID(i,id_level_t::HAD_medium_noiso_v5));
        bool elsIsoCut = (elMiniRelIsoCMS3_EA(i, gconf.ea_version) < 0.1);
        // Store 'good' electrons
        if (elsPtCut && elsEtaCut && elsIDCut && elsIsoCut) {
            goodLeptonIdxs.push_back(i);
            goodLeptonIDs.push_back(-11*els_charge().at(i));
        }
    } // END Loop over electrons ---------------------------

    // START Loop over muons -------------------------------
    for (unsigned int i = 0; i < mus_p4().size(); i++) {
        // Define cuts for readability
        bool musPtCut = (mus_p4().at(i).pt() > 20);
        bool musEtaCut = (abs(mus_p4().at(i).eta()) < 2.4);
        bool musIDCut = (isMediumMuonPOG(i));
        bool musIsoCut = (muMiniRelIsoCMS3_EA(i, gconf.ea_version) < 0.2);
        // Store 'good' muons
        if (musPtCut && musEtaCut && musIDCut && musIsoCut) {
            goodLeptonIdxs.push_back(i); 
            goodLeptonIDs.push_back(-13*mus_charge().at(i)); 
        }
    } // END Loop over muons -------------------------------

    // Find best (highest Pt) lepton, should usually only be one to choose from
    int bestLepton = 0;
    for (unsigned int i = 0; i < goodLeptonIdxs.size(); i++) {
        int tl = goodLeptonIdxs.at(i);
        int bl = goodLeptonIdxs.at(bestLepton);
        LorentzVector thisLepton_p4 = (abs(goodLeptonIDs.at(i)) == 11) ? els_p4().at(tl) : mus_p4().at(tl);
        float thisLepton_pt = thisLepton_p4.pt();
        LorentzVector bestLepton_p4 = (abs(goodLeptonIDs.at(bestLepton)) == 11) ? els_p4().at(bl) : mus_p4().at(bl);
        float bestLepton_pt = bestLepton_p4.pt();
        if (thisLepton_pt > bestLepton_pt) {
            bestLepton = i;
        }
    }

    /* --> Fill Photons <-- */
    LorentzVector bestPhoton_p4; // Declare here, since used in following conditionals
    if (goodPhotons.size() > 0) {
        // Best photon
        int bestPhoton_i = goodPhotons.at(bestPhoton);
        bestPhoton_p4 = photons_p4().at(bestPhoton_i);
        // Best photon kinematics
        recoGamma_pt = bestPhoton_p4.pt();
        recoGamma_eta = bestPhoton_p4.eta();
        recoGamma_phi = bestPhoton_p4.phi();
        // Best photon relative isolation
        recoGamma_relIso = (photons_recoChargedHadronIso().at(bestPhoton_i))/(recoGamma_pt);
        // Best photon gen-reco match
        if (!evt_isRealData()) {
            genRecoGamma_isMatch = (bestMatch != -1) ? 1 : 0;
        }
        // Minimum gen photon-parton dR
        if (minPhotonParton_dR != 999.) {
            minGammaParton_dR = minPhotonParton_dR;
        }
    }

    /* --> Fill Leptons <-- */
    LorentzVector bestLepton_p4; // Declare here, since used in following conditionals
    if (goodLeptonIdxs.size() > 0) {
        // Number of 'good' leptons
        recoWLepton_nLep = goodLeptonIdxs.size();
        // Best lepton
        recoWLepton_id = goodLeptonIDs.at(bestLepton);
        int bestLep_i = goodLeptonIdxs.at(bestLepton);
        // Best lepton kinematics
        bestLepton_p4 = (abs(recoWLepton_id) == 11) ? els_p4().at(bestLep_i) : mus_p4().at(bestLep_i);
        recoWLepton_pt = bestLepton_p4.pt(); 
        recoWLepton_eta = bestLepton_p4.eta(); 
        recoWLepton_phi = bestLepton_p4.phi(); 
        // dR(photon, lepton)
        recoGammaWLepton_dR = dR(recoGamma_phi, recoWLepton_phi, recoGamma_eta, recoWLepton_eta);
    }

    /* --> Fill Best Phi Candidate <-- */
    LorentzVector bestPhi_p4; // Declare here, since used in following conditionals
    if (mesonCandsFromK_mass.size() > 0) {
        // Best K+, K-
        int bestKp_i = mesonCands_posHadronIdx.at(bestPhiCand);
        int bestKm_i = mesonCands_negHadronIdx.at(bestPhiCand);
        LorentzVector bestKp_p4 = pfcands_p4().at(bestKp_i);
        LorentzVector bestKm_p4 = pfcands_p4().at(bestKm_i);
        bestKp_p4.SetE(pow(bestKp_p4.P2()+pow(.493677,2), 0.5)); // Assumed by tracker to be a
        bestKm_p4.SetE(pow(bestKm_p4.P2()+pow(.493677,2), 0.5)); // Pi, so must manually set E 
        // Best Phi
        bestPhi_p4 = bestKp_p4+bestKm_p4;
        // Best Phi cand mass
        recoPhi_mass = mesonCandsFromK_mass.at(bestPhiCand);
        recoPhi_pt = bestPhi_p4.pt();
        recoPhi_eta = bestPhi_p4.eta();
        recoPhi_phi = bestPhi_p4.phi();
        // Best K+, K- kinematics
        recoKp_pt = bestKp_p4.pt();
        recoKp_phi = bestKp_p4.phi();
        recoKp_eta = bestKp_p4.eta();
        recoKm_pt = bestKm_p4.pt();
        recoKm_phi = bestKm_p4.phi();
        recoKm_eta = bestKm_p4.eta();
        recoKpKm_dR = dR(recoKp_phi, recoKm_phi, recoKp_eta, recoKm_eta);
        // Best Phi, K+, K- relative isolation
        float bestKp_iso = pfcands_trackIso().at(bestKp_i);
        float bestKm_iso = pfcands_trackIso().at(bestKm_i);
        bestKp_iso = ((pfcands_dz().at(bestKm_i) >= 0 && pfcands_dz().at(bestKm_i) < 0.1) || pfcands_fromPV().at(bestKm_i) > 1) ? bestKp_iso - bestKm_p4.pt() : bestKp_iso;
        bestKm_iso = ((pfcands_dz().at(bestKp_i) >= 0 && pfcands_dz().at(bestKp_i) < 0.1) || pfcands_fromPV().at(bestKp_i) > 1) ? bestKm_iso - bestKp_p4.pt() : bestKm_iso;
        recoPhi_relIso = max(bestKp_iso, bestKm_iso)/recoPhi_pt;
        recoKp_relIso = bestKp_iso/recoKp_pt;
        recoKm_relIso = bestKm_iso/recoKm_pt;
        // Best Higgs cand mass
        if (goodPhotons.size() > 0) {
            recoHiggs_mass = (bestKp_p4+bestKm_p4+bestPhoton_p4).M();
            recoPhiGamma_dR = dR(recoPhi_phi, recoGamma_phi, recoPhi_eta, recoGamma_eta);
        }
    }

    /* --> Fill Best Rho Candidate <-- */
    LorentzVector bestRho_p4; // Declare here, since used in following conditionals
    if (mesonCandsFromPi_mass.size() > 0) {
        // Best Pi+, Pi-
        int bestPip_i = mesonCands_posHadronIdx.at(bestRhoCand);
        int bestPim_i = mesonCands_negHadronIdx.at(bestRhoCand);
        LorentzVector bestPip_p4 = pfcands_p4().at(bestPip_i);
        LorentzVector bestPim_p4 = pfcands_p4().at(bestPim_i);
        // Best Rho
        bestRho_p4 = bestPip_p4+bestPim_p4;
        // Best Rho cand mass
        recoRho_mass = mesonCandsFromPi_mass.at(bestRhoCand);
        recoRho_pt = bestRho_p4.pt();
        recoRho_eta = bestRho_p4.eta();
        recoRho_phi = bestRho_p4.phi();
        // Best Pi+, Pi- kinematics
        recoPip_pt = bestPip_p4.pt();
        recoPim_phi = bestPip_p4.phi();
        recoPip_eta = bestPip_p4.eta();
        recoPim_pt = bestPim_p4.pt();
        recoPip_phi = bestPim_p4.phi();
        recoPim_eta = bestPim_p4.eta();
        recoPipPim_dR = dR(recoPip_phi, recoPim_phi, recoPip_eta, recoPim_eta);
        // Best Rho, Pi+, Pi- relative isolation
        float bestPip_iso = pfcands_trackIso().at(bestPip_i);
        float bestPim_iso = pfcands_trackIso().at(bestPim_i);
        bestPip_iso = ((pfcands_dz().at(bestPim_i) >= 0 && pfcands_dz().at(bestPim_i) < 0.1) || pfcands_fromPV().at(bestPim_i) > 1) ? bestPip_iso - bestPim_p4.pt() : bestPip_iso;
        bestPim_iso = ((pfcands_dz().at(bestPip_i) >= 0 && pfcands_dz().at(bestPip_i) < 0.1) || pfcands_fromPV().at(bestPip_i) > 1) ? bestPim_iso - bestPip_p4.pt() : bestPim_iso;
        recoRho_relIso = max(bestPip_iso, bestPim_iso)/recoRho_pt;
        recoPip_relIso = bestPip_iso/recoPip_pt;
        recoPim_relIso = bestPim_iso/recoPim_pt;
        // Best Higgs cand mass
        if (goodPhotons.size() > 0) {
            recoHiggs_mass = (bestPip_p4+bestPim_p4+bestPhoton_p4).M();
            recoRhoGamma_dR = dR(recoRho_phi, recoGamma_phi, recoRho_eta, recoGamma_eta);
        }
    }

    /* --> Fill Magic Angles <-- */
    if (goodPhotons.size() > 0 && goodLeptonIdxs.size() > 0 && mesonCandsFromK_mass.size() > 0) {
        // Output from MELA
        MagicAngles mAngles;
        // Variables for MELA
        TLorentzVector lep_p4(bestLepton_p4.Px(),bestLepton_p4.Py(),bestLepton_p4.Pz(),bestLepton_p4.E());
        TLorentzVector mes_p4(bestPhi_p4.Px(),bestPhi_p4.Py(),bestPhi_p4.Pz(),bestPhi_p4.E());
        TLorentzVector gam_p4(bestPhoton_p4.Px(),bestPhoton_p4.Py(),bestPhoton_p4.Pz(),bestPhoton_p4.E());
        // Get magic angles
        mAngles = getAngles(met_pt, met_phi, recoWLepton_id, lep_p4, mes_p4, gam_p4);
        // Fill branches
        recoMagAng_cosThetaStar = mAngles.angles[0];
        recoMagAng_cosTheta1 = mAngles.angles[1];
        recoMagAng_cosTheta2 = mAngles.angles[2];
        recoMagAng_Phi = mAngles.angles[3];
        recoMagAng_Phi1 = mAngles.angles[4];
        recoMagAng_m1 = mAngles.angles[5];
        recoMagAng_m2 = mAngles.angles[6];
    }

    /* --> Fill nCands <-- */
    recoMeson_nCands = mesonCandsFromK_mass.size(); // number of cands from K, pi are the same

    /* --> Fill Triggers <-- */
    // Single muon
    HLT_singleMu = passHLTTriggerPattern("HLT_IsoMu17_eta2p1_v") ||
                   passHLTTriggerPattern("HLT_IsoMu20_v") || passHLTTriggerPattern("HLT_IsoMu20_eta2p1_v") ||
                   passHLTTriggerPattern("HLT_IsoTkMu20_v") || passHLTTriggerPattern("HLT_IsoTkMu20_eta2p1_v") ||
                   passHLTTriggerPattern("HLT_IsoMu24_v") || passHLTTriggerPattern("HLT_IsoTkMu24_v") || passHLTTriggerPattern("HLT_IsoMu24_eta2p1_v") ||
                   passHLTTriggerPattern("HLT_IsoMu27_v") || passHLTTriggerPattern("HLT_IsoTkMu27_v");
    // Single electron
    HLT_singleEl = passHLTTriggerPattern("HLT_Ele27_eta2p1_WPTight_Gsf_v") || // 2016
                   passHLTTriggerPattern("HLT_Ele32_eta2p1_WPTight_Gsf_v") || // 2016
                   passHLTTriggerPattern("HLT_Ele27_WPTight_Gsf_v") ||
                   passHLTTriggerPattern("HLT_Ele32_WPTight_Gsf_v") ||  // 2017,18
                   passHLTTriggerPattern("HLT_Ele35_WPTight_Gsf_v") ||  // 2017,18
                   passHLTTriggerPattern("HLT_Ele38_WPTight_Gsf_v") ||  // 2017,18
                   passHLTTriggerPattern("HLT_Ele28_eta2p1_WPTight_Gsf_HT150_v"); // 2016

    /* --> Fill Filters <-- */
    // Define variable flags
    bool Flag_ecalBadCalib = (config.year == 2016) ? true : filt_ecalBadCalibFilterUpdate(); // 2017, 2018
    bool Flag_eeBadSc = (config.isData) ? filt_eeBadSc() : true; // Data only
    // Check flags
    passFilters = filt_goodVertices() && filt_globalSuperTightHalo2016() && filt_hbheNoise() && 
                  filt_hbheNoiseIso() && filt_ecalTP() && filt_BadPFMuonFilter() &&
                  Flag_ecalBadCalib && Flag_eeBadSc;

    /* --> Fill HEM <-- */
    isHEM = 0;
    if (config.year == 2018) {
        if ((config.isData && run >= 319077) || (!config.isData && event % 1961 < 1286)) {
            bool photon_inRegion = (goodPhotons.size() > 0 && (recoGamma_eta >= -4.7 && recoGamma_eta <= -1.4) && (recoGamma_phi >= -1.6 && recoGamma_phi <= -0.8));
            bool electron_inRegion = (abs(recoWLepton_id) == abs(11) && (recoWLepton_eta >= -4.7 && recoWLepton_eta <= -1.4) && (recoWLepton_phi >= -1.6 && recoWLepton_phi <= -0.8));
            if (photon_inRegion || electron_inRegion) isHEM = 1;
        }
    }

    return;
}

/**
 * Fill relevant branches with gen-reco information
 * @params: none
 * @return: none
 */
void BabyTree::FillGenRecoBranches() {
    
    if (genGamma_pt != -999 && recoGamma_pt != -999) {
        genRecoGamma_dR = dR(genGamma_phi, recoGamma_phi, genGamma_eta, recoGamma_eta);
    }
    if (genHiggsMeson_id == 333 && recoPhi_mass != -999) {
        genRecoPhi_dR = dR(genHiggsMeson_phi, recoPhi_phi, genHiggsMeson_eta, recoPhi_eta);
    }
    // else if (genHiggsMeson_id == 113 && recoHiggsMeson_id == 113) {
    //     genRecoHiggsRho_dR = dR(genHiggsMeson_phi, recoHiggsMeson_phi, genHiggsMeson_eta, recoHiggsMeson_eta);
    // }

    return;
}
