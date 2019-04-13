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

// CMS3
#include "/home/users/jguiang/projects/mt2/MT2Analysis/CORE/CMS3.h"
#include "/home/users/jguiang/projects/mt2/MT2Analysis/CORE/ElectronSelections.h"
#include "/home/users/jguiang/projects/mt2/MT2Analysis/CORE/MuonSelections.h"
#include "/home/users/jguiang/projects/mt2/MT2Analysis/CORE/PhotonSelections.h"
#include "/home/users/jguiang/projects/mt2/MT2Analysis/CORE/IsolationTools.h"
#include "/home/users/jguiang/projects/mt2/MT2Analysis/CORE/TriggerSelections.h"

// Header
#include "mcTree.h"

using namespace std;
using namespace tas;

mcTree::mcTree() {
    /* --> TTree Setup <-- */
    t = new TTree("tree", "tree");
    /* --> Meta Branches Setup <-- */
    // Event info
    b_run = t->Branch("run", &run, "run/I");
    b_lumi = t->Branch("lumi", &lumi, "lumi/I");
    b_event = t->Branch("event", &event, "event/I");
    // Gen-Reco dR
    b_genRecoGamma_dR = t->Branch("genRecoGamma_dR", &genRecoGamma_dR, "genRecoGamma_dR/F");
    b_genRecoPhi_dR = t->Branch("genRecoPhi_dR", &genRecoPhi_dR, "genRecoPhi_dR/F");
    b_genRecoRho_dR = t->Branch("genRecoRho_dR", &genRecoRho_dR, "genRecoRho_dR/F");
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
    b_recoPhi_iso = t->Branch("recoPhi_iso", &recoPhi_iso, "recoPhi_iso/F");
    // Reco K+, K-
    b_recoKm_pt = t->Branch("recoKm_pt", &recoKm_pt, "recoKm_pt/F");
    b_recoKm_eta = t->Branch("recoKm_eta", &recoKm_eta, "recoKm_eta/F");
    b_recoKm_phi = t->Branch("recoKm_phi", &recoKm_phi, "recoKm_phi/F");
    b_recoKm_iso = t->Branch("recoKm_iso", &recoKm_iso, "recoKm_iso/F");
    b_recoKp_pt = t->Branch("recoKp_pt", &recoKp_pt, "recoKp_pt/F");
    b_recoKp_eta = t->Branch("recoKp_eta", &recoKp_eta, "recoKp_eta/F");
    b_recoKp_phi = t->Branch("recoKp_phi", &recoKp_phi, "recoKp_phi/F");
    b_recoKp_iso = t->Branch("recoKp_iso", &recoKp_iso, "recoKp_iso/F");
    b_recoKpKm_dR = t->Branch("recoKpKm_dR", &recoKpKm_dR, "recoKpKm_dR/F");
    // Reco Rho
    b_recoRho_mass = t->Branch("recoRho_mass", &recoRho_mass, "recoRho_mass/F");
    b_recoRho_pt = t->Branch("recoRho_pt", &recoRho_pt, "recoRho_pt/F");
    b_recoRho_eta = t->Branch("recoRho_eta", &recoRho_eta, "recoRho_eta/F");
    b_recoRho_phi = t->Branch("recoRho_phi", &recoRho_phi, "recoRho_phi/F");
    b_recoRho_iso = t->Branch("recoRho_iso", &recoRho_iso, "recoRho_iso/F");
    // Reco Pi+, Pi-
    b_recoPim_pt = t->Branch("recoPim_pt", &recoPim_pt, "recoPim_pt/F");
    b_recoPim_eta = t->Branch("recoPim_eta", &recoPim_eta, "recoPim_eta/F");
    b_recoPim_phi = t->Branch("recoPim_phi", &recoPim_phi, "recoPim_phi/F");
    b_recoPim_iso = t->Branch("recoPim_iso", &recoPim_iso, "recoPim_iso/F");
    b_recoPip_pt = t->Branch("recoPip_pt", &recoPip_pt, "recoPip_pt/F");
    b_recoPip_eta = t->Branch("recoPip_eta", &recoPip_eta, "recoPip_eta/F");
    b_recoPip_phi = t->Branch("recoPip_phi", &recoPip_phi, "recoPip_phi/F");
    b_recoPip_iso = t->Branch("recoPip_iso", &recoPip_iso, "recoPip_iso/F");
    b_recoPipPim_dR = t->Branch("recoPipPim_dR", &recoPipPim_dR, "recoPipPim_dR/F");
    // Reco Photons
    b_recoGamma_pt = t->Branch("recoGamma_pt", &recoGamma_pt, "recoGamma_pt/F");
    b_recoGamma_phi = t->Branch("recoGamma_phi", &recoGamma_phi, "recoGamma_phi/F");
    b_recoGamma_eta = t->Branch("recoGamma_eta", &recoGamma_eta, "recoGamma_eta/F");
    b_recoGamma_iso = t->Branch("recoGamma_iso", &recoGamma_iso, "recoGamma_iso/F");
    // Reco Leptons
    b_recoWLepton_id = t->Branch("recoWLepton_id", &recoWLepton_id, "recoWLepton_id/I");
    b_recoWLepton_pt = t->Branch("recoWLepton_pt", &recoWLepton_pt, "recoWLepton_pt/F");
    b_recoWLepton_eta = t->Branch("recoWLepton_eta", &recoWLepton_eta, "recoWLepton_eta/F");
    b_recoWLepton_phi = t->Branch("recoWLepton_phi", &recoWLepton_phi, "recoWLepton_phi/F");
    b_recoWLepton_nLep = t->Branch("recoWLepton_nLep", &recoWLepton_nLep, "recoWLepton_nLep/I");
}

void mcTree::Reset() {
    // Meta
    run = -999;
    lumi = -999;
    event = -999;
    genRecoGamma_dR = -999;
    genRecoPhi_dR = -999;
    genRecoRho_dR = -999;
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
    recoPhi_iso = -999;
    recoKm_pt = -999;
    recoKm_eta = -999;
    recoKm_phi = -999;
    recoKm_iso = -999;
    recoKp_pt = -999;
    recoKp_eta = -999;
    recoKp_phi = -999;
    recoKp_iso = -999;
    recoKpKm_dR = -999;
    recoRho_mass = -999;
    recoRho_pt = -999;
    recoRho_eta = -999;
    recoRho_phi = -999;
    recoRho_iso = -999;
    recoPim_pt = -999; 
    recoPim_eta = -999;
    recoPim_phi = -999;
    recoPim_iso = -999;
    recoPip_pt = -999;
    recoPip_eta = -999;
    recoPip_phi = -999;
    recoPip_iso = -999;
    recoPipPim_dR = -999;
    recoGamma_pt = -999;
    recoGamma_phi = -999;
    recoGamma_eta = -999;
    recoGamma_iso = -999;
    recoWLepton_id = -999;
    recoWLepton_pt = -999;
    recoWLepton_eta = -999;
    recoWLepton_phi = -999;
    recoWLepton_nLep = -999;

    return;
}

float mcTree::dR(float phi1, float phi2, float eta1, float eta2) {

    float dphi = abs(phi2 - phi1);
    if (dphi > M_PI){
        dphi = 2*M_PI - dphi;
    }
    return sqrt(pow((dphi), 2) + pow((eta2 - eta1), 2));
}

void mcTree::FillGenBranches() {

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
            genHiggs_pt = genps_p4()[i].pt();
            genHiggs_eta = genps_p4()[i].eta();
            genHiggs_phi = genps_p4()[i].phi(); 
        }
        if (thisID == 24 || thisID == -24) {
            // Fill W boson branches
            genW_pt = genps_p4()[i].pt();
            genW_eta = genps_p4()[i].eta();
            genW_phi = genps_p4()[i].phi(); 
        }
        // Add particle id, idx, p4 to relevant systems
        if (thisMotherID == 24 || thisMotherID == -24) {
            if (thisID == -11 || thisID == 11) {
                mothers[W_to_ElNu] = thisMotherIdx;
                daughters[W_to_ElNu].push_back(i);
            }
            if (thisID == -13 || thisID == 13) {
                mothers[W_to_MuNu] = thisMotherIdx;
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
            if ((thisID == 321 || thisID == -321) && genps_id_mother()[thisMotherIdx] == 25) {
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

    // Retrieve products filled in gen-level loop, fill tree branches
    if (daughters[W_to_MuNu].size() == 1 || daughters[W_to_ElNu].size() == 1) {
        decay = (daughters[W_to_MuNu].size() == 1) ? W_to_MuNu : W_to_ElNu;
        p4_sum -= p4_sum; // Reset p4 sum
        for (unsigned int j = 0; j < daughters[decay].size(); j++) {
            // Get daughter index
            int idx = daughters[decay].at(j);
            // Get ID
            int id = genps_id().at(idx);
            // Skip neutrinos
            if (id == 14 || id == 12) {
                continue;
            }
            else {
                // Fill leptons from W branches
                genWLepton_id = id;
                genWLepton_pt = genps_p4()[idx].pt();
                genWLepton_eta = genps_p4()[idx].eta();
                genWLepton_phi = genps_p4()[idx].phi();
                p4_sum += genps_p4()[idx];
            }
        }
        genW_mass = p4_sum.M();
    }

    if (daughters[H_to_PhiGamma].size() == 2 || daughters[H_to_RhoGamma].size() == 2) {
        decay = (daughters[H_to_PhiGamma].size() == 2) ? H_to_PhiGamma : H_to_RhoGamma;
        p4_sum -= p4_sum; // Reset p4 sum
        for (unsigned int j = 0; j < daughters[decay].size(); j++) {
            // Get daughter index
            int idx = daughters[decay].at(j);
            // Get ID
            int id = genps_id().at(idx);
            // Fill photon branches
            if (id == 22) {
                genGamma_pt = genps_p4()[idx].pt();
                genGamma_eta = genps_p4()[idx].eta();
                genGamma_phi = genps_p4()[idx].phi();
            }
            else {
                // Fill mesons from Higgs branches
                genHiggsMeson_id = id;
                genHiggsMeson_pt = genps_p4()[idx].pt();
                genHiggsMeson_eta = genps_p4()[idx].eta();
                genHiggsMeson_phi = genps_p4()[idx].phi();
                p4_sum += genps_p4()[idx];
            }
        }
        genHiggsMesonGamma_dR = dR(genHiggsMeson_phi, genGamma_phi, genHiggsMeson_eta, genGamma_eta);
        genHiggs_mass = p4_sum.M();
    }

    if (daughters[Phi_to_KK].size() == 2) {
        decay = Phi_to_KK;
        p4_sum -= p4_sum; // Reset p4 sum
        for (unsigned int j = 0; j < daughters[decay].size(); j++) {
            // Get daughter index
            int idx = daughters[decay].at(j);
            // Get ID
            int id = genps_id().at(idx);
            // Fill Kaon from Phi branches
            if (id == 321) {
                genKp_pt = genps_p4()[idx].pt();
                genKp_phi = genps_p4()[idx].phi();
                genKp_eta = genps_p4()[idx].eta();
            }
            else if (id == -321) {
                genKm_pt = genps_p4()[idx].pt();
                genKm_phi = genps_p4()[idx].phi();
                genKm_eta = genps_p4()[idx].eta();
            }
            p4_sum += genps_p4()[idx];
        }
        genKpKm_dR = dR(genKp_phi, genKm_phi, genKp_eta, genKm_eta);
        genHiggsMeson_mass = p4_sum.M();
    }

    return;
}

void mcTree::FillRecoBranches() {

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

    /* --> Leptons <-- */
    vector<int> goodLeptonIdxs;
    vector<int> goodLeptonIDs;
    // START Loop over electrons ---------------------------
    for (unsigned int i = 0; i < els_p4().size(); i++) {
        // Define cuts for readability
        bool elsPtCut = (els_p4().at(i).pt() > 20);
        bool elsEtaCut = (els_p4().at(i).eta() < 2.4);
        bool elsIDCut = (isMediumElectronPOGfall17noIso_v2(i));
        bool elsIsoCut = (elMiniRelIsoCMS3_EA(i, 3) < 0.1);
        // Store 'good' electrons
        if (elsPtCut && elsEtaCut && elsIDCut && elsIsoCut) {
            goodLeptonIdxs.push_back(i);
            goodLeptonIDs.push_back(11);
        }
    } // END Loop over electrons ---------------------------

    // START Loop over muons -------------------------------
    for (unsigned int i = 0; i < mus_p4().size(); i++) {
        // Define cuts for readability
        bool musPtCut = (mus_p4().at(i).pt() > 20);
        bool musEtaCut = (mus_p4().at(i).eta() < 2.4);
        bool musIDCut = (isMediumMuonPOG(i));
        bool musIsoCut = (muMiniRelIsoCMS3_EA(i, 3) < 0.2);
        // Store 'good' muons
        if (musPtCut && musEtaCut && musIDCut && musIsoCut) {
            goodLeptonIdxs.push_back(i); 
            goodLeptonIDs.push_back(13); 
        }
    } // END Loop over muons -------------------------------

    // Find best (highest Pt) lepton, should usually only be one to choose from
    int bestLepton = 0;
    if (goodLeptonIdxs.size() > 1) {
        for (unsigned int i = 0; i < goodLeptonIdxs.size(); i++) {
            int tl = goodLeptonIdxs.at(i);
            int bl = goodLeptonIdxs.at(bestLepton);
            LorentzVector thisLepton_p4 = (goodLeptonIDs.at(i) == 11) ? els_p4().at(tl) : mus_p4().at(tl);
            LorentzVector bestLepton_p4 = (goodLeptonIDs.at(bestLepton) == 11) ? els_p4().at(bl) : mus_p4().at(bl);
            float thisLepton_pt = thisLepton_p4.pt();
            float bestLepton_pt = bestLepton_p4.pt();
            if (thisLepton_pt > bestLepton_pt) {
                bestLepton = i;
            }
        }
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
        // bool photonLooseIDCut  = isLoosePhoton(i, analysis_t::HAD, 3);
        bool photonTightIDCut  = isTightPhoton(i, analysis_t::HAD, 3);
        bool photonIsoCut = (photons_recoChargedHadronIso().at(i)/photon_pt < 0.06);
        // Store 'good' photons
        if (photonPtCut && photonEtaCut && photonTightIDCut && photonIsoCut) goodPhotons.push_back(i);
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

    // Retrieve products filled in loops, fill tree branches

    /* --> Photons <-- */
    LorentzVector bestPhoton_p4; // Declare here, since used in following conditionals
    if (goodPhotons.size() > 0) {
        // Best Photon
        int bestPhoton_i = goodPhotons.at(bestPhoton);
        bestPhoton_p4 = photons_p4().at(bestPhoton_i);
        // Best Photon isolation
        recoGamma_iso = photons_recoChargedHadronIso().at(bestPhoton_i);
        // Best Photon kinematics
        recoGamma_pt = bestPhoton_p4.pt();
        recoGamma_eta = bestPhoton_p4.eta();
        recoGamma_phi = bestPhoton_p4.phi();
    }

    /* --> Leptons <-- */
    if (goodLeptonIdxs.size() > 0) {
        // Number of 'good' leptons
        recoWLepton_nLep = goodLeptonIdxs.size();
        // Best lepton
        recoWLepton_id = goodLeptonIDs.at(bestLepton);
        int bestLep_i = goodLeptonIdxs.at(bestLepton);
        // Best lepton kinematics
        LorentzVector bestLepton_p4 = (recoWLepton_id == 11) ? els_p4().at(bestLep_i) : mus_p4().at(bestLep_i);
        recoWLepton_pt = bestLepton_p4.pt(); 
        recoWLepton_eta = bestLepton_p4.eta(); 
        recoWLepton_phi = bestLepton_p4.phi(); 
    }

    /* --> Best Phi Candidate <-- */
    if (mesonCandsFromK_mass.size() > 0) {
        // Best K+, K-
        int bestKp_i = mesonCands_posHadronIdx.at(bestPhiCand);
        int bestKm_i = mesonCands_negHadronIdx.at(bestPhiCand);
        LorentzVector bestKp_p4 = pfcands_p4().at(bestKp_i);
        LorentzVector bestKm_p4 = pfcands_p4().at(bestKm_i);
        bestKp_p4.SetE(pow(bestKp_p4.P2()+pow(.493677,2), 0.5)); // Assumed by tracker to be a
        bestKm_p4.SetE(pow(bestKm_p4.P2()+pow(.493677,2), 0.5)); // Pi, so must manually set E 
        // Best Phi
        LorentzVector bestPhi_p4 = bestKp_p4+bestKm_p4;
        // Best Phi, K+, K- isolation
        float bestKp_iso = pfcands_trackIso().at(bestKp_i);
        float bestKm_iso = pfcands_trackIso().at(bestKm_i);
        bestKp_iso = ((pfcands_dz().at(bestKm_i) >= 0 && pfcands_dz().at(bestKm_i) < 0.1) || pfcands_fromPV().at(bestKm_i) > 1) ? bestKp_iso - bestKm_p4.pt() : bestKp_iso;
        bestKm_iso = ((pfcands_dz().at(bestKp_i) >= 0 && pfcands_dz().at(bestKp_i) < 0.1) || pfcands_fromPV().at(bestKp_i) > 1) ? bestKm_iso - bestKp_p4.pt() : bestKm_iso;
        recoPhi_iso = max(bestKp_iso, bestKm_iso);
        recoKp_iso = bestKp_iso;
        recoKm_iso = bestKm_iso;
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
        // Best Higgs cand mass
        if (goodPhotons.size() > 0) {
            recoHiggs_mass = (bestKp_p4+bestKm_p4+bestPhoton_p4).M();
        }
    }

    /* --> Best Rho Candidate <-- */
    if (mesonCandsFromPi_mass.size() > 0) {
        // Best Pi+, Pi-
        int bestPip_i = mesonCands_posHadronIdx.at(bestRhoCand);
        int bestPim_i = mesonCands_negHadronIdx.at(bestRhoCand);
        LorentzVector bestPip_p4 = pfcands_p4().at(bestPip_i);
        LorentzVector bestPim_p4 = pfcands_p4().at(bestPim_i);
        // Best Rho
        LorentzVector bestRho_p4 = bestPip_p4+bestPim_p4;
        // Best Rho, Pi+, Pi- isolation
        float bestPip_iso = pfcands_trackIso().at(bestPip_i);
        float bestPim_iso = pfcands_trackIso().at(bestPim_i);
        bestPip_iso = ((pfcands_dz().at(bestPim_i) >= 0 && pfcands_dz().at(bestPim_i) < 0.1) || pfcands_fromPV().at(bestPim_i) > 1) ? bestPip_iso - bestPim_p4.pt() : bestPip_iso;
        bestPim_iso = ((pfcands_dz().at(bestPip_i) >= 0 && pfcands_dz().at(bestPip_i) < 0.1) || pfcands_fromPV().at(bestPip_i) > 1) ? bestPim_iso - bestPip_p4.pt() : bestPim_iso;
        recoRho_iso = max(bestPip_iso, bestPim_iso);
        recoPip_iso = bestPip_iso;
        recoPim_iso = bestPim_iso;
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
        // Best Higgs cand mass
        if (goodPhotons.size() > 0) {
            recoHiggs_mass = (bestPip_p4+bestPim_p4+bestPhoton_p4).M();
        }
    }

    // Save number of meson candidates
    recoMeson_nCands = mesonCandsFromK_mass.size(); // number of cands from K, pi are the same

    return;
}

void mcTree::FillGenRecoBranches() {
    
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
