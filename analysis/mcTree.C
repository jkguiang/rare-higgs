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
#include "/home/users/jguiang/projects/mt2/MT2Analysis/CORE/IsolationTools.h"
#include "/home/users/jguiang/projects/mt2/MT2Analysis/CORE/TriggerSelections.h"

// Header
#include "mcTree.h"

using namespace std;
using namespace tas;

mcTree::mcTree() {
    /* --> TTree Setup <-- */
    t = new TTree("tree", "tree");
    /* --> Gen Branches Setup <-- */
    b_genW_pt = t->Branch("genW_pt", &genW_pt, "genW_pt/F");
    b_genW_eta = t->Branch("genW_eta", &genW_eta, "genW_eta/F");
    b_genW_phi = t->Branch("genW_phi", &genW_phi, "genW_phi/F");
    b_genW_mass = t->Branch("genW_mass", &genW_mass, "genW_mass/F");
    b_genWLepton_id = t->Branch("genWLepton_id", &genWLepton_id, "genWLepton_id/F");
    b_genWLepton_pt = t->Branch("genWLepton_pt", &genWLepton_pt, "genWLepton_pt/F");
    b_genWLepton_eta = t->Branch("genWLepton_eta", &genWLepton_eta, "genWLepton_eta/F");
    b_genWLepton_phi = t->Branch("genWLepton_phi", &genWLepton_phi, "genWLepton_phi/F");
    b_genHiggs_pt = t->Branch("genHiggs_pt", &genHiggs_pt, "genHiggs_pt/F");
    b_genHiggs_eta = t->Branch("genHiggs_eta", &genHiggs_eta, "genHiggs_eta/F");
    b_genHiggs_phi = t->Branch("genHiggs_phi", &genHiggs_phi, "genHiggs_phi/F");
    b_genHiggs_mass = t->Branch("genHiggs_mass", &genHiggs_mass, "genHiggs_mass/F");
    b_genHiggsMeson_id = t->Branch("genHiggsMeson_id", &genHiggsMeson_id, "genHiggsMeson_id/F");
    b_genHiggsMeson_pt = t->Branch("genHiggsMeson_pt", &genHiggsMeson_pt, "genHiggsMeson_pt/F");
    b_genHiggsMeson_eta = t->Branch("genHiggsMeson_eta", &genHiggsMeson_eta, "genHiggsMeson_eta/F");
    b_genHiggsMeson_phi = t->Branch("genHiggsMeson_phi", &genHiggsMeson_phi, "genHiggsMeson_phi/F");
    b_genHiggsMeson_mass = t->Branch("genHiggsMeson_mass", &genHiggsMeson_mass, "genHiggsMeson_mass/F");
    b_genHiggsMesonGamma_dR = t->Branch("genHiggsMesonGamma_dR", &genHiggsMesonGamma_dR, "genHiggsMesonGamma_dR/F");
    b_genGamma_pt = t->Branch("genGamma_pt", &genGamma_pt, "genGamma_pt/F");
    b_genGamma_phi = t->Branch("genGamma_phi", &genGamma_phi, "genGamma_phi/F");
    b_genGamma_eta = t->Branch("genGamma_eta", &genGamma_eta, "genGamma_eta/F");
    b_genKm_pt = t->Branch("genKm_pt", &genKm_pt, "genKm_pt/F");
    b_genKm_phi = t->Branch("genKm_phi", &genKm_phi, "genKm_phi/F");
    b_genKm_eta = t->Branch("genKm_eta", &genKm_eta, "genKm_eta/F");
    b_genKp_pt = t->Branch("genKp_pt", &genKp_pt, "genKp_pt/F");
    b_genKp_phi = t->Branch("genKp_phi", &genKp_phi, "genKp_phi/F");
    b_genKp_eta = t->Branch("genKp_eta", &genKp_eta, "genKp_eta/F");
    b_genKpKm_dR = t->Branch("genKpKm_dR", &genKpKm_dR, "genKpKm_dR/F");
    /* --> Reco Branches Setup <-- */
    b_recoPhi_mass = t->Branch("recoPhi_mass", &recoPhi_mass, "recoPhi_mass/F");
    b_recoKm_pt = t->Branch("recoKm_pt", &recoKm_pt, "recoKm_pt/F");
    b_recoKm_phi = t->Branch("recoKm_phi", &recoKm_phi, "recoKm_phi/F");
    b_recoKm_eta = t->Branch("recoKm_eta", &recoKm_eta, "recoKm_eta/F");
    b_recoKp_pt = t->Branch("recoKp_pt", &recoKp_pt, "recoKp_pt/F");
    b_recoKp_phi = t->Branch("recoKp_phi", &recoKp_phi, "recoKp_phi/F");
    b_recoKp_eta = t->Branch("recoKp_eta", &recoKp_eta, "recoKp_eta/F");
    b_recoKpKm_dR = t->Branch("recoKpKm_dR", &recoKpKm_dR, "recoKpKm_dR/F");
}

void mcTree::Reset() {
    // Gen
    genW_pt = 0.0;
    genW_eta = 0.0;
    genW_phi = 0.0;
    genW_mass = 0.0;
    genWLepton_id = 0.0;
    genWLepton_pt = 0.0;
    genWLepton_eta = 0.0;
    genWLepton_phi = 0.0;
    genHiggs_pt = 0.0;
    genHiggs_eta = 0.0;
    genHiggs_phi = 0.0;
    genHiggs_mass = 0.0;
    genHiggsMeson_id = 0.0;
    genHiggsMeson_pt = 0.0;
    genHiggsMeson_eta = 0.0;
    genHiggsMeson_phi = 0.0;
    genHiggsMeson_mass = 0.0;
    genHiggsMesonGamma_dR = 0.0;
    genGamma_pt = 0.0;
    genGamma_phi = 0.0;
    genGamma_eta = 0.0;
    genKm_pt = 0.0;
    genKm_phi = 0.0;
    genKm_eta = 0.0;
    genKp_pt = 0.0;
    genKp_phi = 0.0;
    genKp_eta = 0.0;
    genKpKm_dR = 0.0;
    // Reco
    recoPhi_mass = 0.0;
    recoKm_pt = 0.0;
    recoKm_phi = 0.0;
    recoKm_eta = 0.0;
    recoKp_pt = 0.0;
    recoKp_phi = 0.0;
    recoKp_eta = 0.0;
    recoKpKm_dR = 0.0;
    return;
}

double mcTree::dR(float phi1, float phi2, float eta1, float eta2) {

    double dphi = abs(phi2 - phi1);
    if (dphi > M_PI){
        dphi = 2*M_PI - dphi;
    }
    return sqrt(pow((dphi), 2) + pow((eta2 - eta1), 2));
}

void mcTree::FillRecoBranches() {
    // Store good leptons
    vector<int> goodLeptons;

    /* --> PF Candidates <-- */
    // Store K+, K- indexes
    vector<int> KpIdxs;
    vector<int> KmIdxs;
    // START Loop over pfcands -----------------------------
    for (unsigned int i = 0; i < pfcands_p4().size(); i++) {
        // Find Kaons
        int thisID = pfcands_particleId().at(i);
        if (abs(thisID) == 211) {
            if (thisID > 0) KpIdxs.push_back(i);
            else if (thisID < 0) KmIdxs.push_back(i);
        }
    } // END Loop over pfcands -----------------------------

    // START Loop over K+ mesons ---------------------------
    vector<float> phiCands_mass;
    vector<int> phiCands_KpIdx;
    vector<int> phiCands_KmIdx;
    for (unsigned int i = 0; i < KpIdxs.size(); i++) {
        int Kp_i = KpIdxs.at(i);
        float Kp_pt = pfcands_p4()[Kp_i].pt();
        // Only consider high-pt K+ mesons
        if (Kp_pt < 10) continue;
        // START Loop over K- mesons -----------------------
        for (unsigned int j = 0; j < KmIdxs.size(); j++) {
            int Km_j = KmIdxs.at(j);
            float Km_pt = pfcands_p4()[Km_j].pt();
            // Only consider high-pt K- mesons
            if (Km_pt < 10) continue;
            // Get Kaon pair p4
            LorentzVector Kp_p4 = pfcands_p4()[Kp_i];
            LorentzVector Km_p4 = pfcands_p4()[Km_j];
            // Get Kaon pair eta, phi
            float Kp_phi = Kp_p4.phi();
            float Km_phi = Km_p4.phi();
            float Kp_eta = Kp_p4.eta();
            float Km_eta = Km_p4.eta();
            // Get dR between kaons
            float KpKm_dR = dR(Kp_phi, Km_phi, Kp_eta, Km_eta);
            if (KpKm_dR < 0.1) {
                float KpKm_mass = (Kp_p4 + Km_p4).M();    
                phiCands_mass.push_back(KpKm_mass);
                phiCands_KpIdx.push_back(Kp_i);
                phiCands_KmIdx.push_back(Km_j);
            }
        } // END Loop over K- mesons -----------------------
    } // END Loop over K+ mesons ---------------------------

    // Find best Phi Candidate
    int bestPhiCand = 0;
    float bestPhiMassDiff = 10.0;
    float phiMass = 1.019;
    if (phiCands_KpIdx.size() > 0) {
        // START Loop over Phi Candidates ----------------------
        for (unsigned int i = 0; i < phiCands_KpIdx.size(); i++) {
            if (i%3 == 0) {
                float thisDiff = abs(phiMass - phiCands_mass.at(i));
                if (thisDiff < bestPhiMassDiff) {
                    bestPhiMassDiff = thisDiff;
                    bestPhiCand = i;
                }
            }
        } // END Loop over Phi Candidates ----------------------
    }

    /* --> Electrons <-- */
    // START Loop over electrons ---------------------------
    for (unsigned int i = 0; i < els_p4().size(); i++) {
        // Define cuts for readability
        bool elsPtCut = (els_p4().at(i).pt() > 20);
        bool elsIDCut = (isMediumElectronPOGfall17noIso_v2(i));
        bool elsIsoCut = (elMiniRelIsoCMS3_EA(i, 3) < 0.1);
        // Store 'good' electrons
        if (elsPtCut && elsIDCut && elsIsoCut) {
            goodLeptons.push_back(i);
        }
    } // END Loop over electrons ---------------------------

    /* --> Muons <-- */
    // START Loop over muons -------------------------------
    for (unsigned int i = 0; i < mus_p4().size(); i++) {
        // Define cuts for readability
        bool musIDCut = (isMediumMuonPOG(i));
        bool musIsoCut = (muMiniRelIsoCMS3_EA(i, 3) < 0.2);
        // Store 'good' muons
        if (musIDCut && musIsoCut) {
            goodLeptons.push_back(i);
        }
    } // END Loop over muons -------------------------------

    /* --> Photons <-- */
    // START Loop over photons -----------------------------
    for (unsigned int l = 0; l < evt_nphotons(); l++) {
    } // END Loop over photons -----------------------------

    // Retrieve products filled in loops, fill tree branches
    if (phiCands_KpIdx.size() > 0 && abs(bestPhiMassDiff/phiMass) < 0.34) {
        // Best K+, K-
        LorentzVector bestKp_p4 = pfcands_p4().at(phiCands_KpIdx.at(bestPhiCand));
        LorentzVector bestKm_p4 = pfcands_p4().at(phiCands_KmIdx.at(bestPhiCand));
        // Best Phi cand mass
        recoPhi_mass = phiCands_mass.at(bestPhiCand);
        // Best K+, K- kinematics
        recoKp_pt = bestKp_p4.pt();
        recoKp_phi = bestKp_p4.phi();
        recoKp_eta = bestKp_p4.eta();
        recoKm_pt = bestKm_p4.pt();
        recoKm_phi = bestKm_p4.phi();
        recoKm_eta = bestKm_p4.eta();
        recoKpKm_dR = dR(recoKp_phi, recoKm_phi, recoKp_eta, recoKm_eta);
    }

    return;
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
        p4_sum -= p4_sum;
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
        p4_sum -= p4_sum;
        for (unsigned int j = 0; j < daughters[decay].size(); j++) {
            // Get daughter index
            int idx = daughters[decay].at(j);
            // Get ID
            int id = genps_id().at(idx);
            // Fill photon branches
            if (daughters[decay].at(j) == 22) {
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
        p4_sum -= p4_sum;
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
