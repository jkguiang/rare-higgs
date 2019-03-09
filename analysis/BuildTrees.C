// -*- C++ -*-
// Usage:
// > root -b -q doAll.C

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

// CMS3
#include "/home/users/jguiang/projects/mt2/MT2Analysis/CORE/CMS3.h"
#include "/home/users/jguiang/projects/mt2/MT2Analysis/CORE/MuonSelections.h"
#include "/home/users/jguiang/projects/mt2/MT2Analysis/CORE/TriggerSelections.h"

// Custom
#include "Decay.h"
#include "GenTree.h"

using namespace std;
using namespace tas;

double dR(float phi1, float phi2, float eta1, float eta2) {

    double dphi = abs(phi2 - phi1);
    if (dphi > M_PI){
        dphi = 2*M_PI - dphi;
    }
    return sqrt(pow((dphi), 2) + pow((eta2 - eta1), 2));
}

int BuildTrees(TChain* chain, char sample_name[], bool verbose = false, bool fast = true, int nEvents = -1, string skimFilePrefix = "test") {

    // Benchmark
    TBenchmark *bmark = new TBenchmark();
    bmark->Start("benchmark");

    // Iniitialize TFile
    TFile* f = new TFile(sample_name, "RECREATE");

    // Get Trees
    GenTree* gt = new GenTree();
    TTree* gtree = gt->t;

    /* --> Decay Setup <-- */
    // Vector for communicating product IDs to Decay object
    vector<int> prods;
    // Higgs to Phi, Gamma
    prods = {333, 22};
    Decay* H_to_PhiGamma = new Decay(25, prods);
    // Higgs to Rho, Gamma
    prods = {113, 22};
    Decay* H_to_RhoGamma = new Decay(25, prods);
    // Higgs to K+, K-, Gamma
    prods = {321, -321, 22};
    Decay* H_to_KKGamma = new Decay(25, prods);
    // Phi to K+, K-
    prods = {321, -321};
    Decay* Phi_to_KK = new Decay(333, prods);
    // W to electron, electron neutrino
    prods = {-11};
    Decay* W_to_ElNu = new Decay(24, prods);
    // W to muon, muon neutrino
    prods = {-13};
    Decay* W_to_MuNu = new Decay(24, prods);

    // Loop over events to Analyze
    unsigned int nEventsTotal = 0;
    unsigned int nEventsChain = chain->GetEntries();
    if (nEvents >= 0) nEventsChain = nEvents;
    TObjArray *listOfFiles = chain->GetListOfFiles();
    TIter fileIter(listOfFiles);
    TFile *currentFile = 0;

    // File Loop
    while ( (currentFile = (TFile*)fileIter.Next()) ) {

        // Get File Content
        TFile file(currentFile->GetTitle());
        TTree *tree = (TTree*)file.Get("Events");
        if (fast) TTreeCache::SetLearnEntries(10);
        if (fast) tree->SetCacheSize(128*1024*1024);
        cms3.Init(tree);

        // Loop over Events in current file
        if (nEventsTotal >= nEventsChain) continue;
        unsigned int nEventsTree = tree->GetEntriesFast();
        for (unsigned int event = 0; event < nEventsTree; ++event) {


            // Get Event Content
            if (nEventsTotal >= nEventsChain) continue;
            if (fast) tree->LoadTree(event);
            cms3.GetEntry(event);
            ++nEventsTotal;

            // Progress
            CMS3::progress( nEventsTotal, nEventsChain );

            /* --> Gen-level Analysis Code <-- */

            // Scale Factor
            double sf = 1.0;

            // Reset branch values
            gt->Reset();

            // Reset decays
            H_to_PhiGamma->Reset();
            H_to_RhoGamma->Reset();
            H_to_KKGamma->Reset();
            Phi_to_KK->Reset();
            W_to_ElNu->Reset();
            W_to_MuNu->Reset();

            // START Loop over gen-level data ------------
            for (unsigned int i = 0; i < genps_id().size(); i++) {
                // Get PDG ID of daugher, mother
                int thisID = genps_id().at(i);
                int thisMotherID = genps_id_mother().at(i);
                // Veto duplicates
                if (!genps_isLastCopy().at(i)) continue;
                if (!genps_isLastCopy()[genps_idx_mother().at(i)]) continue;
                // Veto weird photons
                if (thisID == 22 && !genps_isHardProcess().at(i)) continue;
                // Fill branches for mother particles
                if (thisID == 25) {
                    // Fill Higgs boson branches
                    gt->genHiggs_pt = genps_p4()[i].pt();
                    gt->genHiggs_eta = genps_p4()[i].eta();
                    gt->genHiggs_phi = genps_p4()[i].phi(); 
                }
                if (thisID == 24 || thisID == -24) {
                    // Fill W boson branches
                    gt->genW_pt = genps_p4()[i].pt();
                    gt->genW_eta = genps_p4()[i].eta();
                    gt->genW_phi = genps_p4()[i].phi(); 
                }
                // Add particle id, idx, p4 to relevant systems
                if (thisMotherID == 24 || thisMotherID == -24) {
                    if (thisID == -11 || thisID == 11) {
                        W_to_ElNu->Add(thisID, i, genps_p4()[i]);
                    }
                    if (thisID == -13 || thisID == 13) {
                        W_to_MuNu->Add(thisID, i, genps_p4()[i]);
                    }
                }
                if (thisMotherID == 25) {
                    if (thisID == 22 || thisID == 333) {
                        H_to_PhiGamma->Add(thisID, i, genps_p4()[i]);
                    }
                    if (thisID == 22 || thisID == 113) {
                        H_to_RhoGamma->Add(thisID, i, genps_p4()[i]);
                    }
                    if (thisID == 22) {
                        H_to_KKGamma->Add(thisID, i, genps_p4()[i]);
                    }
                }
                if (thisMotherID == 333) {
                    if ((thisID == 321 || thisID == -321) && genps_id_mother()[genps_idx_mother().at(i)] == 25) {
                        Phi_to_KK->Add(thisID, i, genps_p4()[i]);
                        H_to_KKGamma->Add(thisID, i, genps_p4()[i]);
                    }
                }
            } // END Loop over gen-level data ------------

            // Retrieve products filled in gen-level loop, fill tree branches
            vector<int> daughterIDs;
            vector<int> daughterIdxs;

            if (W_to_MuNu->products->Tag() == "24_-13" || W_to_MuNu->products->Tag() == "24_13") {
                // Daughters
                daughterIDs = W_to_MuNu->products->daughterIDs;
                daughterIdxs = W_to_MuNu->products->daughterIdxs;
                for (unsigned int j = 0; j < daughterIDs.size(); j++) {
                    // Get daughter index
                    int k = daughterIdxs.at(j);
                    // Skip neutrinos
                    if (daughterIDs.at(j) == 14) {
                        continue;
                    }
                    else {
                        // Fill leptons from W branches
                        gt->genWLepton_id = daughterIDs.at(j);
                        gt->genWLepton_pt = genps_p4()[k].pt();
                        gt->genWLepton_eta = genps_p4()[k].eta();
                        gt->genWLepton_phi = genps_p4()[k].phi();
                    }
                }
            }

            if (W_to_ElNu->products->Tag() == "24_-11" || W_to_ElNu->products->Tag() == "24_11") {
                // Daughters
                daughterIDs = W_to_ElNu->products->daughterIDs;
                daughterIdxs = W_to_ElNu->products->daughterIdxs;
                for (unsigned int j = 0; j < daughterIDs.size(); j++) {
                    // Get daughter index
                    int k = daughterIdxs.at(j);
                    // Skip neutrinos
                    if (daughterIDs.at(j) == 12) {
                        continue;
                    }
                    else {
                        // Fill leptons from W branches
                        gt->genWLepton_id = daughterIDs.at(j);
                        gt->genWLepton_pt = genps_p4()[k].pt();
                        gt->genWLepton_eta = genps_p4()[k].eta();
                        gt->genWLepton_phi = genps_p4()[k].phi();
                    }
                }
            }

            if (H_to_PhiGamma->Full()) {
                // Daughters
                daughterIDs = H_to_PhiGamma->products->daughterIDs;
                daughterIdxs = H_to_PhiGamma->products->daughterIdxs;
                gt->genHiggs_mass = H_to_PhiGamma->products->Mass();
                for (unsigned int j = 0; j < daughterIDs.size(); j++) {
                    // Get daughter index
                    int k = daughterIdxs.at(j);
                    // Fill photon branches
                    if (daughterIDs.at(j) == 22) {
                        gt->genGamma_pt = genps_p4()[k].pt();
                        gt->genGamma_eta = genps_p4()[k].eta();
                        gt->genGamma_phi = genps_p4()[k].phi();
                    }
                    else {
                        // Fill mesons from Higgs branches
                        gt->genHiggsMeson_id = daughterIDs.at(j);
                        gt->genHiggsMeson_pt = genps_p4()[k].pt();
                        gt->genHiggsMeson_eta = genps_p4()[k].eta();
                        gt->genHiggsMeson_phi = genps_p4()[k].phi();
                    }
                }
                gt->genHiggsMesonGamma_dR = dR(gt->genHiggsMeson_phi, gt->genGamma_phi, gt->genHiggsMeson_eta, gt->genGamma_eta);
            }

            if (H_to_RhoGamma->Full()) {
                // Daughters
                daughterIDs = H_to_RhoGamma->products->daughterIDs;
                daughterIdxs = H_to_RhoGamma->products->daughterIdxs;
                gt->genHiggs_mass = H_to_RhoGamma->products->Mass();
                for (unsigned int j = 0; j < daughterIDs.size(); j++) {
                    // Get daughter index
                    int k = daughterIdxs.at(j);
                    // Fill photon branches
                    if (daughterIDs.at(j) == 22) {
                        gt->genGamma_pt = genps_p4()[k].pt();
                        gt->genGamma_eta = genps_p4()[k].eta();
                        gt->genGamma_phi = genps_p4()[k].phi();
                    }
                    else {
                        // Fill mesons from Higgs branches
                        gt->genHiggsMeson_id = daughterIDs.at(j);
                        gt->genHiggsMeson_pt = genps_p4()[k].pt();
                        gt->genHiggsMeson_eta = genps_p4()[k].eta();
                        gt->genHiggsMeson_phi = genps_p4()[k].phi();
                    }
                }
                gt->genHiggsMesonGamma_dR = dR(gt->genHiggsMeson_phi, gt->genGamma_phi, gt->genHiggsMeson_eta, gt->genGamma_eta);
            }

            if (Phi_to_KK->Full()) {
                // Daughters
                daughterIDs = Phi_to_KK->products->daughterIDs;
                daughterIdxs = Phi_to_KK->products->daughterIdxs;
                for (unsigned int j = 0; j < daughterIDs.size(); j++) {
                    // Get daughter index
                    int k = daughterIdxs.at(j);
                    // Fill Kaon from Phi branches
                    if (daughterIDs.at(j) == 321) {
                        gt->genKp_pt = genps_p4()[k].pt();
                        gt->genKp_phi = genps_p4()[k].phi();
                        gt->genKp_eta = genps_p4()[k].eta();
                    }
                    else if (daughterIDs.at(j) == -321) {
                        gt->genKm_pt = genps_p4()[k].pt();
                        gt->genKm_phi = genps_p4()[k].phi();
                        gt->genKm_eta = genps_p4()[k].eta();
                    }
                }
                gt->genKpKm_dR = dR(gt->genKp_phi, gt->genKm_phi, gt->genKp_eta, gt->genKm_eta);
            }

            // Fill tree
            gtree->Fill();

            /* --> END Gen-level Analysis Code <-- */
            
        }
  
        // Clean Up
        delete tree;
        file.Close();
    }
    if (nEventsChain != nEventsTotal) {
        cout << Form( "ERROR: number of events from files (%d) is not equal to total number of events (%d)", nEventsChain, nEventsTotal ) << endl;
    }
  
    // Write trees
    f->cd();
    gtree->Write();
    f->Close();

    // return
    bmark->Stop("benchmark");
    cout << endl;
    cout << nEventsTotal << " Events Processed" << endl;
    cout << "------------------------------" << endl;
    cout << "CPU  Time: " << Form( "%.01f", bmark->GetCpuTime("benchmark")  ) << endl;
    cout << "Real Time: " << Form( "%.01f", bmark->GetRealTime("benchmark") ) << endl;
    cout << endl;
    delete bmark;
    return 0;
}
