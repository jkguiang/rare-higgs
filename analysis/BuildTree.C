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

using namespace std;
using namespace tas;

int BuildTree(TChain* chain, char sample_name[], bool verbose = false, bool fast = true, int nEvents = -1, string skimFilePrefix = "test") {

    // Benchmark
    TBenchmark *bmark = new TBenchmark();
    bmark->Start("benchmark");

    // Iniitialize TFile
    TFile* f = new TFile(sample_name, "RECREATE");

    /* --> TTree Setup <-- */
    TTree* t = new TTree("tree", "tree");
    // W Boson
    float genW_pt;
    TBranch* b_genW_pt = t->Branch("genW_pt", &genW_pt, "genW_pt/F");
    float genW_eta;
    TBranch* b_genW_eta = t->Branch("genW_eta", &genW_eta, "genW_eta/F");
    float genW_phi;                                                             
    TBranch* b_genW_phi = t->Branch("genW_phi", &genW_phi, "genW_phi/F");
    // Lepton from W Boson
    float genWLepton_id;
    TBranch* b_genWLepton_id = t->Branch("genWLepton_id", &genWLepton_id, "genWLepton_id/F");
    float genWLepton_pt;
    TBranch* b_genWLepton_pt = t->Branch("genWLepton_pt", &genWLepton_pt, "genWLepton_pt/F");
    float genWLepton_eta;
    TBranch* b_genWLepton_eta = t->Branch("genWLepton_eta", &genWLepton_eta, "genWLepton_eta/F");
    float genWLepton_phi;                                                             
    TBranch* b_genWLepton_phi = t->Branch("genWLepton_phi", &genWLepton_phi, "genWLepton_phi/F");
    // Higgs Boson
    float genHiggs_pt;
    TBranch* b_genHiggs_pt = t->Branch("genHiggs_pt", &genHiggs_pt, "genHiggs_pt/F");
    float genHiggs_eta;
    TBranch* b_genHiggs_eta = t->Branch("genHiggs_eta", &genHiggs_eta, "genHiggs_eta/F");
    float genHiggs_phi;                                                             
    TBranch* b_genHiggs_phi = t->Branch("genHiggs_phi", &genHiggs_phi, "genHiggs_phi/F");
    float genHiggs_mass;
    TBranch* b_genHiggs_mass = t->Branch("genHiggs_mass", &genHiggs_mass, "genHiggs_mass/F");
    // Mesons from Higgs
    float genHiggsMeson_id;
    TBranch* b_genHiggsMeson_id = t->Branch("genHiggsMeson_id", &genHiggsMeson_id, "genHiggsMeson_id/F");
    float genHiggsMeson_pt;
    TBranch* b_genHiggsMeson_pt = t->Branch("genHiggsMeson_pt", &genHiggsMeson_pt, "genHiggsMeson_pt/F");
    float genHiggsMeson_eta;
    TBranch* b_genHiggsMeson_eta = t->Branch("genHiggsMeson_eta", &genHiggsMeson_eta, "genHiggsMeson_eta/F");
    float genHiggsMeson_phi;
    TBranch* b_genHiggsMeson_phi = t->Branch("genHiggsMeson_phi", &genHiggsMeson_phi, "genHiggsMeson_phi/F");
    // Gamma from Higgs
    float genGamma_pt;
    TBranch* b_genGamma_pt = t->Branch("genGamma_pt", &genGamma_pt, "genGamma_pt/F");
    float genGamma_phi;
    TBranch* b_genGamma_phi = t->Branch("genGamma_phi", &genGamma_phi, "genGamma_phi/F");
    float genGamma_eta;
    TBranch* b_genGamma_eta = t->Branch("genGamma_eta", &genGamma_eta, "genGamma_eta/F");
    // K- from Phi
    float genKm_pt;
    TBranch* b_genKm_pt = t->Branch("genKm_pt", &genKm_pt, "genKm_pt/F");
    float genKm_phi;
    TBranch* b_genKm_phi = t->Branch("genKm_phi", &genKm_phi, "genKm_phi/F");
    float genKm_eta;
    TBranch* b_genKm_eta = t->Branch("genKm_eta", &genKm_eta, "genKm_eta/F");
    // K+ from Phi
    float genKp_pt;
    TBranch* b_genKp_pt = t->Branch("genKp_pt", &genKp_pt, "genKp_pt/F");
    float genKp_phi;
    TBranch* b_genKp_phi = t->Branch("genKp_phi", &genKp_phi, "genKp_phi/F");
    float genKp_eta;
    TBranch* b_genKp_eta = t->Branch("genKp_eta", &genKp_eta, "genKp_eta/F");

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

            // Reset branch values
            genW_pt = 0.0;
            genW_eta = 0.0;
            genW_phi = 0.0;                                                             
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
            genGamma_pt = 0.0;
            genGamma_phi = 0.0;
            genGamma_eta = 0.0;
            genKm_pt = 0.0;
            genKm_phi = 0.0;
            genKm_eta = 0.0;
            genKp_pt = 0.0;
            genKp_phi = 0.0;
            genKp_eta = 0.0;

            // Get Event Content
            if (nEventsTotal >= nEventsChain) continue;
            if (fast) tree->LoadTree(event);
            cms3.GetEntry(event);
            ++nEventsTotal;

            // Progress
            CMS3::progress( nEventsTotal, nEventsChain );

            /* --> Analysis Code <-- */

            // Scale Factor
            double sf = 1.0;

            // Reset decays
            H_to_PhiGamma->Reset();
            H_to_RhoGamma->Reset();
            H_to_KKGamma->Reset();
            Phi_to_KK->Reset();
            W_to_ElNu->Reset();
            W_to_MuNu->Reset();

            // START Loop over gen-level data ------------
            for (unsigned int i = 0; i < genps_id().size(); i++) {
                // Veto duplicates
                if (!genps_isLastCopy().at(i)) continue;
                if (!genps_isLastCopy()[genps_idx_mother().at(i)]) continue;
                // Get PDG ID of daugher, mother
                int thisID = genps_id().at(i);
                int thisMotherID = genps_id_mother().at(i);
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
                    cout << thisID  << " from " << thisMotherID << "[index = " << genps_idx_mother().at(i) << "]\n";
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
                    if (thisID == 321 || thisID == -321) {
                        Phi_to_KK->Add(thisID, i, genps_p4()[i]);
                        H_to_KKGamma->Add(thisID, i, genps_p4()[i]);
                    }
                }

                // Particle data printout
                if (verbose) {
                    cout << "-------------- PARTICLE DATA --------------\n";
                    cout << "Particle: " << thisID << "\n";
                    cout << "-------------------------------------------\n";
                }

                if ( (H_to_PhiGamma->Full() && H_to_KKGamma->Full() && Phi_to_KK->Full()) || H_to_RhoGamma->Full() ) break;
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
                        genWLepton_id = daughterIDs.at(j);
                        genWLepton_pt = genps_p4()[k].pt();
                        genWLepton_eta = genps_p4()[k].eta();
                        genWLepton_phi = genps_p4()[k].phi();
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
                        genWLepton_id = daughterIDs.at(j);
                        genWLepton_pt = genps_p4()[k].pt();
                        genWLepton_eta = genps_p4()[k].eta();
                        genWLepton_phi = genps_p4()[k].phi();
                    }
                }
            }

            if (H_to_PhiGamma->Full()) {
                // Daughters
                daughterIDs = H_to_PhiGamma->products->daughterIDs;
                daughterIdxs = H_to_PhiGamma->products->daughterIdxs;
                genHiggs_mass = H_to_PhiGamma->products->Mass();
                for (unsigned int j = 0; j < daughterIDs.size(); j++) {
                    // Get daughter index
                    int k = daughterIdxs.at(j);
                    // Fill photon branches
                    if (daughterIDs.at(j) == 22) {
                        genGamma_pt = genps_p4()[k].pt();
                        genGamma_eta = genps_p4()[k].eta();
                        genGamma_phi = genps_p4()[k].phi();
                    }
                    else {
                        // Fill mesons from Higgs branches
                        genHiggsMeson_id = daughterIDs.at(j);
                        genHiggsMeson_pt = genps_p4()[k].pt();
                        genHiggsMeson_eta = genps_p4()[k].eta();
                        genHiggsMeson_phi = genps_p4()[k].phi();
                    }
                }
            }

            if (H_to_RhoGamma->Full()) {
                // Daughters
                daughterIDs = H_to_RhoGamma->products->daughterIDs;
                daughterIdxs = H_to_RhoGamma->products->daughterIdxs;
                genHiggs_mass = H_to_RhoGamma->products->Mass();
                for (unsigned int j = 0; j < daughterIDs.size(); j++) {
                    // Get daughter index
                    int k = daughterIdxs.at(j);
                    // Fill photon branches
                    if (daughterIDs.at(j) == 22) {
                        genGamma_pt = genps_p4()[k].pt();
                        genGamma_eta = genps_p4()[k].eta();
                        genGamma_phi = genps_p4()[k].phi();
                    }
                    else {
                        // Fill mesons from Higgs branches
                        genHiggsMeson_id = daughterIDs.at(j);
                        genHiggsMeson_pt = genps_p4()[k].pt();
                        genHiggsMeson_eta = genps_p4()[k].eta();
                        genHiggsMeson_phi = genps_p4()[k].phi();
                    }
                }
            }

            if (H_to_KKGamma->Full()) {
                // Daughters
                daughterIDs = H_to_KKGamma->products->daughterIDs;
                daughterIdxs = H_to_KKGamma->products->daughterIdxs;
                for (unsigned int j = 0; j < daughterIDs.size(); j++) {
                    // Get daughter index
                    int k = daughterIdxs.at(j);
                    // Fill Kaon from Phi branches
                    if (daughterIDs.at(j) == 321) {
                        genKp_pt = genps_p4()[k].pt();
                        genKp_phi = genps_p4()[k].phi();
                        genKp_eta = genps_p4()[k].eta();
                    }
                    else if (daughterIDs.at(j) == -321) {
                        genKm_pt = genps_p4()[k].pt();
                        genKm_phi = genps_p4()[k].phi();
                        genKm_eta = genps_p4()[k].eta();
                    }
                }
            }

            // Fill tree
            t->Fill();
            
        }
  
        // Clean Up
        delete tree;
        file.Close();
    }
    if (nEventsChain != nEventsTotal) {
        cout << Form( "ERROR: number of events from files (%d) is not equal to total number of events (%d)", nEventsChain, nEventsTotal ) << endl;
    }
  
    // Write Histograms
    f->cd();
    t->Write();
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
