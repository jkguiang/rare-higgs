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

int SanityCheck(TChain* chain, char sample_name[], bool verbose = false, bool fast = true, int nEvents = -1, string skimFilePrefix = "test") {

    // Benchmark
    TBenchmark *bmark = new TBenchmark();
    bmark->Start("benchmark");

    // Directory Setup
    TDirectory *rootdir = gDirectory->GetDirectory("Rint:");

    vector<int> prods;
    // From Higgs
    prods = {333, 22};
    Decay* H_to_PhiGamma = new Decay(25, prods);
    TH1F* H_to_PhiGamma_mass = new TH1F("H_to_PhiGamma_mass", "", 100,120,130);
    prods = {113, 22};
    Decay* H_to_RhoGamma = new Decay(25, prods);
    TH1F* H_to_RhoGamma_mass = new TH1F("H_to_RhoGamma_mass", "", 100,120,130);
    prods = {321, -321, 22};
    Decay* H_to_RhoGamma = new Decay(25, prods);
    TH1F* H_to_KKGamma_mass = new TH1F("H_to_KKGamma_mass", "", 100,120,130);
    // From Phi
    prods = {321, -321};
    Decay* Phi_to_KK = new Decay(333, prods);
    TH1F* Phi_to_KK_mass = new TH1F("Phi_to_KK_mass", "", 100,0,2);

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

            /* --> Analysis Code <-- */

            // Scale1fb
            double sf = 1.0;

            // Add new systems to decays
            H_to_PhiGamma->NewSystem();
            H_to_RhoGamma->NewSystem();
            H_to_KKGamma->NewSystem();
            Phi_to_KK->NewSystem();

            // START Loop over gen-level data ------------
            for (unsigned int i = 0; i < genps_id().size(); i++) {
                // Get PDG ID of daugher, mother
                int thisID = genps_id().at(i);
                int thisMotherID = genps_id_mother().at(i);
                // Reject uninteresting daughters/mothers
                if (thisMotherID != 25 && thisMotherID != 333) continue;
                if (thisID != 333 && thisID != 22 && thisID != 113 && thisID != 321 && thisID != -321) continue;
                if (!genps_isLastCopy().at(i)) continue;
                // Add particle to relevant systems
                if (thisMotherID == 25) {
                    if (thisID == 22 || thisID == 333) {
                        H_to_PhiGamma->Add(thisID, genps_p4()[i]);
                    }
                    if (thisID == 22 || thisID == 113) {
                        H_to_RhoGamma->Add(thisID, genps_p4()[i]);
                    }
                    if (thisID == 22 || thisID == 321 || thisID == -321) {
                        H_to_KKGamma->Add(thisID, genps_p4()[i]);
                    }
                }
                if (thisMotherID == 333) {
                    if (thisID == 321 || thisID == -321) {
                        Phi_to_KK->Add(thisID, genps_p4()[i]);
                    }
                }

                // Particle data printout
                if (verbose) {
                    cout << "-------------- PARTICLE DATA --------------\n";
                    cout << "Particle: " << thisID << "\n";
                    cout << "-------------------------------------------\n";
                }

                if (H_to_PhiGamma->Full() && H_to_RhoGamma->Full() && H_to_KKGamma->Full() && Phi_to_KK->Full()) break;
            } // END Loop over gen-level data ------------

            // Retrieve systems filled in gen-level loop
            System* H_to_PhiGamma_sys = H_to_PhiGamma->GetSystem();
            System* H_to_RhoGamma_sys = H_to_RhoGamma->GetSystem();
            System* H_to_KKGamma_sys = H_to_KKGamma->GetSystem();
            System* Phi_to_KK_sys = Phi_to_KK->GetSystem();
            // Fill Histograms
            H_to_PhiGamma_mass->Fill(H_to_PhiGamma_sys->Mass(), sf);
            H_to_RhoGamma_mass->Fill(H_to_RhoGamma_sys->Mass(), sf);
            H_to_KKGamma_mass->Fill(H_to_KKGamma_sys->Mass(), sf);
            Phi_to_KK_mass->Fill(Phi_to_KK_sys->Mass(), sf);
            
        }
  
        // Clean Up
        delete tree;
        file.Close();
    }
    if (nEventsChain != nEventsTotal) {
        cout << Form( "ERROR: number of events from files (%d) is not equal to total number of events (%d)", nEventsChain, nEventsTotal ) << endl;
    }
  
    // Write Histograms
    TFile* f = new TFile(sample_name, "RECREATE");

    H_to_PhiGamma_mass->Write();
    H_to_RhoGamma_mass->Write();
    H_to_KKGamma_mass->Write();
    Phi_to_KK_mass->Write();
    
    f->Close();

    // Clear histograms  
    H_to_PhiGamma_mass->Delete();
    H_to_RhoGamma_mass->Delete();
    H_to_KKGamma_mass->Delete();
    Phi_to_KK_mass->Delete();

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
