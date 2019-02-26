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

using namespace std;
using namespace tas;

struct System {
    // Member Variables
    int motherIdx;
    int motherID;
    int tag;
    LorentzVector momenta;
    // Overload Constructor
    System(int, int);
    // Destructor
    ~System();
    // Methods
    void Add(int, LorentzVector);
    double Mass();
    string Tag();
};

System::System(int new_motherIdx, int new_motherID) {
    motherIdx = new_motherIdx;
    motherID = new_motherID;
    tag = 0;
}

System::~System() {}

void System::Add(int daughterID, LorentzVector new_p4) {
    momenta += new_p4;
    tag += daughterID;
    return;
}

double System::Mass() {
    return momenta.M();
}

string System::Tag() {
    return to_string(motherID)+"_"+to_string(tag);
}



int SanityCheck(TChain* chain, char sample_name[], bool fast = true, int nEvents = -1, string skimFilePrefix = "test") {

    // Benchmark
    TBenchmark *bmark = new TBenchmark();
    bmark->Start("benchmark");

    // Directory Setup
    TDirectory *rootdir = gDirectory->GetDirectory("Rint:");

    // From Higgs
    TH1F* H_to_PhiGamma_mass = new TH1F("H_to_PhiGamma_mass", "", 100,120,130);
    string H_to_PhiGamma_tag = to_string(25)+to_string(333+22);
    TH1F* H_to_RhoGamma_mass = new TH1F("H_to_RhoGamma_mass", "", 100,120,130);
    string H_to_RhoGamma_tag = to_string(25)+to_string(113+22);
    TH1F* H_to_KKGamma_mass = new TH1F("H_to_KKGamma_mass", "", 100,120,130);
    string H_to_KKGamma_tag = to_string(25)+to_string(130+310+22);
    // From Phi
    TH1F* Phi_to_KK_mass = new TH1F("Phi_to_KK_mass", "", 100,0,2);
    string H_to_KKGamma_tag = to_string(333)+to_string(130+310);

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
            bool verbose = false;

            // Systems
            vector<int> motherIdxs;
            vector<System*> systems;

            // Loop over gen-level data
            for (unsigned int i = 0; i < genps_id().size(); i++) {
                // Get PDG ID of daugher, mother
                int thisID = genps_id().at(i);
                int thisMotherID = genps_id_mother().at(i);
                // Reject uninteresting daughters/mothers
                if (thisMotherID != 25 && thisMotherID != 333) continue;
                if (thisID != 333 && thisID != 22 && thisID != 113 && thisID != 130 && thisID != 310) continue;
                // Get index of mother
                int thisMotherIdx = genps_idx_mother().at(i);
                // Check if system already initialized
                int masterIdx = -1;
                auto iter = find(motherIdxs.begin(), motherIdxs.end(), thisMotherIdx);
                if (iter == motherIdxs.end()) {
                    System* newSystem = new System(thisMotherIdx, thisMotherID);
                    systems.push_back(newSystem);
                    motherIdxs.push_back(thisMotherIdx);
                    masterIdx = (motherIdxs.size()-1);
                }
                else {
                    masterIdx = distance(motherIdxs.begin(), iter);
                }
                // Get current system
                System* curSystem = systems.at(masterIdx); 
                // Add particle to system
                curSystem->Add(thisID, genps_p4()[i]);
                // Particle data printout
                if (verbose) {
                    cout << "-------------- PARTICLE DATA --------------\n";
                    cout << "Particle: " << thisID << "\n";
                    cout << "System: " << curSystem->motherID << "(PDG), " << curSystem->motherIdx << "(Index)\n";
                    cout << "pt: " << curSystem->momenta << "\n";
                    cout << "tag: " << curSystem->Tag() << "\n";
                    cout << "mass: " << curSystem->Mass() << "\n";
                    cout << "-------------------------------------------\n";
                }
                // Get unique tag for system
                string curTag = curSystem->Tag();
                // Fill histogram if is system complete
                if (curTag == H_to_PhiGamma_tag) {
                    H_to_PhiGamma_mass->Fill(curSystem->Mass(), sf);
                    delete curSystem;
                }
                else if (curTag == H_to_RhoGamma_tag) {
                    H_to_RhoGamma_mass->Fill(curSystem->Mass(), sf);
                    delete curSystem;
                }
                else if (curTag == Phi_to_KK_tag) {
                    Phi_to_KK_mass->Fill(curSystem->Mass(), sf);
                    delete curSystem;
                }
            }
            
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
