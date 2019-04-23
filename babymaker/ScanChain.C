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

// CMS3
#include "CORE/CMS3.h"
#include "CORE/ElectronSelections.h"
#include "CORE/MuonSelections.h"
#include "CORE/IsolationTools.h"
#include "CORE/TriggerSelections.h"
#include "CORE/Tools/datasetinfo/getDatasetInfo.h"

// Custom
#include "RPGCORE/mcTree.h"

// Namespaces
using namespace std;
using namespace tas;

// Global Variables
DatasetInfoFromFile datasetInfoFromFile;

int BuildMCTree(TChain* chain, TString outName, TString sampleName, bool verbose = false, bool fast = true, int nEvents = -1, string skimFilePrefix = "test") {

    // Benchmark
    TBenchmark *bmark = new TBenchmark();
    bmark->Start("benchmark");

    // Initialize TFile
    TFile* f = new TFile(outName, "RECREATE");

    // Initialize TTree
    mcTree* mct = new mcTree();
    TTree* mctree = mct->t;

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

            /* --> Start Analysis Code <-- */

            // Scale Factor
            double sf = 1.0;

            // Missing Et
            double met = evt_pfmet();

            // Is Signal
            bool isSignal = sampleName.Contains("WH_HtoRhoGammaPhiGamma");

            // Reset branch values
            mct->Reset();

            // Fill metadata branches
            mct->run = evt_run();
            mct->lumi = evt_lumiBlock();
            mct->event = evt_event();
            if (isSignal) {
                // No sf for signal
                mct->scale1fb = 1;
            }
            else {
                const string datasetName = cms3.evt_dataset().at(0).Data();
                TString cms3_version = cms3.evt_CMS3tag().at(0);
                mct->scale1fb = datasetInfoFromFile.getScale1fbFromFile(datasetName, cms3_version.Data());
            }
            // Fill gen branches
            if (isSignal) mct->FillGenBranches();
            // Fill reco branches
            mct->FillRecoBranches();
            // Fill gen-reco dR branches
            mct->FillGenRecoBranches();

            // Fill tree
            mctree->Fill();

            /* --> END Analysis Code <-- */
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
    mctree->Write();
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
