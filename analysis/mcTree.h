#ifndef MCTREE_H
#define MCTREE_H

class mcTree {

    /* Meta TTree Branches */
    TBranch* b_run;
    TBranch* b_lumi;
    TBranch* b_event;

    /* --> Gen TTree Branches <-- */
    // W Boson
    TBranch* b_genW_pt;
    TBranch* b_genW_eta;
    TBranch* b_genW_phi;
    TBranch* b_genW_mass;
    // Lepton from W Boson
    TBranch* b_genWLepton_id;
    TBranch* b_genWLepton_pt;
    TBranch* b_genWLepton_eta;
    TBranch* b_genWLepton_phi;
    // Higgs Boson
    TBranch* b_genHiggs_pt;
    TBranch* b_genHiggs_eta;
    TBranch* b_genHiggs_phi;
    TBranch* b_genHiggs_mass;
    // Mesons from Higgs
    TBranch* b_genHiggsMeson_id;
    TBranch* b_genHiggsMeson_pt;
    TBranch* b_genHiggsMeson_eta;
    TBranch* b_genHiggsMeson_phi;
    TBranch* b_genHiggsMeson_mass;
    TBranch* b_genHiggsMesonGamma_dR;
    // Gamma from Higgs
    TBranch* b_genGamma_pt;
    TBranch* b_genGamma_phi;
    TBranch* b_genGamma_eta;
    // K- from Phi
    TBranch* b_genKm_pt;
    TBranch* b_genKm_phi;
    TBranch* b_genKm_eta;
    // K+ from Phi
    TBranch* b_genKp_pt;
    TBranch* b_genKp_phi;
    TBranch* b_genKp_eta;
    // dR between Kaons
    TBranch* b_genKpKm_dR;

    /* --> Reco TTree Branches <-- */
    // Reconstructed Higgs
    TBranch* b_recoHiggs_mass;
    // Reconstructed Phi
    TBranch* b_recoPhi_mass;
    TBranch* b_recoPhi_iso;
    // K- from Phi
    TBranch* b_recoKm_pt;
    TBranch* b_recoKm_phi;
    TBranch* b_recoKm_eta;
    // K+ from Phi
    TBranch* b_recoKp_pt;
    TBranch* b_recoKp_phi;
    TBranch* b_recoKp_eta;
    // dR between Kaons
    TBranch* b_recoKpKm_dR;
    // Reconstructed Rho
    TBranch* b_recoRho_mass;
    TBranch* b_recoRho_iso;
    // Pi- from Rho
    TBranch* b_recoPim_pt; 
    TBranch* b_recoPim_phi;
    TBranch* b_recoPim_eta;
    // Pi+ from Rho
    TBranch* b_recoPip_pt;
    TBranch* b_recoPip_phi;
    TBranch* b_recoPip_eta;
    // dR between Pions
    TBranch* b_recoPipPim_dR;
    // Photons
    TBranch* b_recoGamma_pt;
    TBranch* b_recoGamma_phi;
    TBranch* b_recoGamma_eta;
    TBranch* b_recoGamma_iso;

    public:
        // TTree
        TTree* t;

        /* --> Meta TTree Branch Values */
        int run;
        int lumi;
        int event;

        /* --> Gen TTree Branch Values <-- */
        // W Boson
        float genW_pt;
        float genW_eta;
        float genW_phi;                                                             
        float genW_mass;
        // Lepton from W Boson
        float genWLepton_id;
        float genWLepton_pt;
        float genWLepton_eta;
        float genWLepton_phi;                                                             
        // Higgs Boson
        float genHiggs_pt;
        float genHiggs_eta;
        float genHiggs_phi;                                                             
        float genHiggs_mass;
        // Mesons from Higgs
        float genHiggsMeson_id;
        float genHiggsMeson_pt;
        float genHiggsMeson_eta;
        float genHiggsMeson_phi;
        float genHiggsMeson_mass;
        float genHiggsMesonGamma_dR;
        // Gamma from Higgs
        float genGamma_pt;
        float genGamma_phi;
        float genGamma_eta;
        // K- from Phi
        float genKm_pt;
        float genKm_phi;
        float genKm_eta;
        // K+ from Phi
        float genKp_pt;
        float genKp_phi;
        float genKp_eta;
        // dR between Kaons
        float genKpKm_dR;

        /* --> Reco TTree Branch Values <-- */
        // Reconstructed Higgs
        float recoHiggs_mass;
        // Reconstructed Phi
        float recoPhi_mass;
        float recoPhi_iso;
        // K- from Phi
        float recoKm_pt;
        float recoKm_phi;
        float recoKm_eta;
        // K+ from Phi
        float recoKp_pt;
        float recoKp_phi;
        float recoKp_eta;
        // dR between Kaons
        float recoKpKm_dR;
        // Reconstructed Rho
        float recoRho_mass;
        float recoRho_iso;
        // Pi- from Rho
        float recoPim_pt; 
        float recoPim_phi;
        float recoPim_eta;
        // Pi+ from Rho
        float recoPip_pt;
        float recoPip_phi;
        float recoPip_eta;
        // dR between Pions
        float recoPipPim_dR;
        // Photons
        float recoGamma_pt;
        float recoGamma_phi;
        float recoGamma_eta;
        float recoGamma_iso;

        /* --> Methods <-- */
        // Constructor
        mcTree();
        // Reset Vars
        void Reset();
        double dR(float, float, float, float);
        void FillRecoBranches();
        void FillGenBranches();
};

#endif
