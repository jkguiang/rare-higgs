class GenTree {

    /* --> TTree Branches <-- */
    // W Boson
    TBranch* b_genW_pt;
    TBranch* b_genW_eta;
    TBranch* b_genW_phi;
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

    public:
        // TTree
        TTree* t;
        /* --> TTree Branch Values <-- */
        // W Boson
        float genW_pt;
        float genW_eta;
        float genW_phi;                                                             
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

        /* --> Methods <-- */
        // Constructor
        GenTree();
        // Reset Vars
        void Reset();
};

GenTree::GenTree() {
    /* --> TTree Setup <-- */
    t = new TTree("tree", "tree");
    b_genW_pt = t->Branch("genW_pt", &genW_pt, "genW_pt/F");
    b_genW_eta = t->Branch("genW_eta", &genW_eta, "genW_eta/F");
    b_genW_phi = t->Branch("genW_phi", &genW_phi, "genW_phi/F");
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
}

void GenTree::Reset() {
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
    return;
}
