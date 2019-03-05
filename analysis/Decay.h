// C++
#include <iostream>
#include <vector>
#include <string>

// CMS3
#include "/home/users/jguiang/projects/mt2/MT2Analysis/CORE/CMS3.h"

using namespace std;
using namespace tas;

struct Products {
    // --> Member Variables <--
    // Mother identification (for tag)
    int motherID;
    // Products identification
    vector<int> daughterIDs;
    vector<int> daughterIdxs;
    // Physical properties
    LorentzVector p4;
    // Overload Constructor
    Products(int);
    // Destructor
    ~Products();
    // Methods
    void Add(int, int, LorentzVector);
    float Mass();
    string Tag();
};

Products::Products(int new_motherID) {
    motherID = new_motherID;
}

Products::~Products() {}

void Products::Add(int new_daughterID, int new_daughterIdx, LorentzVector new_p4) {
    p4 += new_p4;
    daughterIDs.push_back(new_daughterID);
    daughterIdxs.push_back(new_daughterIdx);
    return;
}

float Products::Mass() {
    return p4.M();
}

string Products::Tag() {
    if (daughterIDs.size() == 0) return to_string(motherID)+"_EMPTY";
    int sum = 0;
    for (unsigned int i = 0; i < daughterIDs.size(); i++) {
        sum += daughterIDs.at(i);
    }
    return to_string(motherID)+"_"+to_string(sum);
}

struct Decay {
    // --> Member Variables <--
    // Decay identification
    int motherID;
    int motherIdx;
    vector<int> daughterIDs;
    Products* products;
    string tag;

    // Overload Constructor
    Decay(int, vector<int>, int);
    // Destructor
    ~Decay();
    // Methods
    void SetMotherIdx(int);
    string Tag();
    void Reset();
    void Add(int, int, LorentzVector);
    bool Full();
};

Decay::Decay(int new_motherID, vector<int> new_daughterIDs, int new_motherIdx = -1) {
    // Save mother/daughter PDG IDs
    motherID = new_motherID;
    motherIdx = new_motherIdx;
    daughterIDs = new_daughterIDs;
    // Make unique tag for system
    tag = Tag();
    // Initialize products
    products = new Products(new_motherID);
}

Decay::~Decay() {}

void Decay::SetMotherIdx(int new_motherIdx) {
    motherIdx = new_motherIdx;
    return;
}

string Decay::Tag() {
    if (daughterIDs.size() == 0) return to_string(motherID)+"_EMPTY";
    int sum = 0;
    for (unsigned int i = 0; i < daughterIDs.size(); i++) {
        sum += daughterIDs.at(i);
    }
    return to_string(motherID)+"_"+to_string(sum);
}

void Decay::Reset() {
    delete products;
    products  = new Products(motherID);
    return;
}

void Decay::Add(int new_daughterID, int new_daughterIdx, LorentzVector new_p4) {
    products->Add(new_daughterID, new_daughterIdx, new_p4);
    return;
}

bool Decay::Full() {
    return (products->Tag() == tag);
}
