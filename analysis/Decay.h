// C++
#include <iostream>
#include <vector>
#include <string>

// CMS3
#include "/home/users/jguiang/projects/mt2/MT2Analysis/CORE/CMS3.h"

using namespace std;
using namespace tas;

struct System {
    // Member Variables
    int motherID;
    vector<int> daughterIDs;
    LorentzVector p4;
    // Overload Constructor
    System(int);
    // Destructor
    ~System();
    // Methods
    void Add(int, LorentzVector);
    double Mass();
    string Tag();
};

System::System(int new_motherID) {
    motherID = new_motherID;
}

System::~System() {}

void System::Add(int new_daughterID, LorentzVector new_p4) {
    p4 += new_p4;
    daughterIDs.push_back(new_daughterID);
    return;
}

double System::Mass() {
    return p4.M();
}

string System::Tag() {
    if (daughterIDs.size() == 0) return to_string(motherID)+"EMPTY";
    int daughterSum = 0;
    for (unsigned int i = 0; i < daughterIDs.size(); i++) {
        daughterSum += daughterIDs.at(i);
    }
    return to_string(motherID)+"_"+to_string(daughterSum);
}

struct Decay {
    // Member Variables
    int motherID;
    vector<int> daughterIDs;
    vector<System*> systems;
    string tag;

    // Overload Constructor
    Decay(int, vector<int>);
    // Destructor
    ~Decay();
    // Methods
    void NewSystem();
    vector<System*> GetSystems();
    System* GetSystem(int);
    void Add(int, LorentzVector, int);
    void Pop();
    string MakeTag();
    bool Full(int);
};

Decay::Decay(int new_motherID, vector<int> new_daughterIDs) {
    motherID = new_motherID;
    daughterIDs = new_daughterIDs;
    tag = MakeTag();
}

Decay::~Decay() {}

string Decay::MakeTag() {
    if (daughterIDs.size() == 0) return to_string(motherID)+"EMPTY";
    int daughterSum = 0;
    for (unsigned int i = 0; i < daughterIDs.size(); i++) {
        daughterSum += daughterIDs.at(i);
    }
    return to_string(motherID)+"_"+to_string(daughterSum);
}

void Decay::NewSystem() {
    System* sys = new System(motherID);
    systems.push_back(sys);
    return;
}

vector<System*> Decay::GetSystems() {
    return systems;
}

System* Decay::GetSystem(int index = -1) {
    if (index == -1) {
        index = systems.size()-1;
    }
    return systems.at(index);
}

void Decay::Add(int new_daughterID, LorentzVector new_p4, int index = -1) {
    if (index == -1) {
        index = systems.size()-1;
    }
    System* sys = systems.at(index);
    sys->Add(new_daughterID, new_p4);
    return;
}

void Decay::Pop() {
    System* sys = systems.at(systems.size()-1);
    systems.pop_back();
    delete sys;
    return;
}

bool Decay::Full(int index = -1) {
    if (index == -1) {
        index = systems.size()-1;
    }
    return (systems.at(index)->Tag() == tag);
}
