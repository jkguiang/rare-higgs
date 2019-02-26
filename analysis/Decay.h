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
    daugherIDs.push_back(new_daughterID);
    return;
}

double System::Mass() {
    return p4.M();
}

string System::Tag() {
    string tag = to_string(motherID);
    for (unsigned int i = 0; i < daughterIDs.size(); i++) {
        tag += to_string(daughter.at(i));
    }
    return tag;
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
    void Add(int, LorentzVector, int);
    string MakeTag();
    bool Full(int);
    vector<System*> GetSystems();
    System* GetSystem(int);
};

Decay::Decay(int new_motherID, vector<int> new_daughterIDs) {
    motherID = new_motherID;
    daugherIDs = new_daughterIDs;
    tag = MakeTag();
}

Decay::~Decay() {}

string Decay::MakeTag() {
    int daughterSum = 0;
    for (unsigned int i = 0; i < daughterIDs.size(); i++) {
        daughterSum += daughter.at(i);
    }
    return to_string(motherID)+to_string(daughterSum);
}

void Decay::NewSystem() {
    System* sys = new System(motherID);
    systems.push_back(sys);
    return;
}

void Decay::Add(int new_daughterID, LorentzVector new_p4, int index = -1) {
    if (index == -1) {
        index = systems.size()-1;
    }
    System* sys = systems.at(index);
    sys->Add(new_daughterID, new_p4);
    return;
}

bool Decay::Full(int index = -1) {
    if (index == -1) {
        index = systems.size()-1;
    }
    return (systems.at(index)->Tag() == tag);
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
