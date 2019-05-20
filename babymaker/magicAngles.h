#include <vector>
#include "TH1F.h"
#include "TH2F.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TLorentzRotation.h"

/* Disclaimer:
 * This file is comprised of code written by 12 contributers, foremost
 * among them being Ulascan Sarica. The original repo containing their
 * work can be found here:
 * https://github.com/cms-analysis/HiggsAnalysis-ZZMatrixElement
 *
 * This standalone file was compiled by Nick Amin
 */

enum FermionMassRemoval{
    NoRemoval,
    ConserveDifermionMass,
    MomentumToEnergy,
    nFermionMassRemovalSchemes
};

bool forbidMassiveLeptons = true;
bool forbidMassiveJets = true;
FermionMassRemoval LeptonMassScheme = ConserveDifermionMass;
FermionMassRemoval JetMassScheme = ConserveDifermionMass;

void scaleMomentumToEnergy(const TLorentzVector& massiveJet, TLorentzVector& masslessJet, double mass){
  double energy, p3, newp3, ratio;
  energy = massiveJet.T();
  p3 = massiveJet.P();
  newp3 = sqrt(max(pow(energy, 2)-mass*fabs(mass), 0.));
  ratio = (p3>0. ? (newp3/p3) : 1.);
  masslessJet.SetXYZT(massiveJet.X()*ratio, massiveJet.Y()*ratio, massiveJet.Z()*ratio, energy);
}

bool isALepton(const int id){
  if (std::abs(id)==11 || std::abs(id)==13 || std::abs(id)==15) return true;
  else return false;
}
bool isAQuark(const int id){
  if (abs(id)<=6 && abs(id)>0) return true;
  else return false;
}
bool isANeutrino(const int id){
  if (std::abs(id)==12 || std::abs(id)==14 || std::abs(id)==16) return true;
  else return false;
}
bool isAGluon(const int id){
  if (std::abs(id)==21) return true;
  else return false;
}
bool isAnUnknownJet(const int id){
  if (id==0) return true;
  else return false;
}
bool isAJet(const int id){
  if (isAnUnknownJet(id) || isAQuark(id) || isAGluon(id)) return true;
  else return false;
}


void constrainedRemovePairMass(TLorentzVector& p1, TLorentzVector& p2, double m1, double m2){
  TLorentzVector nullFourVector(0, 0, 0, 0);
  TLorentzVector p1hat, p2hat;
  if (p1==nullFourVector || p2==nullFourVector) return;

  /***** shiftMass in C++ *****/
  TLorentzVector p12=p1+p2;
  TLorentzVector diffp2p1=p2-p1;
  double p1sq = p1.M2();
  double p2sq = p2.M2();
  double p1p2 = p1.Dot(p2);
  double m1sq = m1*fabs(m1);
  double m2sq = m2*fabs(m2);
  double p12sq = p12.M2();

  TLorentzVector avec=(p1sq*p2 - p2sq*p1 + p1p2*diffp2p1);
  double a = avec.M2();
  double b = (p12sq + m2sq - m1sq) * (pow(p1p2, 2) - p1sq*p2sq);
  double c = pow((p12sq + m2sq - m1sq), 2)*p1sq/4. - pow((p1sq + p1p2), 2)*m2sq;
  double eta =  (-b - sqrt(fabs(b*b -4.*a*c)))/(2.*a);
  double xi = (p12sq + m2sq - m1sq - 2.*eta*(p2sq + p1p2))/(2.*(p1sq + p1p2));

  p1hat = (1.-xi)*p1 + (1.-eta)*p2;
  p2hat = xi*p1 + eta*p2;
  p1=p1hat;
  p2=p2hat;
}

pair<TLorentzVector, TLorentzVector> removeMassFromPair(
  TLorentzVector const& jet1, int const& jet1Id,
  TLorentzVector const& jet2, int const& jet2Id,
  double m1=0, double m2=0
  ){
  TLorentzVector nullFourVector(0, 0, 0, 0);
  TLorentzVector jet1massless(0, 0, 0, 0), jet2massless(0, 0, 0, 0);

  if (forbidMassiveJets && (isAJet(jet1Id) || isAJet(jet2Id))){
    if (JetMassScheme==FermionMassRemoval::NoRemoval){
      jet1massless=jet1;
      jet2massless=jet2;
    }
    else if (jet1==nullFourVector || jet2==nullFourVector || jet1==jet2 || JetMassScheme==MomentumToEnergy){
      scaleMomentumToEnergy(jet1, jet1massless, m1);
      scaleMomentumToEnergy(jet2, jet2massless, m2);
    }
    else if (JetMassScheme==ConserveDifermionMass){
      jet1massless=jet1;
      jet2massless=jet2;
      constrainedRemovePairMass(jet1massless, jet2massless, m1, m2);
    }
    else{
      jet1massless=jet1;
      jet2massless=jet2;
    }
  }
  else if (forbidMassiveLeptons && (isALepton(jet1Id) || isANeutrino(jet1Id) || isALepton(jet2Id) || isANeutrino(jet2Id))){
    if (LeptonMassScheme==FermionMassRemoval::NoRemoval){
      jet1massless=jet1;
      jet2massless=jet2;
    }
    else if (jet1==nullFourVector || jet2==nullFourVector || jet1==jet2 || LeptonMassScheme==MomentumToEnergy){
      scaleMomentumToEnergy(jet1, jet1massless, m1);
      scaleMomentumToEnergy(jet2, jet2massless, m2);
    }
    else if (LeptonMassScheme==ConserveDifermionMass){
      jet1massless=jet1;
      jet2massless=jet2;
      constrainedRemovePairMass(jet1massless, jet2massless, m1, m2);
    }
    else{
      jet1massless=jet1;
      jet2massless=jet2;
    }
  }
  else{
    jet1massless=jet1;
    jet2massless=jet2;
  }

  pair<TLorentzVector, TLorentzVector> result(jet1massless, jet2massless);
  return result;
}

void computeAngles(
  float& costhetastar,
  float& costheta1,
  float& costheta2,
  float& Phi,
  float& Phi1,
  TLorentzVector p4M11, int Z1_lept1Id,
  TLorentzVector p4M12, int Z1_lept2Id,
  TLorentzVector p4M21, int Z2_lept1Id,
  TLorentzVector p4M22, int Z2_lept2Id,
  bool verbose=false
  ){
  TLorentzVector nullFourVector(0, 0, 0, 0);
  if (p4M12==nullFourVector || p4M22==nullFourVector){
      pair<TLorentzVector, TLorentzVector> f13Pair = removeMassFromPair(p4M11, Z1_lept1Id, p4M21, Z2_lept1Id);
    p4M11 = f13Pair.first;
    p4M21 = f13Pair.second;
  }
  else if (p4M11==nullFourVector || p4M21==nullFourVector){
      pair<TLorentzVector, TLorentzVector> f24Pair = removeMassFromPair(p4M12, Z1_lept2Id, p4M22, Z2_lept2Id);
    p4M12 = f24Pair.first;
    p4M22 = f24Pair.second;
  }
  else{
      pair<TLorentzVector, TLorentzVector> f12Pair = removeMassFromPair(p4M11, Z1_lept1Id, p4M12, Z1_lept2Id);
      pair<TLorentzVector, TLorentzVector> f34Pair = removeMassFromPair(p4M21, Z2_lept1Id, p4M22, Z2_lept2Id);
    p4M11 = f12Pair.first;
    p4M12 = f12Pair.second;
    p4M21 = f34Pair.first;
    p4M22 = f34Pair.second;
  }

  //build Z 4-vectors
  TLorentzVector p4Z1 = p4M11 + p4M12;
  TLorentzVector p4Z2 = p4M21 + p4M22;

  // Sort Z1 leptons so that:
  if (
    Z1_lept2Id!=-9000
    &&
    (
    (Z1_lept1Id*Z1_lept2Id<0 && Z1_lept1Id<0) // for OS pairs: lep1 must be the negative one
    ||
    ((Z1_lept1Id*Z1_lept2Id>0 || (Z1_lept1Id==0 && Z1_lept2Id==0)) && p4M11.Phi()<=p4M12.Phi()) //for SS pairs: use random deterministic convention
    )
    ) swap(p4M11, p4M12);
  // Same for Z2 leptons
  if (
    Z2_lept2Id!=-9000
    &&
    (
    (Z2_lept1Id*Z2_lept2Id<0 && Z2_lept1Id<0)
    ||
    ((Z2_lept1Id*Z2_lept2Id>0 || (Z2_lept1Id==0 && Z2_lept2Id==0)) && p4M21.Phi()<=p4M22.Phi())
    )
    ) swap(p4M21, p4M22);

  // BEGIN THE CALCULATION

  // build H 4-vectors
  TLorentzVector p4H = p4Z1 + p4Z2;

  // -----------------------------------

  //// costhetastar
  TVector3 boostX = -(p4H.BoostVector());
  TLorentzVector thep4Z1inXFrame(p4Z1);
  TLorentzVector thep4Z2inXFrame(p4Z2);
  thep4Z1inXFrame.Boost(boostX);
  thep4Z2inXFrame.Boost(boostX);
  TVector3 theZ1X_p3 = TVector3(thep4Z1inXFrame.X(), thep4Z1inXFrame.Y(), thep4Z1inXFrame.Z());
  TVector3 theZ2X_p3 = TVector3(thep4Z2inXFrame.X(), thep4Z2inXFrame.Y(), thep4Z2inXFrame.Z());
  costhetastar = theZ1X_p3.CosTheta();

  TVector3 boostV1(0, 0, 0);
  TVector3 boostV2(0, 0, 0);
  //// --------------------------- costheta1
  if (!(fabs(Z1_lept1Id)==21 || fabs(Z1_lept1Id)==22 || fabs(Z1_lept2Id)==21 || fabs(Z1_lept2Id)==22)){
    boostV1 = -(p4Z1.BoostVector());
    if (boostV1.Mag()>=1. && verbose) {
        std::cout << "Warning: computeAngles: Z1 boost with beta=1, scaling down" << std::endl;
      boostV1*=0.9999/boostV1.Mag();
    }
    TLorentzVector p4M11_BV1(p4M11);
    TLorentzVector p4M12_BV1(p4M12);
    TLorentzVector p4M21_BV1(p4M21);
    TLorentzVector p4M22_BV1(p4M22);
    p4M11_BV1.Boost(boostV1);
    p4M12_BV1.Boost(boostV1);
    p4M21_BV1.Boost(boostV1);
    p4M22_BV1.Boost(boostV1);

    TLorentzVector p4V2_BV1 = p4M21_BV1 + p4M22_BV1;
    //// costheta1
    costheta1 = -p4V2_BV1.Vect().Unit().Dot(p4M11_BV1.Vect().Unit());
  }
  else costheta1 = 0;

  //// --------------------------- costheta2
  if (!(fabs(Z2_lept1Id)==21 || fabs(Z2_lept1Id)==22 || fabs(Z2_lept2Id)==21 || fabs(Z2_lept2Id)==22)){
    boostV2 = -(p4Z2.BoostVector());
    if (boostV2.Mag()>=1. && verbose) {
        std::cout << "Warning: computeAngles: Z2 boost with beta=1, scaling down" << std::endl;
      boostV2*=0.9999/boostV2.Mag();
    }
    TLorentzVector p4M11_BV2(p4M11);
    TLorentzVector p4M12_BV2(p4M12);
    TLorentzVector p4M21_BV2(p4M21);
    TLorentzVector p4M22_BV2(p4M22);
    p4M11_BV2.Boost(boostV2);
    p4M12_BV2.Boost(boostV2);
    p4M21_BV2.Boost(boostV2);
    p4M22_BV2.Boost(boostV2);

    TLorentzVector p4V1_BV2 = p4M11_BV2 + p4M12_BV2;
    //// costheta2
    costheta2 = -p4V1_BV2.Vect().Unit().Dot(p4M21_BV2.Vect().Unit());
  }
  else costheta2 = 0;

  //// --------------------------- Phi and Phi1 (old phistar1 - azimuthal production angle)
  TLorentzVector p4M11_BX(p4M11);
  TLorentzVector p4M12_BX(p4M12);
  TLorentzVector p4M21_BX(p4M21);
  TLorentzVector p4M22_BX(p4M22);

  p4M11_BX.Boost(boostX);
  p4M12_BX.Boost(boostX);
  p4M21_BX.Boost(boostX);
  p4M22_BX.Boost(boostX);
  TLorentzVector p4V1_BX = p4M11_BX + p4M12_BX;

  TVector3 beamAxis(0, 0, 1);
  TVector3 p3V1_BX = p4V1_BX.Vect().Unit();
  TVector3 normal1_BX = (p4M11_BX.Vect().Cross(p4M12_BX.Vect())).Unit();
  TVector3 normal2_BX = (p4M21_BX.Vect().Cross(p4M22_BX.Vect())).Unit();
  TVector3 normalSC_BX = (beamAxis.Cross(p3V1_BX)).Unit();


  //// Phi
  float tmpSgnPhi = p3V1_BX.Dot(normal1_BX.Cross(normal2_BX));
  float sgnPhi = 0;
  if (fabs(tmpSgnPhi)>0.) sgnPhi = tmpSgnPhi/fabs(tmpSgnPhi);
  float dot_BX12 = normal1_BX.Dot(normal2_BX);
  if (fabs(dot_BX12)>=1.) dot_BX12 *= 1./fabs(dot_BX12);
  Phi = sgnPhi * acos(-1.*dot_BX12);


  //// Phi1
  float tmpSgnPhi1 = p3V1_BX.Dot(normal1_BX.Cross(normalSC_BX));
  float sgnPhi1 = 0;
  if (fabs(tmpSgnPhi1)>0.) sgnPhi1 = tmpSgnPhi1/fabs(tmpSgnPhi1);
  float dot_BX1SC = normal1_BX.Dot(normalSC_BX);
  if (fabs(dot_BX1SC)>=1.) dot_BX1SC *= 1./fabs(dot_BX1SC);
  Phi1 = sgnPhi1 * acos(dot_BX1SC);

  if ((isnan(costhetastar) || isnan(costheta1) || isnan(costheta2) || isnan(Phi) || isnan(Phi1)) && (verbose)){
    std::cout << "WARNING: NaN in computeAngles: "
    << costhetastar << " "
    << costheta1  << " "
    << costheta2  << " "
    << Phi  << " "
    << Phi1  << " " << std::endl;
    std::cout << "   boostV1: " <<boostV1.Pt() << " " << boostV1.Eta() << " " << boostV1.Phi() << " " << boostV1.Mag() << std::endl;
    std::cout << "   boostV2: " <<boostV2.Pt() << " " << boostV2.Eta() << " " << boostV2.Phi() << " " << boostV2.Mag() << std::endl;
  }
}

void computeVHAngles(
  float& costhetastar,
  float& costheta1,
  float& costheta2,
  float& Phi,
  float& Phi1,
  float& m1,
  float& m2,
  TLorentzVector p4M11, int Z1_lept1Id,
  TLorentzVector p4M12, int Z1_lept2Id,
  TLorentzVector p4M21, int Z2_lept1Id,
  TLorentzVector p4M22, int Z2_lept2Id,
  TLorentzVector jet1, int jet1Id,
  TLorentzVector jet2, int jet2Id,
  TLorentzVector* injet1, int injet1Id, // Gen. partons in lab frame
  TLorentzVector* injet2, int injet2Id
  ){
  TLorentzVector nullFourVector(0, 0, 0, 0);
  if (p4M12==nullFourVector || p4M22==nullFourVector){
    pair<TLorentzVector, TLorentzVector> f13Pair = removeMassFromPair(p4M11, Z1_lept1Id, p4M21, Z2_lept1Id);
    p4M11 = f13Pair.first;
    p4M21 = f13Pair.second;
  }
  else if (p4M11==nullFourVector || p4M21==nullFourVector){
    pair<TLorentzVector, TLorentzVector> f24Pair = removeMassFromPair(p4M12, Z1_lept2Id, p4M22, Z2_lept2Id);
    p4M12 = f24Pair.first;
    p4M22 = f24Pair.second;
  }
  else{
    pair<TLorentzVector, TLorentzVector> f12Pair = removeMassFromPair(p4M11, Z1_lept1Id, p4M12, Z1_lept2Id);
    pair<TLorentzVector, TLorentzVector> f34Pair = removeMassFromPair(p4M21, Z2_lept1Id, p4M22, Z2_lept2Id);
    p4M11 = f12Pair.first;
    p4M12 = f12Pair.second;
    p4M21 = f34Pair.first;
    p4M22 = f34Pair.second;
  }

  TLorentzVector jet1massless, jet2massless;
  pair<TLorentzVector, TLorentzVector> jetPair = removeMassFromPair(jet1, jet1Id, jet2, jet2Id);
  jet1massless = jetPair.first;
  jet2massless = jetPair.second;

  // Build Z 4-vectors
  TLorentzVector p4Z1 = p4M11 + p4M12;
  TLorentzVector p4Z2 = p4M21 + p4M22;
  TLorentzVector pH = p4Z1 + p4Z2;

  // Apply convention for outgoing particles
  if (
    (jet1Id*jet2Id<0 && jet1Id<0) // for OS pairs: jet1 must be the particle
    ||
    (jet1Id*jet2Id>0 && jet1massless.Phi()<=jet2massless.Phi()) // for SS pairs: use random deterministic convention
    ){
    swap(jet1massless, jet2massless);
    swap(jet1Id, jet2Id);
  }

  //Find the incoming partons - first boost so the pT(HJJ) = 0, then boost away the pz.
  //This preserves the z direction.  Then assume that the partons come in this z direction.
  //This is exactly correct at LO (since pT=0 anyway).
  //Then associate the one going forwards with jet1 and the one going backwards with jet2
  TLorentzRotation movingframe;
  TLorentzVector pHJJ = pH+jet1massless+jet2massless;
  TLorentzVector pHJJ_perp(pHJJ.X(), pHJJ.Y(), 0, pHJJ.T());
  movingframe.Boost(-pHJJ_perp.BoostVector());
  pHJJ.Boost(-pHJJ_perp.BoostVector());
  movingframe.Boost(-pHJJ.BoostVector());
  pHJJ.Boost(-pHJJ.BoostVector());   //make sure to boost HJJ AFTER boosting movingframe

  TLorentzVector P1(0, 0, -pHJJ.T()/2, pHJJ.T()/2);
  TLorentzVector P2(0, 0, pHJJ.T()/2, pHJJ.T()/2);
  // Transform incoming partons back to the original frame
  P1.Transform(movingframe.Inverse());
  P2.Transform(movingframe.Inverse());
  pHJJ.Transform(movingframe.Inverse());
  // movingframe and HJJ_T will not be used anymore
  // Handle gen. partons if they are available
  if (injet1 && injet2 && fabs((*injet1+*injet2).P()-pHJJ.P())<pHJJ.P()*1e-4){
    P1=*injet1;
    P2=*injet2;
    // Apply convention for incoming (!) particles
    if (
      (injet1Id*injet2Id<0 && injet1Id>0) // for OS pairs: parton 2 must be the particle
      ||
      (injet1Id*injet2Id>0 && P1.Z()>=P2.Z()) //for SS pairs: use random deterministic convention
      ){
      swap(P1, P2);
      swap(injet1Id, injet2Id);
    }
  }

  // Rotate every vector such that Z1 - Z2 axis is the "beam axis" analogue of decay
  TLorentzRotation ZZframe;
  TVector3 beamAxis(0, 0, 1);
  if (p4Z1==nullFourVector || p4Z2==nullFourVector){
    TVector3 pNewAxis = (p4Z2-p4Z1).Vect().Unit(); // Let Z2 be in the z direction so that once the direction of H is reversed, Z1 is in the z direction
    if (pNewAxis != nullFourVector.Vect()){
      TVector3 pNewAxisPerp = pNewAxis.Cross(beamAxis);
      ZZframe.Rotate(acos(pNewAxis.Dot(beamAxis)), pNewAxisPerp);

      P1.Transform(ZZframe);
      P2.Transform(ZZframe);
      jet1massless = -jet1massless;
      jet2massless = -jet2massless;
      jet1massless.Transform(ZZframe);
      jet2massless.Transform(ZZframe);
      jet1massless = -jet1massless;
      jet2massless = -jet2massless;
    }
  }
  else{
    ZZframe.Boost(-pH.BoostVector());
    p4Z1.Boost(-pH.BoostVector());
    p4Z2.Boost(-pH.BoostVector());
    TVector3 pNewAxis = (p4Z2-p4Z1).Vect().Unit(); // Let Z2 be in the z direction so that once the direction of H is reversed, Z1 is in the z direction
    TVector3 pNewAxisPerp = pNewAxis.Cross(beamAxis);
    ZZframe.Rotate(acos(pNewAxis.Dot(beamAxis)), pNewAxisPerp);
    P1.Transform(ZZframe);
    P2.Transform(ZZframe);
    jet1massless = -jet1massless;
    jet2massless = -jet2massless;
    jet1massless.Transform(ZZframe);
    jet2massless.Transform(ZZframe);
    jet1massless = -jet1massless;
    jet2massless = -jet2massless;
  }

  computeAngles(
    costhetastar,
    costheta1,
    costheta2,
    Phi,
    Phi1,
    -P1, 23, // Id is 23 to avoid an attempt to remove quark mass
    -P2, 0, // Id is 0 to avoid swapping
    jet1massless, 23,
    jet2massless, 0
  );
  m1 = (P1 + P2).M();
  m2 = (jet1massless + jet2massless).M();
}

/* Disclaimer (continued):
 * The code below was written by Jonathan Guiang
 */

struct MagicAngles {
    float angles[7];
};

MagicAngles getAngles(float met_pt, float met_phi, int lep_id, TLorentzVector lep_p4, TLorentzVector mes_p4, TLorentzVector gam_p4, float nu_pz=-999) {
    // Known neutrino 3-momentum components
    float nu_px = met_pt*cos(met_phi);
    float nu_py = met_pt*sin(met_phi);
    // Best lepton 3-momentum
    float lep_px = lep_p4.Px();
    float lep_py = lep_p4.Py();
    float lep_pz = lep_p4.Pz();
    float E_l = pow(pow(lep_px, 2)+pow(lep_py, 2)+pow(lep_pz, 2), 0.5); // m_l ~ 0
    // Mass of W-boson
    float m = 80.379;
    // Solve for nu_pz if not provided
    if (nu_pz < -998) {
        // Terms of quadratic equation solutions for nu_pz
        float a = ( pow(E_l,2) - pow(lep_pz,2) );
        float b = ( -lep_pz*pow(m,2) - 2*lep_px*lep_pz*nu_px - 2*lep_py*lep_pz*nu_py );
        float c = ( pow(E_l,2)*pow(nu_py,2) - pow(lep_py,2)*pow(nu_py,2) - pow(m,2)*lep_py*nu_py + pow(E_l,2)*pow(nu_px,2) 
                    - pow(lep_px,2)*pow(nu_px,2) - pow(m,2)*lep_px*nu_px -2*lep_px*lep_py*nu_px*nu_py - pow(m,4)/4 );
        float disc = ( pow(b, 2) - 4*a*c );
        if (disc < 0) disc = 0.0;
        // Determine nu_pz
        float nu_pz = 0.0;
        float quad_p = (-b+pow(disc, 0.5))/(2*a);
        float quad_m = (-b-pow(disc, 0.5))/(2*a);
        nu_pz = (abs(quad_p) < abs(quad_m)) ? quad_p : quad_m;
    }
    // MELA output params
    float costhetastar = -999;
    float costheta1 = -999;
    float costheta2 = -999;
    float Phi = -999;
    float Phi1 = -999;
    float m1 = -999;
    float m2 = -999;
    // MELA input params
    TLorentzVector nu_p4(nu_px, nu_py, nu_pz, pow(pow(nu_px, 2)+pow(nu_py, 2)+pow(nu_pz, 2), 0.5)); // m_v ~ 0
    TLorentzVector dummy_p4(0,0,0,0);
    TLorentzVector jet1 = (lep_id > 0) ? lep_p4 : nu_p4;
    int jet1Id = (lep_id > 0) ? lep_id : -lep_id+1;
    TLorentzVector jet2 = (lep_id > 0) ? nu_p4 : lep_p4;
    int jet2Id = (lep_id > 0) ? -lep_id-1 : lep_id;
    // Get angles
    computeVHAngles(
            costhetastar,
            costheta1,
            costheta2,
            Phi,
            Phi1,
            m1,
            m2,
            mes_p4, 23,
            dummy_p4, -9000,
            gam_p4, 23,
            dummy_p4, -9000,
            jet1, jet1Id,
            jet2, jet2Id,
            nullptr, -9000,
            nullptr, -9000
            );

    // Fill return struct
    MagicAngles mAngles;
    mAngles.angles[0] = costhetastar;
    mAngles.angles[1] = costheta1;
    mAngles.angles[2] = costheta2;
    mAngles.angles[3] = Phi;
    mAngles.angles[4] = Phi1;
    mAngles.angles[5] = m1;
    mAngles.angles[6] = m2;

    return mAngles;
}
