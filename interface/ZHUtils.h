#ifndef ZHUtils_H
#define ZHUtils_H

/** \class ZHUtils
 *  No description available.
 *
 *  \author R.-J. Wang 
 */

#include <iostream>
#include <TString.h>
#include <map>

#include <TSystem.h>
#include <TROOT.h>
#include <TMath.h>
#include <TLorentzVector.h>
#include <TRandom3.h>
#include <utility>
#include <vector>

#include "TFile.h"
#include "TH1.h"
#include "TGraph.h"
#include "TH2.h"

#include "Math/LorentzVector.h"
#include "TVector2.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "CMGTools/HiggsAna2l2v/interface/ZZ2l2nuPhysicsEvent.h"

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector;
typedef std::vector<LorentzVector> LorentzVectorCollection;


class ZHUtils
{
 public:
//namespace ZHUtils{

  ZHUtils(const edm::ParameterSet &runProcess);
  double Collins_Soper(const LorentzVector& lepton1, const LorentzVector& lepton2); 
  float weightNLOEWKzh(float pt);
  float weightNLOEWKzz(float pt);
  float weightNLOEWKwz(float pt);
  //float weightNLOEWKwplusz(float pt);
  //float weightNLOEWKwminuz(float pt); 

  double GetNLOZHWeight(const PhysicsEvent_t& phys);
  double GetNLOZZWeight(const PhysicsEvent_t& phys);
  double GetNLOWZWeight(const PhysicsEvent_t& phys);


  //void ZHEventHandler(const edm::ParameterSet &runProcess);
  std::map<TString,float> getWeights(double ValtoWeight, TString wgtName);
  double get2DWeights(double Val_x, double Val_y, TString wgtName, TString cat);


  // prompt rate and fake rate
  void get_frFile(const edm::ParameterSet &runProcess);
  double promptRate(int pdgid, double pt, double abseta);
  double fakeRate(int pdgid, double pt, double abseta, TString key);

  double getN_PFweight(int TL_type, LorentzVector lep1, int id1, LorentzVector lep2, int id2, TString key);
  double getN_FPweight(int TL_type, LorentzVector lep1, int id1, LorentzVector lep2, int id2, TString key);
  double getN_FFweight(int TL_type, LorentzVector lep1, int id1, LorentzVector lep2, int id2, TString key);

  ~ZHUtils();

 private:
  std::map<TString, std::map<TString,TGraph *> > wgtsH_;
  std::map<TString, std::map<TString,TH2F *> > wgts2DH_;

  std::map<TString, TH1F *> fakerate1DH_;

};
#endif

