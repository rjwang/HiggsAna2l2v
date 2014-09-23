/*
 * Wrapper for common operations on a gamma event
 * Get weights/mass shapes from file
 * Analyze event and assign trigger categories, weights and massive candidates
 * $Date: 2013/11/02 21:28:45 $
 * $Revision: 1.0$
 * \author Ren-Jie Wang
 */

#ifndef ZHphotonjetsEventHandler_H
#define ZHphotonjetsEventHandler_H

#include "TString.h"
#include "TSystem.h"
#include "TFile.h"
#include "TH1.h"
#include "TGraph.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "CMGTools/HiggsAna2l2v/interface/ZZ2l2nuPhysicsEvent.h"

struct BosonPlusJets_t {
    Float_t en,qt,eta, weight, mindrbj, minmbj,njets,nvtx,cat;
};


class ZHphotonjetsEventHandler {
public:

    ZHphotonjetsEventHandler(const edm::ParameterSet &runProcess);

    TTree *initSummary() {
        t_=new TTree("data","boson+jets summary");
        t_->Branch("cat",     &evSummary_.cat,      "cat/F");
        t_->Branch("qt",      &evSummary_.qt,       "qt/F");
        t_->Branch("eta",      &evSummary_.eta,       "eta/F");
        t_->Branch("en",      &evSummary_.en,       "en/F");
        t_->Branch("weight",  &evSummary_.weight,   "weight/F");
        t_->Branch("mindrbj", &evSummary_.mindrbj,  "mindrbj/F");
        t_->Branch("minmbj",  &evSummary_.minmbj,   "minmbj/F");
        t_->Branch("njets",   &evSummary_.njets,    "njets/F");
        t_->Branch("nvtx",    &evSummary_.nvtx,     "nvtx/F");
        return t_;
    }

    //to be called event by event
    bool isGood(PhysicsEvent_t &phys,bool is2011);

    //getters
    float triggerThr() {
        return triggerThr_;
    }
    float triggerWeight() {
        return triggerWgt_;
    }
    bool isGood() {
        return isGoodEvent_;
    }
    LorentzVector massiveGamma(TString channel) {
        return massiveGamma_[channel];
    }
    float dphimassiveGamma(TString channel) {
        return dphimassiveGamma_[channel];
    }
    float dphi2lepGamma(TString channel) {
        return dphi2lepGamma_[channel];
    }    
    float coslzmassiveGamma(TString channel) {
        return coslzmassiveGamma_[channel];
    }
    std::map<TString,LorentzVector> getMassiveGamma() {
        return massiveGamma_;
    }
    float getWeight(TString channel) {
        return evWeights_[channel];
    }

    std::map<TString,float> getWeights(ZZ2l2nuSummary_t &ev, PhysicsEvent_t &phys, TString evCategoryLabel, double _metphi_);
    std::map<TString,double> getQtFitWeights(PhysicsEvent_t &phys, TString evCategoryLabel);
    std::map<TString,double> get2011QtFitWeights(PhysicsEvent_t &phys, TString evCategoryLabel);
    double getQtFit(double *x, double *par);

    ~ZHphotonjetsEventHandler();


    //these can be accessed after call isGood(physics)
    bool isGoodEvent_;
    float triggerThr_,triggerWgt_;
    std::map<TString,LorentzVector> massiveGamma_;
    std::map<TString,float> dphimassiveGamma_;
    std::map<TString,float> coslzmassiveGamma_;
    std::map<TString,float> evWeights_;
    std::map<TString,double> QtFitWeights_;
    std::map<TString,float> dphi2lepGamma_;

    BosonPlusJets_t evSummary_;


private:

    TTree *t_;

    bool isMC_;
    /*std::map<TString, std::map<TString,TGraph *> > wgtsH_;*/
    std::map<TString, std::map<TString, TH1 *> > WgtsH_;
    std::map<TString, TH1 *> zmassH_;
    std::map<TString, TH1 *> dphi2leptonH_;
    std::map<TString, TH1 *> qCosl1Z_CSH_;
    std::map<TString, TGraph *> qGammaDataH_;
    std::map<TString, TH2 *> dphi2lepvsqtH_;

    static double par0jets_mmqt[4];
    static double par0jets_eeqt[4];
    static double par0jets_gqt[4];
    static double par0jets_mcgqt[4];
    static double par1jets_mmqt[4];
    static double par1jets_eeqt[4];
    static double par1jets_gqt[4];
    static double par1jets_mcgqt[4];
    
    static double par0jets_mmqt_2011[4];
    static double par0jets_eeqt_2011[4];
    static double par0jets_gqt_2011[4]; //not use
    static double par0jets_mcgqt_2011[4];
    static double par1jets_mmqt_2011[4];
    static double par1jets_eeqt_2011[4];
    static double par1jets_gqt_2011[4];//not use
    static double par1jets_mcgqt_2011[4];

};


// Nov.18 0jets
/*
double ZHphotonjetsEventHandler::par0jets_mmqt[4] = {173.528,-6.77964,19.9762,2.64192};
double ZHphotonjetsEventHandler::par0jets_eeqt[4] = {144.236,-7.32295,20,3.04513};
double ZHphotonjetsEventHandler::par0jets_gqt[4] = {394.111,-6.07104,-2.67577,27.0425};
double ZHphotonjetsEventHandler::par0jets_mcgqt[4] = {399.479,-5.90795,11.0101,12.2024};

double ZHphotonjetsEventHandler::par1jets_mmqt[4] = {542.596,-4.44669,50.4502,16.101};
double ZHphotonjetsEventHandler::par1jets_eeqt[4] = {515.063,-4.28677,36.9853,17.1138};
double ZHphotonjetsEventHandler::par1jets_gqt[4] = {831.514,-5.15587,5.2404,24.3906};
double ZHphotonjetsEventHandler::par1jets_mcgqt[4] = {748.628,-6.06733,20.0528,60.1906};
*/


//Mono-Z, July 09, 2014 
//combine fit and ratio as weight
//
/*
double ZHphotonjetsEventHandler::par0jets_mmqt[4] = {230.801,-4.64064,-1.59745,59.1334};
double ZHphotonjetsEventHandler::par0jets_eeqt[4] = {214.02,-5.10411,20,1.26653};
double ZHphotonjetsEventHandler::par0jets_gqt[4] = {432.336,-5.12048,3.53363e+09,5.97489};
double ZHphotonjetsEventHandler::par0jets_mcgqt[4] = {420.258,-5.12143,20,3.54621};

double ZHphotonjetsEventHandler::par1jets_mmqt[4] = {468.825,-5.38642,20,3.35};
double ZHphotonjetsEventHandler::par1jets_eeqt[4] = {449.804,-5.26356,20,3.55473};
double ZHphotonjetsEventHandler::par1jets_gqt[4] = {659.43,-7.50765,82.0926,62.5261};
double ZHphotonjetsEventHandler::par1jets_mcgqt[4] = {698.496,-7.07511,53.8569,65.1831};
*/

// Mono-Z, Aug 04, 2014
// VV backgrounds are subtracted
/*
double ZHphotonjetsEventHandler::par0jets_mmqt[4] = {221.221,-4.43262,-1.35128,97.8416};
double ZHphotonjetsEventHandler::par0jets_eeqt[4] = {206.381,-5.26251,20,1.17532};
double ZHphotonjetsEventHandler::par0jets_gqt[4] = {432.336,-5.12048,3.53363e+09,5.97489};
double ZHphotonjetsEventHandler::par0jets_mcgqt[4] = {420.258,-5.12143,20,3.54621};

double ZHphotonjetsEventHandler::par1jets_mmqt[4] = {443.872,-7.56699,212.922,49.1785};
double ZHphotonjetsEventHandler::par1jets_eeqt[4] = {449.418,-5.26432,15,2.80591};
double ZHphotonjetsEventHandler::par1jets_gqt[4] = {659.43,-7.50765,82.0926,62.5261};
double ZHphotonjetsEventHandler::par1jets_mcgqt[4] = {698.496,-7.07511,53.8569,65.1831};
*/

// Mono-Z, Aug 05, 2014
// VV backgrounds are subtracted
/*
double ZHphotonjetsEventHandler::par0jets_mmqt[4] = {221.221,-4.43262,-1.35128,97.8416};
double ZHphotonjetsEventHandler::par0jets_eeqt[4] = {192.719,-5.74894,20,0.959371};
double ZHphotonjetsEventHandler::par0jets_gqt[4] = {432.336,-5.12048,3.53363e+09,5.97489};
double ZHphotonjetsEventHandler::par0jets_mcgqt[4] = {420.258,-5.12143,20,3.54621};

double ZHphotonjetsEventHandler::par1jets_mmqt[4] = {464.793,-5.51274,15,7.84662};
double ZHphotonjetsEventHandler::par1jets_eeqt[4] = {448.923,-5.27413,15,3.35496};
double ZHphotonjetsEventHandler::par1jets_gqt[4] = {659.43,-7.50765,82.0926,62.5261};
double ZHphotonjetsEventHandler::par1jets_mcgqt[4] = {698.496,-7.07511,53.8569,65.1831};
*/

// Mono-Z, Aug 17, 2014
// VV backgrounds are subtracted, Z mass window 10GeV
double ZHphotonjetsEventHandler::par0jets_mmqt[4] = {201.73,-5.87057,10,0.989512};
double ZHphotonjetsEventHandler::par0jets_eeqt[4] = {191.875,-5.73781,20,0.148266};
double ZHphotonjetsEventHandler::par0jets_gqt[4] = {432.336,-5.12048,3.53363e+09,5.97489};
double ZHphotonjetsEventHandler::par0jets_mcgqt[4] = {420.258,-5.12143,20,3.54621};

double ZHphotonjetsEventHandler::par1jets_mmqt[4] = {461.541,-5.50603,15,3.16613};
double ZHphotonjetsEventHandler::par1jets_eeqt[4] = {445.851,-5.25844,15,3.0958};
double ZHphotonjetsEventHandler::par1jets_gqt[4] = {659.43,-7.50765,82.0926,62.5261};
double ZHphotonjetsEventHandler::par1jets_mcgqt[4] = {698.496,-7.07511,53.8569,65.1831};

//####################################

//7TeV, Nov20
double ZHphotonjetsEventHandler::par0jets_mmqt_2011[4] = {130.696,-7.2158,20,0.587308};
double ZHphotonjetsEventHandler::par0jets_eeqt_2011[4] = {121.074,-7.10006,20,0.346779};
double ZHphotonjetsEventHandler::par0jets_gqt_2011[4] = {2,-1,20,1};//not use
double ZHphotonjetsEventHandler::par0jets_mcgqt_2011[4] = {413.693,-5.68578,-2.11051,47.7309};

double ZHphotonjetsEventHandler::par1jets_mmqt_2011[4] = {391.674,-4.3801,56.069,14.6251};
double ZHphotonjetsEventHandler::par1jets_eeqt_2011[4] = {363.327,-4.32413,44.137,17.1511};
double ZHphotonjetsEventHandler::par1jets_gqt_2011[4] = {643.536,-5.09531,222258,8.05072};//not use
double ZHphotonjetsEventHandler::par1jets_mcgqt_2011[4] = {799.36,-5.56018,4.84617,56.6004};


#endif
