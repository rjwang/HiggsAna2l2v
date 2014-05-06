/** \class ZHUtils
 *
 *  \author R.-J. Wang
 */


#include "CMGTools/HiggsAna2l2v/interface/ZHUtils.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "TVector2.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TRandom.h"

using namespace std;

//namespace ZHUtils {

ZHUtils::ZHUtils(const edm::ParameterSet &runProcess)
{
    std::vector<std::string> ZjetsWeightsFiles =  runProcess.getParameter<std::vector<std::string> >("weightsFile");
    std::vector<TString> cats;
    cats.push_back("eeeq0jets");
    cats.push_back("eeeq1jets");
    cats.push_back("mumueq0jets");
    cats.push_back("mumueq1jets");

    for(size_t ifile=0; ifile<ZjetsWeightsFiles.size(); ifile++) {
        TString ZjetsWeightsFile(ZjetsWeightsFiles[ifile].c_str());
        gSystem->ExpandPathName(ZjetsWeightsFile);
        TFile *fwgt=TFile::Open(ZjetsWeightsFile);

        std::map<TString, TGraph*> iWgtsH;
        std::map<TString, TH2F*> iWgts2DH;
        TString wgtName;
        if(fwgt) {
            cout << "[ZHUtils] retrieving weights from: " << ZjetsWeightsFile << endl;

            if(ZjetsWeightsFile.Contains("met_met")) wgtName="met_met";
            else if(ZjetsWeightsFile.Contains("MllvsPt")) wgtName="MllvsPt";

            //event weights
            for(size_t ic=0; ic<cats.size(); ic++) {
                TString key = cats[ic] + "_" + wgtName;
                cout << "key: " << key << endl;
                if(wgtName=="met_met") {
                    TGraph *h = (TGraph *) fwgt->Get(key);
                    if(h!=0) iWgtsH[key] = h;
                } else if(wgtName=="MllvsPt") {
                    TH2F *h = (TH2F *) fwgt->Get(key);
                    h->SetDirectory(0); //THIS IS IMPORTANT FOR TH1 Weight File!!!
                    if(h!=0) {
                        iWgts2DH[key] = h;
                        cout << "Saving hists." << endl;
                    }
                }
            }
            wgtsH_[wgtName]=iWgtsH;
            wgts2DH_[wgtName]=iWgts2DH;

        }
        fwgt->Close();
        cout << "[ZHUtils] close weights file: " << ZjetsWeightsFile << endl;

    }
}


double ZHUtils::Collins_Soper(const LorentzVector& lepton1, const LorentzVector& lepton2)
{

    TLorentzVector LVlep1(lepton1.Px(), lepton1.Py(), lepton1.Pz(), lepton1.E());
    TLorentzVector LVlep2(lepton2.Px(), lepton2.Py(), lepton2.Pz(), lepton2.E());
    TLorentzVector LVZ = LVlep1+LVlep2;

    //do transformation to Collins-Soper frame (1 rotation, 2 boosts)

    // 1st transormation - rotate to ptZ direction
    double zrot = -LVZ.Phi();
    LVlep1.RotateZ(zrot);
    LVlep2.RotateZ(zrot);

    // 2nd transformation - boost in z direction
    double beta_boostz = -LVZ.Pz()/LVZ.E();
    LVlep1.Boost(0.,0.,beta_boostz);
    LVlep2.Boost(0.,0.,beta_boostz);


    // 3rd transformation: boost in transverse direction (x-prime)
    double beta_boostx = -(LVlep1.Px()+LVlep2.Px())/(LVlep1.E()+LVlep2.E());
    LVlep1.Boost(beta_boostx,0.,0.);
    LVlep2.Boost(beta_boostx,0.,0.);

    // compute cos(theta*) in Colin-Soper frame
    double cos_theta_CS = LVlep1.CosTheta();
    if (LVZ.Pz() < 0) {
        cos_theta_CS *= -1.;
    }

    return cos_theta_CS;
}


float ZHUtils::weightNLOEWKzh(float pt)
{
    if(pt < 50) return 1;
    return 0.94-(0.2-0.068)/400.*(pt-50.);
}


float ZHUtils::weightNLOEWKzz(float pt)
{
    if(pt < 50) return 1;
    return 1.+(-0.071*pt+0.55)/100.;
}

float ZHUtils::weightNLOEWKwz(float pt)
{
    if(pt < 50) return 1;
    return 1.+(-0.037*pt+1.9)/100.;
}



double ZHUtils::GetNLOZHWeight(const PhysicsEvent_t& phys)
{
    double weight = 1.0;

    double pt_dilep = (phys.genleptons[0]+phys.genleptons[1]).pt();
    weight = weightNLOEWKzh(pt_dilep);
    return weight;
}



double ZHUtils::GetNLOZZWeight(const PhysicsEvent_t& phys)
{
    double weight = 1.0;
    double pt_dilep = (phys.genleptons[0]+phys.genleptons[1]).pt();
    double pt_zvv = (phys.genneutrinos[0]+phys.genneutrinos[1]).pt();
    double trailing_pt = pt_dilep < pt_zvv ? pt_dilep : pt_zvv;
    weight = weightNLOEWKzz(trailing_pt);
    return weight;
}


double ZHUtils::GetNLOWZWeight(const PhysicsEvent_t& phys)
{
    double weight = 1.0;
    int ngenleps = phys.genleptons.size();
    int Nneutrinos = phys.genneutrinos.size();
    if(ngenleps==3 && Nneutrinos==1) {
        int id_lep0 = fabs(phys.genleptons[0].id);
        int id_lep1 = fabs(phys.genleptons[1].id);
        int id_lep2 = fabs(phys.genleptons[2].id);
        int id_Lep0 = phys.genleptons[0].id;
        int id_Lep1 = phys.genleptons[1].id;
        int id_Lep2 = phys.genleptons[2].id;
        bool hassameids = (id_lep0==id_lep1)&&(id_lep1==id_lep2);

        double pt_zll(-1.);
        double pt_wlv(-1);
        if(!hassameids) {
            if(id_lep0==id_lep1) {
                pt_zll = (phys.genleptons[0]+phys.genleptons[1]).pt();
                pt_wlv = (phys.genleptons[2]+phys.genneutrinos[0]).pt();
            } else if(id_lep0==id_lep2) {
                pt_zll = (phys.genleptons[0]+phys.genleptons[2]).pt();
                pt_wlv = (phys.genleptons[1]+phys.genneutrinos[0]).pt();
            } else if(id_lep1==id_lep2) {
                pt_zll = (phys.genleptons[1]+phys.genleptons[2]).pt();
                pt_wlv = (phys.genleptons[0]+phys.genneutrinos[0]).pt();
            }
        } else { // if 3 leptons has same ids

            //find the one different lepton
            int sum_lepid = id_Lep0+id_Lep1+id_Lep2;
            if(id_Lep0 == -sum_lepid) {
                double mass1 = fabs((phys.genleptons[0]+phys.genleptons[1]).mass()-91.);
                double mass2 = fabs((phys.genleptons[0]+phys.genleptons[2]).mass()-91.);
                if(mass1 < mass2) {
                    pt_zll = (phys.genleptons[0]+phys.genleptons[1]).pt();
                    pt_wlv = (phys.genleptons[2]+phys.genneutrinos[0]).pt();
                } else {
                    pt_zll = (phys.genleptons[0]+phys.genleptons[2]).pt();
                    pt_wlv = (phys.genleptons[1]+phys.genneutrinos[0]).pt();
                }
            } else if(id_Lep1 == -sum_lepid) {
                double mass1 = fabs((phys.genleptons[1]+phys.genleptons[0]).mass()-91.);
                double mass2 = fabs((phys.genleptons[1]+phys.genleptons[2]).mass()-91.);
                if(mass1 < mass2) {
                    pt_zll = (phys.genleptons[1]+phys.genleptons[0]).pt();
                    pt_wlv = (phys.genleptons[2]+phys.genneutrinos[0]).pt();
                } else {
                    pt_zll = (phys.genleptons[1]+phys.genleptons[2]).pt();
                    pt_wlv = (phys.genleptons[0]+phys.genneutrinos[0]).pt();
                }
            } else if(id_Lep2 == -sum_lepid) {
                double mass1 = fabs((phys.genleptons[2]+phys.genleptons[0]).mass()-91.);
                double mass2 = fabs((phys.genleptons[2]+phys.genleptons[1]).mass()-91.);
                if(mass1 < mass2) {
                    pt_zll = (phys.genleptons[2]+phys.genleptons[0]).pt();
                    pt_wlv = (phys.genleptons[1]+phys.genneutrinos[0]).pt();
                } else {
                    pt_zll = (phys.genleptons[2]+phys.genleptons[1]).pt();
                    pt_wlv = (phys.genleptons[0]+phys.genneutrinos[0]).pt();
                }
            }
        }//if 3 leptons has same ids

        double trailing_pt = pt_zll < pt_wlv ? pt_zll : pt_wlv;
        weight = weightNLOEWKwz(trailing_pt);
    } //ngenleps==3 && Nneutrinos==1

    return weight;
}


std::map<TString,float> ZHUtils::getWeights(double ValtoWeight, TString wgtName)
{
    std::vector<TString> cats;
    cats.push_back("eeeq0jets");
    cats.push_back("eeeq1jets");
    cats.push_back("mumueq0jets");
    cats.push_back("mumueq1jets");


    std::map<TString, float> evWeights_;
    float weight(1.0);
    for(std::map<TString, std::map<TString,TGraph *> >::iterator wIt = wgtsH_.begin(); wIt != wgtsH_.end(); wIt++) {
        for(size_t ic=0; ic<cats.size(); ic++) {
            weight=1.0;
            TString key = cats[ic] + "_" + wgtName;
            if(wIt->second.find(key) == wIt->second.end()) {
                cout << key << " not found" << endl;
                continue;
            }
            TGraph *h = wIt->second[key];
            if(wgtName == "met_met") weight *= h->Eval(ValtoWeight);
            evWeights_[cats[ic]] = weight;
        }
    }

    return evWeights_;
}




double ZHUtils::get2DWeights(double Val_x, double Val_y, TString wgtName, TString cat)
{

    float weight(1.0);
    for(std::map<TString, std::map<TString,TH2F *> >::iterator wIt = wgts2DH_.begin(); wIt != wgts2DH_.end(); wIt++) {
        weight=1.0;
        TString key = cat + "_" + wgtName;
        if(wIt->second.find(key) == wIt->second.end()) {
            //cout << key << " not found" << endl;
            continue;
        }
        if(wgtName == "MllvsPt") {
            TH2F *h = wIt->second[key];
            if(h==0) cout << "cannot find hist!" << endl;
            int binx = h->GetXaxis()->FindBin(Val_x);
            int biny = h->GetYaxis()->FindBin(Val_y);
            weight *= h->GetBinContent(binx,biny);
        }
    }

    return weight;
}



double ZHUtils::promptRate(int pdgid, double pt, double abseta)
{
    double prompt_rate = 1.0;

    if(pdgid==11) { //electron

        if(abseta<1.0) {
            if(pt<25) 	   prompt_rate = 0.8046;
            else if(pt<30) prompt_rate = 0.8476;
            else if(pt<35) prompt_rate = 0.8857;
            else if(pt<40) prompt_rate = 0.9173;
            else if(pt<45) prompt_rate = 0.9322;
            else if(pt<60) prompt_rate = 0.9418;
            else if(pt<100) prompt_rate = 0.9498;

        } else if(abseta<1.48) {
            if(pt<25)      prompt_rate = 0.7642;
            else if(pt<30) prompt_rate = 0.8116;
            else if(pt<35) prompt_rate = 0.8474;
            else if(pt<40) prompt_rate = 0.8791;
            else if(pt<45) prompt_rate = 0.8932;
            else if(pt<60) prompt_rate = 0.9109;
            else if(pt<100) prompt_rate = 0.9283;

        } else if(abseta<2.0) {
            if(pt<25)      prompt_rate = 0.7406;
            else if(pt<30) prompt_rate = 0.7882;
            else if(pt<35) prompt_rate = 0.8178;
            else if(pt<40) prompt_rate = 0.8554;
            else if(pt<45) prompt_rate = 0.8774;
            else if(pt<60) prompt_rate = 0.8975;
            else if(pt<100) prompt_rate = 0.9273;

        } else {
            if(pt<25)      prompt_rate = 0.7255;
            else if(pt<30) prompt_rate = 0.7770;
            else if(pt<35) prompt_rate = 0.8191;
            else if(pt<40) prompt_rate = 0.8528;
            else if(pt<45) prompt_rate = 0.8790;
            else if(pt<60) prompt_rate = 0.8989;
            else if(pt<100) prompt_rate = 0.9252;
        }

    } else if(pdgid==13) { //muon

        if(abseta<1.0) {
            if(pt<25)      prompt_rate = 0.8929;
            else if(pt<30) prompt_rate = 0.9292;
            else if(pt<35) prompt_rate = 0.9556;
            else if(pt<40) prompt_rate = 0.9732;
            else if(pt<45) prompt_rate = 0.9842;
            else if(pt<60) prompt_rate = 0.9884;
            else if(pt<100) prompt_rate = 0.9904;

        } else if(abseta<1.48) {
            if(pt<25)      prompt_rate = 0.9259;
            else if(pt<30) prompt_rate = 0.9477;
            else if(pt<35) prompt_rate = 0.9651;
            else if(pt<40) prompt_rate = 0.9788;
            else if(pt<45) prompt_rate = 0.9886;
            else if(pt<60) prompt_rate = 0.9915;
            else if(pt<100) prompt_rate = 0.9904;

        } else if(abseta<2.0) {
            if(pt<25)      prompt_rate = 0.9412;
            else if(pt<30) prompt_rate = 0.9593;
            else if(pt<35) prompt_rate = 0.9708;
            else if(pt<40) prompt_rate = 0.9805;
            else if(pt<45) prompt_rate = 0.9870;
            else if(pt<60) prompt_rate = 0.9883;
            else if(pt<100) prompt_rate = 0.9858;

        } else {
            if(pt<25)      prompt_rate = 0.9260;
            else if(pt<30) prompt_rate = 0.9506;
            else if(pt<35) prompt_rate = 0.9664;
            else if(pt<40) prompt_rate = 0.9750;
            else if(pt<45) prompt_rate = 0.9818;
            else if(pt<60) prompt_rate = 0.9838;
            else if(pt<100) prompt_rate = 0.9824;
        }

    }

    return prompt_rate;
}


/*
   float weightNLOEWKzz(float pt)
{
    if(pt < 50) return 1;
    double p_0 = -1.28637;
    double p_1 = -0.0967502;
    double p_2 = 5.38621e-05;
    double p_3 = -1.54371e-08;
    double a = p_0+p_1*pt+p_2*pt*pt+p_3*pow(pt,3);
    return 1.+a/100.;
}


float weightNLOEWKwplusz(float pt)
{
    if(pt < 50) return 1;
    double p_0 = 1.31664;
    double p_1 = -0.0517186;
    double p_2 = -7.28883e-06;
    double p_3 = 2.5116e-08;
    double a = p_0+p_1*pt+p_2*pt*pt+p_3*pow(pt,3);
    return 1.+a/100.;
}

float weightNLOEWKwminuz(float pt)
{
    if(pt < 50) return 1;
    double p_0 = 1.54202;
    double p_1 = -0.0497632;
    double p_2 = -1.26063e-05;
    double p_3 = 2.81725e-08;
    double a = p_0+p_1*pt+p_2*pt*pt+p_3*pow(pt,3);
    return 1.+a/100.;
}
*/

//
ZHUtils::~ZHUtils()
{
}
