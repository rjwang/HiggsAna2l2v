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


void ZHUtils::get_frFile(const edm::ParameterSet &runProcess)
{
    std::vector<std::string> FakeRateFiles = runProcess.getParameter<std::vector<std::string> >("fakeRateFile");
    for(size_t ifile=0; ifile<FakeRateFiles.size(); ifile++) {

        TString frFile(FakeRateFiles[ifile].c_str());
        gSystem->ExpandPathName(frFile);
        TFile *fr_File=TFile::Open(frFile);

        if(fr_File) {
            cout << "[ZHUtils] retrieving fake rates from: " << frFile << endl;
            std::vector<TString> varNames(1,"");
            varNames.push_back("_mtup");
            varNames.push_back("_mtdown");
            varNames.push_back("_metup");
            varNames.push_back("_metdown");
            varNames.push_back("_jptup");
            varNames.push_back("_jptdown");
            varNames.push_back("_dphiup");
            varNames.push_back("_dphidown");
            varNames.push_back("_ewkup");
            varNames.push_back("_ewkdown");
            size_t nvarsToInclude=varNames.size();

            std::vector<TString> tagNames;
            tagNames.push_back("ele");
            tagNames.push_back("mu");

            for(size_t itag=0; itag<tagNames.size(); itag++) {
                for(size_t ivar=0; ivar<nvarsToInclude; ivar++) {
                    TString key = tagNames[itag]+"FakePt_syst"+varNames[ivar];
                    cout << "key: " << key << endl;
                    TH1F *h = (TH1F *) fr_File->Get(key);
                    h->SetDirectory(0); //THIS IS IMPORTANT FOR TH1 Weight File!!!
                    fakerate1DH_[key] = h;
                } // ivar
            } // itag

        } //fr_File
        fr_File->Close();
        cout << "[ZHUtils] close file: " << frFile << endl;
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

// from Michael Brodski <brodski@physik.rwth-aachen.de>
// Oct 30, 2014
double ZHUtils::ZZQCDNLOkfactor(float gen_z_pt /*generator MET*/)
{
    double kfactor = 0.;
    double ZZ_kFactors[]= {0.0,
                           0.94326764345169067,
                           1.0035196542739868,
                           1.1098935604095459,
                           1.2901004552841187,
                           1.1905548572540283,
                           1.3430054187774658,
                           1.1620302200317383,
                           1.3074016571044922,
                           1.2638260126113892,
                           1.2879478931427002,
                           1.196311354637146,
                           1.2901802062988281,
                           1.2089201211929321,
                           1.1979459524154663,
                           1.110659122467041,
                           1.1369823217391968,
                           1.1524343490600586,
                           1.1048606634140015,
                           1.1464763879776001,
                           1.093919038772583,
                           1.0052757263183594,
                           0.98668587207794189,
                           1.0232065916061401,
                           1.0102578401565552,
                           0.99546754360198975,
                           1.0504124164581299,
                           0.89631515741348267,
                           0.88606923818588257,
                           1.0205018520355225,
                          };

    //5-GeV binning, divide by 5
    if(gen_z_pt<150.) {
        int binnumber = (int)gen_z_pt/5;
        kfactor = ZZ_kFactors[binnumber];
    } else { //fitted linear for the k-factor
        kfactor= 0.990098 - 0.000781857*gen_z_pt;
    }

    return kfactor;
}

double ZHUtils::GetQCDNLOZZweight(const PhysicsEvent_t& phys)
{
    double weight = 1.0;
    double pt_zvv = (phys.genneutrinos[0]+phys.genneutrinos[1]).pt();
    weight = ZZQCDNLOkfactor(pt_zvv);
    return weight;
}

double ZHUtils::WZQCDNLOkfactor(float gen_nu_pt /*generator MET*/)
{
    double kfactor = 0.;
    //reweight bin-for-bin below 200Gev, then use an extrapolated linear function
    double WZ_kFactors[]= {0.0,
                           1.9187283515930176,
                           1.851645827293396,
                           1.8116303682327271,
                           1.8450329303741455,
                           1.8722785711288452,
                           2.0393874645233154,
                           1.9451465606689453,
                           1.7974948883056641,
                           1.5585530996322632,
                           1.4086459875106812,
                           1.3027514219284058,
                           1.2167214155197144,
                           1.1457395553588867,
                           1.1231216192245483,
                           1.0461404323577881,
                           1.0161194801330566,
                           0.9654882550239563,
                           0.92649304866790771,
                           0.91419690847396851,
                           0.86213141679763794,
                           0.85179990530014038,
                           0.84018367528915405,
                           0.80480223894119263,
                           0.8126557469367981,
                           0.81011945009231567,
                           0.75017786026000977,
                           0.76554232835769653,
                           0.71625840663909912,
                           0.74609845876693726,
                           0.74994057416915894,
                           0.7460666298866272,
                           0.74263542890548706,
                           0.72279000282287598,
                           0.74116051197052002,
                           0.74059844017028809,
                           0.73111110925674438,
                           0.7415471076965332,
                           0.71208411455154419,
                           0.70913481712341309,
                           0.73165619373321533
                          };
    //5-GeV binning, divide by 5
    if(gen_nu_pt<200.) {
        int binnumber = (int)gen_nu_pt/5;
        kfactor = WZ_kFactors[binnumber];
    } else { //fitted linear for the k-factor
        kfactor= 0.763757 - 0.000217464*gen_nu_pt;
    }

    return kfactor;
}

double ZHUtils::GetQCDNLOWZweight(const PhysicsEvent_t& phys)
{
    double weight = 1.0;
    int Nneutrinos = phys.genneutrinos.size();
    if( Nneutrinos==1) {
        double vpt = phys.genneutrinos[0].pt();
        weight = WZQCDNLOkfactor(vpt);
    }
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
            else 	   prompt_rate = 0.9498;

        } else if(abseta<1.48) {
            if(pt<25)      prompt_rate = 0.7642;
            else if(pt<30) prompt_rate = 0.8116;
            else if(pt<35) prompt_rate = 0.8474;
            else if(pt<40) prompt_rate = 0.8791;
            else if(pt<45) prompt_rate = 0.8932;
            else if(pt<60) prompt_rate = 0.9109;
            else 	   prompt_rate = 0.9283;

        } else if(abseta<2.0) {
            if(pt<25)      prompt_rate = 0.7406;
            else if(pt<30) prompt_rate = 0.7882;
            else if(pt<35) prompt_rate = 0.8178;
            else if(pt<40) prompt_rate = 0.8554;
            else if(pt<45) prompt_rate = 0.8774;
            else if(pt<60) prompt_rate = 0.8975;
            else 	   prompt_rate = 0.9273;

        } else {
            if(pt<25)      prompt_rate = 0.7255;
            else if(pt<30) prompt_rate = 0.7770;
            else if(pt<35) prompt_rate = 0.8191;
            else if(pt<40) prompt_rate = 0.8528;
            else if(pt<45) prompt_rate = 0.8790;
            else if(pt<60) prompt_rate = 0.8989;
            else 	   prompt_rate = 0.9252;
        }

    } else if(pdgid==13) { //muon

        if(abseta<1.0) {
            if(pt<25)      prompt_rate = 0.8929;
            else if(pt<30) prompt_rate = 0.9292;
            else if(pt<35) prompt_rate = 0.9556;
            else if(pt<40) prompt_rate = 0.9732;
            else if(pt<45) prompt_rate = 0.9842;
            else if(pt<60) prompt_rate = 0.9884;
            else 	   prompt_rate = 0.9904;

        } else if(abseta<1.48) {
            if(pt<25)      prompt_rate = 0.9259;
            else if(pt<30) prompt_rate = 0.9477;
            else if(pt<35) prompt_rate = 0.9651;
            else if(pt<40) prompt_rate = 0.9788;
            else if(pt<45) prompt_rate = 0.9886;
            else if(pt<60) prompt_rate = 0.9915;
            else 	   prompt_rate = 0.9904;

        } else if(abseta<2.0) {
            if(pt<25)      prompt_rate = 0.9412;
            else if(pt<30) prompt_rate = 0.9593;
            else if(pt<35) prompt_rate = 0.9708;
            else if(pt<40) prompt_rate = 0.9805;
            else if(pt<45) prompt_rate = 0.9870;
            else if(pt<60) prompt_rate = 0.9883;
            else 	   prompt_rate = 0.9858;

        } else {
            if(pt<25)      prompt_rate = 0.9260;
            else if(pt<30) prompt_rate = 0.9506;
            else if(pt<35) prompt_rate = 0.9664;
            else if(pt<40) prompt_rate = 0.9750;
            else if(pt<45) prompt_rate = 0.9818;
            else if(pt<60) prompt_rate = 0.9838;
            else 	   prompt_rate = 0.9824;
        }

    }

    return prompt_rate;
}




double ZHUtils::fakeRate(int pdgid, double pt, double abseta, TString key)
{
    double fake_rate = 0.;
    TString tag;
    if(pdgid==11) tag="ele"+key;
    else if(pdgid==13) tag="mu"+key;

    TH1F* fr_h = fakerate1DH_[tag];
    if(fr_h==0) cout << "cannot find hist: " << tag << endl;
    int binx = fr_h->GetXaxis()->FindBin(pt);
    fake_rate = fr_h->GetBinContent(binx);

    /*
    if(pdgid==11) { //electron

        if(pt<25)      fake_rate = 0.106236;
        else if(pt<30) fake_rate = 0.104462;
        else if(pt<35) fake_rate = 0.101852;
        else if(pt<40) fake_rate = 0.103804;
        else if(pt<45) fake_rate = 0.106015;
        else if(pt<50) fake_rate = 0.111559;
        else if(pt<60) fake_rate = 0.111799;
        else 	       fake_rate = 0.110794;

    } else if(pdgid==13) { //muon

        if(pt<25)      fake_rate = 0.0803843;
        else if(pt<30) fake_rate = 0.0864601;
        else if(pt<35) fake_rate = 0.085329;
        else if(pt<40) fake_rate = 0.101627;
        else if(pt<45) fake_rate = 0.0939324;
        else if(pt<50) fake_rate = 0.0651487;
        else if(pt<60) fake_rate = 0.0523345;
        else           fake_rate = 0.101818;
    }
    */

    return fake_rate;

}


double ZHUtils::getN_PFweight(int TL_type, LorentzVector lep1, int id1, LorentzVector lep2, int id2, TString key)
{
    double p_1 = promptRate( id1, lep1.pt(), fabs(lep1.eta()) );
    double p_2 = promptRate( id2, lep2.pt(), fabs(lep2.eta()) );
    double f_1 = fakeRate( id1, lep1.pt(), fabs(lep1.eta()), key);
    double f_2 = fakeRate( id2, lep2.pt(), fabs(lep2.eta()), key);

    double weight = 1./(p_1 - f_1);
    weight /= (p_2 - f_2);
    //TL_type: 00(LL) 01(TL) 10(LT) 11(TT)
    if(TL_type == 0) weight *= -f_1*p_2;
    if(TL_type == 1) weight *= (1.-f_1)*p_2;
    if(TL_type == 2) weight *= f_1*(1.-p_2);
    if(TL_type == 3) weight *= -(1.-f_1)*(1.-p_2);

    return weight;
}


double ZHUtils::getN_FPweight(int TL_type, LorentzVector lep1, int id1, LorentzVector lep2, int id2, TString key)
{
    double p_1 = promptRate( id1, lep1.pt(), fabs(lep1.eta()) );
    double p_2 = promptRate( id2, lep2.pt(), fabs(lep2.eta()) );
    double f_1 = fakeRate( id1, lep1.pt(), fabs(lep1.eta()), key);
    double f_2 = fakeRate( id2, lep2.pt(), fabs(lep2.eta()), key);

    double weight = 1./(p_1 - f_1);
    weight /= (p_2 - f_2);
    //TL_type: 00(LL) 01(TL) 10(LT) 11(TT)
    if(TL_type == 0) weight *= -p_1*f_2;
    if(TL_type == 1) weight *= (1.-p_1)*f_2;
    if(TL_type == 2) weight *= p_1*(1.-f_2);
    if(TL_type == 3) weight *= -(1.-p_1)*(1.-f_2);

    return weight;
}

double ZHUtils::getN_FFweight(int TL_type, LorentzVector lep1, int id1, LorentzVector lep2, int id2, TString key)
{
    double p_1 = promptRate( id1, lep1.pt(), fabs(lep1.eta()) );
    double p_2 = promptRate( id2, lep2.pt(), fabs(lep2.eta()) );
    double f_1 = fakeRate( id1, lep1.pt(), fabs(lep1.eta()), key);
    double f_2 = fakeRate( id2, lep2.pt(), fabs(lep2.eta()), key);

    double weight = 1./(p_1 - f_1);
    weight /= (p_2 - f_2);
    //TL_type: 00(LL) 01(TL) 10(LT) 11(TT)
    if(TL_type == 0) weight *= p_1*p_2;
    if(TL_type == 1) weight *= -(1.-p_1)*p_2;
    if(TL_type == 2) weight *= -p_1*(1.-p_2);
    if(TL_type == 3) weight *= (1.-p_1)*(1.-p_2);

    return weight;
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
