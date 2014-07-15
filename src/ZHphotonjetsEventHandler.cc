#include "CMGTools/HiggsAna2l2v/interface/ZHphotonjetsEventHandler.h"

using namespace std;

//
ZHphotonjetsEventHandler::ZHphotonjetsEventHandler(const edm::ParameterSet &runProcess)
    : isGoodEvent_(false)
{
    //trigger thresholds to consider
    isMC_ = runProcess.getParameter<bool>("isMC");

    //open file and retrieve weights + mass shapes
    std::vector<std::string> gammaPtWeightsFiles =  runProcess.getParameter<std::vector<std::string> >("weightsFile");
    TString cats[]   =  {"eq0jets","eq1jets"};//,"lesq1jets"};
    TString dilCats[] = {"ee","mumu","ll"};
    for(size_t ifile=0; ifile<gammaPtWeightsFiles.size(); ifile++) {
        TString gammaPtWeightsFile(gammaPtWeightsFiles[ifile].c_str());
        gSystem->ExpandPathName(gammaPtWeightsFile);
        TFile *fwgt=TFile::Open(gammaPtWeightsFile);

        if(fwgt) {
            cout << "[ZHphotonjetsEventHandler] retrieving weights from: " << gammaPtWeightsFile << endl;
            TString wgtName("qt");
            //TString wgtType( isMC_ ? "mcfitwgts" : "datafitwgts");
            TString wgtType( isMC_ ? "mcwgts" : "datawgts");

            if (gammaPtWeightsFile.Contains("nvtx")) {
                wgtName="nvtx";
                wgtType=( isMC_ ? "mcwgts" : "datawgts");
            }

            if (gammaPtWeightsFile.Contains("metPhi")) {
                wgtName="metPhi";
                wgtType=( isMC_ ? "mcwgts" : "datawgts");
            }


            //event weights
            //std::map<TString, TGraph*> iWgtsH;
            std::map<TString, TH1*> iWgtsH;
            for(size_t ic=0; ic<sizeof(cats)/sizeof(TString); ic++) {
                for(size_t id=0; id<sizeof(dilCats)/sizeof(TString); id++) {
                    TString key = dilCats[id] + cats[ic] + "_" + wgtName + "_" + wgtType;
                    TH1 *h = (TH1 *) fwgt->Get(key);
                    //TGraph *h = (TGraph *) fwgt->Get(key);
                    if(h!=0) {
                        cout << "FOUND HIST: " << key << endl;//"\t";
                        h->SetDirectory(0); //THIS IS IMPORTANT FOR TH1 Weight File!!!
                        key = dilCats[id] + cats[ic];
                        iWgtsH[key] = h;
                        //cout << "Nbins: " << iWgtsH[key]->GetXaxis()->GetNbins() << endl;
                    } //else {
                    //  cout << "CANNOT FIND HIST: " << key << endl;
                    //}

                    key = dilCats[id];

                    //mass shape (take it from the first file with weights)
                    TString hname = key+"_qmass";
                    TH1 *massh = (TH1 *) fwgt->Get(hname);
                    if(massh!=0) {
			cout << "FOUND HIST: " << hname << endl;
                        massh->SetDirectory(0);
                        zmassH_[key]= massh;
                    }

                    //azimuthal angle between lepton pairs
                    TString hname_dphi = key+"_qdphi2lep";
                    TH1 *dphi2lep = (TH1 *) fwgt->Get(hname_dphi);
                    if(dphi2lep!=0) {
			cout << "FOUND HIST: " << hname_dphi << endl;
                        dphi2lep->SetDirectory(0);
                        dphi2leptonH_[key] = dphi2lep;
                    }

                    //azimuthal angle between lepton pairs
                    TString hname_dphivsqt = key+cats[ic]+"_qdphi2lepvsqt";
                    TH2 *dphi2lepvsqt = (TH2 *) fwgt->Get(hname_dphivsqt);
                    if(dphi2lepvsqt!=0) {
			cout << "FOUND HIST: " << hname_dphivsqt << endl;
                        dphi2lepvsqt->SetDirectory(0);
                        dphi2lepvsqtH_[key] = dphi2lepvsqt;
                    }

                    //Collins-Soper angle
                    TString hname_qCoslZ = key+"_qCosl1Z_CS";
                    TH1 *qCosl1Z_CS = (TH1 *) fwgt->Get(hname_qCoslZ);
                    if(qCosl1Z_CS!=0) {
			cout << "FOUND HIST: " << hname_qCoslZ << endl;
                        qCosl1Z_CS->SetDirectory(0);
                        qCosl1Z_CSH_[key] = qCosl1Z_CS;
                    }


                    //Collins-Soper angle
                    TString hname_qgammaData = cats[ic]+"_gammadata";
                    TGraph *qgammaData_h = (TGraph *) fwgt->Get(hname_qgammaData);
                    if(qgammaData_h!=0) {
                        cout << "FOUND HIST: " << qgammaData_h << endl;
                        //qgammaData_h->SetDirect0ory(0);
                        qGammaDataH_[cats[ic]] = qgammaData_h;
                    }


                }
            }
            fwgt->Close();
            WgtsH_[wgtName]=iWgtsH;
        }
    }

    if(WgtsH_.size()) std::cout << "[ZHphotonjetsEventHandler] gamma spectrum will be reweighted using distributions found in "  << gammaPtWeightsFiles.size() << " files" << std::endl;

    //cout << "WgtsH_.size(): " << WgtsH_.size() << endl;
}




// apply trigger weights for 8TeV, can we have one for 7TeV?
bool ZHphotonjetsEventHandler::isGood(PhysicsEvent_t &phys, bool is2011)
{
    //reset
    isGoodEvent_=false;
    massiveGamma_.clear();
    evWeights_.clear();
    triggerWgt_=0;

    //check if it is a gamma event
    if(phys.gammas.size()==0) return isGoodEvent_;
    if(phys.cat<22) return isGoodEvent_;

    float pt = phys.gammas[0].pt();
    if(!is2011) {
        if(pt<39) {
            if( !(phys.gammaTriggerWord & 0x1) ) return isGoodEvent_;
            triggerThr_=22;
            triggerWgt_=phys.gammaPrescale[0];
        } else if(pt>=39  && pt<55) {
            if( !( (phys.gammaTriggerWord>>1) & 0x1) ) return isGoodEvent_;
            triggerThr_=36;
            triggerWgt_=phys.gammaPrescale[1];
        } else if(pt>=55  && pt<82) {
            if( !( (phys.gammaTriggerWord>>2) & 0x1) ) return isGoodEvent_;
            triggerThr_=50;
            triggerWgt_=phys.gammaPrescale[2];
        } else if(pt>=82  && pt<100) {
            if( !( (phys.gammaTriggerWord>>3) & 0x1) ) return isGoodEvent_;
            triggerThr_=75;
            triggerWgt_=phys.gammaPrescale[3];
        } else if(pt>=100) {
            if( !( (phys.gammaTriggerWord>>4) & 0x1) ) return isGoodEvent_;
            triggerThr_=90;
            triggerWgt_=phys.gammaPrescale[4];
        }
    } else {
        triggerThr_ =( phys.cat-22)/1000;
        triggerWgt_=1;
    }

    //all done here
    isGoodEvent_=true;
    return isGoodEvent_;
}


//RJ
std::map<TString,float> ZHphotonjetsEventHandler::getWeights(ZZ2l2nuSummary_t &ev, PhysicsEvent_t &phys, TString evCategoryLabel, double _metphi_)
{
    //loop over categories
    LorentzVector gamma=phys.gammas[0];
    TString dilCats[]= {"ee","mumu","ll"};


    for(size_t id=0; id<sizeof(dilCats)/sizeof(TString); id++) {

        //the key to search for
        TString key = dilCats[id];//+evCategoryLabel;

        //generate a massive gamma (0 if in non-weighting mode)
        float mass(0);
        if(zmassH_.find(key)!=zmassH_.end()) {
            if(zmassH_[key]->Integral())
                //while(fabs(mass-91)>15)
                mass = zmassH_[key]->GetRandom();
            //cout << "Mass: " << mass << endl;
        }
        //cout << "MASS: " << mass<< endl;
        massiveGamma_[dilCats[id]]=LorentzVector(gamma.px(),gamma.py(),gamma.pz(),sqrt(pow(mass,2)+pow(gamma.energy(),2)));

        //generate azimuthal angle between lepton pairs in photon+jets sample
        float lepsdphi(0);
        if(dphi2leptonH_.find(key)!=dphi2leptonH_.end()) {
            if(dphi2leptonH_[key]->Integral())
                lepsdphi = dphi2leptonH_[key]->GetRandom();
        }
        dphimassiveGamma_[dilCats[id]]=lepsdphi;

	//new dphi(l1,l2)
        float dphi2leps(0);
        if(dphi2lepvsqtH_.find(key)!=dphi2lepvsqtH_.end()) {
	    int ptbin_ = dphi2lepvsqtH_[key]->GetXaxis()->FindBin(gamma.pt());
	    TH1* tmp_H = dphi2lepvsqtH_[key]->ProjectionY("",ptbin_,ptbin_);
            if(tmp_H->Integral()){
                dphi2leps = tmp_H->GetRandom();
	    }
        }
        dphi2lepGamma_[dilCats[id]] = dphi2leps;


        //generate Collins-Soper angle in photon+jets sample
        float coslz_cs(0);
        if(qCosl1Z_CSH_.find(key)!=qCosl1Z_CSH_.end()) {
            if(qCosl1Z_CSH_[key]->Integral())
                coslz_cs = qCosl1Z_CSH_[key]->GetRandom();
        }
        coslzmassiveGamma_[dilCats[id]]=coslz_cs;


        //get event weight (will be 0 by default if we're running in weighting mode)
        key = dilCats[id]+evCategoryLabel;
        double weight(1.0);

        //for(std::map<TString, std::map<TString,TGraph *> >::iterator wIt = WgtsH_.begin(); wIt != WgtsH_.end(); wIt++)
        for(std::map<TString, std::map<TString,TH1 *> >::iterator wIt = WgtsH_.begin(); wIt != WgtsH_.end(); wIt++) {

            if(wIt->second.find(key) == wIt->second.end()) {
                //cout << key << " not found" << endl;
                continue;
            }
            //float iweight(1.0);
            //TGraph *h = wIt->second[key];
            TH1 *h = wIt->second[key];

            //if(wIt->first == "qt") weight *= h->GetBinContent( h->FindBin(gamma.pt()) );
            //if(wIt->first == "nvtx") weight *= h->GetBinContent( h->FindBin(ev.nvtx) );
            //if(wIt->first == "metPhi") weight *= h->GetBinContent( h->FindBin(_metphi_) );

	    if(wIt->first == "qt") {
		double qt_weight = h->GetBinContent( h->FindBin(gamma.pt()) );
		weight *= qt_weight>0 ? qt_weight : 1.0;
	    }
	    if(wIt->first == "nvtx"){
		double nvtx_weight = h->GetBinContent( h->FindBin(ev.nvtx) );
		weight *= nvtx_weight>0 ? nvtx_weight : 1.0;
	    }
	    if(wIt->first == "metPhi"){
		double metPhi_weight = h->GetBinContent( h->FindBin(_metphi_) );
		weight *= metPhi_weight >0 ? metPhi_weight : 1.0;
	    }

            //cout << wIt->first << "-" << key << " qt: " << gamma.pt() << ": "<< weight << endl;
            //if(weight<0) weight=0;
        }
        evWeights_[dilCats[id]]=weight;
    }

    return evWeights_;
}


std::map<TString,double> ZHphotonjetsEventHandler::getQtFitWeights(PhysicsEvent_t &phys, TString evCategoryLabel)
{
    //loop over categories
    LorentzVector gamma=phys.gammas[0];
    TString dilCats[]= {"ee","mumu","ll"};

    for(size_t id=0; id<sizeof(dilCats)/sizeof(TString); id++) {
        //get event weight (will be 0 by default if we're running in weighting mode)
        TString key = dilCats[id]+evCategoryLabel;
        double weight(1.0);

        double *Zjets_pars;
        double *Gjets_pars;
        //get all parameters, call ZHphotonjetsEventHandler::getQtFit
        if(isMC_) {
            if(key.Contains("eeeq0jets")) {
                Zjets_pars = &par0jets_eeqt[0];
                Gjets_pars = &par0jets_mcgqt[0];
            } else if(key.Contains("eeeq1jets")) {
                Zjets_pars = &par1jets_eeqt[0];
                Gjets_pars = &par1jets_mcgqt[0];
            } else if(key.Contains("mumueq0jets")) {
                Zjets_pars = &par0jets_mmqt[0];
                Gjets_pars = &par0jets_mcgqt[0];
            } else if(key.Contains("mumueq1jets")) {
                Zjets_pars = &par1jets_mmqt[0];
                Gjets_pars = &par1jets_mcgqt[0];
            }
        } else {
            if(key.Contains("eeeq0jets")) {
                Zjets_pars = &par0jets_eeqt[0];
                Gjets_pars = &par0jets_gqt[0];
            } else if(key.Contains("eeeq1jets")) {
                Zjets_pars = &par1jets_eeqt[0];
                Gjets_pars = &par1jets_gqt[0];
            } else if(key.Contains("mumueq0jets")) {
                Zjets_pars = &par0jets_mmqt[0];
                Gjets_pars = &par0jets_gqt[0];
            } else if(key.Contains("mumueq1jets")) {
                Zjets_pars = &par1jets_mmqt[0];
                Gjets_pars = &par1jets_gqt[0];
            }
        }


        if(!key.Contains("ll") && (key.Contains("eq0jets") || key.Contains("eq1jets")) ) {
            //double pt[1] = {gamma.pt()};
            double pt[] = {gamma.pt()};
            double *pT = &pt[0];
            double iweight = getQtFit(pT,Zjets_pars)/getQtFit(pT,Gjets_pars);
            weight *= iweight;

            //for debug
            /*
                        double *_par = Zjets_pars;
                        cout << key << ": "
                             << _par[0] << "\t"
                             << _par[1] << "\t"
                             << _par[2] << "\t"
                             << _par[3] << "\t"
                             << "Pt: " << pT[0] << "\t"
                             << "weights: " << weight << endl;
            */
        }

        QtFitWeights_[dilCats[id]]=weight;

    }

    return QtFitWeights_;

}





std::map<TString,double> ZHphotonjetsEventHandler::get2011QtFitWeights(PhysicsEvent_t &phys, TString evCategoryLabel)
{
    //loop over categories
    LorentzVector gamma=phys.gammas[0];
    TString dilCats[]= {"ee","mumu","ll"};

    for(size_t id=0; id<sizeof(dilCats)/sizeof(TString); id++) {
        //get event weight (will be 0 by default if we're running in weighting mode)
        TString key = dilCats[id]+evCategoryLabel;
        double weight(1.0);

        double *Zjets_pars;
        double *Gjets_pars;
        //get all parameters, call ZHphotonjetsEventHandler::getQtFit
        if(isMC_) {
            if(key.Contains("eeeq0jets")) {
                Zjets_pars = &par0jets_eeqt_2011[0];
                Gjets_pars = &par0jets_mcgqt_2011[0];
            } else if(key.Contains("eeeq1jets")) {
                Zjets_pars = &par1jets_eeqt_2011[0];
                Gjets_pars = &par1jets_mcgqt_2011[0];
            } else if(key.Contains("mumueq0jets")) {
                Zjets_pars = &par0jets_mmqt_2011[0];
                Gjets_pars = &par0jets_mcgqt_2011[0];
            } else if(key.Contains("mumueq1jets")) {
                Zjets_pars = &par1jets_mmqt_2011[0];
                Gjets_pars = &par1jets_mcgqt_2011[0];
            }
        } else {
            if(key.Contains("eeeq0jets")) {
                Zjets_pars = &par0jets_eeqt_2011[0];
            } else if(key.Contains("eeeq1jets")) {
                Zjets_pars = &par1jets_eeqt_2011[0];
            } else if(key.Contains("mumueq0jets")) {
                Zjets_pars = &par0jets_mmqt_2011[0];
            } else if(key.Contains("mumueq1jets")) {
                Zjets_pars = &par1jets_mmqt_2011[0];
            }
        }


        if(!key.Contains("ll") && (key.Contains("eq0jets") || key.Contains("eq1jets")) ) {
            //double pt[1] = {gamma.pt()};
            double pt[] = {gamma.pt()};
            double *pT = &pt[0];
            double iweight;
	    if(isMC_) iweight = getQtFit(pT,Zjets_pars)/getQtFit(pT,Gjets_pars);
	    else iweight = getQtFit(pT,Zjets_pars)/qGammaDataH_[evCategoryLabel]->Eval(gamma.pt()); 

            weight *= iweight;

            //for debug
            /*
                        double *_par = Zjets_pars;
                        cout << key << ": "
                             << _par[0] << "\t"
                             << _par[1] << "\t"
                             << _par[2] << "\t"
                             << _par[3] << "\t"
                             << "Pt: " << pT[0] << "\t"
                             << "weights: " << weight << endl;
            */
        }

        QtFitWeights_[dilCats[id]]=weight;

    }

    return QtFitWeights_;

}



double ZHphotonjetsEventHandler::getQtFit(double *x, double *par)
{
    return pow(x[0]/par[0],par[1])*exp(-x[0]/par[0])/(1+par[2]*exp(-x[0]/par[3]));
}



//
ZHphotonjetsEventHandler::~ZHphotonjetsEventHandler()
{
}
