#include <iostream>
#include <boost/shared_ptr.hpp>

#include "EgammaAnalysis/ElectronTools/interface/EGammaCutBasedEleId.h"

#include "CMGTools/HiggsAna2l2v/interface/ZZ2l2nuSummaryHandler.h"
#include "CMGTools/HiggsAna2l2v/interface/ZZ2l2nuPhysicsEvent.h"
#include "CMGTools/HiggsAna2l2v/interface/METUtils.h"
#include "CMGTools/HiggsAna2l2v/interface/ZHUtils.h"
#include "CMGTools/HiggsAna2l2v/interface/setStyle.h"
#include "CMGTools/HiggsAna2l2v/interface/plotter.h"
#include "CMGTools/HiggsAna2l2v/interface/ObjectFilters.h"
#include "CMGTools/HiggsAna2l2v/interface/SmartSelectionMonitor.h"
#include "CMGTools/HiggsAna2l2v/interface/TMVAUtils.h"
#include "CMGTools/HiggsAna2l2v/interface/MacroUtils.h"
#include "CMGTools/HiggsAna2l2v/interface/EventCategory.h"
#include "CMGTools/HiggsAna2l2v/interface/LeptonEfficiencySF.h"

#include "CondFormats/JetMETObjects/interface/JetResolution.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

#include "FWCore/FWLite/interface/AutoLibraryLoader.h"
#include "FWCore/PythonParameterSet/interface/MakeParameterSets.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"

#include "CMGTools/HiggsAna2l2v/src/MuScleFitCorrector_v4_1/MuScleFitCorrector.h"

#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TEventList.h"
#include "TROOT.h"
#include "TMath.h"

using namespace std;

int main(int argc, char* argv[])
{
    //##############################################
    //########    GLOBAL INITIALIZATION     ########
    //##############################################

    // check arguments
    if(argc<2) {
        std::cout << "Usage : " << argv[0] << " parameters_cfg.py" << std::endl;
        exit(0);
    }

    // load framework libraries
    gSystem->Load( "libFWCoreFWLite" );
    AutoLibraryLoader::enable();

    // configure the process
    const edm::ParameterSet &runProcess = edm::readPSetsFrom(argv[1])->getParameter<edm::ParameterSet>("runProcess");

    bool use2011Id = runProcess.getParameter<bool>("is2011");
    bool useCHS(true);
    bool useJERsmearing(true);
    bool useJetsOnlyInTracker(false);
    bool usePUsubJetId(true);
    // invert ajets <-> jets ...
    // 2011:  ajets -> non CHS,  jets -> CHS
    // 2012:  ajets -> CHS,      jets -> non CHS
    if(use2011Id) useCHS = !useCHS;
    cout << "Note: will apply " << (use2011Id ? 2011 : 2012) << " version of the id's" << endl;

    bool isMC       = runProcess.getParameter<bool>("isMC");
    int mctruthmode = runProcess.getParameter<int>("mctruthmode");

    TString url=runProcess.getParameter<std::string>("input");
    TString outFileUrl(gSystem->BaseName(url));
    outFileUrl.ReplaceAll(".root","");
    if(mctruthmode!=0) {
        outFileUrl += "_filt";
        outFileUrl += mctruthmode;
    }
    TString outdir=runProcess.getParameter<std::string>("outdir");
    TString outUrl( outdir );
    gSystem->Exec("mkdir -p " + outUrl);
    int fType(0);
    if(url.Contains("DoubleEle")) fType=EE;
    if(url.Contains("DoubleMu"))  fType=MUMU;
    if(url.Contains("MuEG"))      fType=EMU;
    if(url.Contains("SingleMu"))  fType=MUMU;
    if(url.Contains("SingleEle")) fType=EE;
    bool isSingleMuPD(!isMC && url.Contains("SingleMu"));
    bool isSingleElePD(!isMC && url.Contains("SingleEle"));
    bool isMC_ZZ  = isMC && ( string(url.Data()).find("TeV_ZZ")  != string::npos);
    bool isMC_WZ  = isMC && ( string(url.Data()).find("TeV_WZ")  != string::npos);

    // print out event information
    TString outTxtUrl_full= outUrl + "/" + outFileUrl + "_FullList.txt";
    FILE* outTxtFile_full = NULL;
    outTxtFile_full = fopen(outTxtUrl_full.Data(), "w");
    printf("TextFile URL = %s\n",outTxtUrl_full.Data());

    TString outTxtUrl_final= outUrl + "/" + outFileUrl + "_FinalList.txt";
    FILE* outTxtFile_final = NULL;
    outTxtFile_final = fopen(outTxtUrl_final.Data(), "w");
    printf("TextFile URL = %s\n",outTxtUrl_final.Data());

    fprintf(outTxtFile_full,"run lumi event passId,passIdAndIso,passZmass,passZpt,pass3dLeptonVeto,passBveto,passLMetVeto,passdphiZllmetCut,passMet,passBalanceCut,passMTcut evCat \n");
    //fprintf(outTxtFile_final,"run lumi event passId,passIdAndIso,passZmass,passZpt,pass3dLeptonVeto,passBveto,passLMetVeto,passdphiZllmetCut,passMet,passBalanceCut,passMTcut evCat \n");


    //tree info
    int evStart     = runProcess.getParameter<int>("evStart");
    int evEnd       = runProcess.getParameter<int>("evEnd");
    TString dirname = runProcess.getParameter<std::string>("dirName");

    //jet energy scale uncertainties
    TString uncFile = runProcess.getParameter<std::string>("jesUncFileName");
    gSystem->ExpandPathName(uncFile);
    JetCorrectionUncertainty jecUnc(uncFile.Data());

    //ZHinvisible reweighting file input
    ZHUtils myZHUtils(runProcess);

    //systematics
    bool runSystematics                        = runProcess.getParameter<bool>("runSystematics");
    std::vector<TString> varNames(1,"");
    if(runSystematics){
      	cout << "Systematics will be computed for this analysis" << endl;
      	varNames.push_back("_jerup");    varNames.push_back("_jerdown");
      	varNames.push_back("_jesup");    varNames.push_back("_jesdown");
      	varNames.push_back("_umetup");   varNames.push_back("_umetdown");
      	varNames.push_back("_lesup");    varNames.push_back("_lesdown");
      	varNames.push_back("_puup");     varNames.push_back("_pudown");
      	varNames.push_back("_btagup");   varNames.push_back("_btagdown");
	// will need to add ZZ and WZ shape uncertainty
      	//if(isMC_ZZ)             { varNames.push_back("_zzptup");   varNames.push_back("_zzptdown");     }
      	//if(isMC_WZ)             { varNames.push_back("_wzptup");   varNames.push_back("_wzptdown");     }
    }
    size_t nvarsToInclude=varNames.size();


    // Muon scale/resolution corrections
    TString fitParametersFile = "/afs/cern.ch/work/r/rewang/HiggsZZd/HZZ/CMSSW_5_3_11/src/CMGTools/HiggsAna2l2v/src/MuScleFitCorrector_v4_1/";
    if(use2011Id) {
        if(isMC) fitParametersFile += "MuScleFit_2011_MC_44X.txt";
        else     fitParametersFile += "MuScleFit_2011_DATA_44X.txt";
    } else {
        if(isMC) fitParametersFile += "MuScleFit_2012_MC_53X.txt"; // CHANGE!!!
        else {
            if( url.Contains("2012D") ) fitParametersFile += "MuScleFit_2012D_DATA_53X.txt";
            else                        fitParametersFile += "MuScleFit_2012ABC_DATA_53X.txt";
        }
    }
    MuScleFitCorrector *corrector_ = new MuScleFitCorrector(fitParametersFile);




    //##############################################
    //########    INITIATING HISTOGRAMS     ########
    //##############################################
    SmartSelectionMonitor mon;

    // pileup control
    mon.addHistogram( new TH1F( "nvtx",";Vertices;Events",50,0,50) );
    mon.addHistogram( new TH1F( "nvtxraw",";Vertices;Events",50,0,50) );


    double VTXaxis[29] = {0,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,50};
    mon.addHistogram( new TH1F( "nvtx_dy",";Vertices;Events",50,0,50) );
    mon.addHistogram( new TH1F( "nvtx_dy_rebin",";Vertices;Events",28,VTXaxis) );


    TH1F *h=(TH1F*) mon.addHistogram( new TH1F ("eventflow", ";;Events", 12,0,12) );
    h->GetXaxis()->SetBinLabel(1,"Trigger");
    h->GetXaxis()->SetBinLabel(2,"#geq 2 leptons");
    h->GetXaxis()->SetBinLabel(3,"#geq 2 iso leptons");

    //for MC normalization (to 1/pb)
    TH1F* Hcutflow  = (TH1F*) mon.addHistogram(  new TH1F ("cutflow"    , "cutflow"    ,6,0,6) ) ;

    Double_t qtaxis[100];
    for(size_t i=0; i<40; i++)  qtaxis[i]=2.5*i;       //0-97.5
    for(size_t i=0; i<20; i++)  qtaxis[40+i]=100+5*i;  //100-195
    for(size_t i=0; i<15; i++)  qtaxis[60+i]=200+10*i; //200-340
    for(size_t i=0; i<25; i++)  qtaxis[75+i]=350+25*i; //350-976
    double Qtaxis[15]= {0,50,55,60,65,70,75,80,85,90,100,125,150,200,1000};

    // Z+jets control with photon template
    mon.addHistogram( new TH1D( "qt_rebin",        ";#it{p}_{T}^{#gamma} [GeV];Events",99,qtaxis));
    mon.addHistogram( new TH1D( "qt_rebin_new",  ";#it{q}_{T} [GeV];Events",14,Qtaxis));
    mon.addHistogram( new TH1F( "qt",";#it{p}_{T}^{#gamma} [GeV];Events",1000,0,1000));
    mon.addHistogram( new TH1F( "qmass", ";#it{m}_{ll} [GeV];Events", 15,76,106) );

    //reweigting pt of DY MC samples
    mon.addHistogram( new TH2F( "MllvsPt", "; #it{p}_{T}^{ll} [GeV]; #it{m}_{ll} [GeV]; Events", 50,0,500, 100,40,250) );

    //for testing ABCD method
    mon.addHistogram( new TH2F( "dPhivsMET_ABCD", "; E_{T}^{miss} [GeV]; #Delta#phi(Z,E_{T}^{miss}) [rad]; Events", 500,0,500, 50,0,TMath::Pi()) );
    mon.addHistogram( new TH2F( "dPhivsMET_A", "; E_{T}^{miss} [GeV]; #Delta#phi(Z,E_{T}^{miss}) [rad]; Events", 500,0,500, 50,0,TMath::Pi()) );
    mon.addHistogram( new TH2F( "dPhivsMET_B", "; E_{T}^{miss} [GeV]; #Delta#phi(Z,E_{T}^{miss}) [rad]; Events", 500,0,500, 50,0,TMath::Pi()) );
    mon.addHistogram( new TH2F( "dPhivsMET_C", "; E_{T}^{miss} [GeV]; #Delta#phi(Z,E_{T}^{miss}) [rad]; Events", 500,0,500, 50,0,TMath::Pi()) );
    mon.addHistogram( new TH2F( "dPhivsMET_D", "; E_{T}^{miss} [GeV]; #Delta#phi(Z,E_{T}^{miss}) [rad]; Events", 500,0,500, 50,0,TMath::Pi()) );

    mon.addHistogram( new TH2F( "dPhivsBal_ABCD", "; E_{T}^{miss}/#it{q}_{T}; #Delta#phi(Z,E_{T}^{miss}) [rad]; Events", 200,0,5, 50,0,TMath::Pi()) );
    mon.addHistogram( new TH2F( "dPhivsBal_A", "; E_{T}^{miss}/#it{q}_{T}; #Delta#phi(Z,E_{T}^{miss}) [rad]; Events", 200,0,5, 50,0,TMath::Pi()) );
    mon.addHistogram( new TH2F( "dPhivsBal_B", "; E_{T}^{miss}/#it{q}_{T}; #Delta#phi(Z,E_{T}^{miss}) [rad]; Events", 200,0,5, 50,0,TMath::Pi()) );
    mon.addHistogram( new TH2F( "dPhivsBal_C", "; E_{T}^{miss}/#it{q}_{T}; #Delta#phi(Z,E_{T}^{miss}) [rad]; Events", 200,0,5, 50,0,TMath::Pi()) );
    mon.addHistogram( new TH2F( "dPhivsBal_D", "; E_{T}^{miss}/#it{q}_{T}; #Delta#phi(Z,E_{T}^{miss}) [rad]; Events", 200,0,5, 50,0,TMath::Pi()) );


    double zpt_bins[] = {0., 20., 40., 60., 80., 100., 200., 400., 800.};
    const int n_zpt_bins = sizeof(zpt_bins)/sizeof(double) - 1;
    mon.addHistogram( new TH1F( "zpt", 	     ";#it{p}_{T}^{ll} [GeV];Events", 50,0,500) );
    mon.addHistogram( new TH1F( "zpt_rebin", ";#it{p}_{T}^{ll} [GeV];Events", n_zpt_bins, zpt_bins) );
    mon.addHistogram( new TH1F( "zmass_raw", ";#it{m}_{ll} [GeV];Events", 100,40,250) );
    mon.addHistogram( new TH1F( "zmass",     ";#it{m}_{ll} [GeV];Events", 100,40,250) );

    // for HZZd MC truth
    mon.addHistogram( new TH1F( "Zdpt", ";Z_{d} #it{p}_{T} [GeV];Events", 100,0,100) );
    //Collins Soper Frame
    mon.addHistogram( new TH1F( "CoslepZ_CS",";cos(#theta*(l,Z));Events", 80,-1.,1.) );


    h = (TH1F*) mon.addHistogram( new TH1F( "nleptons", ";Lepton multiplicity;Events", 3,2,4) );
    for(int ibin=1; ibin<=h->GetXaxis()->GetNbins(); ibin++) {
        TString label("");
        if(ibin==h->GetXaxis()->GetNbins()) label +="#geq";
        else                                label +="=";
        label += (ibin+1);
        h->GetXaxis()->SetBinLabel(ibin,label);
    }


    h=(TH1F *)mon.addHistogram( new TH1F("npfjets",  ";Jet multiplicity (#it{p}_{T}>30 GeV);Events",5,0,5) );
    for(int ibin=1; ibin<=h->GetXaxis()->GetNbins(); ibin++) {
        TString label("");
        if(ibin==h->GetXaxis()->GetNbins()) label +="#geq";
        else                                label +="=";
        label += (ibin-1);
        h->GetXaxis()->SetBinLabel(ibin,label);
    }


    h=(TH1F *)mon.addHistogram( new TH1F("npfjetsbtags",    ";b-tag multiplicity;Events",5,0,5) );
    for(int ibin=1; ibin<=h->GetXaxis()->GetNbins(); ibin++) {
        TString label("");
        if(ibin==h->GetXaxis()->GetNbins()) label +="#geq";
        else                                label +="=";
        label += (ibin-1);
        h->GetXaxis()->SetBinLabel(ibin,label);
    }


    mon.addHistogram( new TH1D( "balancedif",   ";|E_{T}^{miss}-#it{q}_{T}|/#it{q}_{T};Events", 5,0,1.0) );
    mon.addHistogram( new TH1D( "balance",	";E_{T}^{miss}/#it{q}_{T};Events", 25,0,5) );
    mon.addHistogram( new TH1F( "met_metSB",    ";E_{T}^{miss} [GeV];Events", 50,0,500) );
    mon.addHistogram( new TH1F( "met_redMetSB", ";Reduced E_{T}^{miss} [GeV];Events", 50,0,300) );

    double newMETBin[21]= {0,10,20,30,40,50,60,70,80,90,100,120,140,160,180,200,250,300,350,400,500};
    mon.addHistogram( new TH1F( "met_met",      ";E_{T}^{miss} [GeV];Events", 20, newMETBin));
    mon.addHistogram( new TH1F( "met_redMet",   ";Reduced E_{T}^{miss} [GeV];Events", 50,0,500) );
    double redMet_bins[] = {0., 20., 40., 60., 80., 100., 150., 400., 800.};
    const int n_redMet_bins = sizeof(redMet_bins)/sizeof(double) - 1;
    mon.addHistogram( new TH1F( "met_redMet_rebin",";Reduced E_{T}^{miss} [GeV];Events", n_redMet_bins,redMet_bins) );


    // final distributions
    mon.addHistogram( new TH1F( "mt_final", 		";#it{m}_{T} [GeV];Events", 8,200,1000) );
    mon.addHistogram( new TH1F( "new_mt_final", 	";#it{m}_{T} [GeV];Events", 10,0,1000) );
    mon.addHistogram( new TH1F( "zpt_final",		";#it{p}_{T}^{ll} [GeV];Events", 25,0,500) );
    mon.addHistogram( new TH1F( "zpt_rebin_final",	";#it{p}_{T}^{ll} [GeV];Events", n_zpt_bins, zpt_bins) );
    mon.addHistogram( new TH1F( "met_met_final", 	";E_{T}^{miss} [GeV];Events", 25,0,250) );
    mon.addHistogram( new TH1F( "met_redMet_final",	";Reduced E_{T}^{miss} [GeV];Events", 50,0,500) );
    mon.addHistogram( new TH1F( "met_redMetL_final",	";Longitudinal Reduced E_{T}^{miss} [GeV];Events", 25,-100,400) );
    mon.addHistogram( new TH1F( "met_redMet_rebin_final",";Reduced E_{T}^{miss} [GeV];Events", n_redMet_bins,redMet_bins) );


    //NLO EWK Correction study
    //mon.addHistogram( new TH1F( "v_pt",";#it{p}_{T} [GeV];Events", 500,0,500) );
    //mon.addHistogram( new TH1F( "v_pt_ewkcorr",";#it{p}_{T} EWK Corr [GeV];Events", 500,0,500) );



    //##############################################
    //######## STUFF FOR CUTS OPTIMIZATION  ########
    //##############################################

    /*
          //run full optimization
           std::vector<double> optim_Cuts1_met;
           std::vector<double> optim_Cuts1_balance;
           std::vector<double> optim_Cuts1_dphi;
           std::vector<double> optim_Cuts1_zmass;
           for(double met=90;met<140;met+=10) {
             for(double balance=0.05;balance<=0.4;balance+=0.05) {
               for(double dphi=2.5;dphi<=3.0;dphi+=0.1) {
                 for(double zm=10;zm<=20;zm+=2.5)
        	 	{
        	   		optim_Cuts1_met	.push_back(met);
        	  		optim_Cuts1_balance	.push_back(balance);
        	   		optim_Cuts1_dphi	.push_back(dphi);
        	   		optim_Cuts1_zmass	.push_back(zm);
        	  	}
               }
             }
           }
    */

    //optimization
    std::vector<double> optim_Cuts1_met;
    std::vector<double> optim_Cuts1_balance;
    std::vector<double> optim_Cuts1_dphi;
    std::vector<double> optim_Cuts1_zmass;

    //fast run fixed limit
    for(double met=90; met<=120; met+=10) {
        for(double balance=0.2; balance<=0.25; balance+=0.05) {
            for(double dphi=2.7; dphi<=2.7; dphi+=0.1) {
                for(double zm=15; zm<=15; zm+=2.5) {
                    //for(double zm=0.; zm<=3.14; zm+=0.1) {
                    optim_Cuts1_met     .push_back(met);
                    optim_Cuts1_balance .push_back(balance);
                    optim_Cuts1_dphi    .push_back(dphi);
                    optim_Cuts1_zmass   .push_back(zm);
                }
            }
        }
    }


    //make it as a TProfile so hadd does not change the value
    TProfile* Hoptim_cuts1_met      =  (TProfile*) mon.addHistogram( new TProfile ("optim_cut1_met",    ";cut index;met",     optim_Cuts1_met.size(),0,optim_Cuts1_met.size()) ) ;
    TProfile* Hoptim_cuts1_balance  =  (TProfile*) mon.addHistogram( new TProfile ("optim_cut1_balance",";cut index;dphi",    optim_Cuts1_met.size(),0,optim_Cuts1_met.size()) ) ;
    TProfile* Hoptim_cuts1_dphi     =  (TProfile*) mon.addHistogram( new TProfile ("optim_cut1_dphi",    ";cut index;dphi",  optim_Cuts1_met.size(),0,optim_Cuts1_met.size()) ) ;
    TProfile* Hoptim_cuts1_zmass    =  (TProfile*) mon.addHistogram( new TProfile ("optim_cut1_zm",     ";cut index;zmass",   optim_Cuts1_met.size(),0,optim_Cuts1_met.size()) ) ;
    for(unsigned int index=0; index<optim_Cuts1_met.size(); index++) {
        Hoptim_cuts1_met        ->Fill(index, optim_Cuts1_met[index]);
        Hoptim_cuts1_balance	->Fill(index, optim_Cuts1_balance[index]);
        Hoptim_cuts1_dphi		->Fill(index, optim_Cuts1_dphi[index]);
        Hoptim_cuts1_zmass		->Fill(index, optim_Cuts1_zmass[index]);
    }

    TH1F* Hoptim_systs  =  (TH1F*) mon.addHistogram( new TH1F ("optim_systs"    , ";syst;", nvarsToInclude,0,nvarsToInclude) ) ;

    double BalanceXaxis[4], DphiLLXaxis[4], CosThetaXaxis[4];
    BalanceXaxis[0] = 0.8;
    BalanceXaxis[1] = 0.9;
    BalanceXaxis[2] = 1.1;
    BalanceXaxis[3] = 1.2;
    DphiLLXaxis[0] = 0.;
    DphiLLXaxis[1] = 0.6;
    DphiLLXaxis[2] = 1.4;
    DphiLLXaxis[3] = TMath::Pi();
    CosThetaXaxis[0] = -1.;
    CosThetaXaxis[1] = -0.5;
    CosThetaXaxis[2] = 0.;
    CosThetaXaxis[3] = 1.;

    const int nBinMVA = 11;
    double xbins[nBinMVA+1] = {0, 250, 300, 350, 400, 450, 500, 550, 600, 700, 800, 1200};

    
    // non-resonant background control
    std::vector<TString> allshapesVars;
    allshapesVars.push_back("redMet_shapes");
    allshapesVars.push_back("redMet_rebin_shapes");
    allshapesVars.push_back("mt_shapes");
    allshapesVars.push_back("new_mt_shapes");
    allshapesVars.push_back("dphi_shapes");
    allshapesVars.push_back("met_shapes");
    allshapesVars.push_back("met_rebin_shapes");
    allshapesVars.push_back("zpt_shapes");
    allshapesVars.push_back("zpt_rebin_shapes");
    allshapesVars.push_back("balance_mt_shapes");
    allshapesVars.push_back("dphiLL_mt_shapes");
    allshapesVars.push_back("coslZ_mt_shapes");


    for(size_t ivar=0; ivar<nvarsToInclude; ivar++) {
        Hoptim_systs->GetXaxis()->SetBinLabel(ivar+1, varNames[ivar]);

        //1-D shapes for limit setting
        mon.addHistogram( new TH2F (TString("redMet_shapes")+varNames[ivar],";cut index;red-MET [GeV];#Events (/5GeV)",optim_Cuts1_met.size(),0,optim_Cuts1_met.size(),80,0,800) );
        mon.addHistogram( new TH2F (TString("redMet_rebin_shapes")+varNames[ivar],";cut index;red-MET [GeV];#Events (/5GeV)",optim_Cuts1_met.size(),0,optim_Cuts1_met.size(), n_redMet_bins,redMet_bins) );
        mon.addHistogram( new TH2F (TString("mt_shapes")+varNames[ivar],";cut index; #it{m}_{T} [GeV];#Events (/100GeV)",optim_Cuts1_met.size(),0,optim_Cuts1_met.size(),12,0,1200) );
        mon.addHistogram( new TH2F (TString("new_mt_shapes")+varNames[ivar],";cut index; #it{m}_{T} [GeV];#Events (/100GeV)",optim_Cuts1_met.size(),0,optim_Cuts1_met.size(), nBinMVA, xbins) );//12,0,1200) );
        mon.addHistogram( new TH2F (TString("dphi_shapes")+varNames[ivar],";cut index; #Delta#phi(Z,E_{T}^{miss}) [rad];#events",optim_Cuts1_met.size(),0,optim_Cuts1_met.size(),5,2.5,TMath::Pi()) );
        mon.addHistogram( new TH2F (TString("met_shapes")+varNames[ivar],";cut index;MET [GeV];#events (/5GeV)",optim_Cuts1_met.size(),0,optim_Cuts1_met.size(), 160,0,800) );
        mon.addHistogram( new TH2F (TString("met_rebin_shapes")+varNames[ivar],";cut index;MET [GeV];#events (/5GeV)",optim_Cuts1_met.size(),0,optim_Cuts1_met.size(), n_redMet_bins,redMet_bins) );
        mon.addHistogram( new TH2F (TString("zpt_shapes")+varNames[ivar],";cut index;Z #it{p}_{T} [GeV];#events (/5GeV)",optim_Cuts1_met.size(),0,optim_Cuts1_met.size(), 160,0,800) );
        mon.addHistogram( new TH2F (TString("zpt_rebin_shapes")+varNames[ivar],";cut index;Z #it{p}_{T} [GeV];#events (/5GeV)",optim_Cuts1_met.size(),0,optim_Cuts1_met.size(), n_zpt_bins, zpt_bins) );

        //2-D shapes for limit setting
        TH2F *hh=(TH2F *) mon.addHistogram( new TH2F (TString("balance_mt_shapes")+varNames[ivar],";cut index; #it{m}_{T} [GeV];#Events (/100GeV)",optim_Cuts1_met.size(),0,optim_Cuts1_met.size(),36,0,3600) );
        for(int j=1; j<=36; j++) hh->GetYaxis()->SetBinLabel(j,"");
        hh->GetYaxis()->SetBinLabel(1 ,"0.8<E_{T}^{miss}/#it{q}_{T}<0.9");
        hh->GetYaxis()->SetBinLabel(13,"0.9<E_{T}^{miss}/#it{q}_{T}<1.1");
        hh->GetYaxis()->SetBinLabel(25,"1.1<E_{T}^{miss}/#it{q}_{T}<1.2");
        hh->GetYaxis()->SetBinLabel(5 ,"0<#it{m}_{T}<1200");
        hh->GetYaxis()->SetBinLabel(17,"0<#it{m}_{T}<1200");
        hh->GetYaxis()->SetBinLabel(29,"0<#it{m}_{T}<1200");

        double dphimtshapes_xbins[34];
        for(int i=0; i<34; i++) {
            if(i<=11) dphimtshapes_xbins[i] = xbins[i];
            else if(i<=22) dphimtshapes_xbins[i] = xbins[i-11]+1200;
            else if(i<=33) dphimtshapes_xbins[i] = xbins[i-22]+2400;
        }

        hh=(TH2F *) mon.addHistogram( new TH2F (TString("dphiLL_mt_shapes")+varNames[ivar],";cut index; #it{m}_{T} [GeV];#Events (/100GeV)",optim_Cuts1_met.size(),0,optim_Cuts1_met.size(),33,dphimtshapes_xbins) ); //36,0,3600) );
        for(int j=1; j<=33; j++) hh->GetYaxis()->SetBinLabel(j,"");
        hh->GetYaxis()->SetBinLabel(1 ,"0<#Delta#phi_{ll}<0.6");
        hh->GetYaxis()->SetBinLabel(13,"0.6<#Delta#phi_{ll}<1.4");
        hh->GetYaxis()->SetBinLabel(25,"1.4<#Delta#phi_{ll}<#pi");
        hh->GetYaxis()->SetBinLabel(5 ,"0<#it{m}_{T}<1200");
        hh->GetYaxis()->SetBinLabel(17,"0<#it{m}_{T}<1200");
        hh->GetYaxis()->SetBinLabel(29,"0<#it{m}_{T}<1200");


        hh=(TH2F *) mon.addHistogram( new TH2F (TString("coslZ_mt_shapes")+varNames[ivar],";cut index; #it{m}_{T} [GeV];#Events (/100GeV)",optim_Cuts1_met.size(),0,optim_Cuts1_met.size(),36,0,3600) );
        for(int j=1; j<=36; j++) hh->GetYaxis()->SetBinLabel(j,"");
        hh->GetYaxis()->SetBinLabel(1 ,"-1.0<cos#theta*<-0.5");
        hh->GetYaxis()->SetBinLabel(13,"-0.5<cos#theta*<0.0");
        hh->GetYaxis()->SetBinLabel(25,"0.0<cos#theta*<1.0");
        hh->GetYaxis()->SetBinLabel(5 ,"0<#it{m}_{T}<1200");
        hh->GetYaxis()->SetBinLabel(17,"0<#it{m}_{T}<1200");
        hh->GetYaxis()->SetBinLabel(29,"0<#it{m}_{T}<1200");


        
        // non-resonant background control
        for(size_t j=0; j<allshapesVars.size(); j++)
        {
            TH2F *h2=(TH2F *) mon.addHistogram( new TH2F (allshapesVars[j]+"_NRBctrl"+varNames[ivar],";cut index;Selection region;Events",optim_Cuts1_met.size(),0,optim_Cuts1_met.size(),6,0,6) );
            h2->GetYaxis()->SetBinLabel(1,"M_{in}^{ll}/=0 b-tags");
            h2->GetYaxis()->SetBinLabel(2,"M_{out}^{ll}/=0 b-tags");
            h2->GetYaxis()->SetBinLabel(3,"M_{out+}^{ll}/=0 b-tags");
            h2->GetYaxis()->SetBinLabel(4,"M_{in}^{ll}/#geq 1 b-tag");
            h2->GetYaxis()->SetBinLabel(5,"M_{out}^{ll}/#geq 1 b-tag");
            h2->GetYaxis()->SetBinLabel(6,"M_{out+}^{ll}/#geq 1 b-tag");
        }

    }




    //##############################################
    //######## GET READY FOR THE EVENT LOOP ########
    //##############################################

    //open the file and get events tree
    ZZ2l2nuSummaryHandler evSummaryHandler;
    TFile *file = TFile::Open(url);
    printf("Looping on %s\n",url.Data());
    if(file==0) return -1;
    if(file->IsZombie()) return -1;
    if( !evSummaryHandler.attachToTree( (TTree *)file->Get(dirname) ) ) {
        file->Close();
        return -1;
    }


    //check run range to compute scale factor (if not all entries are used)
    const Int_t totalEntries= evSummaryHandler.getEntries();
    float rescaleFactor( evEnd>0 ?  float(totalEntries)/float(evEnd-evStart) : -1 );
    if(evEnd<0 || evEnd>evSummaryHandler.getEntries() ) evEnd=totalEntries;
    if(evStart > evEnd ) {
        file->Close();
        return -1;
    }

    //MC normalization (to 1/pb)
    float cnorm=1.0;
    if(isMC) {
        TH1F* cutflowH = (TH1F *) file->Get("dataAnalyzer/llvv/cutflow");
        if(cutflowH) cnorm=cutflowH->GetBinContent(1);
        if(rescaleFactor>0) cnorm /= rescaleFactor;
        printf("cnorm = %f\n",cnorm);
    }
    Hcutflow->SetBinContent(1,cnorm);


    //pileup weighting: based on vtx for now...
    std::vector<double> dataPileupDistributionDouble = runProcess.getParameter< std::vector<double> >("datapileup");
    std::vector<float> dataPileupDistribution;
    for(unsigned int i=0; i<dataPileupDistributionDouble.size(); i++) {
        dataPileupDistribution.push_back(dataPileupDistributionDouble[i]);
    }
    std::vector<float> mcPileupDistribution;
    bool useObservedPU(true);
    //bool useObservedPU(use2011Id);
    if(!use2011Id && url.Contains("toZZto2L")) useObservedPU=true;
    if(isMC) {
        TString puDist("dataAnalyzer/llvv/pileuptrue");
        if(useObservedPU) puDist="dataAnalyzer/llvv/pileup";
        TH1F* histo = (TH1F *) file->Get(puDist);
        if(!histo)std::cout<<"pileup histogram is null!!!\n";
        for(int i=1; i<=histo->GetNbinsX(); i++) {
            mcPileupDistribution.push_back(histo->GetBinContent(i));
        }
        delete histo;
        if(dataPileupDistribution.size()==0) dataPileupDistribution=mcPileupDistribution;
    }
    while(mcPileupDistribution.size()<dataPileupDistribution.size())  mcPileupDistribution.push_back(0.0);
    while(mcPileupDistribution.size()>dataPileupDistribution.size())dataPileupDistribution.push_back(0.0);

    gROOT->cd();  //THIS LINE IS NEEDED TO MAKE SURE THAT HISTOGRAM INTERNALLY PRODUCED IN LumiReWeighting ARE NOT DESTROYED WHEN CLOSING THE FILE
    edm::LumiReWeighting *LumiWeights=0;
    PuShifter_t PuShifters;
    if(isMC) {
        LumiWeights= new edm::LumiReWeighting(mcPileupDistribution,dataPileupDistribution);
        PuShifters=getPUshifters(dataPileupDistribution,0.05);
    }

    //event Categorizer
    //EventCategory eventCategoryInst(4); //jet(0,>=1)+vbf binning
    EventCategory eventCategoryInst(1);   //jet(0,1,>=2) binning


    LeptonEfficiencySF lsf(use2011Id ? 2011:2012);

    //##############################################################################################################################
    //################################################           EVENT LOOP         ################################################
    //##############################################################################################################################

    //loop on all the events
    printf("Progressing Bar     :0%%       20%%       40%%       60%%       80%%       100%%\n");
    printf("Scanning the ntuple :");
    int treeStep = (evEnd-evStart)/50;
    if(treeStep==0)treeStep=1;
    DuplicatesChecker duplicatesChecker;
    int nDuplicates(0);
    for( int iev=evStart; iev<evEnd; iev++) {
        if((iev-evStart)%treeStep==0) {
            printf(".");
            fflush(stdout);
        }

        //##############################################   EVENT LOOP STARTS   ##############################################
        //load the event content from tree
        evSummaryHandler.getEntry(iev);
        ZZ2l2nuSummary_t &ev=evSummaryHandler.getEvent();
        if(!isMC && duplicatesChecker.isDuplicate( ev.run, ev.lumi, ev.event) ) {
            nDuplicates++;
            continue;
        }
        PhysicsEvent_t phys=getPhysicsEventFrom(ev);

        //event category
        bool isSameFlavor(ev.cat==MUMU || ev.cat==EE);
        TString tag_cat;
        switch(ev.cat) {
        case MUMU :
            tag_cat = "mumu";
            break;
        case EE   :
            tag_cat = "ee";
            break;
        case EMU  :
            tag_cat = "emu";
            break;
        default   :
            continue;
        }
        //      if(isMC && mctruthmode==1 && !isDYToLL(ev.mccat) && !isZZ2l2nu(ev.mccat) ) continue;
        if(isMC && mctruthmode==1 && !isDYToLL(ev.mccat) ) continue;
        if(isMC && mctruthmode==2 && !isDYToTauTau(ev.mccat) ) continue;

        //require compatibilitiy of the event with the PD
        bool hasTrigger(false);
        bool hasEEtrigger = ev.triggerType & 0x1;
        bool hasMMtrigger = (ev.triggerType >> 1 ) & 0x1;
        bool hasEMtrigger = (ev.triggerType >> 2 ) & 0x1;
        bool hasMtrigger  = (ev.triggerType >> 3 ) & 0x1;
	bool hasEtrigger  = (ev.triggerType >> 4 ) & 0x1;
        if(!isMC) {
            if(ev.cat!=fType) continue;

            if(ev.cat==EE   && !(hasEEtrigger||hasEtrigger) ) continue;
            if(ev.cat==MUMU && !(hasMMtrigger||hasMtrigger) ) continue;
            //if(ev.cat==EMU  && !hasEMtrigger) continue;
	    if(ev.cat==EMU  && !hasEMtrigger && !(hasEtrigger && hasMtrigger) ) continue;

            //this is a safety veto for the single mu PD
            if(isSingleMuPD) {
                if(!hasMtrigger) continue;
                if(hasMtrigger && hasMMtrigger) continue;
            }

            //this is a safety veto for the single Ele PD
            if(isSingleElePD) {
                if(!hasEtrigger) continue;
                if(hasEtrigger && hasEEtrigger) continue;
            }

            hasTrigger=true;
        } else {
            if(ev.cat==EE   && (hasEEtrigger || hasEtrigger) ) hasTrigger=true;
            if(ev.cat==MUMU && (hasMMtrigger || hasMtrigger) ) hasTrigger=true;
            if(ev.cat==EMU  && (hasEMtrigger || (hasEtrigger && hasMtrigger)) ) hasTrigger=true;
            if(use2011Id) hasTrigger = true; //for 2011
        }

        //prepare the tag's vectors for histo filling
        std::vector<TString> tags(1,"all");

        //pileup weight
        float weight = 1.0;
        double TotalWeight_plus = 1.0;
        double TotalWeight_minus = 1.0;
        if(isMC) {
            weight            = LumiWeights->weight(useObservedPU ? ev.ngenITpu : ev.ngenTruepu);
            TotalWeight_plus  = PuShifters[PUUP]->Eval(useObservedPU ? ev.ngenITpu : ev.ngenTruepu);
            TotalWeight_minus = PuShifters[PUDOWN]->Eval(useObservedPU ? ev.ngenITpu : ev.ngenTruepu);
        }
        Hcutflow->Fill(1,1);
        Hcutflow->Fill(2,weight);
        Hcutflow->Fill(3,weight*TotalWeight_minus);
        Hcutflow->Fill(4,weight*TotalWeight_plus);


        //#################################
        //######## EWK Correction  ########
        //#################################
        bool isZH(false);
        bool isWZ(false);
        bool isZZ(false);
        bool isZjets(false);
        bool isMC_HZZd(false);
        if(isMC) { // need to put this an external function
            TString input_File = url.Data();
            if(input_File.Contains("TeV_ZH")) isZH=true;
            if(input_File.Contains("TeV_WZ")) isWZ=true;
            if(input_File.Contains("TeV_ZZ")) isZZ=true;
            if(input_File.Contains("TeV_DYJetsToLL") || input_File.Contains("TeV_ZG")) isZjets=true;
	    if(input_File.Contains("TeV_HZZd")) isMC_HZZd=true;


            if(isZH) {
                double pt_zll = (phys.genleptons[0]+phys.genleptons[1]).pt();
                float ewk_nloweightsZH = myZHUtils.weightNLOEWKsignal(pt_zll);
                //mon.fillHisto("v_pt",tags,pt_zll,weight);
                weight *= ewk_nloweightsZH;
                //mon.fillHisto("v_pt_ewkcorr",tags,pt_zll,weight);
            }


	    if(isMC_HZZd){
	     	double pt_vv = (phys.genneutrinos[0]+phys.genneutrinos[1]).pt();
	  	mon.fillHisto("Zdpt",tags, pt_vv, weight);	
	    }

            if(isZZ) {
                double pt_zll = (phys.genleptons[0]+phys.genleptons[1]).pt();
                double pt_zvv = (phys.genneutrinos[0]+phys.genneutrinos[1]).pt();
                double trailing_zzpt = pt_zvv;
                if(pt_zll<pt_zvv) trailing_zzpt = pt_zll;
                double ewk_nloweightsZZ = myZHUtils.weightNLOEWKzz(trailing_zzpt);
                //mon.fillHisto("v_pt",tags,trailing_zzpt,weight);
                weight *= ewk_nloweightsZZ;
                //fprintf(outTxtFile,"ZZ:%f \n",ewk_nloweightsZZ);
                //mon.fillHisto("v_pt_ewkcorr",tags,trailing_zzpt,weight);
            }

            if(isWZ) {
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

                    double trailing_pt = pt_wlv;
                    if(pt_zll<pt_wlv) trailing_pt = pt_zll;
                    double ewk_nloweightsWZ = myZHUtils.weightNLOEWKwz(trailing_pt);
                    //mon.fillHisto("v_pt",tags,trailing_pt,weight);
                    weight *= ewk_nloweightsWZ;
                    //mon.fillHisto("v_pt_ewkcorr",tags,trailing_pt,weight);

                } //ngenleps==3 && Nneutrinos==1
            } //isWZ
        }



        //#####################################
        //######## Objects Selection  #########
        //#####################################

        //
        //MET variables
        //
        LorentzVector rawMetP4=phys.met[2];
        if(use2011Id) rawMetP4=phys.met[0]; //V3 Test
        //LorentzVector fullTypeIMetP4=phys.met[0];
        //LorentzVector mvaMetP4=phys.met[7];

        //apply JER base corrections to jets (and compute associated variations on the MET variable)
        // std PF
        std::vector<PhysicsObjectJetCollection> variedAJets;
        LorentzVectorCollection zvvs;
        PhysicsObjectJetCollection &recoJets = ( useCHS ? phys.ajets : phys.jets);
        METUtils::computeVariation(recoJets, phys.leptons, rawMetP4, variedAJets, zvvs, &jecUnc);
        if(!useJERsmearing) zvvs[0] = rawMetP4; // this stinks a bit...

        //
        // LEPTON ANALYSIS
        //
        LorentzVector lep1=phys.leptons[0];
        LorentzVector lep2=phys.leptons[1];
        //Muon scale and resolution corrections
        if(fabs(phys.leptons[0].id)==13) {
            TLorentzVector tmpLep1(lep1.Px(), lep1.Py(), lep1.Pz(), lep1.E());
            int tmpCh1 = phys.leptons[0].id < 0.0 ? -1 : 1; // "id" has same sign as charge here (opposite to "genid")
            corrector_->applyPtCorrection(tmpLep1, tmpCh1);
            if( isMC && (!use2011Id) ) corrector_->applyPtSmearing(tmpLep1, tmpCh1);
            lep1 = LorentzVector(tmpLep1.Px(), tmpLep1.Py(), tmpLep1.Pz(), tmpLep1.E());
        }
        if(fabs(phys.leptons[1].id)==13) {
            TLorentzVector tmpLep2(lep2.Px(), lep2.Py(), lep2.Pz(), lep2.E());
            int tmpCh2 = phys.leptons[0].id < 0.0 ? -1 : 1; // "id" has same sign as charge here (opposite to "genid")
            corrector_->applyPtCorrection(tmpLep2, tmpCh2);
            if( isMC && (!use2011Id) ) corrector_->applyPtSmearing(tmpLep2, tmpCh2);
            lep2 = LorentzVector(tmpLep2.Px(), tmpLep2.Py(), tmpLep2.Pz(), tmpLep2.E());
        }
        LorentzVector zll(lep1+lep2);
        double dphi2l = fabs(deltaPhi(lep1.phi(),lep2.phi()));
        //bool passdphi2l(true);
        //bool passdphi2l(cos(dphi2l)>0);

        //two opposite sign leptons
        int  id1=phys.leptons[0].id;
        int  id2=phys.leptons[1].id;
        bool passOppositeSign((id1*id2)<0);

        bool passId(true);
        passId &= passOppositeSign;
        bool passIdAndIso(true);
        passIdAndIso &= passOppositeSign;
        bool passZmass(fabs(zll.mass()-91)<15); //RJ changed

        bool isZsideBand( (zll.mass()>40 && zll.mass()<70) || (zll.mass()>110 && zll.mass()<200));
        bool isZsideBandPlus( (zll.mass()>110 && zll.mass()<200));
        bool passZpt(zll.pt()>50); //for Gamma+Jets, can also be used for control plots

        //check alternative selections for the dilepton
        double llScaleFactor(1.0),llTriggerEfficiency(1.0);
        LorentzVector genZP4(0,0,0,0); // for checks on Sherpa ZZ
        int genmatchid[2] = {-1, -1};
        double genmatchdr[2] = {0.1, 0.1};
        for(int ilep=0; ilep<2; ilep++) {
            TString lepStr( fabs(phys.leptons[ilep].id)==13 ? "mu" : "e");

            //generator level matching
            //  int matchid(0);
            //LorentzVector genP4(0,0,0,0);
            for(size_t igl=0; igl<phys.genleptons.size(); igl++) {
                //if(deltaR(phys.genleptons[igl],phys.leptons[ilep])>0.1) continue;
                if(ilep==1 && int(igl)==genmatchid[0]) continue;
                if(deltaR(phys.genleptons[igl],phys.leptons[ilep])<genmatchdr[ilep]) {
                    genmatchdr[ilep] = deltaR(phys.genleptons[igl],phys.leptons[ilep]);
                    genmatchid[ilep] = igl;
                }
            }
            if(genmatchid[0]>-1 && genmatchid[1]>-1) {
                //genP4=phys.genleptons[igl];
                //matchid=phys.genleptons[igl].id;
                genZP4 = phys.genleptons[0] + phys.genleptons[1];
            }

            //id and isolation
            int lpid=phys.leptons[ilep].pid;
            float relIso2011    = phys.leptons[ilep].relIsoRho(ev.rho);
            float relIso = (lepStr=="mu") ?
                           phys.leptons[ilep].pfRelIsoDbeta() :
                           phys.leptons[ilep].ePFRelIsoCorrected2012(ev.rho,ev.en_sceta[lpid]); //RENJIE
            std::vector<int> passIds;
            std::map<int,bool> passIsos;
            bool hasGoodId(false), isIso(false);
            if(fabs(phys.leptons[ilep].id)==13) {
                if( hasObjectId(ev.mn_idbits[lpid], MID_LOOSE) )    {
                    passIds.push_back(0);
                    passIsos[0]=(relIso<0.2);
                }
                if( hasObjectId(ev.mn_idbits[lpid], MID_TIGHT) )    {
                    passIds.push_back(1);
                    passIsos[1]=(relIso<0.2);
                    if(!use2011Id) {
                        hasGoodId=true;
                        isIso=passIsos[0];
                    }
                }
                if( hasObjectId(ev.mn_idbits[lpid], MID_VBTF2011) ) {
                    passIds.push_back(2);
                    passIsos[2]=(relIso2011<0.15);
                    if(use2011Id) {
                        hasGoodId=true;
                        isIso=passIsos[2];
                    }
                }
                if( hasObjectId(ev.mn_idbits[lpid], MID_SOFT) )     {
                    passIds.push_back(3);
                    passIsos[3]=true;
                }


                llScaleFactor *= lsf.getLeptonEfficiency(phys.leptons[ilep].pt(),fabs(phys.leptons[ilep].eta()),13).first;

                if(use2011Id) {
                    try {
                        //llScaleFactor *= muonScaleFactor(phys.leptons[ilep].pt(),fabs(phys.leptons[ilep].eta()));
                        llTriggerEfficiency *= 1.0;//muonTriggerEfficiency(phys.leptons[ilep].pt(),fabs(phys.leptons[ilep].eta()));
                    } catch(std::string &e) {
                    }
                } else {
                    //llScaleFactor *= 1;
                    llTriggerEfficiency *= 1.0; //muonTriggerEfficiency(phys.leptons[ilep].pt(),fabs(phys.leptons[ilep].eta()),2012);
                }
            } else {
                int wps[]= {EgammaCutBasedEleId::LOOSE,EgammaCutBasedEleId::MEDIUM, EID_VBTF2011, EgammaCutBasedEleId::VETO};
                llScaleFactor *= lsf.getLeptonEfficiency(phys.leptons[ilep].pt(),fabs(phys.leptons[ilep].eta()),11).first;
                for(int iwp=0; iwp<4; iwp++) {
                    if(iwp==2 && hasObjectId(ev.en_idbits[lpid], EID_VBTF2011)) {
                        passIds.push_back(2);
                        passIsos[2]=(relIso2011<0.10);
                        if(use2011Id) {
                            hasGoodId=true;
                            isIso=passIsos[2];
                            try {
                                //llScaleFactor *= electronScaleFactor(phys.leptons[ilep].pt(),fabs(phys.leptons[ilep].eta()));
                                llTriggerEfficiency *= 1.0;//electronTriggerEfficiency(phys.leptons[ilep].pt(),fabs(phys.leptons[ilep].eta()));
                            } catch(std::string &e) {
                            }
                        }
                    } else {
                        bool passWp = EgammaCutBasedEleId::PassWP(EgammaCutBasedEleId::WorkingPoint(wps[iwp]),
                                      (fabs(phys.leptons[ilep].eta())<1.4442),
                                      phys.leptons[ilep].pt(), phys.leptons[ilep].eta(),
                                      ev.en_detain[lpid],  ev.en_dphiin[lpid], ev.en_sihih[lpid], ev.en_hoe[lpid],
                                      ev.en_ooemoop[lpid], phys.leptons[ilep].d0, phys.leptons[ilep].dZ,
                                      0., 0., 0.,
                                      !hasObjectId(ev.en_idbits[lpid], EID_CONVERSIONVETO),0,ev.rho);
                        if(passWp) {
                            passIds.push_back(iwp);
                            passIsos[iwp]=(relIso<0.15);
                            if(wps[iwp]==EgammaCutBasedEleId::MEDIUM && !use2011Id) {
                                hasGoodId=true;
                                isIso=passIsos[iwp];
                            }
                        }
                        if(!use2011Id) {
                            //llScaleFactor *= 1;
                            llTriggerEfficiency *= 1.0; //electronTriggerEfficiency(phys.leptons[ilep].pt(),fabs(phys.leptons[ilep].eta()),2012);
                        }
                    }
                }
            }
            if(!hasGoodId)  {
                passId &=false;    //RJ
                passIdAndIso &=false;
            }
            if(!isIso) passIdAndIso &=false; //RJ

            //        fill control histograms (constrained to the Z mass)
            if(passZmass && isSameFlavor) {
                if(hasGoodId) {
                    mon.fillHisto(lepStr+"reliso",     tags, use2011Id? relIso2011 : relIso,   weight);
                }
            }
        }

        tags.push_back(tag_cat);
        if(tag_cat=="mumu" || tag_cat=="ee") tags.push_back("ll");


        //
        // 3rd LEPTON ANALYSIS
        //  1_ for lepton veto
        //  2_ for WZ Control Sample
        //
        bool pass3dLeptonVeto(true);
        //bool has3dLepton(false);
        //bool WZ3Lepsame(false);
        int nextraleptons(0), nextraleptons_WZ(0);
        std::vector<int> nextraleptonsid;
        std::vector<LorentzVector> extraLeptonsP4, extraLeptonsP4_WZ;
        for(size_t ilep=2; ilep<phys.leptons.size(); ilep++) {
            //lepton type
            bool isGood(false), isGood_WZ(false);
            int lpid=phys.leptons[ilep].pid;
            if(fabs(phys.leptons[ilep].id)==13) {
                if(!use2011Id) {
                    isGood = (hasObjectId(ev.mn_idbits[lpid], MID_LOOSE) && phys.leptons[ilep].pfRelIsoDbeta()<0.2);
                    isGood |= (hasObjectId(ev.mn_idbits[lpid], MID_SOFT) && phys.leptons[ilep].pt()>3);
                    isGood_WZ = (hasObjectId(ev.mn_idbits[lpid], MID_LOOSE) && phys.leptons[ilep].pfRelIsoDbeta()<0.2 && phys.leptons[ilep].pt()>20);
                } else {
                    isGood = (hasObjectId(ev.mn_idbits[lpid], MID_VBTF2011) && phys.leptons[ilep].relIsoRho(ev.rho)<0.15 && phys.leptons[ilep].pt()>10);
                    isGood |= (hasObjectId(ev.mn_idbits[lpid], MID_SOFT2011) && phys.leptons[ilep].pt()>3);
                    isGood_WZ = (hasObjectId(ev.mn_idbits[lpid], MID_VBTF2011) && phys.leptons[ilep].relIsoRho(ev.rho)<0.15 && phys.leptons[ilep].pt()>20);
                }
            } else {
                if(!use2011Id) {
                    isGood = ( hasObjectId(ev.en_idbits[lpid],EID_VETO) && phys.leptons[ilep].ePFRelIsoCorrected2012(ev.rho,ev.en_sceta[lpid])<0.15 && phys.leptons[ilep].pt()>10);
                    isGood_WZ=(hasObjectId(ev.en_idbits[lpid],EID_VETO) && phys.leptons[ilep].ePFRelIsoCorrected2012(ev.rho,ev.en_sceta[lpid])<0.15 && phys.leptons[ilep].pt()>20);
                } else {
                    isGood = ( hasObjectId(ev.en_idbits[lpid],EID_VBTF2011) && phys.leptons[ilep].relIsoRho(ev.rho)<0.1 && phys.leptons[ilep].pt()>10);
                    isGood_WZ=(hasObjectId(ev.en_idbits[lpid],EID_VBTF2011) && phys.leptons[ilep].relIsoRho(ev.rho)<0.1 && phys.leptons[ilep].pt()>20);
                }
            }
            nextraleptons += isGood;
            nextraleptons_WZ += isGood_WZ;
            if(isGood_WZ) nextraleptonsid.push_back(phys.leptons[ilep].id);
            //WZ3Lepsame = ( fabs(id1)==fabs(id2) ) && ( fabs(phys.leptons[ilep].id) == fabs(id2) ); //SAME
            //if(fabs(id1)==15 || fabs(id2)==15 || fabs(phys.leptons[ilep].id)==15) WZ3Lepsame=false; //TEST
            //Oct12. test lepton veto uncertainty

            //if(!isGood) continue;
            //Muon scale and resolution corrections
            // probably unnecessary
            LorentzVector tmpLep = phys.leptons[ilep];
            if(fabs(phys.leptons[ilep].id)==13) {
                TLorentzVector tmpTLep(tmpLep.Px(), tmpLep.Py(), tmpLep.Pz(), tmpLep.E());
                int tmpCh = phys.leptons[ilep].id < 0.0 ? -1 : 1; // "id" has same sign as charge here (opposite to "genid")
                corrector_->applyPtCorrection(tmpTLep, tmpCh);
                if( isMC && (!use2011Id) ) corrector_->applyPtSmearing(tmpTLep, tmpCh);
                tmpLep = LorentzVector(tmpTLep.Px(), tmpTLep.Py(), tmpTLep.Pz(), tmpTLep.E());
            }
            if(isGood) extraLeptonsP4.push_back(tmpLep);
            if(isGood_WZ) extraLeptonsP4_WZ.push_back(tmpLep);
            //extraLeptonsP4.push_back( phys.leptons[ilep] );
        }
        pass3dLeptonVeto=(nextraleptons==0);
        //has3dLepton=(nextraleptons_WZ>0);

        //
        //STD PF JET ANALYSIS
        //
        bool passJetveto(true);
        bool passBveto(true);
        bool passRedMet(true);
        bool passDphijmet(true);
        bool passBalanceCut(true);
        PhysicsObjectJetCollection &aJets = ( useJERsmearing ? variedAJets[0] : recoJets );
        PhysicsObjectJetCollection aGoodIdJets;
        LorentzVector aClusteredMetP4(zll);
        aClusteredMetP4 *= -1;
        LorentzVector Recoil(zvvs[0]);
        Recoil *= -1;
        Recoil -= zll;
        int nABtags(0),nAJetsGood30(0),nAJetsGood15(0), nCSVMtags(0), nCSVTtags(0);
        float mindphijmet(999999.),mindphijmet15(999999.);
        for(size_t ijet=0; ijet<aJets.size(); ijet++) {
            if(aJets[ijet].pt()<15) continue;
            nAJetsGood15++;

            float idphijmet( fabs(deltaPhi(aJets[ijet].phi(),zvvs[0].phi()) ) );
            if(aJets[ijet].pt()>15) if(idphijmet<mindphijmet15)  mindphijmet15=idphijmet;
            if(aJets[ijet].pt()>30) if(idphijmet<mindphijmet)  mindphijmet=idphijmet;

            bool isGoodJet = hasObjectId(aJets[ijet].pid,JETID_LOOSE);
            if(usePUsubJetId) isGoodJet = hasObjectId(aJets[ijet].pid,JETID_CUTBASED_LOOSE);
            if(!isGoodJet) continue;


            aClusteredMetP4 -= aJets[ijet];

            if(useJetsOnlyInTracker && fabs(aJets[ijet].eta())>2.5) continue;

            aGoodIdJets.push_back(aJets[ijet]);
            if(aJets[ijet].pt()>30)nAJetsGood30++;

            if(aJets[ijet].pt()>20 /*&& fabs(aJets[ijet].eta())<2.5*/)  nABtags += (aJets[ijet].btag6>0.244);
            if(aJets[ijet].pt()>20 /*&& fabs(aJets[ijet].eta())<2.5*/)  nCSVMtags += (aJets[ijet].btag6>0.679);
            if(aJets[ijet].pt()>20 /*&& fabs(aJets[ijet].eta())<2.5*/)  nCSVTtags += (aJets[ijet].btag6>0.898);
        }
        //passJetveto=(nAJetsGood30==0);
        passBveto=(nABtags==0);
        //bool passMBveto=(nCSVMtags==0);
        //bool passTBveto=(nCSVTtags==0);
        //passDphijmet=(mindphijmet>0.5);
        //if(nAJetsGood30==0) passDphijmet=(mindphijmet15>0.5);
        //passBalanceCut=(zvvs[0].pt()/zll.pt()>0.8 && zvvs[0].pt()/zll.pt()<1.2);
        passBalanceCut=(zvvs[0].pt()/zll.pt()>0.75 && zvvs[0].pt()/zll.pt()<1.25); //Guillelmo's method
        //bool passBalanceCutWWCtrl=(zvvs[0].pt()/zll.pt()>0.4 && zvvs[0].pt()/zll.pt()<1.8);


        //ad-hoc cut for obvious correlations between MET and a lepton
        //double dphil1met=fabs(deltaPhi(lep1.phi(),zvvs[0].phi()));
        //double dphil2met=fabs(deltaPhi(lep2.phi(),zvvs[0].phi()));
        bool passLMetVeto(true);
        //if(!use2011Id && zvvs[0].pt()>60 && min(dphil1met,dphil2met)<0.2) passLMetVeto=false;

        //other mets
        METUtils::stRedMET aRedMetOut;
        LorentzVector null(0,0,0,0);
        LorentzVector aRedMet=METUtils::redMET(METUtils::INDEPENDENTLYMINIMIZED, zll, 0, null, 0, aClusteredMetP4, zvvs[0], true, &aRedMetOut);
        double aRedMetL=aRedMetOut.redMET_l;
        //double aRedMetT=aRedMetOut.redMET_t;
        TVector2 aClusteredMet2(aClusteredMetP4.px(),aClusteredMetP4.py());
        //double clusteredMetL=aRedMetOut.a_l*aClusteredMet2;
        //double clusteredMetT=aRedMetOut.a_t*aClusteredMet2;
        passRedMet=(zvvs[0].pt()>120); //ReducedMET: aRedMet.pt(); PFMET: zvvs[0].pt()
        //bool passRedMetWWCtrl=(aRedMet.pt()>65);

        TVector2 RecoilMet2(Recoil.px(),Recoil.py());
        //double RecoilMetL=aRedMetOut.a_l*RecoilMet2;
        //double RecoilMetT=aRedMetOut.a_t*RecoilMet2;


        // compute dphi(zpt,redMet)
        double dphiZllmet=fabs(deltaPhi(zll.phi(),zvvs[0].phi()));
        //double dphiMet_mvaMet=fabs(deltaPhi(mvaMetP4.phi(),zvvs[0].phi()));
        //double dphiZllredMet=fabs(deltaPhi(zll.phi(),aRedMet.phi()));
        bool passdphiZllmetCut(dphiZllmet>2.7);///2.6);
        //bool passdphiMetmvaMet(dphiMet_mvaMet<0.2);

        //transverse masses
        double aMT=METUtils::transverseMass(zll,zvvs[0],true);
        double aMTmassless=METUtils::transverseMass(zll,zvvs[0],false);
        double balanceDif=fabs(zvvs[0].pt()-zll.pt())/zll.pt();
        TVector2 dil2(zll.px(),zll.py());
        TVector2 met2(zvvs[0].px(),zvvs[0].py());
        double axialMet=dil2*met2;
        axialMet /= -zll.pt();
        //double pfMetCompL = aRedMetOut.a_l*met2;
        //double pfMetCompT = aRedMetOut.a_t*met2;

        //Collins-Soper Frame
        double colin_soper = myZHUtils.Collins_Soper(lep1,lep2);


        //RJ
        bool passMTcut(aMT>220 && aMT<1200);
        //passMTcut &= passdphiMetmvaMet;
        //passMTcut &= passdphi2l;


        TString evCat(tag_cat);
        {
            int eventSubCat  = eventCategoryInst.Get(phys,&aGoodIdJets);
            TString tag_subcat = eventCategoryInst.GetLabel(eventSubCat);
            evCat+=tag_subcat;
        }

        
        fprintf(outTxtFile_full,"%d | %d |  %d | %d%d%d%d%d%d%d%d%d%d%d | %s ",ev.run,ev.lumi,ev.event,
                passId,passIdAndIso,passZmass,passZpt,pass3dLeptonVeto,passBveto,passLMetVeto,passdphiZllmetCut,
                passRedMet,passBalanceCut,passMTcut,evCat.Data());
	for(unsigned int i=0; i<aJets.size(); i++) {fprintf(outTxtFile_full,"| %f",aJets[i].pt());}
		
        fprintf(outTxtFile_full,"| %d",nAJetsGood30);
        fprintf(outTxtFile_full,"\n");
        




        //#########################################################
        //####  RUN PRESELECTION AND CONTROL REGION PLOTS  ########
        //#########################################################

        if(isMC && use2011Id) weight *= llScaleFactor*llTriggerEfficiency;
        if(hasTrigger)                 {
            mon.fillHisto("eventflow",tags,0,weight);
        }
        if(hasTrigger && passId)       {
            mon.fillHisto("eventflow",tags,1,weight);
        }
        if(hasTrigger && passIdAndIso) {
            mon.fillHisto("eventflow",tags,2,weight);
        } else continue;

        //event category
        int eventSubCat  = eventCategoryInst.Get(phys,&aGoodIdJets);
        TString tag_subcat = eventCategoryInst.GetLabel(eventSubCat);

        tags.push_back(tag_cat+tag_subcat);
        if(tag_subcat=="eq0jets" || tag_subcat=="eq1jets") tags.push_back(tag_cat+"lesq1jets");
        if(tag_cat=="mumu" || tag_cat=="ee") {
            tags.push_back("ll"+tag_subcat);
            if(tag_subcat=="eq0jets" || tag_subcat=="eq1jets") tags.push_back("lllesq1jets");
        }



	mon.fillHisto("zmass_raw",       tags, zll.mass(), weight);
	//
        // Reweighting MC DY samples
	//
        mon.fillHisto("MllvsPt", tags, zll.pt(), zll.mass(), weight);

        if(isZjets) {
            double DY_weights = myZHUtils.get2DWeights(zll.pt(), zll.mass(),"MllvsPt",tag_cat+tag_subcat);
            weight *= DY_weights;
        }
	mon.fillHisto("zmass",       tags, zll.mass(), weight);

        // new method estimate Z+jets: ABCD
        if(passZmass && passZpt && pass3dLeptonVeto && passBveto && passLMetVeto) {
            mon.fillHisto("dPhivsMET_ABCD", tags, zvvs[0].pt(), dphiZllmet, weight);
            if(zvvs[0].pt() > 120  && dphiZllmet > 2.7)  mon.fillHisto("dPhivsMET_A", tags, zvvs[0].pt(), dphiZllmet, weight);
            if(zvvs[0].pt() > 120  && dphiZllmet <= 2.7) mon.fillHisto("dPhivsMET_B", tags, zvvs[0].pt(), dphiZllmet, weight);
            if(zvvs[0].pt() <= 120 && dphiZllmet > 2.7)  mon.fillHisto("dPhivsMET_C", tags, zvvs[0].pt(), dphiZllmet, weight);
            if(zvvs[0].pt() <= 120 && dphiZllmet <= 2.7) mon.fillHisto("dPhivsMET_D", tags, zvvs[0].pt(), dphiZllmet, weight);

            double METBal = fabs(1-zvvs[0].pt()/zll.pt());
            mon.fillHisto("dPhivsBal_ABCD", tags, METBal, dphiZllmet, weight);
            if(METBal < 0.25  && dphiZllmet > 2.7) mon.fillHisto("dPhivsBal_A", tags, METBal, dphiZllmet, weight);
            if(METBal < 0.25  && dphiZllmet <= 2.7) mon.fillHisto("dPhivsBal_B", tags, METBal, dphiZllmet, weight);
            if(METBal >= 0.25  && dphiZllmet > 2.7) mon.fillHisto("dPhivsBal_C", tags, METBal, dphiZllmet, weight);
            if(METBal >= 0.25  && dphiZllmet <= 2.7) mon.fillHisto("dPhivsBal_D", tags, METBal, dphiZllmet, weight);
        }



        //##############################################
        //########  Main Event Selection        ########
        //##############################################


        if(passZmass) {
            mon.fillHisto("nvtx"     ,   tags, ev.nvtx,      weight);
            mon.fillHisto("nvtxraw"  ,   tags, ev.nvtx,      1);
            mon.fillHisto("zpt"      ,   tags, zll.pt(),     weight);
            mon.fillHisto("zpt_rebin",   tags, zll.pt(),     weight);
	    mon.fillHisto("CoslepZ_CS",  tags, colin_soper,  weight);

            if(passZpt) {
                //analyze lepton kinematics
                LorentzVector leadingLep(phys.leptons[0].pt()>phys.leptons[1].pt() ? phys.leptons[0]: phys.leptons[1]);
                LorentzVector trailerLep(phys.leptons[0].pt()>phys.leptons[1].pt() ? phys.leptons[1]: phys.leptons[0]);
                mon.fillHisto("nleptons",tags,2+nextraleptons,weight);

                if(pass3dLeptonVeto) {
                    //final jet control
                    mon.fillHisto("npfjets",              tags, nAJetsGood30,weight);

                    // jet veto
                    if(passJetveto) { //set passJetveto always true
                        mon.fillHisto("npfjetsbtags",  tags, nABtags ,weight);

                        //b-veto
                        if(passBveto) {
                            mon.fillHisto("qt", tags, zll.pt(), weight, true);
                            mon.fillHisto("qmass",tags, zll.mass(),weight);
                            mon.fillHisto("qt_rebin", tags, zll.pt(), weight, true);
                            mon.fillHisto("qt_rebin_new", tags, zll.pt(), weight, true);
                            mon.fillHisto("nvtx_dy"     ,   tags, ev.nvtx,      weight);
                            mon.fillHisto("nvtx_dy_rebin"     ,   tags, ev.nvtx,      weight);
                            mon.fillHisto("met_met_phi_dy",   tags,zvvs[0].phi(),weight);

                            if(passDphijmet && passLMetVeto) {
                                //RJ, Data-Driven DY sample
                                mon.fillHisto("met_met",         tags,  zvvs[0].pt(),  weight,true);
                                mon.fillHisto("met_redMet",tags,aRedMet.pt(),weight);
                                mon.fillHisto("met_redMet_rebin",tags,aRedMet.pt(),weight);
                                mon.fillHisto("balance",tags, zvvs[0].pt()/zll.pt(),weight);
                                mon.fillHisto("balancedif",tags, balanceDif,weight);

                                if(passdphiZllmetCut) {
                                    mon.fillHisto("deltaleptonpt",tags, leadingLep.pt()-trailerLep.pt()    ,weight);
                                    mon.fillHisto("deltazpt",tags, zll.pt()-zvvs[0].pt(),weight);

                                    if(passRedMet) {

                                        if(passBalanceCut) {   //RJ

                                            if(passMTcut) { //RJ

                                                // Final distributions
                                                mon.fillHisto("mt_final",          tags, aMT,          weight);
                                                mon.fillHisto("new_mt_final",	   tags, aMTmassless,  weight);
                                                mon.fillHisto("zpt_final",         tags, zll.pt(),     weight);
                                                mon.fillHisto("zpt_rebin_final",   tags, zll.pt(),     weight);
                                                mon.fillHisto("met_met_final",     tags, zvvs[0].pt(), weight);
                                                mon.fillHisto("met_redMet_final",  tags, aRedMet.pt(), weight);
                                                mon.fillHisto("met_redMet_rebin_final", tags, aRedMet.pt(), weight);
                                                mon.fillHisto("met_redMetL_final", tags, aRedMetL,     weight);

                                                TString evCat(tag_cat);
                                                {
                                                    int eventSubCat  = eventCategoryInst.Get(phys,&aGoodIdJets);
                                                    TString tag_subcat = eventCategoryInst.GetLabel(eventSubCat);
                                                    evCat+=tag_subcat;
                                                }

                                                if(!evCat.Contains("geq2jets")) {
                                                    //fprintf(outTxtFile_final,"%d | %d |  %d | %d%d%d%d%d%d%d%d%d%d%d | %s \n",ev.run,ev.lumi,ev.event,
                                                    //      passId,passIdAndIso,passZmass,passZpt,pass3dLeptonVeto,passBveto,passLMetVeto,passdphiZllmetCut,
                                                    //      passRedMet,passBalanceCut,passMTcut, evCat.Data());

                                                    /*
                                                                            fprintf(outTxtFile_final,"%d | %d |  %d | %d%d%d%d%d%d%d%d%d%d%d | %s ",ev.run,ev.lumi,ev.event,
                                                                                    passId,passIdAndIso,passZmass,passZpt,pass3dLeptonVeto,passBveto,passLMetVeto,passdphiZllmetCut,
                                                                                    passRedMet,passBalanceCut,passMTcut,evCat.Data());

                                                                            if(aJets.size()>0 && zvvs.size()>0) {
                                                                                fprintf(outTxtFile_final,"| %f | %f | %f | %f",aJets[0].pt(),aJets[0].eta(), aJets[0].btag2, zvvs[0].pt());
                                                                            }
                                                                            fprintf(outTxtFile_final,"\n");
                                                    */
                                                }


                                            } //end MT cut
                                        }//end passBalanceCut
                                    } // end passRedMet
                                }//end passdphiZllmetCut
                            }//end passLMetVeto,
                        }//end passBveto
                    }//end passJetveto
                }//3lept
            }//end passZpt
        }//end passZmass



        //RENJIE, Top Control, MET control in the sideband
        bool passSB( ((zll.mass()>40 && zll.mass()<70) || (zll.mass()>110 && zll.mass()<200)) && zll.pt()>55 );
        if(passSB && pass3dLeptonVeto && passDphijmet && !passBveto && passLMetVeto) {
            mon.fillHisto("met_metSB",tags,zvvs[0].pt(),weight);
            mon.fillHisto("met_redMetSB",tags,aRedMet.pt(),weight);
        }


        

        //##############################################################################
        //### HISTOS FOR STATISTICAL ANALYSIS (include systematic variations)
        //##############################################################################

        
        //Fill histogram for posterior optimization, or for control regions
        for(size_t ivar=0; ivar<nvarsToInclude; ivar++) {
            float iweight = weight;                                               //nominal
            if(varNames[ivar]=="_puup")        iweight *=TotalWeight_plus;        //pu up
            if(varNames[ivar]=="_pudown")      iweight *=TotalWeight_minus;       //pu down
            

            //recompute MET/MT if JES/JER was varied
            LorentzVector zvv = zvvs[ivar>8 ? 0 : ivar];
            PhysicsObjectJetCollection &varJets = ( ivar<=4 ? variedAJets[ivar] : ( useJERsmearing ? variedAJets[0] : recoJets ) );
            PhysicsObjectJetCollection tightVarJets;
            LorentzVector clusteredMetP4(zll);
            clusteredMetP4 *= -1;
            bool passLocalJetveto(true);
            bool passLocalBveto(true); //bool passLocalBveto(passBveto);
            bool passLocalDphijmet(true);
            bool passLocalBalanceCut(true);
            float localmindphijmet(999999.),localmindphijmet15(999999.);
            int localNAJetsGood30(0);
            for(size_t ijet=0; ijet<varJets.size(); ijet++) {
                if(varJets[ijet].pt()<15) continue;

                // dphi
                float idphijmet( fabs(deltaPhi(varJets[ijet].phi(),zvv.phi()) ) );
                if(varJets[ijet].pt()>15) if(idphijmet<localmindphijmet15) localmindphijmet15 = idphijmet;
                if(varJets[ijet].pt()>30) if(idphijmet<localmindphijmet)   localmindphijmet   = idphijmet;

                bool isGoodJet = hasObjectId(aJets[ijet].pid,JETID_LOOSE);
                if(usePUsubJetId) isGoodJet = hasObjectId(aJets[ijet].pid,JETID_CUTBASED_LOOSE);
                if(!isGoodJet) continue;

                clusteredMetP4 -= varJets[ijet];

                if(useJetsOnlyInTracker && fabs(varJets[ijet].eta())>2.5) continue;

                tightVarJets.push_back( varJets[ijet] );
                if(varJets[ijet].pt()>30)localNAJetsGood30++;

                if(varJets[ijet].pt()>20 /*&& fabs(varJets[ijet].eta())<2.5*/) {
                    if(ivar==11)      passLocalBveto &= (varJets[ijet].btag2<0.250);
                    else if(ivar==12) passLocalBveto &= (varJets[ijet].btag2<0.240);
                    else              passLocalBveto &= (varJets[ijet].btag2<0.244);
                }
            }
            //passLocalJetveto=(localNAJetsGood30==0);
            //passLocalDphijmet=(localmindphijmet>0.5); //RJ
            //if(localNAJetsGood30==0) passLocalDphijmet=(localmindphijmet15>0.5);
            //passLocalBalanceCut=(zvv.pt()/zll.pt()>0.8 && zvv.pt()/zll.pt()<1.2);

            //double dphil1met=fabs(deltaPhi(lep1.phi(),zvv.phi()));
            //double dphil2met=fabs(deltaPhi(lep2.phi(),zvv.phi()));
            bool passLocalLMetVeto(true);
            //if(!use2011Id && zvv.pt()>60 && min(dphil1met,dphil2met)<0.2) passLMetVeto=false;

            float mt = METUtils::transverseMass(zll,zvv,true);
            LorentzVector nullP4(0,0,0,0);
            LorentzVector redMet = METUtils::redMET(METUtils::INDEPENDENTLYMINIMIZED, zll, 0, nullP4, 0, clusteredMetP4, zvv, true);

            float new_mt = METUtils::transverseMass(zll,zvv,false); //RJ


            // with standard Z mass, Z pt, RedMet (will be variated later)
            bool passPreselection(passZmass && passZpt && pass3dLeptonVeto && passLocalJetveto && passLocalBveto && passLocalDphijmet && passLocalLMetVeto && passLocalBalanceCut && passRedMet);
            //bool passPreselectionMbvetoMzmass(passZpt && pass3dLeptonVeto && passLocalJetveto && passLocalDphijmet && passLocalLMetVeto && passLocalBalanceCut && passRedMet);

            //RJ compute dphi(zpt,redMet) for statistics
            double LocaldphiZllmet=fabs(deltaPhi(zll.phi(),zvv.phi()));
            //bool passLocaldphiZllmetCut(LocaldphiZllmet>2.6);

            //RJ
            bool passLocalMTcut(mt>220 && mt<1200);
            //passLocalMTcut &= passdphiMetmvaMet;
            //passLocalMTcut &= passdphi2l;

            //re-assign the event category if jets were varied
            int eventSubCat  = eventCategoryInst.Get(phys,&tightVarJets);
            TString tag_subcat = eventCategoryInst.GetLabel(eventSubCat);
            tags.clear();

            if(tag_subcat != "geq2jets") {
                //tags.push_back(tag_cat);
                tags.push_back(tag_cat+tag_subcat);
                if(tag_cat=="mumu" || tag_cat=="ee") tags.push_back(string("ll")+tag_subcat);
            }

            //fill shapes
            for(unsigned int index=0; index<optim_Cuts1_met.size(); index++) {

                float minMet=optim_Cuts1_met[index];
                float minBalance=optim_Cuts1_balance[index];
                float minDphi=optim_Cuts1_dphi[index];
                float deltaZ=optim_Cuts1_zmass[index];

                bool passLocalRedMet(zvv.pt()>minMet); //Reduced MET: redMet.pt() ; PFMET: zvv.pt()
                bool passLocalRedMet_RJ(redMet.pt()>65);  //Reduced MET: redMet.pt() ; PFMET: zvv.pt()  for passPreselectionMjvetoMbvetoMzmass
                bool passLocalZmass(fabs(zll.mass()-91)<deltaZ);
                bool passLocalZpt(zll.pt()>50.); //fix Zpt cut
                passLocalBalanceCut=(zvv.pt()/zll.pt()>(1.-minBalance) && zvv.pt()/zll.pt()<(1.+minBalance));
                bool passLocalBalanceCut_RJ=(zvv.pt()/zll.pt()>0.4 && zvv.pt()/zll.pt()<1.8);
                bool passLocaldphiZllmetCut(LocaldphiZllmet>minDphi);
                passPreselection = (passLocalMTcut && passLocaldphiZllmetCut && passLocalZmass && passLocalZpt && pass3dLeptonVeto && passLocalJetveto && passLocalBveto && passLocalDphijmet && passLocalBalanceCut && passLocalRedMet);
                //passPreselectionMbvetoMzmass = (passLocalMTcut && passLocaldphiZllmetCut && passLocalZpt && pass3dLeptonVeto && passLocalJetveto && passLocalDphijmet && passLocalBalanceCut && passLocalRedMet); //not use
                bool passPreselectionMjvetoMbvetoMzmass = ( passLocalRedMet_RJ && passLocalZpt && pass3dLeptonVeto && passLocalDphijmet && passLocalBalanceCut_RJ );
                //bool passPreselectionMjvetoMbvetoMzmass = (passLocalRedMet && passLocalZpt && pass3dLeptonVeto && passLocalDphijmet && passLocalBalanceCut && passLocaldphiZllmetCut);

                if( passPreselection ) {
                    mon.fillHisto(TString("redMet_shapes")+varNames[ivar],tags,index, redMet.pt(),iweight);
                    mon.fillHisto(TString("redMet_rebin_shapes")+varNames[ivar],tags,index, redMet.pt(),iweight);
                    mon.fillHisto(TString("mt_shapes")+varNames[ivar],tags,index,new_mt,iweight);
                    mon.fillHisto(TString("new_mt_shapes")+varNames[ivar],tags,index,new_mt,iweight);
                    /* balance_mt_shapes */
                    float i_balance = zvv.pt()/zll.pt();
                    if(i_balance>=BalanceXaxis[0] && i_balance< BalanceXaxis[1]) mon.fillHisto(TString("balance_mt_shapes")+varNames[ivar],tags,index,new_mt,iweight);
                    if(i_balance>=BalanceXaxis[1] && i_balance< BalanceXaxis[2]) mon.fillHisto(TString("balance_mt_shapes")+varNames[ivar],tags,index,new_mt+1200.,iweight);
                    if(i_balance>=BalanceXaxis[2] && i_balance<=BalanceXaxis[3]) mon.fillHisto(TString("balance_mt_shapes")+varNames[ivar],tags,index,new_mt+2400.,iweight);
                    /* dphiLL_mt_shapes */
                    if(dphi2l>=DphiLLXaxis[0] && dphi2l< DphiLLXaxis[1]) mon.fillHisto(TString("dphiLL_mt_shapes")+varNames[ivar],tags,index,new_mt,iweight);
                    if(dphi2l>=DphiLLXaxis[1] && dphi2l< DphiLLXaxis[2]) mon.fillHisto(TString("dphiLL_mt_shapes")+varNames[ivar],tags,index,new_mt+1200.,iweight);
                    if(dphi2l>=DphiLLXaxis[2] && dphi2l<=DphiLLXaxis[3]) mon.fillHisto(TString("dphiLL_mt_shapes")+varNames[ivar],tags,index,new_mt+2400.,iweight);
                    /* coslZ_mt_shapes */
                    if(colin_soper>=CosThetaXaxis[0] && colin_soper< CosThetaXaxis[1]) mon.fillHisto(TString("coslZ_mt_shapes")+varNames[ivar],tags,index,new_mt,iweight);
                    if(colin_soper>=CosThetaXaxis[1] && colin_soper< CosThetaXaxis[2]) mon.fillHisto(TString("coslZ_mt_shapes")+varNames[ivar],tags,index,new_mt+1200.,iweight);
                    if(colin_soper>=CosThetaXaxis[2] && colin_soper<=CosThetaXaxis[3]) mon.fillHisto(TString("coslZ_mt_shapes")+varNames[ivar],tags,index,new_mt+2400.,iweight);
                    mon.fillHisto(TString("dphi_shapes")+varNames[ivar],tags,index,LocaldphiZllmet,iweight);
                    mon.fillHisto(TString("met_shapes")+varNames[ivar],tags,index,zvv.pt(),iweight);
                    mon.fillHisto(TString("met_rebin_shapes")+varNames[ivar],tags,index,zvv.pt(),iweight);
                    mon.fillHisto(TString("zpt_shapes")+varNames[ivar],tags,index, zll.pt(),iweight);
                    mon.fillHisto(TString("zpt_rebin_shapes")+varNames[ivar],tags,index, zll.pt(),iweight);
                }
                
                
                
                if( passPreselectionMjvetoMbvetoMzmass && passLocalZmass && passLocalJetveto && passLocalBveto ) {
                    for(size_t j=0; j<allshapesVars.size(); j++){
                        mon.fillHisto(allshapesVars[j]+"_NRBctrl"+varNames[ivar],tags,index,0,iweight);
                    }
                }
                if( passPreselectionMjvetoMbvetoMzmass && isZsideBand && passLocalJetveto && passLocalBveto ) {
                    for(size_t j=0; j<allshapesVars.size(); j++){
                        mon.fillHisto(allshapesVars[j]+"_NRBctrl"+varNames[ivar],tags,index,1,iweight);
                    }
                }
                if( passPreselectionMjvetoMbvetoMzmass && isZsideBandPlus && passLocalJetveto && passLocalBveto ) {
                    for(size_t j=0; j<allshapesVars.size(); j++){
                        mon.fillHisto(allshapesVars[j]+"_NRBctrl"+varNames[ivar],tags,index,2,iweight);
                    }
                }

                //RJ, changing eventCategory for NRB control sample, i.e., no jetveto
                //prepare the tag's vectors for histo filling
                std::vector<TString> NRBtags(1,"all");
                NRBtags.clear();
                //NRBtags.push_back(tag_cat);
                NRBtags.push_back(tag_cat+"eq0jets");
                NRBtags.push_back(tag_cat+"eq1jets");
                if(tag_cat=="mumu" || tag_cat=="ee") {
                    NRBtags.push_back(string("ll")+"eq0jets");
                    NRBtags.push_back(string("ll")+"eq1jets");
                }

                if( passPreselectionMjvetoMbvetoMzmass && passLocalZmass && !passLocalBveto ) {
                    for(size_t j=0; j<allshapesVars.size(); j++){
                        mon.fillHisto(allshapesVars[j]+"_NRBctrl"+varNames[ivar],NRBtags,index,3,iweight);
                    }
                }
                if( passPreselectionMjvetoMbvetoMzmass && isZsideBand && !passLocalBveto ) {
                    for(size_t j=0; j<allshapesVars.size(); j++){
                        mon.fillHisto(allshapesVars[j]+"_NRBctrl"+varNames[ivar],NRBtags,index,4,iweight);
                    }
                }
                if( passPreselectionMjvetoMbvetoMzmass && isZsideBandPlus && !passLocalBveto ) {
                    for(size_t j=0; j<allshapesVars.size(); j++){
                        mon.fillHisto(allshapesVars[j]+"_NRBctrl"+varNames[ivar],NRBtags,index,5,iweight);
                    }
                }

            }//all shape variables END
        }//Systematic variation END 
    }//Event loop END

    printf("\n");
    file->Close();

    //##############################################
    //########     SAVING HISTO TO FILE     ########
    //##############################################
    //save control plots to file
    outUrl += "/";
    outUrl += outFileUrl + ".root";
    printf("Results save in %s\n", outUrl.Data());

    //save all to the file
    TFile *ofile=TFile::Open(outUrl, "recreate");
    mon.Write();
    ofile->Close();

    if(outTxtFile_full)fclose(outTxtFile_full);
    if(outTxtFile_final)fclose(outTxtFile_final);
}





