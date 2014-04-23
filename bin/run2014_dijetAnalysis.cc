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

std::vector<LorentzVector> checkDiMuonMass(std::vector<LorentzVector> leps);
std::vector<LorentzVector> checkDiElectronMass(std::vector<LorentzVector> leps);

int main(int argc, char* argv[])
{
    //##################################################################################
    //##########################    GLOBAL INITIALIZATION     ##########################
    //##################################################################################

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
/*
    bool isSingleMuPD(!isMC && url.Contains("SingleMu"));
    bool isSingleElePD(!isMC && url.Contains("SingleEle"));
    bool isMC_ZZ  = isMC && ( string(url.Data()).find("TeV_ZZ")  != string::npos);
    bool isMC_WZ  = isMC && ( string(url.Data()).find("TeV_WZ")  != string::npos);
    bool isMC_ZH  = isMC && ( string(url.Data()).find("TeV_ZH")  != string::npos);
    bool isMC_HZZd= isMC && ( string(url.Data()).find("TeV_HZZd")  != string::npos);
*/
    bool isV0JetsMC(isMC && (url.Contains("DYJetsToLL_50toInf") || url.Contains("WJets")));
/*
    bool isMC_DY  = isMC && ( (string(url.Data()).find("TeV_DYJetsToLL")!= string::npos)
			|| (string(url.Data()).find("TeV_ZG")!= string::npos) );
*/
    // print out event information
/*
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
*/

    //tree info
    int evStart     = runProcess.getParameter<int>("evStart");
    int evEnd       = runProcess.getParameter<int>("evEnd");
    TString dirname = runProcess.getParameter<std::string>("dirName");

    //jet energy scale uncertainties
    TString uncFile = runProcess.getParameter<std::string>("jesUncFileName");
    gSystem->ExpandPathName(uncFile);
    JetCorrectionUncertainty jecUnc(uncFile.Data());

    //ZHinvisible reweighting file input
    //ZHUtils myZHUtils(runProcess);

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
    //size_t nvarsToInclude=varNames.size();




    //##################################################################################
    //##########################    INITIATING HISTOGRAMS     ##########################
    //##################################################################################

    SmartSelectionMonitor mon;


    TH1F *h=(TH1F*) mon.addHistogram( new TH1F ("eventflow", ";;Events", 12,0,12) );
    h->GetXaxis()->SetBinLabel(1,"Trigger");
    h->GetXaxis()->SetBinLabel(2,"#geq 2 leptons");
    h->GetXaxis()->SetBinLabel(3,"#geq 2 iso leptons");

    //for MC normalization (to 1/pb)
    TH1F* Hcutflow  = (TH1F*) mon.addHistogram(  new TH1F ("cutflow"    , "cutflow"    ,6,0,6) ) ;

    //mon.addHistogram( new TH1F( "zpt_raw",       ";#it{p}_{T}^{ll} [GeV];Events", 50,0,500) );
    //mon.addHistogram( new TH1F( "zmass_raw", ";#it{m}_{ll} [GeV];Events", 100,40,250) );
    //mon.addHistogram( new TH1F( "elewmt_raw",          ";#it{m}_{T} [GeV];Events", 50,0,20) );
    //mon.addHistogram( new TH1F( "redmet_raw",        ";Red-E_{T}^{miss} [GeV];Events", 50,0,500) );
    //mon.addHistogram( new TH1F( "met_sig_raw",        ";E_{T}^{miss} Significance;Events", 50,0,100) );
    //mon.addHistogram( new TH1F( "mtless_ZMet_raw",      ";#it{m}_{T}(Z, E_{T}^{miss}) [GeV];Events", 50,0,500) );


    // plots for control plot 
    // all plots
    mon.addHistogram( new TH1F( "pfmet_raw",        ";E_{T}^{miss} [GeV];Events", 50,0,100) );

    // electron channel
    mon.addHistogram( new TH1F( "elepfmet_raw",        ";E_{T}^{miss} [GeV];Events", 50,0,100) );
    mon.addHistogram( new TH1F( "elewmt_raw",          ";#it{m}_{T} [GeV];Events", 50,0,100) );
    mon.addHistogram( new TH1F( "elept_raw",          ";#it{p}_{T} [GeV];Events", 50,0,100) );
    mon.addHistogram( new TH1F( "elewmt_metgep50",          ";#it{m}_{T} [GeV];Events", 50,0,100) );
    mon.addHistogram( new TH1F( "eledphilepj_raw",     ";#Delta#phi(l,j) [rad];Events", 50,0,TMath::Pi()) );
	
    // muon channel 
    mon.addHistogram( new TH1F( "mupfmet_raw",        ";E_{T}^{miss} [GeV];Events", 50,0,100) );
    mon.addHistogram( new TH1F( "muwmt_raw",          ";#it{m}_{T} [GeV];Events", 50,0,100) );
    mon.addHistogram( new TH1F( "mupt_raw",          ";#it{p}_{T} [GeV];Events", 50,0,100) );
    mon.addHistogram( new TH1F( "muwmt_metgep50",          ";#it{m}_{T} [GeV];Events", 50,0,100) );
    mon.addHistogram( new TH1F( "mudphilepj_raw",     ";#Delta#phi(l,j) [rad];Events", 50,0,TMath::Pi()) );
   



    //plots for driving fake rates
    double fakePt[11]= {20,25,30,35,40,45,50,60,70,80,100};
    double mufakePt[9] = {20,25,30,35,40,45,50,60,100};
    double fakeEta[9]= {-2.5,-2.,-1.479,-1.,0.,1.,1.479,2.,2.5};
    
    mon.addHistogram( new TH1F( "eleLooseFakePt",       ";Loose #it{p}_{T}^{e} [GeV];Events", 10,fakePt) );
    mon.addHistogram( new TH1F( "eleTightFakePt",       ";Tight #it{p}_{T}^{e} [GeV];Events", 10,fakePt) );
    mon.addHistogram( new TH1F( "eleLooseFakeEta",       ";Loose #it{p}_{T}^{e} [GeV];Events", 8,fakeEta) );
    mon.addHistogram( new TH1F( "eleTightFakeEta",       ";Tight #it{p}_{T}^{e} [GeV];Events", 8,fakeEta) );

    mon.addHistogram( new TH1F( "eleLooseFake_etabin1",       ";Loose #it{p}_{T}^{e} [GeV];Events", 10,fakePt) );
    mon.addHistogram( new TH1F( "eleTightFake_etabin1",       ";Tight #it{p}_{T}^{e} [GeV];Events", 10,fakePt) );
    mon.addHistogram( new TH1F( "eleLooseFake_etabin2",       ";Loose #it{p}_{T}^{e} [GeV];Events", 10,fakePt) );
    mon.addHistogram( new TH1F( "eleTightFake_etabin2",       ";Tight #it{p}_{T}^{e} [GeV];Events", 10,fakePt) );
    mon.addHistogram( new TH1F( "eleLooseFake_etabin3",       ";Loose #it{p}_{T}^{e} [GeV];Events", 10,fakePt) );
    mon.addHistogram( new TH1F( "eleTightFake_etabin3",       ";Tight #it{p}_{T}^{e} [GeV];Events", 10,fakePt) );
    mon.addHistogram( new TH1F( "eleLooseFake_etabin4",       ";Loose #it{p}_{T}^{e} [GeV];Events", 10,fakePt) );
    mon.addHistogram( new TH1F( "eleTightFake_etabin4",       ";Tight #it{p}_{T}^{e} [GeV];Events", 10,fakePt) );


    mon.addHistogram( new TH1F( "muLooseFakePt",        ";Loose #it{p}_{T}^{#mu} [GeV];Events", 8,mufakePt) );
    mon.addHistogram( new TH1F( "muTightFakePt",        ";Tight #it{p}_{T}^{#mu} [GeV];Events", 8,mufakePt) );
    mon.addHistogram( new TH1F( "muLooseFakeEta",        ";Loose #it{p}_{T}^{#mu} [GeV];Events", 8,fakeEta) );
    mon.addHistogram( new TH1F( "muTightFakeEta",        ";Tight #it{p}_{T}^{#mu} [GeV];Events", 8,fakeEta) );

    mon.addHistogram( new TH1F( "muLooseFake_etabin1",       ";Loose #it{p}_{T}^{#mu} [GeV];Events", 8,mufakePt) );
    mon.addHistogram( new TH1F( "muTightFake_etabin1",       ";Tight #it{p}_{T}^{#mu} [GeV];Events", 8,mufakePt) );
    mon.addHistogram( new TH1F( "muLooseFake_etabin2",       ";Loose #it{p}_{T}^{#mu} [GeV];Events", 8,mufakePt) );
    mon.addHistogram( new TH1F( "muTightFake_etabin2",       ";Tight #it{p}_{T}^{#mu} [GeV];Events", 8,mufakePt) );
    mon.addHistogram( new TH1F( "muLooseFake_etabin3",       ";Loose #it{p}_{T}^{#mu} [GeV];Events", 8,mufakePt) );
    mon.addHistogram( new TH1F( "muTightFake_etabin3",       ";Tight #it{p}_{T}^{#mu} [GeV];Events", 8,mufakePt) );
    mon.addHistogram( new TH1F( "muLooseFake_etabin4",       ";Loose #it{p}_{T}^{#mu} [GeV];Events", 8,mufakePt) );
    mon.addHistogram( new TH1F( "muTightFake_etabin4",       ";Tight #it{p}_{T}^{#mu} [GeV];Events", 8,mufakePt) );







    //##################################################################################
    //#############         GET READY FOR THE EVENT LOOP           #####################
    //##################################################################################

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
        TH1F* cutflowH = (TH1F *) file->Get("dataAnalyzer/dijet/cutflow");
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
    if(isMC) {
        TString puDist("dataAnalyzer/dijet/pileuptrue");
        if(useObservedPU) puDist="dataAnalyzer/dijet/pileup";
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
    //EventCategory eventCategoryInst(1);   //jet(0,1,>=2) binning


    LeptonEfficiencySF lsf(use2011Id ? 2011:2012);

    //####################################################################################################################
    //###########################################           EVENT LOOP         ###########################################
    //####################################################################################################################

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

        if(isV0JetsMC){
            if(ev.mc_nup>5) continue;
        }

        PhysicsEvent_t phys=getPhysicsEventFrom(ev);

        //event category
        //bool isSameFlavor(ev.cat==MUMU || ev.cat==EE);
        TString tag_cat;
/*
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
*/
	//split inclusive DY sample into DYToLL and DYToTauTau
        //if(isMC && mctruthmode==1 && !isDYToLL(ev.mccat) ) continue;
        //if(isMC && mctruthmode==2 && !isDYToTauTau(ev.mccat) ) continue;

        //require compatibilitiy of the event with the PD
        bool hasTrigger(false);
        bool hasMtrigger = ev.triggerType & 0x1;
        bool hasEtrigger = (ev.triggerType >> 1 ) & 0x1;
        if(!isMC) {
            if( !(hasMtrigger||hasEtrigger) ) continue;
            hasTrigger=true;
        } else {
            if( (hasMtrigger || hasEtrigger) ) hasTrigger=true;
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






	// ewk correction for diboson processes
	//if(isMC_ZH) weight *= myZHUtils.GetNLOZHWeight(phys);
	//if(isMC_ZZ) weight *= myZHUtils.GetNLOZZWeight(phys);
        //if(isMC_WZ) weight *= myZHUtils.GetNLOWZWeight(phys);


        //#########################################################################
        //#####################      Objects Selection       ######################
        //#########################################################################

        //
        //MET variables
        //
        LorentzVector rawMetP4=phys.met[1];
        if(use2011Id) rawMetP4=phys.met[0];
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

        std::vector<LorentzVector> looseMuons_raw;
        std::vector<LorentzVector> tightMuons_raw;
        std::vector<LorentzVector> looseElectrons_raw;
        std::vector<LorentzVector> tightElectrons_raw;


/*
        LorentzVector lep1=phys.leptons[0];
        LorentzVector lep2=phys.leptons[1];
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
*/
        //check alternative selections for the dilepton
        double llScaleFactor(1.0),llTriggerEfficiency(1.0);
        //LorentzVector genZP4(0,0,0,0); // for checks on Sherpa ZZ
        //int genmatchid[2] = {-1, -1};
        //double genmatchdr[2] = {0.1, 0.1};



        // loop all leptons begin
        for(size_t ilep=0; ilep<phys.leptons.size(); ilep++) {
            TString lepStr( fabs(phys.leptons[ilep].id)==13 ? "mu" : "e");

            //id and isolation
            int lpid=phys.leptons[ilep].pid;
            //float relIso2011 = phys.leptons[ilep].relIsoRho(ev.rho);
            float relIso = (lepStr=="mu") ?
                           phys.leptons[ilep].pfRelIsoDbeta() :
                           phys.leptons[ilep].ePFRelIsoCorrected2012(ev.rho,ev.en_sceta[lpid]);

            std::vector<int> passIds;
            std::map<int,bool> passIsos;

            //bool hasGoodId(false), isIso(false);

            if(fabs(phys.leptons[ilep].id)==13) {
                if( hasObjectId(ev.mn_idbits[lpid], MID_PF) 
		    && hasObjectId(ev.mn_idbits[lpid], MID_GLOBAL)
		    && phys.leptons[ilep].trkchi2 < 10 
		    && ev.mn_validMuonHits[lpid] > 0 
		    && ev.mn_nMatchedStations[lpid] > 1 
                    && fabs(phys.leptons[ilep].d0)< 2.0/*0.2*/
                    && fabs(phys.leptons[ilep].dZ)<0.5
                    && phys.leptons[ilep].trkValidPixelHits>0
                    && ev.mn_trkLayersWithMeasurement[lpid]>5
		    && relIso< 1.0 /*0.2*/ )
		{
		    looseMuons_raw.push_back(phys.leptons[ilep]);
                }
                if( hasObjectId(ev.mn_idbits[lpid], MID_TIGHT) && relIso<0.2)    {
		    tightMuons_raw.push_back(phys.leptons[ilep]);
                }

                llScaleFactor *= lsf.getLeptonEfficiency(phys.leptons[ilep].pt(),fabs(phys.leptons[ilep].eta()),13).first;

            } else {
                int wps[]= {EgammaCutBasedEleId::LOOSE,EgammaCutBasedEleId::MEDIUM, EgammaCutBasedEleId::VETO};
                llScaleFactor *= lsf.getLeptonEfficiency(phys.leptons[ilep].pt(),fabs(phys.leptons[ilep].eta()),11).first;
                for(int iwp=0; iwp<3; iwp++) {

                        bool passWp = EgammaCutBasedEleId::PassWP(EgammaCutBasedEleId::WorkingPoint(wps[iwp]),
                                      (fabs(phys.leptons[ilep].eta())<1.4442),
                                      phys.leptons[ilep].pt(), phys.leptons[ilep].eta(),
                                      ev.en_detain[lpid],  ev.en_dphiin[lpid], ev.en_sihih[lpid], ev.en_hoe[lpid],
                                      ev.en_ooemoop[lpid], phys.leptons[ilep].d0, phys.leptons[ilep].dZ,
                                      0., 0., 0.,
                                      !hasObjectId(ev.en_idbits[lpid], EID_CONVERSIONVETO),0,ev.rho);
                        if(passWp && iwp==2) {
			    looseElectrons_raw.push_back(phys.leptons[ilep]);
                        }
			if(passWp && iwp==1 && relIso<0.15) {
			    tightElectrons_raw.push_back(phys.leptons[ilep]);
			}
                        if(!use2011Id) {
                            //llScaleFactor *= 1;
                            llTriggerEfficiency *= 1.0; //electronTriggerEfficiency(phys.leptons[ilep].pt(),fabs(phys.leptons[ilep].eta()),2012);
                        }
                }
            }

        } // loop all leptons end




        //
        //STD PF JET ANALYSIS
        //
        //bool passJetveto(true);
        //bool passBveto(true);
        //bool passMetCut(true);
        //bool passDphijmet(true);
        //bool passBalanceCut(true);
        PhysicsObjectJetCollection &aJets = ( useJERsmearing ? variedAJets[0] : recoJets );
        PhysicsObjectJetCollection aGoodIdJets;
        //LorentzVector aClusteredMetP4(zll);
        //aClusteredMetP4 *= -1;
        //LorentzVector Recoil(zvvs[0]);
        //Recoil *= -1;
        //Recoil -= zll;
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


            //aClusteredMetP4 -= aJets[ijet];

            if(useJetsOnlyInTracker && fabs(aJets[ijet].eta())>2.5) continue;

            aGoodIdJets.push_back(aJets[ijet]);
            if(aJets[ijet].pt()>30)nAJetsGood30++;

            if(aJets[ijet].pt()>20 /*&& fabs(aJets[ijet].eta())<2.5*/)  nABtags += (aJets[ijet].btag6>0.244);
            if(aJets[ijet].pt()>20 /*&& fabs(aJets[ijet].eta())<2.5*/)  nCSVMtags += (aJets[ijet].btag6>0.679);
            if(aJets[ijet].pt()>20 /*&& fabs(aJets[ijet].eta())<2.5*/)  nCSVTtags += (aJets[ijet].btag6>0.898);
        }
        //passJetveto=(nAJetsGood30==0);
        //passBveto=(nABtags==0);
        //bool passMBveto=(nCSVMtags==0);
        //bool passTBveto=(nCSVTtags==0);
        //passDphijmet=(mindphijmet>0.5);
        //if(nAJetsGood30==0) passDphijmet=(mindphijmet15>0.5);
        //passBalanceCut=(zvvs[0].pt()/zll.pt()>0.8 && zvvs[0].pt()/zll.pt()<1.2);
        //passBalanceCut=(zvvs[0].pt()/zll.pt()>0.75 && zvvs[0].pt()/zll.pt()<1.25);
	//bool passBalanceCut05=(zvvs[0].pt()/zll.pt()>0.5 && zvvs[0].pt()/zll.pt()<1.5);
	//bool passBalanceCut075=(zvvs[0].pt()/zll.pt()>0.25 && zvvs[0].pt()/zll.pt()<1.75);
        //bool passBalanceCutWWCtrl=(zvvs[0].pt()/zll.pt()>0.4 && zvvs[0].pt()/zll.pt()<1.8);


        //ad-hoc cut for obvious correlations between MET and a lepton
        //double dphil1met=fabs(deltaPhi(lep1.phi(),zvvs[0].phi()));
        //double dphil2met=fabs(deltaPhi(lep2.phi(),zvvs[0].phi()));
        //bool passLMetVeto(true);
        //if(!use2011Id && zvvs[0].pt()>60 && min(dphil1met,dphil2met)<0.2) passLMetVeto=false;

        //other mets
        //METUtils::stRedMET aRedMetOut;
        LorentzVector null(0,0,0,0);
        //LorentzVector aRedMet=METUtils::redMET(METUtils::INDEPENDENTLYMINIMIZED, zll, 0, null, 0, aClusteredMetP4, zvvs[0], true, &aRedMetOut);
        //double aRedMetL=aRedMetOut.redMET_l;
        //double aRedMetT=aRedMetOut.redMET_t;
        //TVector2 aClusteredMet2(aClusteredMetP4.px(),aClusteredMetP4.py());
        //double clusteredMetL=aRedMetOut.a_l*aClusteredMet2;
        //double clusteredMetT=aRedMetOut.a_t*aClusteredMet2;
        //passMetCut=(zvvs[0].pt()>120); //ReducedMET: aRedMet.pt(); PFMET: zvvs[0].pt()
	/*
	bool passMetCut60=(zvvs[0].pt()>60);
	bool passMetCut80=(zvvs[0].pt()>80);
	bool passMetCut100=(zvvs[0].pt()>100);
	bool passMetCut120=(zvvs[0].pt()>120);
	*/
        //bool passMetCutWWCtrl=(aRedMet.pt()>65);

        //TVector2 RecoilMet2(Recoil.px(),Recoil.py());
        //double RecoilMetL=aRedMetOut.a_l*RecoilMet2;
        //double RecoilMetT=aRedMetOut.a_t*RecoilMet2;


        // compute dphi(zpt,redMet)
        //double dphiZllmet=fabs(deltaPhi(zll.phi(),zvvs[0].phi()));
        //double dphiMet_mvaMet=fabs(deltaPhi(mvaMetP4.phi(),zvvs[0].phi()));
        //double dphiZllredMet=fabs(deltaPhi(zll.phi(),aRedMet.phi()));
        //bool passdphiZllmetCut(dphiZllmet>2.7);///2.6);
	//bool passdphiZllmetCut20(dphiZllmet>2.0);
	//bool passdphiZllmetCut24(dphiZllmet>2.4);
        //bool passdphiMetmvaMet(dphiMet_mvaMet<0.2);

/*
        //transverse masses
        double aMT=METUtils::transverseMass(zll,zvvs[0],true);
        double aMTmassless=METUtils::transverseMass(zll,zvvs[0],false);
	double mtless_MetLLp = METUtils::transverseMass(lep1,zvvs[0],false);
	double mtless_MetTLp = METUtils::transverseMass(lep2,zvvs[0],false);
        double balanceDif=fabs(zvvs[0].pt()-zll.pt())/zll.pt();
        TVector2 dil2(zll.px(),zll.py());
        TVector2 met2(zvvs[0].px(),zvvs[0].py());
        double axialMet=dil2*met2;
        axialMet /= -zll.pt();
        //double pfMetCompL = aRedMetOut.a_l*met2;
        //double pfMetCompT = aRedMetOut.a_t*met2;


        bool passMTcut(aMT>220 && aMT<1200);

*/



        //#########################################################
        //####  RUN PRESELECTION AND CONTROL REGION PLOTS  ########
        //#########################################################

        if(isMC && use2011Id) weight *= llScaleFactor*llTriggerEfficiency;

        if(hasTrigger)                 {
            mon.fillHisto("eventflow",tags,0,weight);
        }
	else continue;



	//mon.fillHisto("zmass_raw"			,tags, zll.mass(), weight);
        //mon.fillHisto("zpt_raw"				,tags, zll.pt(),   weight);
	mon.fillHisto("pfmet_raw"			,tags, zvvs[0].pt(), weight);
	//mon.fillHisto("redmet_raw"			,tags, aRedMet.pt(), weight);
	//mon.fillHisto("met_sig_raw"			,tags, phys.met_sig[1], weight);
	//mon.fillHisto("mtless_ZMet_raw"                 ,tags, aMTmassless, weight);
	//mon.fillHisto("nleptons_raw"			,tags, nextraleptons, weight);
	////mon.fillHisto("npfjets_raw"			,tags, nAJetsGood30, weight);
	////mon.fillHisto("npfbjets_raw"			,tags, nABtags, weight);

        //##############################################
        //########  Main Event Selection        ########
        //##############################################


	// investigating dijet ctrl sammple
	// control plots
	LorentzVector LeadingJet;
	double largest_pt=-999;
	for(size_t j=0; j<aGoodIdJets.size(); j++){
		double pt = aGoodIdJets[j].pt();
		if(pt>largest_pt) {
			largest_pt = pt;
			LeadingJet = aGoodIdJets[j];
		}
	}



	if( hasEtrigger && ( isMC || (!isMC && fType==EE) ) && looseElectrons_raw.size()>0 && aGoodIdJets.size()>0)
	{
		mon.fillHisto("elepfmet_raw"                       ,tags, zvvs[0].pt(), weight);
		for(size_t j=0; j<looseElectrons_raw.size(); j++){
			double Wmt = METUtils::transverseMass(looseElectrons_raw[j],zvvs[0],false);
			double dphi=fabs(deltaPhi(looseElectrons_raw[j].phi(),LeadingJet.phi()));
			mon.fillHisto("eledphilepj_raw",tags, dphi, weight);	
			mon.fillHisto("elewmt_raw",tags, Wmt, weight);
			mon.fillHisto("elept_raw",tags, looseElectrons_raw[j].pt(), weight);
			if(zvvs[0].pt()>50) mon.fillHisto("elewmt_metgep50",tags, Wmt, weight);
		}
	}

	if( hasMtrigger && (isMC || (!isMC && fType==MUMU) ) && looseMuons_raw.size()>0 && aGoodIdJets.size()>0) 
	{
		mon.fillHisto("mupfmet_raw"                       ,tags, zvvs[0].pt(), weight);
		for(size_t j=0; j<looseMuons_raw.size(); j++){
			double Wmt = METUtils::transverseMass(looseMuons_raw[j],zvvs[0],false);
			double dphi=fabs(deltaPhi(looseMuons_raw[j].phi(),LeadingJet.phi()));
			mon.fillHisto("mudphilepj_raw",tags, dphi, weight);
			mon.fillHisto("muwmt_raw",tags, Wmt, weight);
			mon.fillHisto("mupt_raw",tags, looseMuons_raw[j].pt(), weight);
			if(zvvs[0].pt()>50) mon.fillHisto("muwmt_metgep50",tags, Wmt, weight);
		}
	} 




	// drive fake rate
	// fixed dijet control sample
	if(zvvs[0].pt()<20)
	{
		std::vector<LorentzVector> looseMuons = checkDiMuonMass(looseMuons_raw);
		std::vector<LorentzVector> tightMuons = checkDiMuonMass(tightMuons_raw); 
		std::vector<LorentzVector> looseElectrons = checkDiElectronMass(looseElectrons_raw);
		std::vector<LorentzVector> tightElectrons = checkDiElectronMass(tightElectrons_raw); 

		if( isMC || (!isMC && fType==EE) ){
			for(size_t j=0; j<looseElectrons.size(); j++){
				if(aGoodIdJets.size()>0){
					double dphi=fabs(deltaPhi(looseElectrons[j].phi(),LeadingJet.phi()));
					if(dphi<1) continue;
				}
				//double Wmt = METUtils::transverseMass(looseElectrons[j],zvvs[0],false);
				//if (Wmt > 20) continue;

				mon.fillHisto("eleLooseFakePt",tags, looseElectrons[j].pt(), weight);
				double eta = looseElectrons[j].eta();
				mon.fillHisto("eleLooseFakeEta",tags, eta, weight);
				if(fabs(eta)<1.0)        mon.fillHisto("eleLooseFake_etabin1",tags, looseElectrons[j].pt(), weight);
				else if(fabs(eta)<1.479) mon.fillHisto("eleLooseFake_etabin2",tags, looseElectrons[j].pt(), weight); 
				else if(fabs(eta)<2.00)  mon.fillHisto("eleLooseFake_etabin3",tags, looseElectrons[j].pt(), weight);
				else 			 mon.fillHisto("eleLooseFake_etabin4",tags, looseElectrons[j].pt(), weight);
			
        		}
			for(size_t j=0; j<tightElectrons.size(); j++){
				if(aGoodIdJets.size()>0){
					double dphi=fabs(deltaPhi(tightElectrons[j].phi(),LeadingJet.phi()));
					if(dphi<1) continue;
				}
				//double Wmt = METUtils::transverseMass(tightElectrons[j],zvvs[0],false);
				//if (Wmt > 20) continue;
				mon.fillHisto("eleTightFakePt",tags, tightElectrons[j].pt(), weight);
                        	double eta = tightElectrons[j].eta();
				mon.fillHisto("eleTightFakeEta",tags, eta, weight);
                        	if(fabs(eta)<1.0) 	 mon.fillHisto("eleTightFake_etabin1",tags, tightElectrons[j].pt(), weight);
                        	else if(fabs(eta)<1.479) mon.fillHisto("eleTightFake_etabin2",tags, tightElectrons[j].pt(), weight);
                        	else if(fabs(eta)<2.00)  mon.fillHisto("eleTightFake_etabin3",tags, tightElectrons[j].pt(), weight);
                        	else 			 mon.fillHisto("eleTightFake_etabin4",tags, tightElectrons[j].pt(), weight);

        		}
		}

		if( isMC || (!isMC && fType==MUMU) ){
			for(size_t j=0; j<looseMuons.size(); j++){
                        	if(aGoodIdJets.size()>0){
                                	double dphi=fabs(deltaPhi(looseMuons[j].phi(),LeadingJet.phi()));
                                	if(dphi<1) continue;
                        	}
				double Wmt = METUtils::transverseMass(looseMuons[j],zvvs[0],false);
				//mon.fillHisto("wmt_raw",tags, Wmt, weight);
				if (Wmt > 20) continue;
					
				double eta = looseMuons[j].eta();
				mon.fillHisto("muLooseFakePt",tags, looseMuons[j].pt(), weight);
				mon.fillHisto("muLooseFakeEta",tags, eta, weight);
				if(fabs(eta)<1.0)        mon.fillHisto("muLooseFake_etabin1",tags, looseMuons[j].pt(), weight);
				else if(fabs(eta)<1.479) mon.fillHisto("muLooseFake_etabin2",tags, looseMuons[j].pt(), weight); 
				else if(fabs(eta)<2.00)  mon.fillHisto("muLooseFake_etabin3",tags, looseMuons[j].pt(), weight);
				else 			 mon.fillHisto("muLooseFake_etabin4",tags, looseMuons[j].pt(), weight);
	        	}
			for(size_t j=0; j<tightMuons.size(); j++){
                        	if(aGoodIdJets.size()>0){
                                	double dphi=fabs(deltaPhi(tightMuons[j].phi(),LeadingJet.phi()));
                                	if(dphi<1) continue;
                        	}
				double Wmt = METUtils::transverseMass(tightMuons[j],zvvs[0],false);
				if (Wmt > 20) continue;
			
				double eta = tightMuons[j].eta();
				mon.fillHisto("muTightFakePt",tags, tightMuons[j].pt(), weight);
				mon.fillHisto("muTightFakeEta",tags, eta, weight);
                        	if(fabs(eta)<1.0)        mon.fillHisto("muTightFake_etabin1",tags, tightMuons[j].pt(), weight);
                        	else if(fabs(eta)<1.479) mon.fillHisto("muTightFake_etabin2",tags, tightMuons[j].pt(), weight);     
                        	else if(fabs(eta)<2.00)  mon.fillHisto("muTightFake_etabin3",tags, tightMuons[j].pt(), weight);
                        	else                     mon.fillHisto("muTightFake_etabin4",tags, tightMuons[j].pt(), weight);
        		}
		}
	}









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
/*
    if(outTxtFile_full)fclose(outTxtFile_full);
    if(outTxtFile_final)fclose(outTxtFile_final);
*/
}



std::vector<LorentzVector> checkDiMuonMass(std::vector<LorentzVector> leps)
{
    if(leps.size() <= 1 ) return leps;
    std::vector<LorentzVector> toSaveLeps;
    for(size_t j=0; j<leps.size(); j++){
	bool iskeep(true);
	for(size_t k=0; k<leps.size(); k++)
       	{
	    if(k==j) continue;
	    LorentzVector dilep(leps[j]+leps[k]);
	    double mass = dilep.mass();
	    if(mass < 20 || (mass > 76 && mass < 106))
	    {
		iskeep = false;
		break;
	    }
        }
	if(iskeep) toSaveLeps.push_back(leps[j]); 
    }

    return toSaveLeps;
}


std::vector<LorentzVector> checkDiElectronMass(std::vector<LorentzVector> leps)
{
    if(leps.size() <= 1 ) return leps;
    std::vector<LorentzVector> toSaveLeps;
    for(size_t j=0; j<leps.size(); j++){
        bool iskeep(true);
        for(size_t k=0; k<leps.size(); k++)
        {
            if(k==j) continue;
            LorentzVector dilep(leps[j]+leps[k]);
            double mass = dilep.mass();
            if(mass < 20 || (mass > 60 && mass < 120))
            {
                iskeep = false;
                break;
            }
        }
        if(iskeep) toSaveLeps.push_back(leps[j]);  
    }

    return toSaveLeps;
}
