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
    bool isSingleMuPD(!isMC && url.Contains("SingleMu"));
    bool isSingleElePD(!isMC && url.Contains("SingleEle"));
    /*
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


    //##################################################################################
    //##########################    INITIATING     TREES      ##########################
    //##################################################################################

    TTree *llvv_tree = new TTree("llvv_tree","Tree for TagandProbe");

    int tag_elepassTight(0);
    int tag_elepassLoose(0);
    int tag_mupassTight(0);
    int tag_mupassLoose(0);
    int probe_elepassTight(0);
    int probe_elepassLoose(0);
    int probe_mupassTight(0);
    int probe_mupassLoose(0);

    float probe_mu_pt;
    float probe_mu_eta;
    float probe_mu_abseta;
    float tag_mu_pt;
    float tag_mu_eta;
    float tag_mu_abseta;

    float mass;

    float probe_ele_pt;
    float probe_ele_eta;
    float probe_ele_abseta;
    float tag_ele_pt;
    float tag_ele_eta;
    float tag_ele_abseta;

    if(fType==EE) {
        llvv_tree->Branch("tag_ele_passTight", &tag_elepassTight ,"tag_ele_passTight/I");
        llvv_tree->Branch("tag_ele_passLoose", &tag_elepassLoose ,"tag_ele_passLoose/I");
        llvv_tree->Branch("probe_ele_passTight", &probe_elepassTight ,"probe_ele_passTight/I");
        llvv_tree->Branch("probe_ele_passLoose", &probe_elepassLoose ,"probe_ele_passLoose/I");

        llvv_tree->Branch("probe_ele_pt", &probe_ele_pt ,"probe_ele_pt/F");
        llvv_tree->Branch("probe_ele_eta", &probe_ele_eta ,"probe_ele_eta/F");
        llvv_tree->Branch("probe_ele_abseta", &probe_ele_abseta ,"probe_ele_abseta/F");
        llvv_tree->Branch("tag_ele_pt", &tag_ele_pt ,"tag_ele_pt/F");
        llvv_tree->Branch("tag_ele_eta", &tag_ele_eta ,"tag_ele_eta/F");
        llvv_tree->Branch("tag_ele_abseta", &tag_ele_abseta ,"tag_ele_abseta/F");
        llvv_tree->Branch("mass", &mass, "mass/F");
    }

    if(fType==MUMU) {
        llvv_tree->Branch("tag_mu_passTight", &tag_mupassTight ,"tag_mu_passTight/I");
        llvv_tree->Branch("tag_mu_passLoose", &tag_mupassLoose ,"tag_mu_passLoose/I");
        llvv_tree->Branch("probe_mu_passTight", &probe_mupassTight ,"probe_mu_passTight/I");
        llvv_tree->Branch("probe_mu_passLoose", &probe_mupassLoose ,"probe_mu_passLoose/I");

        llvv_tree->Branch("probe_mu_pt", &probe_mu_pt ,"probe_mu_pt/F");
        llvv_tree->Branch("probe_mu_eta", &probe_mu_eta ,"probe_mu_eta/F");
        llvv_tree->Branch("probe_mu_abseta", &probe_mu_abseta ,"probe_mu_abseta/F");
        llvv_tree->Branch("tag_mu_pt", &tag_mu_pt ,"tag_mu_pt/F");
        llvv_tree->Branch("tag_mu_eta", &tag_mu_eta ,"tag_mu_eta/F");
	llvv_tree->Branch("tag_mu_abseta", &tag_mu_abseta ,"tag_mu_abseta/F");
        llvv_tree->Branch("mass", &mass, "mass/F");
    }


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
    //mon.addHistogram( new TH1F( "pfmet_raw",        ";E_{T}^{miss} [GeV];Events", 50,0,100) );
    //mon.addHistogram( new TH1F( "elewmt_raw",          ";#it{m}_{T} [GeV];Events", 50,0,20) );
    //mon.addHistogram( new TH1F( "redmet_raw",        ";Red-E_{T}^{miss} [GeV];Events", 50,0,500) );
    //mon.addHistogram( new TH1F( "met_sig_raw",        ";E_{T}^{miss} Significance;Events", 50,0,100) );
    //mon.addHistogram( new TH1F( "mtless_ZMet_raw",      ";#it{m}_{T}(Z, E_{T}^{miss}) [GeV];Events", 50,0,500) );



    double fakePt[11]= {20,25,30,35,40,45,50,60,70,80,100};
    double fakeEta[9]= {-2.5,-2.,-1.479,-1.,0.,1.,1.479,2.,2.5};

    mon.addHistogram( new TH1F( "eleLooseLepEffPt",       ";Loose #it{p}_{T}^{e} [GeV];Events", 10,fakePt) );
    mon.addHistogram( new TH1F( "eleTightLepEffPt",       ";Tight #it{p}_{T}^{e} [GeV];Events", 10,fakePt) );
    mon.addHistogram( new TH1F( "eleLooseLepEffEta",       ";Loose #it{p}_{T}^{e} [GeV];Events", 8,fakeEta) );
    mon.addHistogram( new TH1F( "eleTightLepEffEta",       ";Tight #it{p}_{T}^{e} [GeV];Events", 8,fakeEta) );

    mon.addHistogram( new TH1F( "eleLooseLepEff_etabin1",       ";Loose #it{p}_{T}^{e} [GeV];Events", 10,fakePt) );
    mon.addHistogram( new TH1F( "eleTightLepEff_etabin1",       ";Tight #it{p}_{T}^{e} [GeV];Events", 10,fakePt) );
    mon.addHistogram( new TH1F( "eleLooseLepEff_etabin2",       ";Loose #it{p}_{T}^{e} [GeV];Events", 10,fakePt) );
    mon.addHistogram( new TH1F( "eleTightLepEff_etabin2",       ";Tight #it{p}_{T}^{e} [GeV];Events", 10,fakePt) );
    mon.addHistogram( new TH1F( "eleLooseLepEff_etabin3",       ";Loose #it{p}_{T}^{e} [GeV];Events", 10,fakePt) );
    mon.addHistogram( new TH1F( "eleTightLepEff_etabin3",       ";Tight #it{p}_{T}^{e} [GeV];Events", 10,fakePt) );
    mon.addHistogram( new TH1F( "eleLooseLepEff_etabin4",       ";Loose #it{p}_{T}^{e} [GeV];Events", 10,fakePt) );
    mon.addHistogram( new TH1F( "eleTightLepEff_etabin4",       ";Tight #it{p}_{T}^{e} [GeV];Events", 10,fakePt) );


    mon.addHistogram( new TH1F( "muLooseLepEffPt",        ";Loose #it{p}_{T}^{#mu} [GeV];Events", 10,fakePt) );
    mon.addHistogram( new TH1F( "muTightLepEffPt",        ";Tight #it{p}_{T}^{#mu} [GeV];Events", 10,fakePt) );
    mon.addHistogram( new TH1F( "muLooseLepEffEta",        ";Loose #it{p}_{T}^{#mu} [GeV];Events", 8,fakeEta) );
    mon.addHistogram( new TH1F( "muTightLepEffEta",        ";Tight #it{p}_{T}^{#mu} [GeV];Events", 8,fakeEta) );

    mon.addHistogram( new TH1F( "muLooseLepEff_etabin1",       ";Loose #it{p}_{T}^{#mu} [GeV];Events", 10,fakePt) );
    mon.addHistogram( new TH1F( "muTightLepEff_etabin1",       ";Tight #it{p}_{T}^{#mu} [GeV];Events", 10,fakePt) );
    mon.addHistogram( new TH1F( "muLooseLepEff_etabin2",       ";Loose #it{p}_{T}^{#mu} [GeV];Events", 10,fakePt) );
    mon.addHistogram( new TH1F( "muTightLepEff_etabin2",       ";Tight #it{p}_{T}^{#mu} [GeV];Events", 10,fakePt) );
    mon.addHistogram( new TH1F( "muLooseLepEff_etabin3",       ";Loose #it{p}_{T}^{#mu} [GeV];Events", 10,fakePt) );
    mon.addHistogram( new TH1F( "muTightLepEff_etabin3",       ";Tight #it{p}_{T}^{#mu} [GeV];Events", 10,fakePt) );
    mon.addHistogram( new TH1F( "muLooseLepEff_etabin4",       ";Loose #it{p}_{T}^{#mu} [GeV];Events", 10,fakePt) );
    mon.addHistogram( new TH1F( "muTightLepEff_etabin4",       ";Tight #it{p}_{T}^{#mu} [GeV];Events", 10,fakePt) );







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

        if(isV0JetsMC) {
            if(ev.mc_nup>5) continue;
        }

        PhysicsEvent_t phys=getPhysicsEventFrom(ev);

        //event category
        //bool isSameFlavor(ev.cat==MUMU || ev.cat==EE);
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

        //split inclusive DY sample into DYToLL and DYToTauTau
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

        LorentzVector lep1=phys.leptons[0];
        LorentzVector lep2=phys.leptons[1];
        LorentzVector zll(lep1+lep2);
        //bool passtightZmass(fabs(zll.mass()-91)<5);

        //check alternative selections for the dilepton
        double llScaleFactor(1.0),llTriggerEfficiency(1.0);



        tag_elepassTight=0;
        tag_elepassLoose=0;
        tag_mupassTight=0;
        tag_mupassLoose=0;
        probe_elepassTight=0;
        probe_elepassLoose=0;
        probe_mupassTight=0;
        probe_mupassLoose=0;

        // looping leptons (what I want)  begin
        for(size_t ilep=0; ilep<2; ilep++) {

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
                        && phys.leptons[ilep].trkchi2 < 50/*10*/
                        && ev.mn_validMuonHits[lpid] > 0
                        && ev.mn_nMatchedStations[lpid] > 1
                        && fabs(phys.leptons[ilep].d0)< 100/*0.2*/
                        && fabs(phys.leptons[ilep].dZ)<0.5
                        && phys.leptons[ilep].trkValidPixelHits>0
                        && ev.mn_trkLayersWithMeasurement[lpid]>5
                        //&& relIso< 1.0 /*0.2*/
                  ) {
                    looseMuons_raw.push_back(phys.leptons[ilep]);
                    if(ilep==0) {
                        tag_mupassLoose=1;
                        tag_mu_pt=phys.leptons[ilep].pt();
                        tag_mu_eta=phys.leptons[ilep].eta();
			tag_mu_abseta=fabs(phys.leptons[ilep].eta());
                    }
                    if(ilep==1) {
                        probe_mupassLoose=1;
                        probe_mu_pt=phys.leptons[ilep].pt();
                        probe_mu_eta=phys.leptons[ilep].eta();
                        probe_mu_abseta=fabs(phys.leptons[ilep].eta());
                    }
                    mass = (phys.leptons[0]+phys.leptons[1]).mass();
                }
                if( hasObjectId(ev.mn_idbits[lpid], MID_TIGHT) && relIso<0.2)    {
                    tightMuons_raw.push_back(phys.leptons[ilep]);
                    if(ilep==0) tag_mupassTight=1;
                    if(ilep==1) probe_mupassTight=1;
                }

                llScaleFactor *= lsf.getLeptonEfficiency(phys.leptons[ilep].pt(),fabs(phys.leptons[ilep].eta()),13).first;

            } else {
                int wps[]= {    EgammaCutBasedEleId::LOOSE, // 0
                                EgammaCutBasedEleId::MEDIUM,  // 1
                                EgammaCutBasedEleId::VETO, //2
                                EgammaCutBasedEleId::FakeRateLOOSE //3
                           };
                llScaleFactor *= lsf.getLeptonEfficiency(phys.leptons[ilep].pt(),fabs(phys.leptons[ilep].eta()),11).first;
                for(int iwp=0; iwp<4; iwp++) {

                    bool passWp = EgammaCutBasedEleId::PassWP(EgammaCutBasedEleId::WorkingPoint(wps[iwp]),
                                  (fabs(phys.leptons[ilep].eta())<1.4442),
                                  phys.leptons[ilep].pt(), phys.leptons[ilep].eta(),
                                  ev.en_detain[lpid],  ev.en_dphiin[lpid], ev.en_sihih[lpid], ev.en_hoe[lpid],
                                  ev.en_ooemoop[lpid], phys.leptons[ilep].d0, phys.leptons[ilep].dZ,
                                  0., 0., 0.,
                                  !hasObjectId(ev.en_idbits[lpid], EID_CONVERSIONVETO),0,ev.rho);
                    if(passWp && iwp==3) {
                        looseElectrons_raw.push_back(phys.leptons[ilep]);
                        if(ilep==0) {
                            tag_elepassLoose=1;
                            tag_ele_pt=phys.leptons[ilep].pt();
                            tag_ele_eta=phys.leptons[ilep].eta();
                            tag_ele_abseta=fabs(phys.leptons[ilep].eta());
                        }
                        if(ilep==1) {
                            probe_elepassLoose=1;
                            probe_ele_pt=phys.leptons[ilep].pt();
                            probe_ele_eta=phys.leptons[ilep].eta();
                            probe_ele_abseta=fabs(phys.leptons[ilep].eta());
                        }
                        mass = (phys.leptons[0]+phys.leptons[1]).mass();
                    }
                    if(passWp && iwp==1 && relIso<0.15) {
                        tightElectrons_raw.push_back(phys.leptons[ilep]);
                        if(ilep==0) tag_elepassTight=1;
                        if(ilep==1) probe_elepassTight=1;
                    }
                    if(!use2011Id) {
                        //llScaleFactor *= 1;
                        llTriggerEfficiency *= 1.0; //electronTriggerEfficiency(phys.leptons[ilep].pt(),fabs(phys.leptons[ilep].eta()),2012);
                    }
                }
            }


        } // loop all leptons end


        llvv_tree->Fill();

        // swap the tag and probe leg, increase statistics
        bool doswap(true);
        if(doswap) {
            //electron
            bool temp_tag_elepassLoose = tag_elepassLoose;
            float temp_tag_ele_pt = tag_ele_pt;
            float temp_tag_ele_eta = tag_ele_eta;
            float temp_tag_ele_abseta = tag_ele_abseta;
            bool temp_tag_elepassTight = tag_elepassTight;

            tag_elepassLoose = probe_elepassLoose;
            tag_ele_pt = probe_ele_pt;
            tag_ele_eta = probe_ele_eta;
            tag_ele_abseta = probe_ele_abseta;
            tag_elepassTight = probe_elepassTight;

            probe_elepassLoose = temp_tag_elepassLoose;
            probe_ele_pt = temp_tag_ele_pt;
            probe_ele_eta = temp_tag_ele_eta;
            probe_ele_abseta = temp_tag_ele_abseta;
            probe_elepassTight = temp_tag_elepassTight;

            //muon
            bool temp_tag_mupassLoose = tag_mupassLoose;
            float temp_tag_mu_pt = tag_mu_pt;
            float temp_tag_mu_eta = tag_mu_eta;
            float temp_tag_mu_abseta = tag_mu_abseta;
            bool temp_tag_mupassTight = tag_mupassTight;

            tag_mupassLoose = probe_mupassLoose;
            tag_mu_pt = probe_mu_pt;
            tag_mu_eta = probe_mu_eta;
            tag_mu_abseta = probe_mu_abseta;
            tag_mupassTight = probe_mupassTight;

            probe_mupassLoose = temp_tag_mupassLoose;
            probe_mu_pt = temp_tag_mu_pt;
            probe_mu_eta = temp_tag_mu_eta;
            probe_mu_abseta = temp_tag_mu_abseta;
            probe_mupassTight = temp_tag_mupassTight;


            llvv_tree->Fill();
        }




        //
        //STD PF JET ANALYSIS
        //
        PhysicsObjectJetCollection &aJets = ( useJERsmearing ? variedAJets[0] : recoJets );
        PhysicsObjectJetCollection aGoodIdJets;

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

            if(useJetsOnlyInTracker && fabs(aJets[ijet].eta())>2.5) continue;

            aGoodIdJets.push_back(aJets[ijet]);
            if(aJets[ijet].pt()>30)nAJetsGood30++;

            if(aJets[ijet].pt()>20 /*&& fabs(aJets[ijet].eta())<2.5*/)  nABtags += (aJets[ijet].btag6>0.244);
            if(aJets[ijet].pt()>20 /*&& fabs(aJets[ijet].eta())<2.5*/)  nCSVMtags += (aJets[ijet].btag6>0.679);
            if(aJets[ijet].pt()>20 /*&& fabs(aJets[ijet].eta())<2.5*/)  nCSVTtags += (aJets[ijet].btag6>0.898);
        }




        //#########################################################
        //####  RUN PRESELECTION AND CONTROL REGION PLOTS  ########
        //#########################################################

        if(isMC && use2011Id) weight *= llScaleFactor*llTriggerEfficiency;

        if(hasTrigger)                 {
            mon.fillHisto("eventflow",tags,0,weight);
        } else continue;



        //mon.fillHisto("zmass_raw"			,tags, zll.mass(), weight);
        //mon.fillHisto("zpt_raw"				,tags, zll.pt(),   weight);
        //mon.fillHisto("pfmet_raw"			,tags, zvvs[0].pt(), weight);
        //mon.fillHisto("redmet_raw"			,tags, aRedMet.pt(), weight);
        //mon.fillHisto("met_sig_raw"			,tags, phys.met_sig[1], weight);
        //mon.fillHisto("mtless_ZMet_raw"                 ,tags, aMTmassless, weight);
        //mon.fillHisto("nleptons_raw"			,tags, nextraleptons, weight);
        ////mon.fillHisto("npfjets_raw"			,tags, nAJetsGood30, weight);
        ////mon.fillHisto("npfbjets_raw"			,tags, nABtags, weight);

        //##############################################
        //########  Main Event Selection        ########
        //##############################################



        //electron
        if( ((isMC && ev.cat==EE) || (!isMC && fType==EE)) && tightElectrons_raw.size()>0 ) {
            for(size_t j=0; j<tightElectrons_raw.size(); j++) {

                LorentzVector Tag = tightElectrons_raw[j];
                for(size_t l=0; l<looseElectrons_raw.size(); l++) {

                    // check if exist a ``Loose'' lepton build tag-and-probe Z candidates
                    LorentzVector Probe = looseElectrons_raw[l];
                    if(Tag.eta()==Probe.eta()) continue; // remove the overlap
                    LorentzVector dilep(Tag+Probe);
                    if(fabs(dilep.mass()-91)>5) continue; // remove the non-Z candidates
                    mon.fillHisto("eleLooseLepEffPt",tags, Probe.pt(), weight);
                    mon.fillHisto("eleLooseLepEffEta",tags, Probe.eta(), weight);
                    double eta = Probe.eta();
                    if(fabs(eta)<1.0)        mon.fillHisto("eleLooseLepEff_etabin1",tags, Probe.pt(), weight);
                    else if(fabs(eta)<1.479) mon.fillHisto("eleLooseLepEff_etabin2",tags, Probe.pt(), weight);
                    else if(fabs(eta)<2.00)  mon.fillHisto("eleLooseLepEff_etabin3",tags, Probe.pt(), weight);
                    else			 mon.fillHisto("eleLooseLepEff_etabin4",tags, Probe.pt(), weight);

                    // check this probe lepton is also ``Tight''
                    for(size_t k=0; k<tightElectrons_raw.size(); k++) {
                        LorentzVector TightProbe = tightElectrons_raw[k];
                        if(Probe.eta()!=TightProbe.eta()) continue; // check if the lepton is same or not
                        mon.fillHisto("eleTightLepEffPt",tags, TightProbe.pt(), weight);
                        mon.fillHisto("eleTightLepEffEta",tags, TightProbe.eta(), weight);
                        double eta = TightProbe.eta();
                        if(fabs(eta)<1.0)        mon.fillHisto("eleTightLepEff_etabin1",tags, TightProbe.pt(), weight);
                        else if(fabs(eta)<1.479) mon.fillHisto("eleTightLepEff_etabin2",tags, TightProbe.pt(), weight);
                        else if(fabs(eta)<2.00)  mon.fillHisto("eleTightLepEff_etabin3",tags, TightProbe.pt(), weight);
                        else                     mon.fillHisto("eleTightLepEff_etabin4",tags, TightProbe.pt(), weight);

                    }

                }

            }
        }




        //muon
        if( ((isMC && ev.cat==MUMU) || (!isMC && fType==MUMU)) && tightMuons_raw.size()>0 ) {
            for(size_t j=0; j<tightMuons_raw.size(); j++) {

                LorentzVector Tag = tightMuons_raw[j];
                for(size_t l=0; l<looseMuons_raw.size(); l++) {

                    // check if exist a ``Loose'' lepton build tag-and-probe Z candidates
                    LorentzVector Probe =looseMuons_raw[l];
                    if(Tag.eta()==Probe.eta()) continue; // remove the overlap
                    LorentzVector dilep(Tag+Probe);
                    if(fabs(dilep.mass()-91)>5) continue; // remove the non-Z candidates
                    mon.fillHisto("muLooseLepEffPt",tags, Probe.pt(), weight);
                    mon.fillHisto("muLooseLepEffEta",tags, Probe.eta(), weight);
                    double eta = Probe.eta();
                    if(fabs(eta)<1.0)        mon.fillHisto("muLooseLepEff_etabin1",tags, Probe.pt(), weight);
                    else if(fabs(eta)<1.479) mon.fillHisto("muLooseLepEff_etabin2",tags, Probe.pt(), weight);
                    else if(fabs(eta)<2.00)  mon.fillHisto("muLooseLepEff_etabin3",tags, Probe.pt(), weight);
                    else			 mon.fillHisto("muLooseLepEff_etabin4",tags, Probe.pt(), weight);

                    // check this probe lepton is also ``Tight''
                    for(size_t k=0; k<tightMuons_raw.size(); k++) {
                        LorentzVector TightProbe = tightMuons_raw[k];
                        if(Probe.eta()!=TightProbe.eta()) continue; // check if the lepton is same or not
                        mon.fillHisto("muTightLepEffPt",tags, TightProbe.pt(), weight);
                        mon.fillHisto("muTightLepEffEta",tags, TightProbe.eta(), weight);
                        double eta = TightProbe.eta();
                        if(fabs(eta)<1.0)        mon.fillHisto("muTightLepEff_etabin1",tags, TightProbe.pt(), weight);
                        else if(fabs(eta)<1.479) mon.fillHisto("muTightLepEff_etabin2",tags, TightProbe.pt(), weight);
                        else if(fabs(eta)<2.00)  mon.fillHisto("muTightLepEff_etabin3",tags, TightProbe.pt(), weight);
                        else                     mon.fillHisto("muTightLepEff_etabin4",tags, TightProbe.pt(), weight);

                    }

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

    ofile->mkdir("tpTree");
    ofile->cd("tpTree");
    llvv_tree->Write();
    ofile->Close();
    /*
        if(outTxtFile_full)fclose(outTxtFile_full);
        if(outTxtFile_final)fclose(outTxtFile_final);
    */
}


