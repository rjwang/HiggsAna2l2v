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

    bool isV0JetsMC(isMC && (url.Contains("DYJetsToLL_50toInf") || url.Contains("WJets")));

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
    //bool runSystematics                        = runProcess.getParameter<bool>("runSystematics");

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
    //mon.addHistogram( new TH1F( "pfmet_raw",        ";E_{T}^{miss} [GeV];Events", 50,0,100) );


    // electron channel
    mon.addHistogram( new TH1F( "elepfmet_raw_L",        ";E_{T}^{miss} [GeV];Events", 50,0,100) );
    mon.addHistogram( new TH1F( "elewmt_raw_L",          ";#it{m}_{T}(e, E_{T}^{miss}) [GeV];Events", 50,0,150) );
    mon.addHistogram( new TH1F( "elept_raw_L",           ";#it{p}^{e}_{T} [GeV];Events", 50,20,100) );
    mon.addHistogram( new TH1F( "elewmt_metgep50_L",     ";#it{m}_{T}(e, E_{T}^{miss}) [GeV];Events", 50,0,150) );
    mon.addHistogram( new TH1F( "eledphilepj_raw_L",     ";#Delta#it{phi}(l,j) [rad];Events", 50,0,TMath::Pi()) );
    mon.addHistogram( new TH1F( "eleZmass_raw_L", 	 ";#it{m}_{ee} [GeV];Events", 50,0,300) );
    mon.addHistogram( new TH1F( "eleLJetpt_raw_L",        ";#it{p}_{T}(Leading Jet) [GeV];Events", 50,15,200) );

    mon.addHistogram( new TH1F( "elepfmet_raw_T",        ";E_{T}^{miss} [GeV];Events", 50,0,100) );
    mon.addHistogram( new TH1F( "elewmt_raw_T",          ";#it{m}_{T}(e, E_{T}^{miss}) [GeV];Events", 50,0,150) );
    mon.addHistogram( new TH1F( "elept_raw_T",          ";#it{p}^{e}_{T} [GeV];Events", 50,20,100) );
    mon.addHistogram( new TH1F( "elewmt_metgep50_T",          ";#it{m}_{T}(e, E_{T}^{miss}) [GeV];Events", 50,0,150) );
    mon.addHistogram( new TH1F( "eledphilepj_raw_T",     ";#Delta#it{phi}(l,j) [rad];Events", 50,0,TMath::Pi()) );
    mon.addHistogram( new TH1F( "eleZmass_raw_T",        ";#it{m}_{ee} [GeV];Events", 50,0,300) );
    mon.addHistogram( new TH1F( "eleLJetpt_raw_T",        ";#it{p}_{T}(Leading Jet) [GeV];Events", 50,15,200) );


    // muon channel
    mon.addHistogram( new TH1F( "mupfmet_raw_L",        ";E_{T}^{miss} [GeV];Events", 50,0,100) );
    mon.addHistogram( new TH1F( "muwmt_raw_L",          ";#it{m}_{T}(#mu, E_{T}^{miss}) [GeV];Events", 50,0,150) );
    mon.addHistogram( new TH1F( "mupt_raw_L",          ";#it{p}^{#mu}_{T} [GeV];Events", 50,20,100) );
    mon.addHistogram( new TH1F( "muwmt_metgep50_L",          ";#it{m}_{T}(#mu, E_{T}^{miss}) [GeV];Events", 50,0,150) );
    mon.addHistogram( new TH1F( "mudphilepj_raw_L",     ";#Delta#it{phi}(l,j) [rad];Events", 50,0,TMath::Pi()) );
    mon.addHistogram( new TH1F( "muZmass_raw_L",        ";#it{m}_{#mu#mu} [GeV];Events", 50,0,300) );
    mon.addHistogram( new TH1F( "muLJetpt_raw_L",        ";#it{p}_{T}(Leading Jet) [GeV];Events", 50,15,200) );

    mon.addHistogram( new TH1F( "mupfmet_raw_T",        ";E_{T}^{miss} [GeV];Events", 50,0,100) );
    mon.addHistogram( new TH1F( "muwmt_raw_T",          ";#it{m}_{T}(#mu, E_{T}^{miss}) [GeV];Events", 50,0,150) );
    mon.addHistogram( new TH1F( "mupt_raw_T",          ";#it{p}^{#mu}_{T} [GeV];Events", 50,20,100) );
    mon.addHistogram( new TH1F( "muwmt_metgep50_T",          ";#it{m}_{T}(#mu, E_{T}^{miss}) [GeV];Events", 50,0,150) );
    mon.addHistogram( new TH1F( "mudphilepj_raw_T",     ";#Delta#it{phi}(l,j) [rad];Events", 50,0,TMath::Pi()) );
    mon.addHistogram( new TH1F( "muZmass_raw_T",        ";#it{m}_{#mu#mu} [GeV];Events", 50,0,300) );
    mon.addHistogram( new TH1F( "muLJetpt_raw_T",        ";#it{p}_{T}(Leading Jet) [GeV];Events", 50,15,200) );




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

        if(isV0JetsMC && ev.mc_nup>5) continue;

        /*
        if(isMC){
            int npromptGammas = ((ev.mccat>>28)&0xf);
            if(mctruthmode==22 && npromptGammas>0) continue;
        }
        */

        PhysicsEvent_t phys=getPhysicsEventFrom(ev);

        //event category
        //bool isSameFlavor(ev.cat==MUMU || ev.cat==EE);
        TString tag_cat;

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

        //check alternative selections for the dilepton
        double llScaleFactor(1.0),llTriggerEfficiency(1.0);

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
                        && phys.leptons[ilep].trkchi2 < 50/*10*/
                        && ev.mn_validMuonHits[lpid] > 0
                        && ev.mn_nMatchedStations[lpid] > 1
                        && fabs(phys.leptons[ilep].d0)< 100/*0.2*/
                        && fabs(phys.leptons[ilep].dZ)< 0.5
                        && phys.leptons[ilep].trkValidPixelHits>0
                        && ev.mn_trkLayersWithMeasurement[lpid]>5
                        //&& relIso< 1.0/*0.2*/
                  ) {
                    looseMuons_raw.push_back(phys.leptons[ilep]);
                }
                if( hasObjectId(ev.mn_idbits[lpid], MID_TIGHT) && relIso<0.2) {
                    tightMuons_raw.push_back(phys.leptons[ilep]);
                }

                llScaleFactor *= lsf.getLeptonEfficiency(phys.leptons[ilep].pt(),fabs(phys.leptons[ilep].eta()),13).first;

            } else {
                int wps[]= { 	EgammaCutBasedEleId::LOOSE, // 0
                                EgammaCutBasedEleId::MEDIUM,  // 1
                                EgammaCutBasedEleId::VETO //2
                           };
                llScaleFactor *= lsf.getLeptonEfficiency(phys.leptons[ilep].pt(),fabs(phys.leptons[ilep].eta()),11).first;
                for(int iwp=0; iwp<3; iwp++) {

                    bool passWp = EgammaCutBasedEleId::PassWP(EgammaCutBasedEleId::WorkingPoint(wps[iwp]),
                                  (fabs(phys.leptons[ilep].eta())<1.4442),
                                  phys.leptons[ilep].pt(), phys.leptons[ilep].eta(),
                                  ev.en_detain[lpid],  ev.en_dphiin[lpid], ev.en_sihih[lpid], ev.en_hoe[lpid],
                                  ev.en_ooemoop[lpid], phys.leptons[ilep].d0, phys.leptons[ilep].dZ,
                                  0., 0., 0.,
                                  !hasObjectId(ev.en_idbits[lpid], EID_CONVERSIONVETO),0,ev.rho);
                    if(passWp && iwp==2 /*&& relIso<1.0*/) {
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


            //aClusteredMetP4 -= aJets[ijet];

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


        // investigating dijet ctrl sammple
        // control plots
        LorentzVector LeadingJet;
        double largest_pt=-999;
        for(size_t j=0; j<aGoodIdJets.size(); j++) {
            double pt = aGoodIdJets[j].pt();
            if(pt>largest_pt) {
                largest_pt = pt;
                LeadingJet = aGoodIdJets[j];
            }
        }


        bool passLooseEle(looseElectrons_raw.size()>0 && aGoodIdJets.size()>0 /*&& LeadingJet.pt()>30*/);
        bool passTightEle(tightElectrons_raw.size()>0 && aGoodIdJets.size()>0 /*&& LeadingJet.pt()>30*/);
        bool passLooseMu(looseMuons_raw.size()>0 && aGoodIdJets.size()>0 /*&& LeadingJet.pt()>30*/);
        bool passTightMu(tightMuons_raw.size()>0 && aGoodIdJets.size()>0 /*&& LeadingJet.pt()>30*/);

        hasEtrigger &= (isMC || (!isMC && fType==EE));
        hasMtrigger &= (isMC || (!isMC && fType==MUMU));


        if(hasEtrigger && passLooseEle) {
            mon.fillHisto("elepfmet_raw_L",   tags, zvvs[0].pt(), weight);
            mon.fillHisto("eleLJetpt_raw_L",  tags, LeadingJet.pt(), weight);
            for(size_t j=0; j<looseElectrons_raw.size(); j++) {
                double Wmt = METUtils::transverseMass(looseElectrons_raw[j],zvvs[0],false);
                double dphi=fabs(deltaPhi(looseElectrons_raw[j].phi(),LeadingJet.phi()));
                mon.fillHisto("eledphilepj_raw_L",tags, dphi, weight);
                mon.fillHisto("elewmt_raw_L",tags, Wmt, weight);
                mon.fillHisto("elept_raw_L",tags, looseElectrons_raw[j].pt(), weight);
                if(zvvs[0].pt()>50) mon.fillHisto("elewmt_metgep50_L",tags, Wmt, weight);

                for(size_t k=j+1; k<looseElectrons_raw.size(); k++) {
                    LorentzVector dilep(looseElectrons_raw[j]+looseElectrons_raw[k]);
                    mon.fillHisto("eleZmass_raw_L", tags, dilep.mass(), weight);
                }
            }
        }

        if(hasEtrigger && passTightEle) {
            mon.fillHisto("elepfmet_raw_T",   tags, zvvs[0].pt(), weight);
            mon.fillHisto("eleLJetpt_raw_T",  tags, LeadingJet.pt(), weight);
            for(size_t j=0; j<tightElectrons_raw.size(); j++) {
                double Wmt = METUtils::transverseMass(tightElectrons_raw[j],zvvs[0],false);
                double dphi=fabs(deltaPhi(tightElectrons_raw[j].phi(),LeadingJet.phi()));
                mon.fillHisto("eledphilepj_raw_T",tags, dphi, weight);
                mon.fillHisto("elewmt_raw_T",tags, Wmt, weight);
                mon.fillHisto("elept_raw_T",tags, tightElectrons_raw[j].pt(), weight);
                if(zvvs[0].pt()>50) mon.fillHisto("elewmt_metgep50_T",tags, Wmt, weight);

                for(size_t k=j+1; k<tightElectrons_raw.size(); k++) {
                    LorentzVector dilep(tightElectrons_raw[j]+tightElectrons_raw[k]);
                    mon.fillHisto("eleZmass_raw_T", tags, dilep.mass(), weight);
                }
            }
        }




        if(hasMtrigger && passLooseMu) {
            mon.fillHisto("mupfmet_raw_L",   tags, zvvs[0].pt(), weight);
            mon.fillHisto("muLJetpt_raw_L",  tags, LeadingJet.pt(), weight);
            for(size_t j=0; j<looseMuons_raw.size(); j++) {
                double Wmt = METUtils::transverseMass(looseMuons_raw[j],zvvs[0],false);
                double dphi=fabs(deltaPhi(looseMuons_raw[j].phi(),LeadingJet.phi()));
                mon.fillHisto("mudphilepj_raw_L",tags, dphi, weight);
                mon.fillHisto("muwmt_raw_L",tags, Wmt, weight);
                mon.fillHisto("mupt_raw_L",tags, looseMuons_raw[j].pt(), weight);
                if(zvvs[0].pt()>50) mon.fillHisto("muwmt_metgep50_L",tags, Wmt, weight);

                for(size_t k=j+1; k<looseMuons_raw.size(); k++) {
                    LorentzVector dilep(looseMuons_raw[j]+looseMuons_raw[k]);
                    mon.fillHisto("muZmass_raw_L", tags, dilep.mass(), weight);
                }

            }
        }

        if(hasMtrigger && passTightMu) {
            mon.fillHisto("mupfmet_raw_T",   tags, zvvs[0].pt(), weight);
            mon.fillHisto("muLJetpt_raw_T",  tags, LeadingJet.pt(), weight);
            for(size_t j=0; j<tightMuons_raw.size(); j++) {
                double Wmt = METUtils::transverseMass(tightMuons_raw[j],zvvs[0],false);
                double dphi=fabs(deltaPhi(tightMuons_raw[j].phi(),LeadingJet.phi()));
                mon.fillHisto("mudphilepj_raw_T",tags, dphi, weight);
                mon.fillHisto("muwmt_raw_T",tags, Wmt, weight);
                mon.fillHisto("mupt_raw_T",tags, tightMuons_raw[j].pt(), weight);
                if(zvvs[0].pt()>50) mon.fillHisto("muwmt_metgep50_T",tags, Wmt, weight);

                for(size_t k=j+1; k<tightMuons_raw.size(); k++) {
                    LorentzVector dilep(tightMuons_raw[j]+tightMuons_raw[k]);
                    mon.fillHisto("muZmass_raw_T", tags, dilep.mass(), weight);
                }

            }
        }




        //##############################################
        //########  Main Fake Rate Calculation  ########
        //##############################################




        // drive fake rate
        // fixed dijet control sample
        if(zvvs[0].pt()<20) {
            std::vector<LorentzVector> looseMuons = checkDiMuonMass(looseMuons_raw);
            std::vector<LorentzVector> tightMuons = checkDiMuonMass(tightMuons_raw);
            std::vector<LorentzVector> looseElectrons = checkDiElectronMass(looseElectrons_raw);
            std::vector<LorentzVector> tightElectrons = checkDiElectronMass(tightElectrons_raw);

            if(hasEtrigger) {
                if(passLooseEle) {
                    for(size_t j=0; j<looseElectrons.size(); j++) {

                        double dphi=fabs(deltaPhi(looseElectrons[j].phi(),LeadingJet.phi()));
                        if(dphi<2) continue;

                        //double Wmt = METUtils::transverseMass(looseElectrons[j],zvvs[0],false);
                        //if (Wmt > 20) continue;

                        mon.fillHisto("eleLooseFakePt",tags, looseElectrons[j].pt(), weight);
                        double eta = looseElectrons[j].eta();
                        mon.fillHisto("eleLooseFakeEta",tags, eta, weight);
                        if(fabs(eta)<1.0)        mon.fillHisto("eleLooseFake_etabin1",tags, looseElectrons[j].pt(), weight);
                        else if(fabs(eta)<1.479) mon.fillHisto("eleLooseFake_etabin2",tags, looseElectrons[j].pt(), weight);
                        else if(fabs(eta)<2.00)  mon.fillHisto("eleLooseFake_etabin3",tags, looseElectrons[j].pt(), weight);
                        else                     mon.fillHisto("eleLooseFake_etabin4",tags, looseElectrons[j].pt(), weight);
                    } //loop looseElectrons
                }//passLooseEle

                if(passTightEle) {
                    for(size_t j=0; j<tightElectrons.size(); j++) {
                        double dphi=fabs(deltaPhi(tightElectrons[j].phi(),LeadingJet.phi()));
                        if(dphi<2) continue;
                        //double Wmt = METUtils::transverseMass(tightElectrons[j],zvvs[0],false);
                        //if (Wmt > 20) continue;
                        mon.fillHisto("eleTightFakePt",tags, tightElectrons[j].pt(), weight);
                        double eta = tightElectrons[j].eta();
                        mon.fillHisto("eleTightFakeEta",tags, eta, weight);
                        if(fabs(eta)<1.0) 	 mon.fillHisto("eleTightFake_etabin1",tags, tightElectrons[j].pt(), weight);
                        else if(fabs(eta)<1.479) mon.fillHisto("eleTightFake_etabin2",tags, tightElectrons[j].pt(), weight);
                        else if(fabs(eta)<2.00)  mon.fillHisto("eleTightFake_etabin3",tags, tightElectrons[j].pt(), weight);
                        else                     mon.fillHisto("eleTightFake_etabin4",tags, tightElectrons[j].pt(), weight);
                    }//loop tightElectrons
                } //passTightEle

            } //hasEtrigger


            if(hasMtrigger) {
                if(passLooseMu) {
                    for(size_t j=0; j<looseMuons.size(); j++) {
                        double dphi=fabs(deltaPhi(looseMuons[j].phi(),LeadingJet.phi()));
                        if(dphi<2) continue;
                        double Wmt = METUtils::transverseMass(looseMuons[j],zvvs[0],false);
                        if (Wmt > 20) continue;

                        double eta = looseMuons[j].eta();
                        mon.fillHisto("muLooseFakePt",tags, looseMuons[j].pt(), weight);
                        mon.fillHisto("muLooseFakeEta",tags, eta, weight);
                        if(fabs(eta)<1.0)        mon.fillHisto("muLooseFake_etabin1",tags, looseMuons[j].pt(), weight);
                        else if(fabs(eta)<1.479) mon.fillHisto("muLooseFake_etabin2",tags, looseMuons[j].pt(), weight);
                        else if(fabs(eta)<2.00)  mon.fillHisto("muLooseFake_etabin3",tags, looseMuons[j].pt(), weight);
                        else 			 mon.fillHisto("muLooseFake_etabin4",tags, looseMuons[j].pt(), weight);
                    }
                }//passLooseMu

                if(passTightMu) {
                    for(size_t j=0; j<tightMuons.size(); j++) {
                        double dphi=fabs(deltaPhi(tightMuons[j].phi(),LeadingJet.phi()));
                        if(dphi<2) continue;
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
                } //passTightMu
            }//hasMtrigger

        }//passMET cut









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

}



std::vector<LorentzVector> checkDiMuonMass(std::vector<LorentzVector> leps)
{
    if(leps.size() <= 1 ) return leps;
    std::vector<LorentzVector> toSaveLeps;
    for(size_t j=0; j<leps.size(); j++) {
        bool iskeep(true);
        for(size_t k=0; k<leps.size(); k++) {
            if(k==j) continue;
            LorentzVector dilep(leps[j]+leps[k]);
            double mass = dilep.mass();
            if(mass < 20 || (mass > 76 && mass < 106)) {
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
    for(size_t j=0; j<leps.size(); j++) {
        bool iskeep(true);
        for(size_t k=0; k<leps.size(); k++) {
            if(k==j) continue;
            LorentzVector dilep(leps[j]+leps[k]);
            double mass = dilep.mass();
            if(mass < 20 || (mass > 60 && mass < 120)) {
                iskeep = false;
                break;
            }
        }
        if(iskeep) toSaveLeps.push_back(leps[j]);
    }

    return toSaveLeps;
}
