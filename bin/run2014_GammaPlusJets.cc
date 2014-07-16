#include <iostream>
#include <boost/shared_ptr.hpp>

#include "CMGTools/HiggsAna2l2v/interface/ZZ2l2nuSummaryHandler.h"
#include "CMGTools/HiggsAna2l2v/interface/ZZ2l2nuPhysicsEvent.h"
#include "CMGTools/HiggsAna2l2v/interface/ZHphotonjetsEventHandler.h"
#include "CMGTools/HiggsAna2l2v/interface/METUtils.h"
#include "CMGTools/HiggsAna2l2v/interface/setStyle.h"
#include "CMGTools/HiggsAna2l2v/interface/plotter.h"
#include "CMGTools/HiggsAna2l2v/interface/ObjectFilters.h"
#include "CMGTools/HiggsAna2l2v/interface/SmartSelectionMonitor.h"
#include "CMGTools/HiggsAna2l2v/interface/EventCategory.h"
#include "CMGTools/HiggsAna2l2v/interface/MacroUtils.h"

#include "FWCore/FWLite/interface/AutoLibraryLoader.h"
#include "FWCore/PythonParameterSet/interface/MakeParameterSets.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"

#include "Math/GenVector/Boost.h"

#include "TROOT.h"
#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TRandom2.h"

#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

using namespace std;

//double Vgamma_kFactor = 1.0;//2.79712417075170716e+00;

LorentzVector min(const LorentzVector& a, const LorentzVector& b)
{
    if(a.pt()<=b.pt())return a;
    return b;
}

TString getJetRegion(float eta)
{
    TString reg("TK");
    if(fabs(eta)>2.5)  reg="HEin";
    if(fabs(eta)>2.75) reg="HEout";
    if(fabs(eta)>3)    reg="HF";
    return reg;
}

//
int main(int argc, char* argv[])
{
    // load framework libraries
    gSystem->Load( "libFWCoreFWLite" );
    AutoLibraryLoader::enable();

    //check arguments
    if ( argc < 2 ) {
        std::cout << "Usage : " << argv[0] << " parameters_cfg.py" << std::endl;
        return 0;
    }

    //get configuration
    const edm::ParameterSet &runProcess = edm::readPSetsFrom(argv[1])->getParameter<edm::ParameterSet>("runProcess");

    bool use2011Id = runProcess.getParameter<bool>("is2011");
    cout << "Note: will apply " << (use2011Id ? 2011 : 2012) << " version of the id's" << endl;
    bool useCHS(true);
    bool useJERsmearing(true);
    bool useJetsOnlyInTracker(false);
    bool usePUsubJetId(true);
    // invert ajets <-> jets ...
    // 2011:  ajets -> non CHS,  jets -> CHS
    // 2012:  ajets -> CHS,      jets -> non CHS
    if(use2011Id) useCHS = !useCHS;

    bool isMC = runProcess.getParameter<bool>("isMC");
    int mctruthmode = runProcess.getParameter<int>("mctruthmode");
    TString myTag = runProcess.getParameter<std::string>("tag");
    bool applyQtweights = runProcess.getParameter<bool>("applyQtweights");

    TString dirname = runProcess.getParameter<std::string>("dirName");
    TString uncFile =  runProcess.getParameter<std::string>("jesUncFileName");
    gSystem->ExpandPathName(uncFile);
    JetCorrectionUncertainty jecUnc(uncFile.Data());

    TString url=runProcess.getParameter<std::string>("input");
    TString outFileUrl(gSystem->BaseName(url));
    outFileUrl.ReplaceAll(".root","");
    if(mctruthmode!=0) {
        outFileUrl += "_filt";
        outFileUrl += mctruthmode;
    }
    TString outdir=runProcess.getParameter<std::string>("outdir");

    bool isV0JetsMC(isMC && url.Contains("WJets"));

    TString outTxtUrl= outdir + "/" + outFileUrl + ".txt";
    FILE* outTxtFile = NULL;
    //if(!isMC)
    outTxtFile = fopen(outTxtUrl.Data(), "w");
    printf("TextFile URL = %s\n",outTxtUrl.Data());


    TRandom2 rndGen;
    //event Categorizer
    //  EventCategory eventCategoryInst(0); //inclusive analysis
    EventCategory eventCategoryInst(1); //jet binning
    //  EventCategory eventCategoryInst(2); //vbf binning
    //  EventCategory eventCategoryInst(3); //jet+vbf binning
    //  EventCategory eventCategoryInst(4);  //0,>=1jet, VBF
    ZHphotonjetsEventHandler gammaEvHandler(runProcess);


    //book histograms
    SmartSelectionMonitor mon;
    TH1F* Hcutflow     = (TH1F*) mon.addHistogram(  new TH1F ("cutflow"    , "cutflow"    ,6,0,6) ) ;

    //##############################################
    //######## STUFF FOR CUTS OPTIMIZATION  ########
    //##############################################

    //optimization
    std::vector<double> optim_Cuts1_MET;
    std::vector<double> optim_Cuts1_Balance;
    std::vector<double> optim_Cuts1_DphiZMET;

    bool runOptimization = runProcess.getParameter<bool>("runOptimization");
    if(runOptimization) {
        for(double met=90; met<=100; met+=10) {
            for(double balance=0.2; balance<=0.2; balance+=0.05) {
                for(double dphi=2.7; dphi<=2.7; dphi+=0.1) {
                    optim_Cuts1_MET     .push_back(met);
                    optim_Cuts1_Balance .push_back(balance);
                    optim_Cuts1_DphiZMET.push_back(dphi);
                }
            }
        }
    }
    size_t nOptims = optim_Cuts1_MET.size();


    //make it as a TProfile so hadd does not change the value
    TProfile* Hoptim_cuts1_MET      = (TProfile*) mon.addHistogram( new TProfile ("optim_cut1_MET",";cut index;met",nOptims,0,nOptims) );
    TProfile* Hoptim_cuts1_Balance  = (TProfile*) mon.addHistogram( new TProfile ("optim_cut1_Balance",";cut index;dphi",nOptims,0,nOptims) );
    TProfile* Hoptim_cuts1_DphiZMET = (TProfile*) mon.addHistogram( new TProfile ("optim_cut1_DphiZMET",";cut index;dphi",nOptims,0,nOptims) );

    for(unsigned int index=0; index<nOptims; index++) {
        Hoptim_cuts1_MET        ->Fill(index, optim_Cuts1_MET[index]);
        Hoptim_cuts1_Balance    ->Fill(index, optim_Cuts1_Balance[index]);
        Hoptim_cuts1_DphiZMET   ->Fill(index, optim_Cuts1_DphiZMET[index]);
    }



    //############################
    //1D shapes for limit setting
    //############################
    mon.addHistogram( new TH2F (TString("mt_shapes"),";cut index; #it{m}_{T} [GeV];#Events (/100GeV)",nOptims,0,nOptims,12,0,1200) );

    //############################
    //2D shapes for limit setting
    //############################








    //##############################################
    //######## CONTROL PLOTS FOR SELECTION #########
    //##############################################
    TH1F *h1=(TH1F*) mon.addHistogram( new TH1F ("eventflow_gamma", ";;Events", 5,0,5) );
    h1->GetXaxis()->SetBinLabel(1,"PreSelection");
    h1->GetXaxis()->SetBinLabel(2,"#Delta #phi(Z,E_{T}^{miss})>2.6");
    h1->GetXaxis()->SetBinLabel(3,"Reduced E_{T}^{miss}>110"); //RJ
    h1->GetXaxis()->SetBinLabel(4,"0.8<E_{T}^{miss}/p_{T}^{ll}<1.2");
    h1->GetXaxis()->SetBinLabel(5,"M_{T}>220");

    double qtaxis[100];
    for(size_t i=0; i<40; i++)  qtaxis[i]=2.5*i;       //0-97.5
    for(size_t i=0; i<20; i++)  qtaxis[40+i]=100+5*i;  //100-195
    for(size_t i=0; i<15; i++)  qtaxis[60+i]=200+10*i; //200-340
    for(size_t i=0; i<25; i++)  qtaxis[75+i]=350+25*i; //350-976
    mon.addHistogram( new TH1D( "qt_dataDY",  ";#it{p}_{T}^{#it{#gamma}} [GeV];Events",99,qtaxis));

    double qtaxisType1[65];
    for(size_t i=0; i<40; i++)  qtaxisType1[i]=2.5*i;
    for(size_t i=0; i<20; i++)  qtaxisType1[40+i]=100+5*i;
    for(size_t i=0; i<3; i++)   qtaxisType1[60+i]=200+50*i;
    for(size_t i=0; i<2; i++)   qtaxisType1[63+i]=350+650*i;
    mon.addHistogram( new TH1D( "qtType1_dataDY",        ";#it{p}_{T}^{#it{#gamma}} [GeV];Events",64,qtaxisType1));


    mon.addHistogram( new TH1F( "qt",        ";#it{p}_{T}^{#it{#gamma}} [GeV];Events / (1.0 GeV/c)",1000,0,1000));
    mon.addHistogram( new TH1F( "nvtx_dataDY",      ";Vertices;Events", 50,0,50) );
    double METBins[21]= {0,10,20,30,40,50,60,70,80,90,100,120,140,160,180,200,250,300,350,400,500};
    mon.addHistogram( new TH1F( "pfmet_dataDY",      ";E_{T}^{miss} [GeV];Events", 20, METBins));

    double METBinsType1[11]= {0,10,20,30,40,60,90,130,180,300,500};
    mon.addHistogram( new TH1F( "pfmetType1_dataDY", ";E_{T}^{miss} [GeV];Events", 10, METBinsType1));

    mon.addHistogram( new TH1F( "DphiZMET_dataDY",   ";#Delta#it{#phi}(#it{l^{+}l^{-}},E_{T}^{miss});Events", 50,0,TMath::Pi()) );
    mon.addHistogram( new TH1D( "balancedif_dataDY", ";|E_{T}^{miss}-#it{q}_{T}|/#it{q}_{T};Events", 5,0,1.0) );


    //##############################################
    //###### open the file and get events tree #####
    //##############################################

    //open the file and get events tree
    TFile *file = TFile::Open(url);
    if(file==0) return -1;
    if(file->IsZombie()) return -1;
    ZZ2l2nuSummaryHandler evSummaryHandler;
    if( !evSummaryHandler.attachToTree( (TTree *)file->Get(dirname) ) ) {
        file->Close();
        return -1;
    }
    float cnorm=1.0;
    if(isMC) {
        TH1F* cutflowH = (TH1F *) file->Get("dataAnalyzer/llvv/cutflow");
        if(cutflowH) cnorm=cutflowH->GetBinContent(1);
        printf("cnorm = %f\n",cnorm);
    }
    Hcutflow->SetBinContent(1,cnorm);

    //duplicate checker
    DuplicatesChecker duplicatesChecker;
    int NumberOfDuplicated(0);

    //pileup reweighting
    bool disableJERSmear(false);
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
        if(dataPileupDistribution.size()==0) {
            dataPileupDistribution=mcPileupDistribution;
            disableJERSmear=true;
            cout << "No PU reweighting or JER smearing will be applied" << endl;
        }
    } else mcPileupDistribution=dataPileupDistribution;
    while(mcPileupDistribution.size()<dataPileupDistribution.size())  mcPileupDistribution.push_back(0.0);
    while(mcPileupDistribution.size()>dataPileupDistribution.size())  dataPileupDistribution.push_back(0.0);


    gROOT->cd();  //THIS LINE IS NEEDED TO MAKE SURE THAT HISTOGRAM INTERNALLY PRODUCED IN LumiReWeighting ARE NOT DESTROYED WHEN CLOSING THE FILE
    edm::LumiReWeighting LumiWeights(mcPileupDistribution,dataPileupDistribution);


    //##############################################
    //############## RUN ANALYSIS ##################
    //##############################################
    //run the analysis
    int evStart(0),evEnd(evSummaryHandler.getEntries());
    for( int iev=0; iev<evEnd; iev++) {
        if(iev%1000==0) printf("\r [ %d/100 ] ",int(100*float(iev-evStart)/float(evEnd)));
        evSummaryHandler.getEntry(iev);

        //interpret event
        ZZ2l2nuSummary_t &ev = evSummaryHandler.getEvent();
        if( !isMC && duplicatesChecker.isDuplicate(ev.run,ev.lumi, ev.event)) {
            NumberOfDuplicated++;
            continue;
        }


        if(!use2011Id) {
            if((mctruthmode==22||mctruthmode==111)) {
                if(isMC) {
                    int npromptGammas = ((ev.mccat>>28)&0xf) ;
                    if(mctruthmode==22  && npromptGammas<1) continue;
                    if(mctruthmode==111 && npromptGammas>0) continue;
                }
            } else continue;
        }

        if(isV0JetsMC && ev.mc_nup>5)                          continue;

        //interpret event
        PhysicsEvent_t phys  = getPhysicsEventFrom(ev);
        bool isGammaEvent    = gammaEvHandler.isGood(phys,use2011Id);
        if(!isGammaEvent) continue;


        //event weight
        float weight = 1.0;
        if(isMC)                  {
            weight = LumiWeights.weight( ev.ngenITpu );
        }
        if(!isMC && !use2011Id) 	{
            weight = gammaEvHandler.triggerWeight();
        }

        Hcutflow->Fill(1,1);
        Hcutflow->Fill(2,weight);
        Hcutflow->Fill(3,weight);
        Hcutflow->Fill(4,weight);
        Hcutflow->Fill(5,1);


        //check which event type is required to use (dilepton/or photon)
        std::vector<TString> dilCats;
        dilCats.push_back("ee");
        dilCats.push_back("mumu");
        dilCats.push_back("ll"); // here ll is not ee+mumu, just ee or mumu

        //build the gamma candidate
        float r9               = phys.gammas[0].r9*(isMC ? 1.005 : 1.0);
        //float sietaieta        = phys.gammas[0].sihih;
        bool hasElectronVeto   = phys.gammas[0].hasElectronVeto;
        LorentzVector gamma    = phys.gammas[0];
        //bool isConv=(phys.gammas[0].isConv);

        //
        // EXTRA LEPTONS
        //
        int nextraleptons(0);
        for(size_t ilep=0; ilep<phys.leptons.size(); ilep++) {
            int lpid=phys.leptons[ilep].pid;
            //FIXME: looks like sometimes the gamma comes double counted as a soft-id electron ???
            if(isGammaEvent && fabs(phys.leptons[ilep].id)==11) {
                double dr( deltaR(phys.leptons[ilep],gamma) );
                if(dr<0.1) continue;
            }

            if(fabs(phys.leptons[ilep].id)==13) { //muon
                if(!use2011Id) {
                    if(  (hasObjectId(ev.mn_idbits[lpid], MID_LOOSE) && phys.leptons[ilep].pfRelIsoDbeta()<0.2 && phys.leptons[ilep].pt()>10)
                            || (hasObjectId(ev.mn_idbits[lpid], MID_SOFT) && phys.leptons[ilep].pt()>3) )
                        nextraleptons++;
                }
            } else {
                if(!use2011Id) {
                    if( hasObjectId(ev.en_idbits[lpid],EID_LOOSE) && phys.leptons[ilep].ePFRelIsoCorrected2012(ev.rho,ev.en_sceta[lpid])<0.15
                            && phys.leptons[ilep].pt()>10)
                        nextraleptons++;
                }
            }

        }



        //
        //MET variables
        //
        LorentzVector metP4=phys.met[1];
        if(use2011Id) metP4=phys.met[0];

        //apply JER base corrections to jets (and compute associated variations on the MET variable)
        // std PF
        std::vector<PhysicsObjectJetCollection> variedAJets;
        LorentzVectorCollection zvvs;
        PhysicsObjectJetCollection &recoJets = ( useCHS ? phys.ajets : phys.jets);
        METUtils::computeVariation(recoJets, phys.leptons, metP4, variedAJets, zvvs, &jecUnc);
        if(!useJERsmearing) zvvs[0] = metP4; // this stinks a bit...


        //count the jets
        int nbtags(0),nAJetsGood30(0),nAJetsGood15(0),nClusteredMet15(0);
        double ht(0);
        double mindphijmet(9999.), mindphijmet15(9999.);
        PhysicsObjectJetCollection aGoodIdJets;
        LorentzVector clusteredMet(gamma);
        clusteredMet *= -1;
        LorentzVector Recoil(metP4);
        Recoil *= -1;
        Recoil -= gamma;
        LorentzVector mht(0,0,0,0),unclusteredMet(0,0,0,0);

        PhysicsObjectJetCollection & jetsToUse= (useJERsmearing ? variedAJets[0] : recoJets );
        //if(disableJERSmear) jetsToUse=phys.ajets;
        for(size_t ijet=0; ijet<jetsToUse.size(); ijet++) {
            LorentzVector ijetP4=jetsToUse[ijet];
            if(ijetP4.pt()<15) continue;
            if(deltaR(ijetP4,gamma)<0.2) continue;

            double idphijmet( fabs(deltaPhi(jetsToUse[ijet].phi(),zvvs[0].phi()) ) );
            if( jetsToUse[ijet].pt()>15 ) if(idphijmet<mindphijmet15)  mindphijmet15=idphijmet;
            if( jetsToUse[ijet].pt()>30 ) if(idphijmet<mindphijmet) mindphijmet=idphijmet;

            bool isGoodJet = hasObjectId(jetsToUse[ijet].pid,JETID_LOOSE);
            if(usePUsubJetId)  isGoodJet =hasObjectId(jetsToUse[ijet].pid,JETID_CUTBASED_LOOSE);
            if(!isGoodJet) continue;

            nClusteredMet15++;
            clusteredMet -= ijetP4;
            mht -= ijetP4;
            ht += ijetP4.pt();

            if(useJetsOnlyInTracker && fabs(ijetP4.eta())>2.5) continue;
            nAJetsGood15++;

            aGoodIdJets.push_back(jetsToUse[ijet]);
            if(ijetP4.pt()>30) nAJetsGood30++;

            if( ijetP4.pt()>20 /*&& fabs(ijetP4.eta())<2.5*/) nbtags += (jetsToUse[ijet].btag6>0.244);
        }



        //
        // EVENT SELECTION
        //
        bool passMultiplicityVetoes (nextraleptons==0);
        bool passVgammaCtrl(nextraleptons==1 || nextraleptons==2);
        //gmmma pt threshold should be consistent with dilepton pt //RJ
        bool passKinematics         (gamma.pt()>50. /*&& gamma.pt()<550.*/ && fabs(gamma.eta())<1.4442 && r9<1.0 && r9>0.9 && !hasElectronVeto);
        if(!isMC) passKinematics &= (gamma.pt()>gammaEvHandler.triggerThr());
        bool passBveto              (nbtags==0);
        bool passMETcut(metP4.pt()>90); //Reduced MET: redMet.pt() ; PFMET: metP4.pt()
        double balance=metP4.pt()/gamma.pt();
        double balanceDif=fabs(1.-balance);
        bool passBalanceCut(balance > 0.80 && balance < 1.20);

        double dphiZllmet=fabs(deltaPhi(gamma.phi(),metP4.phi()));
        bool passDphiZMETcut(dphiZllmet>2.7);


        //event category
        int eventSubCat    = eventCategoryInst.Get(phys,&aGoodIdJets);
        TString tag_subcat = eventCategoryInst.GetLabel(eventSubCat);
        std::vector<TString> subcats;
        subcats.push_back("");

        if(tag_subcat=="eq0jets" || tag_subcat=="eq1jets") subcats.push_back("lesq1jets");
        if(tag_subcat != "geq2jets") subcats.push_back(tag_subcat);






        //now do the control plots
        std::map<TString, float> Gweights_; //old mehthod to retrive qt weight,
        if(isGammaEvent) Gweights_ = gammaEvHandler.getWeights(ev,phys,tag_subcat,metP4.phi());
        //if(isGammaEvent && (tag_subcat=="eq0jets" || tag_subcat=="eq1jets")) Gweights_ = gammaEvHandler.getWeights(ev,phys,"lesq1jets",metP4.phi());

        //testing with new qt fit weights
        std::map<TString, double> qtWeights;
        if(isGammaEvent && !use2011Id) qtWeights = gammaEvHandler.getQtFitWeights(phys,tag_subcat);
        if(isGammaEvent && use2011Id)  qtWeights = gammaEvHandler.get2011QtFitWeights(phys,tag_subcat);

        //now do the control plots
        if(!passKinematics) continue;
        if(!passBveto) continue;

        //fill the histograms now
        for(size_t idc=0; idc<dilCats.size(); idc++) {
            LorentzVector iboson(gammaEvHandler.massiveGamma(dilCats[idc]));
            float zmass=iboson.mass();
            //float dphi2lep = gammaEvHandler.dphimassiveGamma(dilCats[idc]);
            float coslz_cs = gammaEvHandler.coslzmassiveGamma(dilCats[idc]);
            //bool passdphi2l(true);
            //bool passdphi2l(cos(dphi2lep)>0);

            //Float_t mt( METUtils::transverseMass(iboson,metP4,true) );
            double mt_massless = METUtils::transverseMass(iboson,metP4,false);
            //bool passMTcut(mt>220 && mt<1200);
            //passMTcut &= passdphi2l;
            float iweight=weight;

	    
            if(isGammaEvent && dilCats[idc] != "ll" && applyQtweights) {

                if(tag_subcat=="eq0jets") {
                    if(gamma.pt()>150)	iweight*=qtWeights[dilCats[idc]];   //fit method
                    else		iweight*=Gweights_[dilCats[idc]];   //ratio method
                }
                if(tag_subcat=="eq1jets") {
                    if(gamma.pt()>200)	iweight*=qtWeights[dilCats[idc]];   //fit method
                    else 		iweight*=Gweights_[dilCats[idc]];   //ratio method
                }
            }
	    

            //if(isGammaEvent && dilCats[idc] == "ll") {
            //    //iweight*=qtWeights["ee"]+qtWeights["mumu"]; //new method
            //    iweight*=Gweights_["ee"]+Gweights_["mumu"];   //new method
            //}



            for(size_t isc=0; isc<subcats.size(); isc++) {
                TString ctf=dilCats[idc]+subcats[isc];

                //Vgamma Contrl
                if(passVgammaCtrl) {
                    mon.fillHisto("met_met_VgammaCtrl",         ctf, metP4.pt(),iweight);
                }


                if(!passMultiplicityVetoes) continue;


                mon.fillHisto("eventflow_gamma", ctf,0,iweight);
                //for reweighting,
                mon.fillHisto("qt",ctf, gamma.pt(),iweight,true);
                mon.fillHisto("qt_dataDY",ctf, gamma.pt(),iweight,true);
                mon.fillHisto("qtType1_dataDY",ctf, gamma.pt(),iweight,true);
                mon.fillHisto("nvtx_dataDY", ctf, ev.nvtx, iweight);
                mon.fillHisto("pfmet_dataDY", ctf, metP4.pt(), iweight, true);
                mon.fillHisto("pfmetType1_dataDY", ctf, metP4.pt(), iweight, true);
                mon.fillHisto("DphiZMET_dataDY",ctf, dphiZllmet, iweight);
                mon.fillHisto("balancedif_dataDY",ctf, balanceDif, iweight);


                if(passDphiZMETcut) {
                    mon.fillHisto("eventflow_gamma", ctf,1,iweight);

                    if(passMETcut) {
                        mon.fillHisto("eventflow_gamma", ctf,2,iweight);

                        if(passBalanceCut) {
                            mon.fillHisto("eventflow_gamma", ctf,3,iweight);

                        } //passBalanceCut
                    }//end passMETcut
                }//passDphiZMETcut


                //##############################
                //### for shape analysis and optimization
                //##############################

                bool passBaseSelection( passKinematics && passMultiplicityVetoes && passBveto );
                for(unsigned int index=0; index<nOptims; index++) {

                    double minMET = optim_Cuts1_MET[index];
                    double minBalance = optim_Cuts1_Balance[index];
                    double minDphi = optim_Cuts1_DphiZMET[index];

                    bool passLocalMETcut(metP4.pt()>minMET);
                    bool passLocalBalanceCut=(balance>(1.-minBalance) && balance<(1.+minBalance));
                    bool passLocalDphiZMETcut(dphiZllmet>minDphi);

                    bool passOptimSelection(passBaseSelection && passLocalMETcut && passLocalBalanceCut && passLocalDphiZMETcut);

                    // fill shapes for limit setting
                    if( passOptimSelection ) {
                        mon.fillHisto(TString("mt_shapes"),ctf, index, mt_massless, iweight);
                    }
                }






            }//end for subcats
        }//end for event
    }

    //all done with the events file
    file->Close();
    cout << "Sample treated as MC? " << isMC << endl;

    //save histograms to file
    TString outUrl( outdir );
    gSystem->Exec("mkdir -p " + outUrl);
    outUrl += "/";
    outUrl += outFileUrl + ".root";
    printf("Results save in %s\n", outUrl.Data());

    TFile *ofile=TFile::Open(outUrl, "recreate");
    mon.Write();
    ofile->Close();


    if(outTxtFile)fclose(outTxtFile);
}


