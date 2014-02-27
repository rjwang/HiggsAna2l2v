#include <string>
#include <vector>
#include <iostream>
#include <list>
#include <iterator>
#include <algorithm>
#include <unordered_map>
#include <sstream>

#include "TROOT.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TChain.h"
#include "TObject.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TLegend.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TTree.h"
#include "TF1.h"
#include "TCutG.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TMultiGraph.h"
#include "TPaveText.h"
#include "THStack.h"
#include "Rtypes.h"
#include "TList.h"

#include "TLatex.h"
#include "TPad.h"


#include "../src/tdrstyle.C"
#include "../src/JSONWrapper.cc"

int cutIndex=-1;
string cutIndexStr="";
double iLumi = 2007;
double iEcm=8;
bool showChi2 = false;
bool showUnc=false;
double baseRelUnc=0.044;
bool noLog=false;
bool isSim=false;
bool isDataBlind=false;
bool do2D  = true;
bool do1D  = true;
bool doTex = true;
bool StoreInFile = true;
bool doPlot = true;
bool splitCanvas = true;
bool onlyCutIndex = false;
bool nosig = false;
string inDir   = "OUTNew/";
string jsonFile = "../../data/beauty-samples.json";
string outDir  = "Img/";
string plotExt = ".png";
string plotExt2 = ".pdf";
string outFile = "plotter.root";
string cutflowhisto = "all_cutflow";
string ctrlSample = "";
int rebin = 1; //RJ
float cutLine1, cutLine2;
float leftcutLine, rightcutLine;


struct stSampleInfo{ double PURescale_up; double PURescale_down; double initialNumberOfEvents;};
std::unordered_map<string, stSampleInfo> sampleInfoMap;
std::unordered_map<string, bool> FileExist;

struct NameAndType{
   std::string name;
   bool type; 
   bool isIndexPlot;
   NameAndType(std::string name_,  bool type_, bool isIndexPlot_){name = name_; type = type_; isIndexPlot = isIndexPlot_;}

   bool operator==(const NameAndType& a){ return a.name == name;}
   bool operator< (const NameAndType& a){ return a.name < name;}
 };



TObject* GetObjectFromPath(TDirectory* File, std::string Path, bool GetACopy=false)
{
   size_t pos = Path.find("/");
   if(pos < 256){
      std::string firstPart = Path.substr(0,pos);
      std::string endPart   = Path.substr(pos+1,Path.length());
      TDirectory* TMP = (TDirectory*)File->Get(firstPart.c_str());
//      if(TMP!=NULL)return GetObjectFromPath(TMP,endPart,GetACopy);
      if(TMP!=NULL){
         TObject* TMP2 =  GetObjectFromPath(TMP,endPart,GetACopy);
//         if(!TMP2)printf("BUG: %s\n",Path.c_str());
         return TMP2;
      }

//      printf("BUG: %s\n",Path.c_str());
      return NULL;
   }else{
      TObject* TMP = File->Get(Path.c_str());
//      if(!TMP)printf("BUG: %s\n",Path.c_str());

      if(GetACopy){
         return TMP->Clone();
      }else{
         return TMP;
      }
   }

}

void GetListOfObject(JSONWrapper::Object& Root, std::string RootDir, std::list<NameAndType>& histlist, TDirectory* dir=NULL, std::string parentPath=""){

  if(parentPath=="" && !dir)
    {
      int dataProcessed = 0;
      int signProcessed = 0;
      int bckgProcessed = 0; 

      std::vector<JSONWrapper::Object> Process = Root["proc"].daughters();
      //loop on all Proc
      for(size_t ip=0; ip<Process.size(); ip++){
          if(Process[ip]["isinvisible"].toBool())continue;
	  bool isData (  Process[ip]["isdata"].toBool()  );
          bool isSign ( !isData &&  Process[ip].isTag("spimpose") && Process[ip]["spimpose"].toBool());
  	  bool isMC   = !isData && !isSign; 
	  string filtExt("");
	  if(Process[ip].isTag("mctruthmode") ) { char buf[255]; sprintf(buf,"_filt%d",(int)Process[ip]["mctruthmode"].toInt()); filtExt += buf; }

          //just to make it faster, only consider the first 3 sample of a same kind
          if(isData){if(dataProcessed>=1){continue;}else{dataProcessed++;}}
          if(isSign){if(signProcessed>=2){continue;}else{signProcessed++;}}
          if(isMC  ){if(bckgProcessed>8) {continue;}else{bckgProcessed++;}}

	  std::vector<JSONWrapper::Object> Samples = (Process[ip])["data"].daughters();
          //to make it faster only consider the first samples 
	  //for(size_t id=0; id<Samples.size()&&id<2; id++){
	  for(size_t id=0; id<Samples.size(); id++){
	      int split = 1;
	      if(Samples[id].isTag("split"))split = Samples[id]["split"].toInt();
	      string segmentExt;
	      if(split>1) { char buf[255]; sprintf(buf,"_%i",0); segmentExt += buf; }
              string FileName = RootDir + (Samples[id])["dtag"].toString() +  (Samples[id].isTag("suffix")?(Samples[id])["suffix"].toString():string("")) + segmentExt + filtExt + ".root";

              printf("Adding all objects from %25s to the list of considered objects\n",  FileName.c_str());
	      TFile* file = new TFile(FileName.c_str());
	      if(file->IsZombie())
		{
		  if(isData) dataProcessed--;
		  if(isSign) signProcessed--;
		  if(isMC) bckgProcessed--;
		}
	      GetListOfObject(Root,RootDir,histlist,(TDirectory*)file,"" );
	      file->Close();
	    }
	}

    //for(std::list<NameAndType>::iterator it= histlist.begin(); it!= histlist.end(); it++){printf("%s\n",it->name.c_str()); }

      return;
   }

   if(dir==NULL)return;

   TList* list = dir->GetListOfKeys();
   for(int i=0;i<list->GetSize();i++){
      TObject* tmp = GetObjectFromPath(dir,list->At(i)->GetName(),false);
      if(tmp->InheritsFrom("TTree")) continue;
//      if(tmp->InheritsFrom("TH1")){
//         histlist.push_back(NameAndType(parentPath+list->At(i)->GetName(), !(tmp->InheritsFrom("TH2") || tmp->InheritsFrom("TH3")) ) );
//      }else if(tmp->InheritsFrom("TDirectory")){
//         GetListOfObject(Root,RootDir,histlist,(TDirectory*)tmp,parentPath+ list->At(i)->GetName()+"/" );
//      }


      if(tmp->InheritsFrom("TDirectory")){
         GetListOfObject(Root,RootDir,histlist,(TDirectory*)tmp,parentPath+ list->At(i)->GetName()+"/" );
      }else if(tmp->InheritsFrom("TH1")){
        bool isTH1 = !(tmp->InheritsFrom("TH2") || tmp->InheritsFrom("TH3"));
        bool hasIndex = string(((TH1*)tmp)->GetXaxis()->GetTitle()).find("cut index")<string::npos;
        if(hasIndex){isTH1=true;}
	histlist.push_back(NameAndType(parentPath+list->At(i)->GetName(), isTH1, hasIndex ) );
      }

      delete tmp;
   }
   

}

void GetInitialNumberOfEvents(JSONWrapper::Object& Root, std::string RootDir, NameAndType HistoProperties){
   std::vector<JSONWrapper::Object> Process = Root["proc"].daughters();
   for(unsigned int i=0;i<Process.size();i++){
      if(!StoreInFile && Process[i]["isinvisible"].toBool())continue;
      string filtExt("");
      if(Process[i].isTag("mctruthmode") ) { char buf[255]; sprintf(buf,"_filt%d",(int)Process[i]["mctruthmode"].toInt()); filtExt += buf; }
      std::vector<JSONWrapper::Object> Samples = (Process[i])["data"].daughters();
      for(unsigned int j=0;j<Samples.size();j++){
	 int split = 1;
	 if(Samples[j].isTag("split"))split = Samples[j]["split"].toInt();

         TH1* tmphist = NULL;
         for(int s=0;s<split;s++){
	   string segmentExt; if(split>1) { char buf[255]; sprintf(buf,"_%i",s); segmentExt += buf; }
            string FileName = RootDir + (Samples[j])["dtag"].toString() + ((Samples[j].isTag("suffix"))?(Samples[j])["suffix"].toString():string("")) + segmentExt + filtExt + ".root";
            TFile* File = new TFile(FileName.c_str());
            bool& fileExist = FileExist[FileName];
            if(!File || File->IsZombie() || !File->IsOpen() || File->TestBit(TFile::kRecovered) ){fileExist=false;  continue; }else{fileExist=true;}

            TH1* tmptmphist = (TH1*) GetObjectFromPath(File,HistoProperties.name); 
	    if(tmptmphist)
	      {
		if(!tmphist){gROOT->cd(); tmphist = (TH1*)tmptmphist->Clone(tmptmphist->GetName());}else{tmphist->Add(tmptmphist);}
		delete tmptmphist;
	      }
            delete File;
         }
         if(!tmphist)continue;
         
	 bool isMC( !Process[i]["isdata"].toBool() && !Process[i]["isdatadriven"].toBool() );

         stSampleInfo& sampleInfo = sampleInfoMap[(Samples[j])["dtag"].toString()];

         double PUCentralnnorm =  1; if(tmphist->GetBinContent(3)>0)PUCentralnnorm = tmphist->GetBinContent(2) / tmphist->GetBinContent(3);
         double PUDownnorm     =  1; if(tmphist->GetBinContent(4)>0)PUDownnorm     = tmphist->GetBinContent(3) / tmphist->GetBinContent(4);
         double PUUpnorm       =  1; if(tmphist->GetBinContent(5)>0)PUUpnorm       = tmphist->GetBinContent(3) / tmphist->GetBinContent(5);
         sampleInfo.PURescale_down = PUDownnorm;
         sampleInfo.PURescale_up   = PUUpnorm;
         if(isMC)printf("PU Renormalization %25s Shift Down --> %6.2f  Central = %6.2f  Up Down --> %6.2f\n",(Samples[j])["dtag"].toString().c_str(),PUDownnorm, PUCentralnnorm, PUUpnorm);	


         double cnorm = 1.0;
         if(tmphist)cnorm = tmphist->GetBinContent(1);
         if(cnorm<=0 || !isMC)cnorm = 1.0;         
         if(cnorm==1 && isMC)printf("is there a problem with %s ? cnorm = %f\n",(Samples[j])["dtag"].toString().c_str(), cnorm);          
         if(!isMC)PUCentralnnorm = 1;

          double VBFMCRescale = tmphist->GetXaxis()->GetNbins()>5 ? tmphist->GetBinContent(6) / tmphist->GetBinContent(2) : 1.0;
	  if(VBFMCRescale!=0)          cnorm *= VBFMCRescale;
          //printf("VBFMCRescale for sample %s is %f\n", (Samples[j])["dtag"].toString().c_str(), VBFMCRescale );
         sampleInfo.initialNumberOfEvents = cnorm / PUCentralnnorm;
         delete tmphist;
      }   
   }
}



void SavingToFile(JSONWrapper::Object& Root, std::string RootDir, NameAndType HistoProperties, TFile* OutputFile){
   std::vector<TObject*> ObjectToDelete;

   std::vector<JSONWrapper::Object> Process = Root["proc"].daughters();
   for(unsigned int i=0;i<Process.size();i++){
      TH1* hist = NULL;
      //TList *tree_list = new TList;;//RJ
      std::vector<JSONWrapper::Object> Samples = (Process[i])["data"].daughters();
      for(unsigned int j=0;j<Samples.size();j++){
         double Weight = 1.0;
         if(!Process[i]["isdata"].toBool() && !Process[i]["isdatadriven"].toBool())Weight*= iLumi;
         if(Samples[j].isTag("xsec")     )Weight*= Samples[j]["xsec"].toDouble();
	 string filtExt("");
	 if(Process[i].isTag("mctruthmode") ) { char buf[255]; sprintf(buf,"_filt%d",(int)Process[i]["mctruthmode"].toInt()); filtExt += buf; }	 

	 std::vector<JSONWrapper::Object> BR = Samples[j]["br"].daughters();for(unsigned int b=0;b<BR.size();b++){Weight*=BR[b].toDouble();}
         stSampleInfo& sampleInfo = sampleInfoMap[(Samples[j])["dtag"].toString()];
         Weight /= sampleInfo.initialNumberOfEvents;

         if(HistoProperties.name.find("puup"  )!=string::npos){Weight *= sampleInfo.PURescale_up;}
         if(HistoProperties.name.find("pudown")!=string::npos){Weight *= sampleInfo.PURescale_down;}

         if(HistoProperties.name.find("optim_cut")!=string::npos){Weight=1.0;}

         int split = 1;
         if(Samples[j].isTag("split"))split = Samples[j]["split"].toInt();
         TH1* tmphist = NULL;
         for(int s=0;s<split;s++){
	   string segmentExt; if(split>1){ char buf[255]; sprintf(buf,"_%i",s); segmentExt += buf;}
	   
            string FileName = RootDir + (Samples[j])["dtag"].toString() + ((Samples[j].isTag("suffix"))?(Samples[j])["suffix"].toString():string("")) +  segmentExt + filtExt + ".root";
            if(!FileExist[FileName])continue;
            TFile* File = new TFile(FileName.c_str());
            if(!File || File->IsZombie() || !File->IsOpen() || File->TestBit(TFile::kRecovered) )continue;
            TH1* tmptmphist = (TH1*) GetObjectFromPath(File,HistoProperties.name); 
            if(!tmptmphist){delete File;continue;}
            if(!tmphist){gROOT->cd(); tmphist = (TH1*)tmptmphist->Clone(tmptmphist->GetName());}else{tmphist->Add(tmptmphist);}

	    //cout << "TEST 1---------------"<<endl;
	    //RJ
	    //TTree *tmptmptree = (TTree*)File->Get("zhinv_tree");
	    //tmptmptree->SetName(Process[i]["tag"].c_str());
	    //tmptmptree->SetWeight(Weight);
	    //cout << "TEST 2---------------"<<endl;
	    //tree_list->Add(tmptmptree);
	    //cout << "TEST 3---------------"<<endl;
            //tmptmptree->Print();
	    //delete tmptmptree; //RJ

            delete tmptmphist;
            delete File;
         }
         if(!tmphist)continue;
         if(!hist){gROOT->cd(); hist = (TH1*)tmphist->Clone(tmphist->GetName());hist->Scale(Weight);}else{hist->Add(tmphist,Weight);}
         delete tmphist;
      }   
      if(!hist)continue;

      string dirName = Process[i]["tag"].c_str();while(dirName.find("/")!=std::string::npos)dirName.replace(dirName.find("/"),1,"-");
      OutputFile->cd();
      //cout << "TEST 4---------------" << endl;
      //string treeName = "tree_"+Process[i]["tag"].c_str();
      ////TTree *zhinv_tree  = TTree::MergeTrees(tree_list);
      //TTree *zhinv_tree     = new TTree(treeName,"");
      //cout << "TEST 4---------------1" << endl;
      ////zhinv_tree->SetName("tree_"+Process[i]["tag"].c_str());
      //zhinv_tree->Write();
      //cout << "TEST 5---------------" << endl;
      TDirectory* subdir = OutputFile->GetDirectory(dirName.c_str());
      if(!subdir || subdir==OutputFile) subdir = OutputFile->mkdir(dirName.c_str());
      subdir->cd();
      hist->Write();
      delete hist;
      //delete zhinv_tree;
   }
}


void SavingTreeToFile(JSONWrapper::Object& Root, std::string RootDir, TFile* OutputFile){
   std::vector<TObject*> ObjectToDelete;
   std::vector<JSONWrapper::Object> Process = Root["proc"].daughters();
   for(unsigned int i=0;i<Process.size();i++){
      TList *tree_list = new TList();//RJ
      std::vector<JSONWrapper::Object> Samples = (Process[i])["data"].daughters();
      for(unsigned int j=0;j<Samples.size();j++){
         double Weight = 1.0;
         if(!Process[i]["isdata"].toBool() && !Process[i]["isdatadriven"].toBool())Weight*= iLumi;
         if(Samples[j].isTag("xsec")     )Weight*= Samples[j]["xsec"].toDouble();
	 string filtExt("");
	 if(Process[i].isTag("mctruthmode") ) { char buf[255]; sprintf(buf,"_filt%d",(int)Process[i]["mctruthmode"].toInt()); filtExt += buf; }	 

	 std::vector<JSONWrapper::Object> BR = Samples[j]["br"].daughters();for(unsigned int b=0;b<BR.size();b++){Weight*=BR[b].toDouble();}
         stSampleInfo& sampleInfo = sampleInfoMap[(Samples[j])["dtag"].toString()];
         Weight /= sampleInfo.initialNumberOfEvents;

         int split = 1;
         if(Samples[j].isTag("split"))split = Samples[j]["split"].toInt();
         for(int s=0;s<split;s++){
	    string segmentExt; if(split>1){ char buf[255]; sprintf(buf,"_%i",s); segmentExt += buf;}
            string FileName = RootDir + (Samples[j])["dtag"].toString() + ((Samples[j].isTag("suffix"))?(Samples[j])["suffix"].toString():string("")) +  segmentExt + filtExt + ".root";
            if(!FileExist[FileName])continue;
            TFile* File = new TFile(FileName.c_str());
            if(!File || File->IsZombie() || !File->IsOpen() || File->TestBit(TFile::kRecovered) )continue;
	    
	    if(TTree *tmptmptree = (TTree*)File->Get("zhinv_tree")){
		cout << "Adding trees ... ... " << FileName.c_str() << endl;
	    	tmptmptree->SetName(Process[i]["tag"].c_str());
	    	tmptmptree->SetWeight(Weight);
	    	tree_list->Add(tmptmptree);
	    }
         } //split
      } //Sample   

      OutputFile->cd();
      TTree *zhinv_tree = TTree::MergeTrees(tree_list);
      zhinv_tree->Write();

      delete zhinv_tree;
      delete tree_list;
   }
}


void fixExtremities(TH1* h,bool addOverflow, bool addUnderflow)
{
  if(h==0) return;

  if(addUnderflow){
      double fbin  = h->GetBinContent(0) + h->GetBinContent(1);
      double fbine = sqrt(h->GetBinError(0)*h->GetBinError(0)
                          + h->GetBinError(1)*h->GetBinError(1));
      h->SetBinContent(1,fbin);
      h->SetBinError(1,fbine);
      h->SetBinContent(0,0);
      h->SetBinError(0,0);
    }
  
  if(addOverflow){  
      int nbins = h->GetNbinsX();
      double fbin  = h->GetBinContent(nbins) + h->GetBinContent(nbins+1);
      double fbine = sqrt(h->GetBinError(nbins)*h->GetBinError(nbins) 
                          + h->GetBinError(nbins+1)*h->GetBinError(nbins+1));
      h->SetBinContent(nbins,fbin);
      h->SetBinError(nbins,fbine);
      h->SetBinContent(nbins+1,0);
      h->SetBinError(nbins+1,0);
    }
}


void Draw2DHistogramSplitCanvas(JSONWrapper::Object& Root, std::string RootDir, NameAndType HistoProperties){
   if(HistoProperties.isIndexPlot && cutIndex<0)return;

   std::string SaveName = "";

   TPaveText* T = new TPaveText(0.2,0.995,0.84,0.95, "NDC");
   T->SetFillColor(0);
   T->SetFillStyle(0);  T->SetLineColor(0);
   T->SetTextAlign(32);
   char Buffer[1024]; 
   if(isSim)	sprintf(Buffer, "CMS simulation, #it{ZH #rightarrow l^{+}l^{-}+#slash{E}_{T}}, #sqrt{s}=%.1f TeV, L=%.1f fb^{-1}", iEcm, iLumi/1000);
   else		sprintf(Buffer, "CMS preliminary, #it{ZH #rightarrow l^{+}l^{-}+#slash{E}_{T}}, #sqrt{s}=%.1f TeV, L=%.1f fb^{-1}", iEcm, iLumi/1000);
   T->AddText(Buffer);

   std::vector<JSONWrapper::Object> Process = Root["proc"].daughters();
   std::vector<TObject*> ObjectToDelete;
   for(unsigned int i=0;i<Process.size();i++){
      if(Process[i]["isinvisible"].toBool())continue;
      string filtExt("");
      if(Process[i].isTag("mctruthmode") ) { char buf[255]; sprintf(buf,"_filt%d",(int)Process[i]["mctruthmode"].toInt()); filtExt += buf; }

      TCanvas* c1 = new TCanvas("c1","c1",500,500);
      c1->SetLogz(true);

      TH1* hist = NULL;
      std::vector<JSONWrapper::Object> Samples = (Process[i])["data"].daughters();
      for(unsigned int j=0;j<Samples.size();j++){
         double Weight = 1.0;
         if(!Process[i]["isdata"].toBool()  && !Process[i]["isdatadriven"].toBool())Weight*= iLumi;
         if(Samples[j].isTag("xsec")     )Weight*= Samples[j]["xsec"].toDouble();
         std::vector<JSONWrapper::Object> BR = Samples[j]["br"].daughters();for(unsigned int b=0;b<BR.size();b++){Weight*=BR[b].toDouble();}
         stSampleInfo& sampleInfo = sampleInfoMap[(Samples[j])["dtag"].toString()];
         Weight /= sampleInfo.initialNumberOfEvents;

         if(HistoProperties.name.find("puup"  )!=string::npos){Weight *= sampleInfo.PURescale_up;}
         if(HistoProperties.name.find("pudown")!=string::npos){Weight *= sampleInfo.PURescale_down;}

         int split = 1;
         if(Samples[j].isTag("split"))split = Samples[j]["split"].toInt();
         TH1* tmphist = NULL;
         for(int s=0;s<split;s++){
	   string segmentExt; if(split>1) { char buf[255]; sprintf(buf,"_%i",s); segmentExt += buf; }
	    
            string FileName = RootDir + (Samples[j])["dtag"].toString() + ((Samples[j].isTag("suffix"))?(Samples[j])["suffix"].toString():string("")) + segmentExt + filtExt + ".root";
            if(!FileExist[FileName])continue;
            TFile* File = new TFile(FileName.c_str());
            if(!File || File->IsZombie() || !File->IsOpen() || File->TestBit(TFile::kRecovered) )continue;
            TH1* tmptmphist = (TH1*) GetObjectFromPath(File,HistoProperties.name); 
            if(!tmptmphist){delete File;continue;}
            if(!tmphist){gROOT->cd(); tmphist = (TH1*)tmptmphist->Clone(tmptmphist->GetName());}else{tmphist->Add(tmptmphist);}
            //if(Process[i]["isdata"].toBool())printf("%s --> %f*%f(%f)\n",(Samples[j])["dtag"].toString().c_str(), tmptmphist->Integral(),Weight, initialNumberOfEvents[(Samples[j])["dtag"].toString()]);
            delete tmptmphist;
            delete File;
         }
         if(!tmphist)continue;
         if(!hist){gROOT->cd(); hist = (TH1*)tmphist->Clone(tmphist->GetName());hist->Scale(Weight);}else{hist->Add(tmphist,Weight);}
         delete tmphist;
      }   
      if(!hist)continue;

      SaveName = hist->GetName();
      ObjectToDelete.push_back(hist);
      hist->SetTitle("");
      hist->SetStats(kFALSE);

      gStyle->SetPalette(1);
      hist->Draw("COLZ");
  
      TPaveText* leg = new TPaveText(0.75,0.35,0.90,0.20, "NDC");
      leg->SetFillColor(0);
      leg->SetFillStyle(0);  leg->SetLineColor(0);
      leg->SetTextAlign(12);
      leg->AddText(Process[i]["tag"].c_str());
      //leg->Draw("same");
      ObjectToDelete.push_back(leg);
//      delete leg;
//      delete hist;

      T->Draw("same");


          //RJ adding ee or mumu channel number!
          TPaveText* T2 = new TPaveText(0.2+0.55,0.93,0.36+0.55,0.83, "NDC");
          T2->SetFillColor(0);
          T2->SetFillStyle(0);  T2->SetLineColor(0);
          T2->SetTextAlign(22);
          char Buffer2[1024];
          TString ssName = SaveName;
          if(SaveName.find("eeeq0jets") != string::npos)        sprintf(Buffer2,"#it{ee channel, 0 jets}");
          if(SaveName.find("eeeq1jets") != string::npos)        sprintf(Buffer2,"#it{ee channel, 1 jets}");
          if(SaveName.find("mumueq0jets") != string::npos)      sprintf(Buffer2,"#it{#mu#mu channel, 0 jets}");
          if(SaveName.find("mumueq1jets") != string::npos)      sprintf(Buffer2,"#it{#mu#mu channel, 1 jets}");
          if(SaveName.find("emueq0jets") != string::npos)       sprintf(Buffer2,"#it{e#mu channel, 0 jets}");
	  if(SaveName.find("emueq1jets") != string::npos)       sprintf(Buffer2,"#it{e#mu channel, 1 jets}");
          T2->AddText(Buffer2);

          if(SaveName.find("eeeq0jets") != string::npos || SaveName.find("eeeq1jets") != string::npos ||
                SaveName.find("mumueq0jets") != string::npos || SaveName.find("mumueq1jets") != string::npos ||
                SaveName.find("emueq0jets") != string::npos || SaveName.find("emueq1jets") != string::npos )
          {
                T2->Draw("same");
          }


      string SavePath = SaveName + "_" + (Process[i])["tag"].toString() + plotExt;
      while(SavePath.find("*")!=std::string::npos)SavePath.replace(SavePath.find("*"),1,"");
      while(SavePath.find("#")!=std::string::npos)SavePath.replace(SavePath.find("#"),1,"");
      while(SavePath.find("{")!=std::string::npos)SavePath.replace(SavePath.find("{"),1,"");
      while(SavePath.find("}")!=std::string::npos)SavePath.replace(SavePath.find("}"),1,"");
      while(SavePath.find("(")!=std::string::npos)SavePath.replace(SavePath.find("("),1,"");
      while(SavePath.find(")")!=std::string::npos)SavePath.replace(SavePath.find(")"),1,"");
      while(SavePath.find("^")!=std::string::npos)SavePath.replace(SavePath.find("^"),1,"");
      while(SavePath.find(" ")!=std::string::npos)SavePath.replace(SavePath.find(" "),1,"");
      while(SavePath.find("/")!=std::string::npos)SavePath.replace(SavePath.find("/"),1,"-");
      if(outDir.size()) SavePath = outDir +"/"+ SavePath;
      system(string(("rm -f ") + SavePath).c_str());
      c1->SaveAs(SavePath.c_str());


      SavePath = SaveName + "_" + (Process[i])["tag"].toString() + plotExt2;
      while(SavePath.find("*")!=std::string::npos)SavePath.replace(SavePath.find("*"),1,"");
      while(SavePath.find("#")!=std::string::npos)SavePath.replace(SavePath.find("#"),1,"");
      while(SavePath.find("{")!=std::string::npos)SavePath.replace(SavePath.find("{"),1,"");
      while(SavePath.find("}")!=std::string::npos)SavePath.replace(SavePath.find("}"),1,"");
      while(SavePath.find("(")!=std::string::npos)SavePath.replace(SavePath.find("("),1,"");
      while(SavePath.find(")")!=std::string::npos)SavePath.replace(SavePath.find(")"),1,"");
      while(SavePath.find("^")!=std::string::npos)SavePath.replace(SavePath.find("^"),1,"");
       while(SavePath.find(" ")!=std::string::npos)SavePath.replace(SavePath.find(" "),1,"");
      while(SavePath.find("/")!=std::string::npos)SavePath.replace(SavePath.find("/"),1,"-");
      if(outDir.size()) SavePath = outDir +"/"+ SavePath;
      system(string(("rm -f ") + SavePath).c_str());
      c1->SaveAs(SavePath.c_str());

      delete c1;
   }

   for(unsigned int d=0;d<ObjectToDelete.size();d++){delete ObjectToDelete[d];}ObjectToDelete.clear();
   delete T;
}


void Draw2DHistogram(JSONWrapper::Object& Root, std::string RootDir, NameAndType HistoProperties){
   if(HistoProperties.isIndexPlot && cutIndex<0)return;

   std::string SaveName = "";

   TPaveText* T = new TPaveText(0.40,0.995,0.85,0.945, "NDC");
   T->SetFillColor(0);
   T->SetFillStyle(0);  T->SetLineColor(0);
   T->SetTextAlign(32);
   char Buffer[1024]; sprintf(Buffer, "CMS preliminary, #sqrt{s}=%.1f TeV, #scale[0.5]{#int} L=%.1f fb^{-1}", iEcm, iLumi/1000);
   T->AddText(Buffer);

   std::vector<JSONWrapper::Object> Process = Root["proc"].daughters();
   int NSampleToDraw = 0;
   for(unsigned int i=0;i<Process.size();i++){
      if(Process[i]["isinvisible"].toBool())continue;
      NSampleToDraw++;
   }
   int CanvasX = 3;
   int CanvasY = ceil(NSampleToDraw/CanvasX);
   TCanvas* c1 = new TCanvas("c1","c1",CanvasX*350,CanvasY*350);
   c1->Divide(CanvasX,CanvasY,0,0);


   std::vector<TObject*> ObjectToDelete;
   for(unsigned int i=0;i<Process.size();i++){
      if(Process[i]["isinvisible"].toBool())continue;

      TVirtualPad* pad = c1->cd(i+1);
      pad->SetLogz(true);
      pad->SetTopMargin(0.0); pad->SetBottomMargin(0.10);  pad->SetRightMargin(0.20);

      TH1* hist = NULL;
      std::vector<JSONWrapper::Object> Samples = (Process[i])["data"].daughters();
      for(unsigned int j=0;j<Samples.size();j++){
         double Weight = 1.0;
         if(!Process[i]["isdata"].toBool()  && !Process[i]["isdatadriven"].toBool())Weight*= iLumi;
	  string filtExt("");
	  if(Process[i].isTag("mctruthmode") ) { char buf[255]; sprintf(buf,"_filt%d",(int)Process[i]["mctruthmode"].toInt()); filtExt += buf; }
         if(Samples[j].isTag("xsec")     )Weight*= Samples[j]["xsec"].toDouble();
         std::vector<JSONWrapper::Object> BR = Samples[j]["br"].daughters();for(unsigned int b=0;b<BR.size();b++){Weight*=BR[b].toDouble();}
         stSampleInfo& sampleInfo = sampleInfoMap[(Samples[j])["dtag"].toString()];
         Weight /= sampleInfo.initialNumberOfEvents;

         if(HistoProperties.name.find("puup"  )!=string::npos){Weight *= sampleInfo.PURescale_up;}
         if(HistoProperties.name.find("pudown")!=string::npos){Weight *= sampleInfo.PURescale_down;}         

         int split = 1;
         if(Samples[j].isTag("split"))split = Samples[j]["split"].toInt();
         TH1* tmphist = NULL;
         for(int s=0;s<split;s++){
	   string segmentExt; if(split>1) { char buf[255]; sprintf(buf,"_%i",s); segmentExt += buf; } 
	    
            string FileName = RootDir + (Samples[j])["dtag"].toString() + ((Samples[j].isTag("suffix"))?(Samples[j])["suffix"].toString():string("")) + segmentExt + filtExt + ".root";
            if(!FileExist[FileName])continue;
            TFile* File = new TFile(FileName.c_str());
            if(!File || File->IsZombie() || !File->IsOpen() || File->TestBit(TFile::kRecovered) )continue;
            TH1* tmptmphist = (TH1*) GetObjectFromPath(File,HistoProperties.name); 
            if(!tmptmphist){delete File;continue;}
            if(!tmphist){gROOT->cd(); tmphist = (TH1*)tmptmphist->Clone(tmptmphist->GetName());}else{tmphist->Add(tmptmphist);}
            //if(Process[i]["isdata"].toBool())printf("%s --> %f\n",(Samples[j])["dtag"].toString().c_str(), tmptmphist->Integral());
            delete tmptmphist;
            delete File;
         }
         if(!tmphist)continue;
         if(!hist){gROOT->cd(); hist = (TH1*)tmphist->Clone(tmphist->GetName());hist->Scale(Weight);}else{hist->Add(tmphist,Weight);}
         delete tmphist;
      }   
      if(!hist)continue;

      SaveName = hist->GetName();
      ObjectToDelete.push_back(hist);
      hist->SetTitle("");
      hist->SetStats(kFALSE);

      hist->Draw("COLZ");
  
      TPaveText* leg = new TPaveText(0.10,0.995,0.30,0.90, "NDC");
      leg->SetFillColor(0);
      leg->SetFillStyle(0);  leg->SetLineColor(0);
      leg->SetTextAlign(12);
      leg->AddText(Process[i]["tag"].c_str());
      leg->Draw("same");
      ObjectToDelete.push_back(leg);
//      delete leg;
//      delete hist;
   }
   c1->cd(0);
   //T->Draw("same");
   string SavePath = SaveName + plotExt;
   while(SavePath.find("*")!=std::string::npos)SavePath.replace(SavePath.find("*"),1,"");
   while(SavePath.find("#")!=std::string::npos)SavePath.replace(SavePath.find("#"),1,"");
   while(SavePath.find("{")!=std::string::npos)SavePath.replace(SavePath.find("{"),1,"");
   while(SavePath.find("}")!=std::string::npos)SavePath.replace(SavePath.find("}"),1,"");
   while(SavePath.find("(")!=std::string::npos)SavePath.replace(SavePath.find("("),1,"");
   while(SavePath.find(")")!=std::string::npos)SavePath.replace(SavePath.find(")"),1,"");
   while(SavePath.find("^")!=std::string::npos)SavePath.replace(SavePath.find("^"),1,"");
   while(SavePath.find("/")!=std::string::npos)SavePath.replace(SavePath.find("/"),1,"-");
   if(outDir.size()) SavePath = outDir +"/"+ SavePath; 
   system(string(("rm -f ") + SavePath).c_str());
   c1->SaveAs(SavePath.c_str());

   SavePath = SaveName + plotExt2;
   while(SavePath.find("*")!=std::string::npos)SavePath.replace(SavePath.find("*"),1,"");
   while(SavePath.find("#")!=std::string::npos)SavePath.replace(SavePath.find("#"),1,"");
   while(SavePath.find("{")!=std::string::npos)SavePath.replace(SavePath.find("{"),1,"");
   while(SavePath.find("}")!=std::string::npos)SavePath.replace(SavePath.find("}"),1,"");
   while(SavePath.find("(")!=std::string::npos)SavePath.replace(SavePath.find("("),1,"");
   while(SavePath.find(")")!=std::string::npos)SavePath.replace(SavePath.find(")"),1,"");
   while(SavePath.find("^")!=std::string::npos)SavePath.replace(SavePath.find("^"),1,"");
   while(SavePath.find("/")!=std::string::npos)SavePath.replace(SavePath.find("/"),1,"-");
   if(outDir.size()) SavePath = outDir +"/"+ SavePath; 
   system(string(("rm -f ") + SavePath).c_str());
   c1->SaveAs(SavePath.c_str());



   for(unsigned int d=0;d<ObjectToDelete.size();d++){delete ObjectToDelete[d];}ObjectToDelete.clear();
   delete c1;
   delete T;
}


void Draw1DHistogram(JSONWrapper::Object& Root, std::string RootDir, NameAndType HistoProperties){
   if(HistoProperties.isIndexPlot && cutIndex<0)return;

   //TCanvas* c1 = new TCanvas("c1","c1",800,600); //change from 800,800, RJ
   TCanvas* c1 = new TCanvas("c1","c1",700,700); 
   TPad* t1 = new TPad();
   if(isDataBlind) t1 = new TPad("t1","t1", 0.0, 0.0, 1.0, 1.0);
   else t1 = new TPad("t1","t1", 0.0, 0.3, 1.0, 1.0);

   t1->Draw();
   t1->cd();
   if(!isDataBlind) t1->SetBottomMargin(0); //RJ
   TString name_denRelUncH;
   if(!noLog) t1->SetLogy(true);
   float maximumFound(noLog);

   //TLegend* legA  = new TLegend(0.845,0.2,0.99,0.99, "NDC"); 
   TLegend* legA  = new TLegend();
   if(!isDataBlind) legA  = new TLegend(0.45,0.7,0.9,0.95, "NDC"); //RJ
   else legA = new TLegend(0.55,0.78,0.9,0.95, "NDC");

   legA->SetNColumns(3); //RJ
   //   TLegend* legA  = new TLegend(0.51,0.93,0.67,0.75, "NDC"); 
   // TLegend* legB  = new TLegend(0.67,0.93,0.83,0.75, "NDC");
   THStack* stack = new THStack("MC","MC");
   TH1 *     mc   = NULL;
   TH1 *     mcPlusRelUnc = NULL;
   TH1 *     mctotalUnc = NULL;
   std::vector<TH1 *> spimpose;
   std::vector<TString> spimposeOpts;
   TH1*     data = NULL;
   std::vector<TObject*> ObjectToDelete;
   std::string SaveName = "";
   std::vector<JSONWrapper::Object> Process = Root["proc"].daughters();
   std::vector<TH1*> fakeStack;
   for(unsigned int i=0;i<Process.size();i++){
      if(Process[i]["isinvisible"].toBool())continue;
      TH1* hist = NULL;
      std::vector<JSONWrapper::Object> Samples = (Process[i])["data"].daughters();
      for(unsigned int j=0;j<Samples.size();j++){
         double Weight = 1.0;
         if(!Process[i]["isdata"].toBool() && !Process[i]["isdatadriven"].toBool() )Weight*= iLumi;
	  string filtExt("");
	  if(Process[i].isTag("mctruthmode") ) { char buf[255]; sprintf(buf,"_filt%d",(int)Process[i]["mctruthmode"].toInt()); filtExt += buf; }
         if(Samples[j].isTag("xsec")     )Weight*= Samples[j]["xsec"].toDouble();
         std::vector<JSONWrapper::Object> BR = Samples[j]["br"].daughters();for(unsigned int b=0;b<BR.size();b++){Weight*=BR[b].toDouble();}
         stSampleInfo& sampleInfo = sampleInfoMap[(Samples[j])["dtag"].toString()];
         Weight /= sampleInfo.initialNumberOfEvents;
	 //	 cout << (Samples[j])["dtag"].toString() << " " << (iLumi/1000)*(1/Weight) <<  " /fb" << endl;
         if(HistoProperties.name.find("puup"  )!=string::npos){Weight *= sampleInfo.PURescale_up  ;}
         if(HistoProperties.name.find("pudown")!=string::npos){Weight *= sampleInfo.PURescale_down;}
         if(HistoProperties.name.find("optim_cut")!=string::npos){Weight=1.0;}

         int split = 1; 
         if(Samples[j].isTag("split"))split = Samples[j]["split"].toInt();
         TH1* tmphist = NULL;
         for(int s=0;s<split;s++){
	   string segmentExt; if(split>1) { char buf[255]; sprintf(buf,"_%i",s); segmentExt += buf; }

	    string FileName = RootDir + (Samples[j])["dtag"].toString() + ((Samples[j].isTag("suffix"))?(Samples[j])["suffix"].toString():string("")) + segmentExt + filtExt + ".root";
            if(!FileExist[FileName]){continue;}
            TFile* File = new TFile(FileName.c_str());
            if(!File || File->IsZombie() || !File->IsOpen() || File->TestBit(TFile::kRecovered) )continue;
            TH1* tmptmphist = NULL; 
            if(HistoProperties.isIndexPlot && cutIndex>=0){
               TH2* tmp2D = (TH2*) GetObjectFromPath(File,HistoProperties.name);
               if(tmp2D){tmptmphist = tmp2D->ProjectionY((string(tmp2D->GetName())+cutIndexStr).c_str(),cutIndex,cutIndex); delete tmp2D;}
//            }else if(HistoProperties.name.find("mt_shape")){
//               TH2* tmp2D = (TH2*) GetObjectFromPath(File,HistoProperties.name);
//               if(tmp2D){tmp2D->GetXaxis()->SetTitle("#events"); tmptmphist = tmp2D->ProjectionX((string(tmp2D->GetName())+cutIndexStr).c_str()); delete tmp2D;}
            }else if(!HistoProperties.isIndexPlot){
               tmptmphist = (TH1*) GetObjectFromPath(File,HistoProperties.name);
	    }
	    if(!tmptmphist){delete File;continue;}
            if(!tmphist){gROOT->cd(); tmphist = (TH1*)tmptmphist->Clone(tmptmphist->GetName());}else{tmphist->Add(tmptmphist);}
            //if(Process[i]["isdata"].toBool())printf("%s --> %f\n",(Samples[j])["dtag"].toString().c_str(), tmptmphist->Integral());
            delete tmptmphist;
            delete File;
         }
         if(!tmphist)continue;
         if(!hist){gROOT->cd(); hist = (TH1*)tmphist->Clone(tmphist->GetName());hist->Scale(Weight);}else{hist->Add(tmphist,Weight);}
         delete tmphist;
      }   
      if(!hist)continue;

      SaveName = hist->GetName();
      TString postfix(""); postfix+=i;
      hist->SetName(SaveName+postfix);
      hist->Rebin(rebin); //RJ
      if(Process[i].isTag("color" ) )hist->SetLineColor  ((int)Process[i][ "color"].toDouble()); else hist->SetLineColor  (1);
      if(Process[i].isTag("color" ) )hist->SetMarkerColor((int)Process[i][ "color"].toDouble()); else hist->SetMarkerColor(1);
      if(Process[i].isTag("color" ) )hist->SetFillColor  ((int)Process[i][ "color"].toDouble()); else hist->SetFillColor  (0);
      if(Process[i].isTag("lcolor") )hist->SetLineColor  ((int)Process[i]["lcolor"].toDouble());
      if(Process[i].isTag("mcolor") )hist->SetMarkerColor((int)Process[i]["mcolor"].toDouble()); 
      if(Process[i].isTag("fcolor") )hist->SetFillColor  ((int)Process[i]["fcolor"].toDouble()); 
      if(Process[i].isTag("lwidth") )hist->SetLineWidth  ((int)Process[i]["lwidth"].toDouble());// else hist->SetLineWidth  (1);
      if(Process[i].isTag("lstyle") )hist->SetLineStyle  ((int)Process[i]["lstyle"].toDouble());// else hist->SetLinStyle  (1);
      if(Process[i].isTag("fill"  ) )hist->SetFillColor  ((int)Process[i]["fill"  ].toDouble());
      if(Process[i].isTag("marker") )hist->SetMarkerStyle((int)Process[i]["marker"].toDouble());// else hist->SetMarkerStyle(1);
      
      fixExtremities(hist,true,true);
      hist->SetTitle("");
      hist->SetStats(kFALSE);
      hist->SetMinimum(2e-2);//5e-1);//2e-2);
      //hist->SetMaximum(1E6);
      hist->SetMaximum(hist->GetBinContent(hist->GetMaximumBin())*1.10);
      TString tSaveName = hist->GetName();
      if(tSaveName.Contains("phi_dy")) hist->SetMaximum(hist->GetBinContent(hist->GetMaximumBin())*15);
      if(tSaveName.Contains("nvtx")) hist->SetMaximum(hist->GetBinContent(hist->GetMaximumBin())*20);
      if(tSaveName.Contains("nleptons")) hist->SetMaximum(hist->GetBinContent(hist->GetMaximumBin())*15);
      if(tSaveName.Contains("npfjets")) hist->SetMaximum(hist->GetBinContent(hist->GetMaximumBin())*20);
      if(tSaveName.Contains("npfjetsbtags")) hist->SetMaximum(hist->GetBinContent(hist->GetMaximumBin())*15);
      

      ObjectToDelete.push_back(hist);
      if(Process[i].isTag("normto")) hist->Scale( Process[i]["normto"].toDouble()/hist->Integral() );
      if((!Process[i].isTag("spimpose") || !Process[i]["spimpose"].toBool()) && !Process[i]["isdata"].toBool()){
         //Add to Stack
	//cout << "------------" << endl;
	//cout << "stack: " << Process[i]["tag"].c_str() << endl;
	if(Process[i]["tag"].c_str() == ctrlSample) 
		stack->Add(hist, "HIST");               
	else 	fakeStack.push_back(hist); //RENJIE

         legA->AddEntry(hist, Process[i]["tag"].c_str(), "F");	 
         if(!mc){mc = (TH1D*)hist->Clone("mc");}else{mc->Add(hist);}
      }
      else if(Process[i].isTag("spimpose") && Process[i]["spimpose"].toBool())
	{
   	  //legB->AddEntry(hist, Process[i]["tag"].c_str(), "L");
	  if(!nosig) legA->AddEntry(hist, Process[i]["tag"].c_str(), Process[i]["isdata"].toBool() ? "P" : "L" );
	  spimposeOpts.push_back( Process[i]["isdata"].toBool() ? "e1" : "hist" );
	  spimpose.push_back(hist);
	  if(maximumFound<hist->GetMaximum()) maximumFound=hist->GetMaximum()*5.0;//1.1;
	}
      else{
	if(Process[i]["isdata"].toBool()){
	  if(!data){
	    data = hist; 
	    legA->AddEntry(hist, Process[i]["tag"].c_str(), "P"); 
	  }
	  else data->Add(hist);
	}
      }
   }
   if(mc &&  maximumFound<mc->GetMaximum()) maximumFound=mc->GetMaximum()*5.0;//1.1;
   if(data && maximumFound<data->GetMaximum()) maximumFound=data->GetMaximum()*5.0;//1.1;

   bool canvasIsFilled(false);

   //RENJIE
   for(unsigned int j=0; j<fakeStack.size(); j++){
	stack->Add(fakeStack[j], "HIST");
   }

   if(stack && stack->GetStack() && stack->GetStack()->GetEntriesFast()>0){
     stack->Draw("");
     
     stack->GetYaxis()->SetTitleOffset(0.9/*0.8*/);
     stack->GetYaxis()->SetLabelSize(0.06);
     stack->GetYaxis()->SetTitleSize(0.06);

     TH1 *hist=(TH1*)stack->GetStack()->At(0);
     if(stack->GetXaxis())
       {
	 if(isDataBlind) stack->GetXaxis()->SetTitle(hist->GetXaxis()->GetTitle());
	 //stack->GetXaxis()->SetTitleSize(0.); 
	 name_denRelUncH = hist->GetXaxis()->GetTitle(); //RJ

	 //Set YTitle, RJ
	 float binsize = hist->GetBinWidth(1);
	 std::ostringstream strs;
	 strs << binsize;
	 TString binSize = strs.str();
	 if(name_denRelUncH.Contains("multiplicity")) stack->GetYaxis()->SetTitle("Events");
	 else{
	 	if(name_denRelUncH.Contains("GeV") && 
			!name_denRelUncH.Contains(">") && !name_denRelUncH.Contains("<")) 
				stack->GetYaxis()->SetTitle("Events / "+binSize+" GeV");
	 	else 			stack->GetYaxis()->SetTitle("Events / "+binSize);
	 }

	 //stack->GetYaxis()->SetTitle(hist->GetYaxis()->GetTitle());
	 stack->SetMinimum(hist->GetMinimum());
	 stack->SetMaximum(maximumFound);
	 if(noLog)
	   {
	     stack->SetMaximum(maximumFound);
	     //stack->GetXaxis()->SetRangeUser(hist->GetMinimum(),maximumFound);
	   }
       }
     ObjectToDelete.push_back(stack);
     canvasIsFilled=true;

     if(showUnc && mc)
       {
	 mcPlusRelUnc = (TH1 *) mc->Clone("totalmcwithunc"); mcPlusRelUnc->SetDirectory(0); 
	 for(int ibin=1; ibin<=mcPlusRelUnc->GetXaxis()->GetNbins(); ibin++)
	   {
	     Double_t error=sqrt(pow(mcPlusRelUnc->GetBinError(ibin),2)+pow(mcPlusRelUnc->GetBinContent(ibin)*baseRelUnc,2));
	     mcPlusRelUnc->SetBinError(ibin,error);
	   }
	 mcPlusRelUnc->SetFillStyle(3427);
	 mcPlusRelUnc->SetFillColor(kGray+1);
	 mcPlusRelUnc->SetMarkerStyle(1);
       }
   }
   if(data){
       if(!isDataBlind) data->Draw(canvasIsFilled ? "E1 same" : "E1");
       canvasIsFilled=true;
   }
   for(size_t ip=0; ip<spimpose.size(); ip++){
     TString opt=spimposeOpts[ip];
     if(!nosig) spimpose[ip]->Draw(opt + (canvasIsFilled ? "same": "") );
     canvasIsFilled=true;
   }

   if(mc)
   {
	mctotalUnc=(TH1D *) mc->Clone("mcrelunc");
       	for(int xbin=1; xbin<=mctotalUnc->GetXaxis()->GetNbins(); xbin++)
          {
          	if(mctotalUnc->GetBinContent(xbin)==0) continue;
           	mctotalUnc->SetBinContent(xbin,mctotalUnc->GetBinContent(xbin));
           	mctotalUnc->SetBinError(xbin,mctotalUnc->GetBinError(xbin));
          }
        mctotalUnc->SetDirectory(0);

        TGraphErrors *mcgr=new TGraphErrors;
        for(int ibin=1; ibin<=mctotalUnc->GetXaxis()->GetNbins(); ibin++)
          {
        	mcgr->SetPoint(ibin-1,mctotalUnc->GetXaxis()->GetBinCenter(ibin),mctotalUnc->GetBinContent(ibin));
                mcgr->SetPointError(ibin-1,mctotalUnc->GetXaxis()->GetBinWidth(ibin)/2,mctotalUnc->GetBinError(ibin));
          }
       	mcgr->SetFillStyle(3254);
       	mcgr->SetFillColor(1);
        mcgr->SetMarkerStyle(1);
       	mcgr->Draw("2 same");
	legA->AddEntry(mcgr,"Stat. Unc.", "F");
   }

   //compare data and MC
   if(showChi2 && data && mc && data->Integral()>0 && mc->Integral()>0)
     {
       TPaveText *pave = new TPaveText(0.6,0.85,0.8,0.9,"NDC");
       pave->SetBorderSize(0);
       pave->SetFillStyle(0);
       pave->SetTextAlign(32);
       pave->SetTextFont(42);
       char buf[100];
       sprintf(buf,"#chi^{2}/ndof : %3.2f", data->Chi2Test(mc,"WWCHI2/NDF") );
       pave->AddText(buf);
       pave->Draw();
     }
   
   TPaveText* T = new TPaveText(0.1,0.995,0.84,0.95, "NDC");
   T->SetFillColor(0);
   T->SetFillStyle(0);  T->SetLineColor(0);
   T->SetTextAlign(22);
   char Buffer[1024]; 
//   if(isSim) sprintf(Buffer, "CMS simulation, ZH#rightarrow ll+invisible, #sqrt{s}=%.1f TeV, #scale[0.5]{#int} L=%.1f fb^{-1}", iEcm, iLumi/1000);
//   else      sprintf(Buffer, "CMS preliminary, ZH#rightarrow ll+invisible, #sqrt{s}=%.1f TeV, #scale[0.5]{#int} L=%.1f fb^{-1}", iEcm, iLumi/1000);
   if(isSim) sprintf(Buffer, "CMS simulation, #it{ZH #rightarrow l^{+}l^{-}+#slash{E}_{T}}, #sqrt{s}=%.1f TeV, L=%.1f fb^{-1}", iEcm, iLumi/1000);
   else      sprintf(Buffer, "CMS preliminary, #it{ZH #rightarrow l^{+}l^{-}+#slash{E}_{T}}, #sqrt{s}=%.1f TeV, L=%.1f fb^{-1}", iEcm, iLumi/1000);

   T->AddText(Buffer);
   T->Draw("same");
   T->SetBorderSize(0);


          //RJ adding ee or mumu channel number!
          TPaveText* T2 = new TPaveText(0.2,0.93,0.36,0.83, "NDC");
          T2->SetFillColor(0);
          T2->SetFillStyle(0);  T2->SetLineColor(0);
          T2->SetTextAlign(22);
          char Buffer2[1024];
          TString ssName = SaveName;
          if(SaveName.find("eeeq0jets") != string::npos)	sprintf(Buffer2,"#it{ee channel, 0 jets}");
	  if(SaveName.find("eeeq1jets") != string::npos)        sprintf(Buffer2,"#it{ee channel, 1 jets}");
          if(SaveName.find("mumueq0jets") != string::npos) 	sprintf(Buffer2,"#it{#mu#mu channel, 0 jets}");
	  if(SaveName.find("mumueq1jets") != string::npos)      sprintf(Buffer2,"#it{#mu#mu channel, 1 jets}");
	  if(SaveName.find("emu") != string::npos)     sprintf(Buffer2,"#it{e#mu channel}");
          T2->AddText(Buffer2);

          if(SaveName.find("eeeq0jets") != string::npos || SaveName.find("eeeq1jets") != string::npos || 
		SaveName.find("mumueq0jets") != string::npos || SaveName.find("mumueq1jets") != string::npos || 
		SaveName.find("emu") != string::npos) 
	  {
	  	T2->Draw("same");
	  }

   //RJ
   TGraph *leftcut;
   if(leftcutLine != 0){
   	const int n = 35;
	double x1[n],y1[n];
	for (int i=0;i<n;i++) {
		x1[i] = leftcutLine;
		double ratio = (double) i/n;
     		y1[i] = maximumFound*ratio;
   	}
   	leftcut = new TGraph(n,x1,y1);
   	leftcut->SetLineColor(6);
	leftcut->SetFillColor(6);
   	leftcut->SetLineWidth(603);
   	leftcut->SetFillStyle(3005);
	leftcut->Draw("same");
   }


   TGraph *rightcut;
   if(rightcutLine != 0){
   	const int n = 35;
	double x1[n],y1[n];
	for (int i=0;i<n;i++) {
		x1[i] = rightcutLine;
		double ratio = (double) i/n;
     		y1[i] = maximumFound*ratio;
   	}
   	rightcut = new TGraph(n,x1,y1);
   	rightcut->SetLineColor(6);
	rightcut->SetFillColor(6);
   	rightcut->SetLineWidth(-603);
   	rightcut->SetFillStyle(3005);
	rightcut->Draw("same");
   }

   //RJ
   TLine *llcut1;
   if (cutLine1 != 0){
	llcut1 = new TLine(cutLine1, 0, cutLine1, maximumFound);
	llcut1->SetLineWidth(2);
	llcut1->SetLineStyle(2);
	llcut1->SetLineColor(kBlack);
	llcut1->Draw("same");
  }

   TLine *llcut2;
   if (cutLine2 != 0){
        llcut2 = new TLine(cutLine2, 0, cutLine2, maximumFound);
        llcut2->SetLineWidth(2);
        llcut2->SetLineStyle(2);
        llcut2->SetLineColor(kBlack);
        llcut2->Draw("same");
  }

   
   
   legA->SetFillColor(0); legA->SetFillStyle(0); legA->SetLineColor(0);
   legA->SetHeader("");
   legA->Draw("same");
   legA->SetTextFont(42);
   //    legB->SetFillColor(0); legB->SetFillStyle(0); legB->SetLineColor(0);
   //    legB->SetHeader("");
   //    legB->Draw("same");
   //   legB->SetTextFont(42);
   
   //
   std::vector<TH1 *> compDists;
   if(data)                   compDists.push_back(data);
   else if(spimpose.size()>0) compDists=spimpose;
   if(mc && compDists.size() && !isDataBlind)
     {
       c1->cd();
       TPad* t2 = new TPad("t2","t2", 0.0, 0.0, 1.0, 0.3);
       t2->Draw();
       t2->cd();
       t2->SetGridy(true);
       t2->SetTopMargin(0);
       t2->SetBottomMargin(0.5);

       //mc stats
       TH1D *denRelUncH=0;
       if(mcPlusRelUnc) denRelUncH=(TH1D *) mcPlusRelUnc->Clone("mcrelunc");
       else             denRelUncH=(TH1D *) mc->Clone("mcrelunc");
       for(int xbin=1; xbin<=denRelUncH->GetXaxis()->GetNbins(); xbin++)
	 {
	   if(denRelUncH->GetBinContent(xbin)==0) continue;
	   Double_t err=denRelUncH->GetBinError(xbin)/denRelUncH->GetBinContent(xbin);
	   denRelUncH->SetBinContent(xbin,1);
	   denRelUncH->SetBinError(xbin,err);
	 }
       //TGraphErrors *denRelUnc=new TGraphErrors(denRelUncH);
       TGraphErrors *denRelUnc=new TGraphErrors(denRelUncH->GetXaxis()->GetNbins());
       int icutgRJ=0;
       double llmin=-999;
       double llmax=999;
       for(int ibin=1; ibin<=denRelUncH->GetXaxis()->GetNbins(); ibin++){
           	//double err=denRelUncH->GetBinError(ibin)/denRelUncH->GetBinContent(ibin);
		if(ibin == 1) llmin=denRelUncH->GetXaxis()->GetBinLowEdge(ibin);
		llmax=denRelUncH->GetXaxis()->GetBinUpEdge(ibin);

		denRelUnc->SetPoint(icutgRJ,denRelUncH->GetXaxis()->GetBinCenter(ibin),1);
		if(denRelUncH->GetBinContent(ibin) < 2e-2) continue;
		denRelUnc->SetPointError(icutgRJ,denRelUncH->GetXaxis()->GetBinWidth(ibin)/2.0, denRelUncH->GetBinError(ibin)/denRelUncH->GetBinContent(ibin));
          	icutgRJ++;
       } //RJ

       denRelUnc->SetLineColor(1);
       denRelUnc->SetFillStyle(3254);
       denRelUnc->SetFillColor(kRed);
       denRelUnc->SetMarkerColor(1);
       denRelUnc->SetMarkerStyle(1);
       denRelUncH->Reset("ICE");       
       denRelUncH->Draw();
       denRelUnc->Draw("2");

       TLine *llbase = new TLine(llmin,1.,llmax,1.);
       llbase->SetLineWidth(1);
       llbase->SetLineStyle(1);
       llbase->SetLineColor(kBlue);
       llbase->Draw("same");


       float yscale = (1.0-0.2)/(0.18-0);       
       denRelUncH->GetYaxis()->SetTitle("Data/MC");
       denRelUncH->SetMinimum(0.09);
       denRelUncH->SetMaximum(1.99);
       //denRelUncH->GetXaxis()->SetTitle("");
       denRelUncH->GetXaxis()->SetTitle(name_denRelUncH); //RJ
       //denRelUncH->SetMinimum(0);
       //denRelUncH->SetMaximum(data->GetBinContent(data->GetMaximumBin())*1.10);
       denRelUncH->GetXaxis()->SetTitleOffset(1.2/*1.3*/);
       denRelUncH->GetXaxis()->SetLabelSize(0.03*yscale);
       denRelUncH->GetXaxis()->SetTitleSize(/*0.03*/0.04*yscale);
       denRelUncH->GetXaxis()->SetTickLength(0.03*yscale);
       denRelUncH->GetYaxis()->SetTitleOffset(0.39/*0.33*/);
       denRelUncH->GetYaxis()->SetNdivisions(5);
       denRelUncH->GetYaxis()->SetLabelSize(0.03*yscale);
       denRelUncH->GetYaxis()->SetTitleSize(0.03*yscale);
       
       //add comparisons
       for(size_t icd=0; icd<compDists.size(); icd++)
	 {
	   TString name("CompHistogram"); name+=icd;
	   TH1D *dataToObsH = (TH1D*)compDists[icd]->Clone(name);
	   dataToObsH->Divide(mc);
	   dataToObsH->Draw("E1 same");
	 }
     }
   else
     {
       //if not comparison resize the canvas
       c1->SetWindowSize(600,400);
       c1->SetCanvasSize(600,400);
       t1->SetPad(0,0,1,1);
       t1->SetBottomMargin(0.1);
       stack->GetXaxis()->SetTitle(name_denRelUncH); //RJ
     }

   c1->Modified();
   c1->Update();
   c1->cd();
   string SavePath = SaveName + plotExt;
   while(SavePath.find("*")!=std::string::npos)SavePath.replace(SavePath.find("*"),1,"");
   while(SavePath.find("#")!=std::string::npos)SavePath.replace(SavePath.find("#"),1,"");
   while(SavePath.find("{")!=std::string::npos)SavePath.replace(SavePath.find("{"),1,"");
   while(SavePath.find("}")!=std::string::npos)SavePath.replace(SavePath.find("}"),1,"");
   while(SavePath.find("(")!=std::string::npos)SavePath.replace(SavePath.find("("),1,"");
   while(SavePath.find(")")!=std::string::npos)SavePath.replace(SavePath.find(")"),1,"");
   while(SavePath.find("^")!=std::string::npos)SavePath.replace(SavePath.find("^"),1,"");
   while(SavePath.find("/")!=std::string::npos)SavePath.replace(SavePath.find("/"),1,"-");
   if(outDir.size()) SavePath = outDir +"/"+ SavePath;
   system(string(("rm -f ") + SavePath).c_str());
   c1->SaveAs(SavePath.c_str());

   SavePath = SaveName + plotExt2;
   while(SavePath.find("*")!=std::string::npos)SavePath.replace(SavePath.find("*"),1,"");
   while(SavePath.find("#")!=std::string::npos)SavePath.replace(SavePath.find("#"),1,"");
   while(SavePath.find("{")!=std::string::npos)SavePath.replace(SavePath.find("{"),1,"");
   while(SavePath.find("}")!=std::string::npos)SavePath.replace(SavePath.find("}"),1,"");
   while(SavePath.find("(")!=std::string::npos)SavePath.replace(SavePath.find("("),1,"");
   while(SavePath.find(")")!=std::string::npos)SavePath.replace(SavePath.find(")"),1,"");
   while(SavePath.find("^")!=std::string::npos)SavePath.replace(SavePath.find("^"),1,"");
   while(SavePath.find("/")!=std::string::npos)SavePath.replace(SavePath.find("/"),1,"-");
   if(outDir.size()) SavePath = outDir +"/"+ SavePath;
   system(string(("rm -f ") + SavePath).c_str());
   c1->SaveAs(SavePath.c_str());





   delete c1;
   for(unsigned int d=0;d<ObjectToDelete.size();d++){delete ObjectToDelete[d];}ObjectToDelete.clear();
   delete legA;
   //   delete legB;
   delete T;
}

std::string toLatexRounded(double value, double error){
   if(value==0.0 && error==0.0)return string("");

   double power = floor(log10(value));
   if(power<=-3)     {power=power+3;
   }else if(power>=2){power=power-2;
   }else             {power=0;}       

   value = value / pow(10,power);
   error = error / pow(10,power);
   int ValueFloating = 1+std::max(-1*log10(error),0.0);
   int ErrorFloating = ValueFloating;

   char tmpchar[255];
   if(power!=0){
      sprintf(tmpchar,"$(%.*f\\pm%.*f)\\times 10^{%g}$",ValueFloating,value,ErrorFloating,error,power);
   }else{
      sprintf(tmpchar,"$%.*f\\pm%.*f$",ValueFloating,value,ErrorFloating,error);
   }
   return string(tmpchar);
}



void ConvertToTex(JSONWrapper::Object& Root, std::string RootDir, NameAndType HistoProperties){
   if(HistoProperties.isIndexPlot && cutIndex<0)return;

   FILE* pFile = NULL;

   std::vector<TObject*> ObjectToDelete;
   TH1* stack = NULL; 
   std::vector<JSONWrapper::Object> Process = Root["proc"].daughters();
   for(unsigned int i=0;i<Process.size();i++){
      TH1* hist = NULL;
      std::vector<JSONWrapper::Object> Samples = (Process[i])["data"].daughters();




      for(unsigned int j=0;j<Samples.size();j++){
         double Weight = 1.0;
         if(!Process[i]["isdata"].toBool()  && !Process[i]["isdatadriven"].toBool())Weight*= iLumi;
	  string filtExt("");
	  if(Process[i].isTag("mctruthmode") ) { char buf[255]; sprintf(buf,"_filt%d",(int)Process[i]["mctruthmode"].toInt()); filtExt += buf; }
         if(Samples[j].isTag("xsec")     )Weight*= Samples[j]["xsec"].toDouble();
         std::vector<JSONWrapper::Object> BR = Samples[j]["br"].daughters();for(unsigned int b=0;b<BR.size();b++){Weight*=BR[b].toDouble();}
         stSampleInfo& sampleInfo = sampleInfoMap[(Samples[j])["dtag"].toString()];
         Weight /= sampleInfo.initialNumberOfEvents;

         if(HistoProperties.name.find("puup"  )!=string::npos){Weight *= sampleInfo.PURescale_up  ;}
         if(HistoProperties.name.find("pudown")!=string::npos){Weight *= sampleInfo.PURescale_down;}


         int split = 1;
         if(Samples[j].isTag("split"))split = Samples[j]["split"].toInt();
         TH1* tmphist = NULL;
         for(int s=0;s<split;s++){
	   string segmentExt; if(split>1) { char buf[255]; sprintf(buf,"_%i",s); segmentExt += buf; }

            string FileName = RootDir + (Samples[j])["dtag"].toString() + ((Samples[j].isTag("suffix"))?(Samples[j])["suffix"].toString():string("")) + segmentExt + filtExt + ".root";
            if(!FileExist[FileName])continue;
            TFile* File = new TFile(FileName.c_str());
            if(!File || File->IsZombie() || !File->IsOpen() || File->TestBit(TFile::kRecovered) )continue;
            TH1* tmptmphist = NULL;
            if(HistoProperties.isIndexPlot && cutIndex>=0){
               TH2* tmp2D = (TH2*) GetObjectFromPath(File,HistoProperties.name);
               if(tmp2D){tmptmphist = tmp2D->ProjectionY("_py",cutIndex,cutIndex); delete tmp2D;}
            }else if(!HistoProperties.isIndexPlot){
               tmptmphist = (TH1*) GetObjectFromPath(File,HistoProperties.name);
            }
            if(!tmptmphist){delete File;continue;}
            if(!tmphist){gROOT->cd(); tmphist = (TH1*)tmptmphist->Clone(tmptmphist->GetName());}else{tmphist->Add(tmptmphist);}
            delete tmptmphist;
            delete File;
         }
         if(!tmphist)continue;
         if(!hist){gROOT->cd(); hist = (TH1*)tmphist->Clone(tmphist->GetName());hist->Scale(Weight);}else{hist->Add(tmphist,Weight);}
         delete tmphist;
      }
      if(!hist)continue;

      if(!pFile){
         string SavePath = string(hist->GetName()) + ".tex";
         while(SavePath.find("*")!=std::string::npos)SavePath.replace(SavePath.find("*"),1,"");
         while(SavePath.find("#")!=std::string::npos)SavePath.replace(SavePath.find("#"),1,"");
         while(SavePath.find("{")!=std::string::npos)SavePath.replace(SavePath.find("{"),1,"");
         while(SavePath.find("}")!=std::string::npos)SavePath.replace(SavePath.find("}"),1,"");
         while(SavePath.find("(")!=std::string::npos)SavePath.replace(SavePath.find("("),1,"");
         while(SavePath.find(")")!=std::string::npos)SavePath.replace(SavePath.find(")"),1,"");
         while(SavePath.find("^")!=std::string::npos)SavePath.replace(SavePath.find("^"),1,"");
         while(SavePath.find("/")!=std::string::npos)SavePath.replace(SavePath.find("/"),1,"-");
         SavePath = outDir + SavePath;
         system(string(("rm -f ") + SavePath).c_str());
         pFile = fopen(SavePath.c_str(), "w");

         fprintf(pFile, "\\begin{table}[htp]\n");
         fprintf(pFile, "\\begin{center}\n");
         fprintf(pFile, "\\caption{}\n");
         fprintf(pFile, "\\label{tab:table}\n");

         string colfmt = "|l";  for(int b=1;b<=hist->GetXaxis()->GetNbins();b++){colfmt = colfmt + "c";} colfmt+="|";
         string colname = "";   for(int b=1;b<=hist->GetXaxis()->GetNbins();b++){
            std::string tempcolname =  hist->GetXaxis()->GetBinLabel(b);
            if(tempcolname.find("_")!=std::string::npos || tempcolname.find("^")!=std::string::npos)tempcolname = string("$") + tempcolname + "$";
            colname = colname + "& " + tempcolname;
         }
         fprintf(pFile, "\\begin{tabular}{ %s } \\hline\n", colfmt.c_str());
         fprintf(pFile, "Process %s \\\\ \\hline\\hline\n", colname.c_str());
      }
      ObjectToDelete.push_back(hist);

      std::string CleanTag = Process[i]["tag"].c_str();
      if(CleanTag.find("#")!=std::string::npos)CleanTag = string("$") + CleanTag + "$";
      while(CleanTag.find("#")!=std::string::npos)CleanTag.replace(CleanTag.find("#"),1,"\\");

      if((!Process[i].isTag("spimpose") || !Process[i]["spimpose"].toBool()) && !Process[i]["isdata"].toBool()){
         //Add to Stack
         if(!stack){stack = (TH1*)hist->Clone("Stack");}else{stack->Add(hist,1.0);}

         char numberastext[2048]; numberastext[0] = '\0';
         for(int b=1;b<=hist->GetXaxis()->GetNbins();b++){sprintf(numberastext,"%s & %s",numberastext, toLatexRounded(hist->GetBinContent(b), hist->GetBinError(b)).c_str());}
         fprintf(pFile, "%s %s \\\\\n",CleanTag.c_str(), numberastext);         
       }else{
          //Add to Canvas   
          if(stack){
            char numberastext[2048]; numberastext[0] = '\0';
            for(int b=1;b<=hist->GetXaxis()->GetNbins();b++){sprintf(numberastext,"%s & %s",numberastext, toLatexRounded(stack->GetBinContent(b), stack->GetBinError(b)).c_str());}
            fprintf(pFile, "%s %s \\\\\n\\hline\n","Total expected", numberastext);
             ObjectToDelete.push_back(stack);
             stack=NULL;
          }

          if(Process[i]["isdata"].toBool()){fprintf(pFile,"\\hline\n");}

          char numberastext[2048]; numberastext[0] = '\0';
          for(int b=1;b<=hist->GetXaxis()->GetNbins();b++){sprintf(numberastext,"%s & %s",numberastext, toLatexRounded(hist->GetBinContent(b), hist->GetBinError(b)).c_str());}
          fprintf(pFile, "%s %s \\\\\n",CleanTag.c_str(), numberastext);
       }

   }

   if(pFile){
      fprintf(pFile,"\\hline\n");
      fprintf(pFile,"\\end{tabular}\n");
      fprintf(pFile,"\\end{center}\n");
      fprintf(pFile,"\\end{table}\n");
      fclose(pFile);
   }
   for(unsigned int d=0;d<ObjectToDelete.size();d++){delete ObjectToDelete[d];}ObjectToDelete.clear();
}


int main(int argc, char* argv[]){
   setTDRStyle();  
   gStyle->SetPadTopMargin   (0.06);
   gStyle->SetPadBottomMargin(0.12);
   //gStyle->SetPadRightMargin (0.16);
   gStyle->SetPadRightMargin (0.06);//RJ
   gStyle->SetPadLeftMargin  (0.14);
   gStyle->SetTitleSize(0.04, "XYZ");
   gStyle->SetTitleXOffset(1.1);
   gStyle->SetTitleYOffset(1.45);
   gStyle->SetPalette(1);
   gStyle->SetNdivisions(505);

   std::vector<string> histoNameMask;
   std::vector<string> histoNameMaskStart;
   std::vector<string> nohistoNameMask;

   for(int i=1;i<argc;i++){
     string arg(argv[i]);
     //printf("--- %i - %s\n",i,argv[i]);

     if(arg.find("--help")!=string::npos){
        printf("--help   --> print this helping text\n");

        printf("--iLumi   --> integrated luminosity to be used for the MC rescale\n");
        printf("--iEcm    --> center of mass energy (TeV) = 8 TeV by default\n");
        printf("--isSim   --> print CMS Simulation instead of the standard title\n");
        printf("--inDir   --> path to the directory containing the .root files to process\n");
        printf("--outDir  --> path of the directory that will contains the output plots and tables\n");
        printf("--outFile --> path of the output summary .root file\n");
        printf("--json    --> containing list of process (and associated style) to process to process\n");
        printf("--only    --> processing only the objects containing the following argument in their name\n");
        printf("--OnlyStartWith  --> processing only the objects starting with the following argument in their name\n");
	printf("--NoPlotwith  --> processing only the objects starting without the following argument in their name\n");
        printf("--index   --> will do the projection on that index for histos of type cutIndex\n");
        printf("--chi2    --> show the data/MC chi^2\n"); 
        printf("--showUnc --> show stat uncertainty (if number is given use it as relative bin by bin uncertainty (e.g. lumi)\n"); 
	printf("--noLog   --> use linear scale\n");
        printf("--no1D   --> Skip processing of 1D objects\n");
        printf("--no2D   --> Skip processing of 2D objects\n");
        printf("--noTex  --> Do not create latex table (when possible)\n");
        printf("--noRoot --> Do not make a summary .root file\n");
        printf("--noPlot --> Do not creates plot files (useful to speedup processing)\n");
	printf("--plotExt --> extension to save\n");
	printf("--cutflow --> name of the histogram with the original number of events (cutflow by default)\n");
        printf("--splitCanvas --> (only for 2D plots) save all the samples in separated pltos\n");
	printf("--rebin	 --> number of bin (default = 1) for rebinning histograms, combine --only option to rebin some specific histograms\n"); //RJ 
	printf("--isDataBlind --> Blind the Data point from 1D Histograms\n"); //RJ
	printf("--addcutLine1 --> Add cut Line 1\n"); //RJ
	printf("--addcutLine2 --> Add cut Line 2\n"); //RJ
	printf("--addleftcutLine --> Add left cut Line \n"); //RJ
	printf("--addrightcutLine --> Add right cut Line \n"); //RJ

        printf("command line example: runPlotter --json ../data/beauty-samples.json --iLumi 2007 --inDir OUT/ --outDir OUT/plots/ --outFile plotter.root --noRoot --noPlot\n");
	return 0;
     }

     if(arg.find("--iLumi"  )!=string::npos && i+1<argc){ sscanf(argv[i+1],"%lf",&iLumi); i++; printf("Lumi = %f\n", iLumi); }
     if(arg.find("--iEcm"  )!=string::npos && i+1<argc){ sscanf(argv[i+1],"%lf",&iEcm); i++; printf("Ecm = %f TeV\n", iEcm); }

     if(arg.find("--inDir"  )!=string::npos && i+1<argc){ inDir    = argv[i+1];  i++;  printf("inDir = %s\n", inDir.c_str());  }
     if(arg.find("--outDir" )!=string::npos && i+1<argc){ outDir   = argv[i+1];  i++;  printf("outDir = %s\n", outDir.c_str());  }
     if(arg.find("--outFile")!=string::npos && i+1<argc){ outFile  = argv[i+1];  i++; printf("output file = %s\n", outFile.c_str()); }
     if(arg.find("--json"   )!=string::npos && i+1<argc){ jsonFile = argv[i+1];  i++;  }
     if(arg.find("--OnlyStartWith"   )!=string::npos && i+1<argc){ histoNameMaskStart.push_back(argv[i+1]); printf("Only process Histo starting with '%s'\n", argv[i+1]); i++;  }
     if(arg.find("--only"   )!=string::npos && i+1<argc)         { histoNameMask.push_back(argv[i+1]); printf("Only process Histo containing '%s'\n", argv[i+1]); i++;  }
     if(arg.find("--NoPlotwith")!=string::npos && i+1<argc)	{ nohistoNameMask.push_back(argv[i+1]); printf("Only process Histo not containing '%s'\n", argv[i+1]); i++;  }
     if(arg.find("--index"  )!=string::npos && i+1<argc)         { sscanf(argv[i+1],"%d",&cutIndex); i++; onlyCutIndex=(cutIndex>=0); printf("index = %i\n", cutIndex);  }
     if(arg.find("--chi2"  )!=string::npos)                      { showChi2 = true;  }
     if(arg.find("--showUnc") != string::npos) { 
       showUnc=true; 
       if(i+1<argc) { 
	 string nextArg(argv[i+1]); 
	 if(nextArg.find("--")==string::npos)
	   {
	     sscanf(argv[i+1],"%lf",&baseRelUnc);
	     i++;
	   }
       }
       printf("Uncertainty band will be included for MC with base relative uncertainty of: %3.2f",baseRelUnc);
     }
     if(arg.find("--isSim")!=string::npos){ isSim = true;    }
     if(arg.find("--isDataBlind")!=string::npos){ isDataBlind = true; printf("isDataBlind\n");} //RJ
     if(arg.find("--nosig")!=string::npos){nosig = true; printf("noSig plot=1\n");}//RJ
     if(arg.find("--noLog")!=string::npos){ noLog = true;    }
     if(arg.find("--no2D"  )!=string::npos){ do2D = false;    }
     if(arg.find("--no1D"  )!=string::npos){ do1D = false;    }
     if(arg.find("--noTex" )!=string::npos){ doTex= false;    }
     if(arg.find("--noRoot")!=string::npos){ StoreInFile = false;    }
     if(arg.find("--noPlot")!=string::npos){ doPlot = false;    }
     if(arg.find("--plotExt" )!=string::npos && i+1<argc){ plotExt   = argv[i+1];  i++;  printf("saving plots as = %s\n", plotExt.c_str());  }
     if(arg.find("--cutflow" )!=string::npos && i+1<argc){ cutflowhisto   = argv[i+1];  i++;  printf("Normalizing from 1st bin in = %s\n", cutflowhisto.c_str());  }
     if(arg.find("--splitCanvas")!=string::npos){ splitCanvas = true;    }
     if(arg.find("--addcutLine1"  )!=string::npos && i+1<argc) { sscanf(argv[i+1],"%f",&cutLine1); i++; printf("cutLine1 = %f\n", cutLine1); } //RJ
     if(arg.find("--addcutLine2"  )!=string::npos && i+1<argc) { sscanf(argv[i+1],"%f",&cutLine2); i++; printf("cutLine2 = %f\n", cutLine2); } //RJ
     if(arg.find("--addleftcutLine"  )!=string::npos && i+1<argc) { sscanf(argv[i+1],"%f",&leftcutLine); i++; printf("leftcutLine = %f\n", leftcutLine); } //RJ
     if(arg.find("--addrightcutLine"  )!=string::npos && i+1<argc) { sscanf(argv[i+1],"%f",&rightcutLine); i++; printf("rightcutLine = %f\n", rightcutLine); } //RJ
     if(arg.find("--ctrlSample" )!=string::npos && i+1<argc){ ctrlSample   = argv[i+1];  i++;  printf("Control Sample = %s\n", ctrlSample.c_str());  }
     if(arg.find("--rebin"  )!=string::npos && i+1<argc) { sscanf(argv[i+1],"%d",&rebin); i++; printf("rebin = %d\n", rebin); } //RJ

   } 
   system( (string("mkdir -p ") + outDir).c_str());
   if(ctrlSample.find("_")!=std::string::npos)ctrlSample.replace(ctrlSample.find("_"),1," ");

   char buf[255];
   sprintf(buf, "_Index%d", cutIndex);
   cutIndexStr = buf;

   JSONWrapper::Object Root(jsonFile, true);
   GetInitialNumberOfEvents(Root,inDir,NameAndType(cutflowhisto,true, false));  //Used to get the rescale factor based on the total number of events geenrated
   std::list<NameAndType> histlist;
   GetListOfObject(Root,inDir,histlist);
   histlist.sort();
   histlist.unique();   


   TFile* OutputFile = NULL;
   if(StoreInFile) OutputFile = new TFile(outFile.c_str(),"RECREATE");
   if(!doPlot){
   	printf("Progressing Bar              :0%%       20%%       40%%       60%%       80%%       100%%\n");
  	printf("                             :");
   }
   int TreeStep = histlist.size()/50;if(TreeStep==0)TreeStep=1;
   string csvFile(outDir +"/histlist.csv");
   system(("echo \"\" > " + csvFile).c_str());

   //SavingTreeToFile(Root,inDir, OutputFile); //RJ, runMVAtree

   int ictr(0);
   for(std::list<NameAndType>::iterator it= histlist.begin(); it!= histlist.end(); it++,ictr++)
     {
       if(!doPlot){ if(ictr%TreeStep==0){printf(".");fflush(stdout);} }
       bool passMasking(false);
       for(unsigned int i=0;i<histoNameMask.size();i++){

	if(histoNameMaskStart.size()==0){
		if(it->name.find(histoNameMask[i])!=std::string::npos)
		passMasking = true;
	} else {
		for(unsigned int j=0;j<histoNameMaskStart.size();j++){

			if(it->name.find(histoNameMask[i])!=std::string::npos && 
				it->name.find(histoNameMaskStart[j])!=std::string::npos) passMasking = true;
		}
	}
       }
       if(histoNameMask.size()==0 && histoNameMaskStart.size()==0)passMasking = true;

	//excluding some hists
       if(passMasking && nohistoNameMask.size()!=0) {
		for(unsigned int j=0;j<nohistoNameMask.size();j++)
		{
			if(it->name.find(nohistoNameMask[j])!=std::string::npos){
				passMasking = false;
				break;}
		}
       }


       if(!passMasking)continue;

       system(("echo \"" + it->name + "\" >> " + csvFile).c_str());
       //RJ, do not convert to Tex files
/*
       if(doTex && (it->name.find("eventflow")!=std::string::npos || it->name.find("evtflow")!=std::string::npos) 
       		&& it->name.find("optim_eventflow")==std::string::npos)
	{
		ConvertToTex(Root,inDir,*it); 
	}
*/

       if(doPlot && do2D  && !it->type){                      if(!splitCanvas){ Draw2DHistogram(Root,inDir,*it);}else{Draw2DHistogramSplitCanvas(Root,inDir,*it);}}
       if(doPlot && do1D  &&  it->type){                                        Draw1DHistogram(Root,inDir,*it); }	

       if(StoreInFile && do2D  && !it->type){                                  	SavingToFile(Root,inDir,*it, OutputFile); }
       if(StoreInFile && do1D  &&  it->type){					SavingToFile(Root,inDir,*it, OutputFile); }

     }


   printf("\n");
   if(StoreInFile) OutputFile->Close();
   
   system(("python ${CMSSW_BASE}/src/CMGTools/HiggsAna2l2v/data/html/generateJSONplotterFromList.py -i " + csvFile + " -o "+outDir+"/plotter.json").c_str());
   //system(("rm " + csvFile).c_str());
   system(("cp ${CMSSW_BASE}/src/CMGTools/HiggsAna2l2v/data/html/index.html " + outDir).c_str());
   printf("You can browse the results using %sindex.html\n",outDir.c_str());
}
