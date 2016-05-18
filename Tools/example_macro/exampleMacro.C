#include "TROOT.h"
#include "TRint.h"
#include "TStyle.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TEfficiency.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "THStack.h"
#include "TLegend.h"
#include "TTree.h"
#include "TBranch.h"
#include "TLorentzVector.h"
#include "Math/Point3D.h"
#include "TChain.h"
#include "Riostream.h"

#include "../src/MuonPogTree.h"
#include "../src/Utils.h"
#include "tdrstyle.C"

#include <cstdlib>
#include <iostream>
#include <algorithm>
#include <sstream>
#include <iostream>
#include <fstream> 
#include <vector>
#include <regex>
#include <map>


#include <boost/filesystem.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <boost/property_tree/exceptions.hpp>
#include <boost/lexical_cast.hpp>

// Helper classes defintion *****
// 1. SampleConfig : configuration class containing sample parameters
// 2. AlgoConfig   : configuration class containing macro logic parameters
// ******************************

namespace muon_pog {
 
  class SampleConfig {

  public :

    // config parameters (public for direct access)

    TString fileName;  
    TString sampleName;  
    Float_t cSection;
    Float_t nEvents;
    Bool_t applyReweighting;
        
    SampleConfig() {};
    
#ifndef __MAKECINT__ // avoid CINT vs boost problems   
    SampleConfig(boost::property_tree::ptree::value_type & vt); 
#endif

    ~SampleConfig() {};

  private:

  };

  class AlgoConfig {

  public :
    
    // config parameters (public for direct access)
    
    Float_t muon_minPt;      
    TString muon_ID;

    Float_t muon_PFIsoCut_Loose;
    Float_t muon_PFIsoCut_Tight;
    Float_t muon_PFNoPUIsoCut_Loose;
    Float_t muon_PFNoPUIsoCut_Tight;
    Float_t muon_TkrIsoCut_Loose;
    Float_t muon_TkrIsoCut_Tight;
    
    AlgoConfig() {};
    
#ifndef __MAKECINT__ // avoid CINT vs boost problems 
    AlgoConfig(boost::property_tree::ptree::value_type & vt); 
#endif

    ~AlgoConfig() {};
    
  private:

  };
  
}

// Helper classes defintion *****
// 1. parseConfig : parse the full cfg file (fill sample configs and algo one)
// ******************************

namespace muon_pog {
  void parseConfig(const std::string configFile, AlgoConfig & algoConfig,
		   std::vector<SampleConfig> & sampleConfigs);  
  void EfficiencyHisto(std::map<TString, std::map<TString, TH1F *> > *histos, TString samname, TString hisnameN, TString hisnameD);
  TGraphErrors *EfficiencyPlot(std::map<TString, std::map<TString, TH1F *> > *histos, std::vector<SampleConfig> sampleConfigs, TString hisname);
}



// The main program *************
// 1. Get configuration file and produces configurations
// 2. Create Plotters and loop on the event to fill them
// 3. Writes results in configurable outuput file
// ******************************

int main(int argc, char* argv[]){
  using namespace muon_pog;


  if (argc != 3) 
    {
      std::cout << "Usage : "
		<< argv[0] << " PATH_TO_CONFIG_FILE PATH_TO_OUTPUT_DIR\n";
      exit(100);
    }

  std::string configFile(argv[1]);
  
  std::cout << "[" << argv[0] << "] Using config file " << configFile << std::endl;

  // Output directory
  TString dirName = argv[2];
  system("mkdir -p " + dirName);
  TFile* outputFile = TFile::Open(dirName + "/results.root","RECREATE"); // CB find a better name for output file  

  // Set it to kTRUE if you do not run interactively
  gROOT->SetBatch(kTRUE); 

  // Initialize Root application
  TRint* app = new TRint("CMS Root Application", &argc, argv);

  setTDRStyle();
 
  AlgoConfig algoConfig;
  std::vector<SampleConfig> sampleConfigs;
  
  parseConfig(configFile,algoConfig,sampleConfigs);

  TString IsoHistos[3] = {"passPFIso","passPFNoPUIso","passTkrIso"};
  TString IsoWP[2]     = {"Loose","Tight"};


  TH1::SetDefaultSumw2(kTRUE);  
  // Kinematic binning
  Float_t pTBin[19]   = {15,20,25,30,35,40,45,50,55,60,70,80,90,100,120,140,160,180,200};
  Float_t etaBin[10]  = {0.0, 0.25, 0.50, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.4};
  Float_t nVtxBin[15] = {0,10,12,14,16,18,20,22,24,26,28,30,32,36,40};

  std::map<TString, std::map<TString, TH1F *> > histos;
  for (auto sampleConfig : sampleConfigs)
    {

      // Global histograms
      histos[sampleConfig.sampleName]["hEvtNVertex"] = new TH1F("hEvtNVertex_" + sampleConfig.sampleName, "hEvtNVertex_" + sampleConfig.sampleName + ";N. Vertex; entries", 30,0,30); 
      histos[sampleConfig.sampleName]["hNGENMuons"] = new TH1F("hNGENMuons_" + sampleConfig.sampleName, "hNGENMuons_" + sampleConfig.sampleName + ";N. of GEN Muons; entries", 5,0,5); 
      histos[sampleConfig.sampleName]["hNMuons"] = new TH1F("hNMuons_" + sampleConfig.sampleName, "hNMuons_" + sampleConfig.sampleName + ";N. of Muons; entries", 5,0,5); 
      // Histos for all muons
      histos[sampleConfig.sampleName]["hNVertex"] = new TH1F("hNVertex_" + sampleConfig.sampleName, "hNVertex_" + sampleConfig.sampleName + ";N. Vertex; entries", 14, nVtxBin); 
      histos[sampleConfig.sampleName]["hpT"] = new TH1F("hpT_" + sampleConfig.sampleName, "hpT_" + sampleConfig.sampleName + ";pT [GeV]; entries", 18, pTBin); 
      histos[sampleConfig.sampleName]["hEta"] = new TH1F("hEta_" + sampleConfig.sampleName, "hEta_" + sampleConfig.sampleName + ";#eta; entries", 9, etaBin); 
      // Histos for PF isolated muons
      for(int h=0; h<3; h++){
	for(int w=0; w<2; w++){
	  histos[sampleConfig.sampleName]["hNVertex_" + IsoHistos[h] + "_" + IsoWP[w]] = new TH1F("hNVertex_" + IsoHistos[h] + "_" + IsoWP[w] + "_" + sampleConfig.sampleName, "hNVertex_" + IsoHistos[h] + "_" + IsoWP[w] + "_" + sampleConfig.sampleName + ";N. Vertex; entries", 14, nVtxBin); 
	  histos[sampleConfig.sampleName]["hpT_" + IsoHistos[h] + "_" + IsoWP[w]] = new TH1F("hpT_" + IsoHistos[h] + "_" + IsoWP[w] + "_" + sampleConfig.sampleName, "hpT_" + IsoHistos[h] + "_" + IsoWP[w] + "_" + sampleConfig.sampleName + ";pT [GeV]; entries", 18, pTBin); 
	  histos[sampleConfig.sampleName]["hEta_" + IsoHistos[h] + "_" + IsoWP[w]] = new TH1F("hEta_" + IsoHistos[h] + "_" + IsoWP[w] + "_" + sampleConfig.sampleName, "hEta_" + IsoHistos[h] + "_" + IsoWP[w] + "_" + sampleConfig.sampleName + ";#eta; entries", 9, etaBin ); 
	  // Efficiencies
	  histos[sampleConfig.sampleName]["hNVertex_" + IsoHistos[h] + "_" + IsoWP[w] + "_Eff"] = new TH1F("hNVertex_" + IsoHistos[h] + "_" + IsoWP[w] + "_Eff" + "_" + sampleConfig.sampleName, "hNVertex_" + IsoHistos[h] + "_" + IsoWP[w] + "_Eff" + "_" + sampleConfig.sampleName + ";N. Vertex; Eff.", 14, nVtxBin); 
	  histos[sampleConfig.sampleName]["hpT_" + IsoHistos[h] + "_" + IsoWP[w] + "_Eff"] = new TH1F("hpT_" + IsoHistos[h] + "_" + IsoWP[w] + "_Eff" + "_" + sampleConfig.sampleName, "hpT_" + IsoHistos[h] + "_" + IsoWP[w] + "_Eff" + "_" + sampleConfig.sampleName + ";pT [GeV]; Eff.", 18, pTBin); 
	  histos[sampleConfig.sampleName]["hEta_" + IsoHistos[h] + "_" + IsoWP[w] + "_Eff"] = new TH1F("hEta_" + IsoHistos[h] + "_" + IsoWP[w] + "_Eff" + "_" + sampleConfig.sampleName, "hEta_" + IsoHistos[h] + "_" + IsoWP[w] + "_Eff" + "_" + sampleConfig.sampleName + ";#eta; Eff.", 9, etaBin ); 
	}
      }
      // Detector isolation
      histos[sampleConfig.sampleName]["hTrackerIso"] = new TH1F("hTrackerIso_" + sampleConfig.sampleName, "htrackerIso_" + sampleConfig.sampleName + ";Tracker Isolation; entries", 100,0,1); 
      histos[sampleConfig.sampleName]["hEMCalIso"] = new TH1F("hEMCalIso_" + sampleConfig.sampleName, "hEMCalIso_" + sampleConfig.sampleName + ";ECAL Isolation; entries", 50,0,10); 
      histos[sampleConfig.sampleName]["hHCalIso"] = new TH1F("hHCalIso_"   + sampleConfig.sampleName, "hHCalIso_" + sampleConfig.sampleName + ";HCAL Isolation; entries", 50,0,10);       
      // -- PF Isolation
      histos[sampleConfig.sampleName]["hPfIso"]   = new TH1F("hPfIso_"   + sampleConfig.sampleName, "hPfIso_" + sampleConfig.sampleName + ";rel PF isolation (dBeta); entries", 100,0,1);
      histos[sampleConfig.sampleName]["hPfIsoCh"] = new TH1F("hPfIsoCh_" + sampleConfig.sampleName, "hPfIsoCh_" + sampleConfig.sampleName + ";PF Ch isolation; entries", 100,0,10); 
      histos[sampleConfig.sampleName]["hPfIsoNe"] = new TH1F("hPfIsoNe_" + sampleConfig.sampleName, "hPfIsoNe_" + sampleConfig.sampleName + ";PF Ne isolation; entries", 100,0,10); 
      histos[sampleConfig.sampleName]["hPfIsoPh"] = new TH1F("hPfIsoPh_" + sampleConfig.sampleName, "hPfIsoPh_" + sampleConfig.sampleName + ";PF Ph isolation; entries", 100,0,10); 
      histos[sampleConfig.sampleName]["hPfIsoPU"] = new TH1F("hPfIsoPU_" + sampleConfig.sampleName, "hPfIsoPU_" + sampleConfig.sampleName + ";PF PU isolation; entries", 100,0,10); 
      // Isolation histograms      
      histos[sampleConfig.sampleName]["hPfIsoCut"]   = new TH1F("hPfIsoCut_" + sampleConfig.sampleName, "hPfIsoCut_" + sampleConfig.sampleName + ";PF isolation CUT; entries", 100,0.0,1.0); 
      histos[sampleConfig.sampleName]["hPfNoPUIsoCut"]= new TH1F("hPfNoPUIsoCut_" + sampleConfig.sampleName, "hPfNoPUIsoCut_" + sampleConfig.sampleName + ";Uncorr. PF isolation CUT; entries", 100,0.0,1.0); 
      histos[sampleConfig.sampleName]["hPfIsoChCut"] = new TH1F("hPfIsoChCut_" + sampleConfig.sampleName, "hPfIsoChCut_" + sampleConfig.sampleName + ";PF Ch isolation CUT; entries", 100,0.0,10); 
      histos[sampleConfig.sampleName]["hPfIsoNeCut"] = new TH1F("hPfIsoNeCut_" + sampleConfig.sampleName, "hPfIsoNeCut_" + sampleConfig.sampleName + ";PF Ne isolation CUT; entries", 100,0.0,10); 
      histos[sampleConfig.sampleName]["hPfIsoPhCut"] = new TH1F("hPfIsoPhCut_" + sampleConfig.sampleName, "hPfIsoPhCut_" + sampleConfig.sampleName + ";PF Ph isolation CUT; entries", 100,0.0,10); 
      histos[sampleConfig.sampleName]["hPfIsoPUCut"] = new TH1F("hPfIsoPUCut_" + sampleConfig.sampleName, "hPfIsoPUCut_" + sampleConfig.sampleName + ";PF PU isolation CUT; entries", 100,0.0,10); 
      histos[sampleConfig.sampleName]["hTrackerIsoCut"] = new TH1F("hTrackerIsoCut_" + sampleConfig.sampleName, "hTrackerIsoCut_" + sampleConfig.sampleName + ";Tracker isolation CUT; entries", 100,0.0,1.0); 
     
    }      
  
  for (auto sampleConfig : sampleConfigs){
      
      TString fileName = sampleConfig.fileName;
      std::cout << "[" << argv[0] << "] Processing file "
		<< fileName.Data() << std::endl;   
      
      // Initialize pointers to summary and full event structure
      
      muon_pog::Event*   ev   = new muon_pog::Event();
      
      // TTree* tree;
      // TBranch* evBranch;
      
      // // Open file, get tree, set branches
      
      // TFile* inputFile = TFile::Open(fileName,"READONLY");
      // tree = (TTree*)inputFile->Get("MUONPOGTREE");
      // if (!tree) inputFile->GetObject("MuonPogTree/MUONPOGTREE",tree);

      // evBranch = tree->GetBranch("event");
      // evBranch->SetAddress(&ev);


      TChain *tree;
      tree = new TChain("MuonPogTree/MUONPOGTREE");

      std::ifstream inList;
      inList.open(Form(fileName,std::ifstream::in));
      while (1) {
	TString inFile;
	inList >> inFile;
	if (!inList.good()) break;
	std::cout << "Loading " << inFile << std::endl;
	tree->Add(inFile);
      }      
      tree->SetBranchAddress( "event",  &ev );
      std::cout << "Total number of files: " << tree->GetNtrees() << std::endl;

      // Watch number of entries
      int nEntries = tree->GetEntries();
      if (sampleConfig.nEvents > 0 && sampleConfig.nEvents < nEntries)
	nEntries = sampleConfig.nEvents;
      std::cout << "[" << argv[0] << "] Number of entries being processed = " << nEntries << std::endl;
      
      int nFilteredEvents = 0;
      
      for (Long64_t iEvent=0; iEvent<nEntries; ++iEvent) 
	{
	  if (tree->LoadTree(iEvent)<0) break;
	  tree->GetEntry(iEvent);
	  
	  if (iEvent % 25000 == 0 )
	    std::cout << "[" << argv[0] << "] processing event : " << iEvent << "\r" << std::flush;
	  
	  float weight = ev->genInfos.size() > 0 ?
	    ev->genInfos[0].genWeight/fabs(ev->genInfos[0].genWeight) : 1.;
	  
	  // -- GEN Particles
	  bool genMuonFromZ = false;
	  std::vector<GenParticle> genMuons;
	  if(sampleConfig.sampleName.Contains("DY")){
	    for (auto & genpart : ev->genParticles){
	      if (std::abs(genpart.pdgId) == 13){
		bool isPromptFinalState = (genpart.flags.at(0) && genpart.status == 1); // 0 -> IsPrompt
		if (isPromptFinalState){
		  for (auto & genmom : genpart.mothers)
		    if (genpart.pt >15 && std::abs(genmom) == 23) genMuonFromZ = true;
		}	      
	      }// if (genMuon)
	      
	      if (genMuonFromZ){
		genMuons.push_back(genpart);
		genMuonFromZ = false;
	      }
	    } // for(genpart)
	  }// if(DY)

	  int nGenMuon = genMuons.size();

	  // -- Vertex
	  ROOT::Math::XYZPoint pVreco;
	  pVreco.SetXYZ(ev->primaryVertex[0], ev->primaryVertex[1], ev->primaryVertex[2]);

	  // -- Muons
	  int EvtGoodMuon = 0;

	  for (auto & muon : ev->muons){
	    
	    if (muon.pt > algoConfig.muon_minPt &&
		hasGoodId(muon, algoConfig.muon_ID)){
	      
	      // Match RECO-GEN Muon
	      if(sampleConfig.sampleName.Contains("DY")){
		bool genrecoMatch = false;
		for (auto & genMuon : genMuons){
		  TVector3 muonV, genMuonV;
		  muonV.SetPtEtaPhi(muon.pt,muon.eta,muon.phi);
		  genMuonV.SetPtEtaPhi(genMuon.pt,genMuon.eta,genMuon.phi);
		  
		  float DeltaMupT = std::abs(muon.pt - genMuon.pt);
		  float DeltaMuR = muonV.DeltaR(genMuonV);
		  
		  if(DeltaMupT < 5 && DeltaMuR < 0.2) genrecoMatch = true;
		}// for(GoodGenMuon)
		
		if (!genrecoMatch) continue;
	      } // if(DY)

	      EvtGoodMuon++;
	      
	      histos[sampleConfig.sampleName]["hNVertex"]->Fill(ev->nVtx, weight);
	      histos[sampleConfig.sampleName]["hpT"]     ->Fill(muon.pt, weight);
	      histos[sampleConfig.sampleName]["hEta"]    ->Fill(muon.eta, weight);	      

	      // Detector Isolation
	      histos[sampleConfig.sampleName]["hTrackerIso"]->Fill(muon.trackerIso/muon.pt, weight);
	      histos[sampleConfig.sampleName]["hEMCalIso"]  ->Fill(muon.EMCalIso, weight);
	      histos[sampleConfig.sampleName]["hHCalIso"]   ->Fill(muon.HCalIso, weight);

	      // PF Isolation and its components
	      histos[sampleConfig.sampleName]["hPfIso"]  ->Fill(muon.isoPflow04, weight);
	      histos[sampleConfig.sampleName]["hPfIsoCh"]->Fill(muon.chargedHadronIso, weight);
	      histos[sampleConfig.sampleName]["hPfIsoNe"]->Fill(muon.neutralHadronIso, weight); 
	      histos[sampleConfig.sampleName]["hPfIsoPh"]->Fill(muon.photonIso, weight); 
	      histos[sampleConfig.sampleName]["hPfIsoPU"]->Fill(muon.chargedHadronIsoPU, weight); 

	      // Efficiencies Vs Kinematic Variables
	      Float_t IsoCutVal[3][2];
	      // passPFIso
	      IsoCutVal[0][0] = algoConfig.muon_PFIsoCut_Loose; // Loose
	      IsoCutVal[0][1] = algoConfig.muon_PFIsoCut_Tight; // Tight
	      // passPFNoPUIso
	      IsoCutVal[1][0] = algoConfig.muon_PFNoPUIsoCut_Loose; // Loose
	      IsoCutVal[1][1] = algoConfig.muon_PFNoPUIsoCut_Tight; // Tight
	      // passTkrIso
	      IsoCutVal[2][0] = algoConfig.muon_TkrIsoCut_Loose; // Loose
	      IsoCutVal[2][1] = algoConfig.muon_TkrIsoCut_Tight; // Tight

	      float PFNoPUIso = (muon.chargedHadronIso + muon.neutralHadronIso + muon.photonIso)/muon.pt;
	      float IsoVal[3] = {muon.isoPflow04, PFNoPUIso, muon.trackerIso/muon.pt};
	      
	      for(int h=0; h<3; h++){
		for(int w=0; w<2; w++){
		  if(IsoVal[h] < IsoCutVal[h][w]){
		    histos[sampleConfig.sampleName]["hNVertex_" + IsoHistos[h] + "_" + IsoWP[w]]->Fill(ev->nVtx, weight);
		    histos[sampleConfig.sampleName]["hpT_"      + IsoHistos[h] + "_" + IsoWP[w]]->Fill(muon.pt,  weight);
		    histos[sampleConfig.sampleName]["hEta_"     + IsoHistos[h] + "_" + IsoWP[w]]->Fill(muon.eta, weight);
		  }
		}
	      }

	      // ROC Curves: Isolation Cut
	      histos[sampleConfig.sampleName]["hPfIsoCut"]     ->Fill(-10, weight); // Underflows to count the total number of muons
	      histos[sampleConfig.sampleName]["hPfNoPUIsoCut"] ->Fill(-10, weight);
	      histos[sampleConfig.sampleName]["hPfIsoChCut"]   ->Fill(-10, weight);
	      histos[sampleConfig.sampleName]["hPfIsoNeCut"]   ->Fill(-10, weight); 
	      histos[sampleConfig.sampleName]["hPfIsoPhCut"]   ->Fill(-10, weight); 
	      histos[sampleConfig.sampleName]["hPfIsoPUCut"]   ->Fill(-10, weight); 
	      histos[sampleConfig.sampleName]["hTrackerIsoCut"]->Fill(-10, weight); 

	      float isocut_max    = 1.0;
	      float isocut_bin    = 0.01;
  	      float isoComcut_max = 10;
	      float isoComcut_bin = 0.1;
	      float isoTkrcut_max = 1.0;
	      float isoTkrcut_bin = 0.01;
 	      for(int step = 0; step < 100; step++){
		if(muon.isoPflow04 < isocut_max) 
		  histos[sampleConfig.sampleName]["hPfIsoCut"]->Fill(isocut_max - (isocut_bin/2.), weight);
		if(PFNoPUIso < isocut_max) 
		  histos[sampleConfig.sampleName]["hPfNoPUIsoCut"]->Fill(isocut_max - (isocut_bin/2.), weight);
		if(muon.chargedHadronIso < isoComcut_max) 
		  histos[sampleConfig.sampleName]["hPfIsoChCut"]->Fill(isoComcut_max - (isoComcut_bin/2.), weight);
		if(muon.neutralHadronIso < isoComcut_max) 
		  histos[sampleConfig.sampleName]["hPfIsoNeCut"]->Fill(isoComcut_max - (isoComcut_bin/2.), weight);
		if(muon.photonIso < isoComcut_max) 
		  histos[sampleConfig.sampleName]["hPfIsoPhCut"]->Fill(isoComcut_max - (isoComcut_bin/2.), weight);
		if(muon.chargedHadronIsoPU < isoComcut_max) 
		  histos[sampleConfig.sampleName]["hPfIsoPUCut"]->Fill(isoComcut_max - (isoComcut_bin/2.), weight);
		if(muon.trackerIso/muon.pt < isocut_max) 
		  histos[sampleConfig.sampleName]["hTrackerIsoCut"]->Fill(isocut_max - (isocut_bin/2.), weight);

		isocut_max    = isocut_max - isocut_bin; // binning
		isoComcut_max = isoComcut_max - isoComcut_bin; 
	      }
	      
	    } // if(goodMuon)
	  }  // for(muons)	  


	  // -- Global Event Variables (filled when event has good muons)
	  if (EvtGoodMuon > 0) {
	    histos[sampleConfig.sampleName]["hNGENMuons"] ->Fill(nGenMuon, weight);
	    histos[sampleConfig.sampleName]["hNMuons"]    ->Fill(EvtGoodMuon, weight);
	    histos[sampleConfig.sampleName]["hEvtNVertex"]->Fill(ev->nVtx, weight);
	  }
	  
	} // for(iEvent)
      
      // Normalization of efficiency histograms to 1
      if (histos[sampleConfig.sampleName]["hNMuons"]->Integral() != 0.0){
	histos[sampleConfig.sampleName]["hPfIsoCut"]    ->Scale(1.0/histos[sampleConfig.sampleName]["hPfIsoCut"]->GetBinContent(0));
	histos[sampleConfig.sampleName]["hPfNoPUIsoCut"]->Scale(1.0/histos[sampleConfig.sampleName]["hPfNoPUIsoCut"]->GetBinContent(0));
	histos[sampleConfig.sampleName]["hPfIsoChCut"]  ->Scale(1.0/histos[sampleConfig.sampleName]["hPfIsoChCut"]->GetBinContent(0));
	histos[sampleConfig.sampleName]["hPfIsoNeCut"]  ->Scale(1.0/histos[sampleConfig.sampleName]["hPfIsoNeCut"]->GetBinContent(0));
	histos[sampleConfig.sampleName]["hPfIsoPhCut"]  ->Scale(1.0/histos[sampleConfig.sampleName]["hPfIsoPhCut"]->GetBinContent(0));
	histos[sampleConfig.sampleName]["hPfIsoPUCut"]  ->Scale(1.0/histos[sampleConfig.sampleName]["hPfIsoPUCut"]->GetBinContent(0));

	histos[sampleConfig.sampleName]["hTrackerIsoCut"]->Scale(1.0/histos[sampleConfig.sampleName]["hTrackerIsoCut"]->GetBinContent(0));
      }

      // Efficiency Plots Vs Kinematic Variables
      
      for(int h=0; h<3; h++){
	for(int w=0; w<2; w++){
	  EfficiencyHisto(&histos, sampleConfig.sampleName, "hNVertex_" + IsoHistos[h] + "_" + IsoWP[w], "hNVertex");
	  EfficiencyHisto(&histos, sampleConfig.sampleName, "hpT_" + IsoHistos[h] + "_" + IsoWP[w],      "hpT");      
	  EfficiencyHisto(&histos, sampleConfig.sampleName, "hEta_" + IsoHistos[h] + "_" + IsoWP[w],     "hEta");      
	}
      }
      
      delete ev;
      tree->Delete();
      std::cout << std::endl;
  } // for(sampleConfig)

  TGraphErrors *gPFIso, *gPFNoPUIso, *gTrackerIso;
  if (sampleConfigs.size() < 2){
    std::cout << "Histograms produced only for " << sampleConfigs[0].fileName << std::endl ; 
    std::cout << "ROC curves have NOT been created!" << std::endl ; 
  }
  else if (sampleConfigs.size() == 2){
    std::cout << "Producing ROC curves for " << std::endl ; 
    std::cout << sampleConfigs[0].fileName << std::endl ; 
    std::cout << sampleConfigs[1].fileName << std::endl ; 

    gPFIso      = EfficiencyPlot(&histos, sampleConfigs, "hPfIsoCut");
    gPFNoPUIso  = EfficiencyPlot(&histos, sampleConfigs, "hPfNoPUIsoCut");
    gTrackerIso = EfficiencyPlot(&histos, sampleConfigs, "hTrackerIsoCut");

  } // if (sampleConfigs.size() == 2)

  outputFile->cd();
  gPFIso->Write();  
  gPFNoPUIso->Write();  
  gTrackerIso->Write();  
  outputFile->Write();
  
  if (!gROOT->IsBatch()) app->Run();

  return 0;

}


muon_pog::AlgoConfig::AlgoConfig(boost::property_tree::ptree::value_type & vt)
{

  try
    {
      muon_minPt     = vt.second.get<Float_t>("muon_minPt");
      muon_ID        = vt.second.get<std::string>("muon_ID");
      muon_PFIsoCut_Loose = vt.second.get<Float_t>("muon_PFIsoCut_Loose");
      muon_PFIsoCut_Tight = vt.second.get<Float_t>("muon_PFIsoCut_Tight");
      muon_PFNoPUIsoCut_Loose = vt.second.get<Float_t>("muon_PFNoPUIsoCut_Loose");
      muon_PFNoPUIsoCut_Tight = vt.second.get<Float_t>("muon_PFNoPUIsoCut_Tight");
      muon_TkrIsoCut_Loose = vt.second.get<Float_t>("muon_TkrIsoCut_Loose");
      muon_TkrIsoCut_Tight = vt.second.get<Float_t>("muon_TkrIsoCut_Tight");
    }

  catch (boost::property_tree::ptree_bad_data bd)
    {
      std::cout << "[AlgoConfig] Can't get data : has error : "
		<< bd.what() << std::endl;
      throw std::runtime_error("Bad INI variables");
    }

}

muon_pog::SampleConfig::SampleConfig(boost::property_tree::ptree::value_type & vt)
{
  
  try
    {
      fileName         = TString(vt.second.get<std::string>("fileName").c_str());
      sampleName       = TString(vt.first.c_str());
      cSection         = vt.second.get<Float_t>("cSection");
      nEvents          = vt.second.get<Float_t>("nEvents"); //CB do we really need this? can't we take nEvents from the file itself?
      applyReweighting = vt.second.get<Bool_t> ("applyReweighting");
    }
  
  catch (boost::property_tree::ptree_bad_data bd)
    {
      std::cout << "[TagAndProbeConfig] Can't get data : has error : "
		<< bd.what() << std::endl;
      throw std::runtime_error("Bad INI variables");
    }
  
}

void muon_pog::parseConfig(const std::string configFile, muon_pog::AlgoConfig & algoConfig,
			   std::vector<muon_pog::SampleConfig> & sampleConfigs)
{

  boost::property_tree::ptree pt;
  
  try
    {
      boost::property_tree::ini_parser::read_ini(configFile, pt);
    }
  catch (boost::property_tree::ini_parser::ini_parser_error iniParseErr)
    {
      std::cout << "[parseConfig] Can't open : " << iniParseErr.filename()
		<< "\n\tin line : " << iniParseErr.line()
		<< "\n\thas error :" << iniParseErr.message()
		<< std::endl;
      throw std::runtime_error("Bad INI parsing");
    }
  
  for( auto vt : pt )
    {
      if (vt.first.find("Algo") != std::string::npos)
	algoConfig = muon_pog::AlgoConfig(vt);
      else
	sampleConfigs.push_back(muon_pog::SampleConfig(vt));
    }
}


/***********************************************
             Efficiency Histograms
**********************************************/
void muon_pog::EfficiencyHisto(std::map<TString, std::map<TString, TH1F *> > *histos, TString samname, TString hisnameN, TString hisnameD){

  TH1F *histo1 = (TH1F*)(*histos)[samname][hisnameN]->Clone();
  TH1F *histo2 = (TH1F*)(*histos)[samname][hisnameD]->Clone();

  for(int i_bin=1; i_bin <= histo1->GetNbinsX() ; i_bin++){
    
    float npass_hist        = 0.0;
    float npass_error_hist2 = 0.0;
    float n_hist            = 0.0;
    float n_error_hist2     = 0.0;
    float eff               = 0.0;
    float eff_error2        = 0.0;
    
    npass_hist = histo1->GetBinContent(i_bin); 
    n_hist     = histo2->GetBinContent(i_bin); 

    npass_error_hist2 = (histo1->GetBinError(i_bin))*(histo1->GetBinError(i_bin)); 
    n_error_hist2     = (histo2->GetBinError(i_bin))*(histo2->GetBinError(i_bin)); 
    
    if(n_hist != 0.0 && npass_hist != 0.0){
      eff        = npass_hist / n_hist;
      eff_error2 = (npass_hist / n_hist)*(npass_hist / n_hist)*
	( npass_error_hist2*(1.0/(npass_hist*npass_hist)) + (n_error_hist2*(1.0/(n_hist*n_hist))) );
    }
    else {
      eff        = 0.0;
      eff_error2 = 0.0;
    }
    
    (*histos)[samname][hisnameN + "_Eff"]->SetBinContent(i_bin, eff);
    (*histos)[samname][hisnameN + "_Eff"]->SetBinError  (i_bin, sqrt(eff_error2));
    
  } // for(i_bin)
  histo1->Delete();
  histo2->Delete();
  
}

TGraphErrors* muon_pog::EfficiencyPlot(std::map<TString, std::map<TString, TH1F *> > *histos, std::vector<SampleConfig> sampleConfigs, TString hisname){
  
  const Int_t n = (*histos)[sampleConfigs[0].sampleName][hisname]->GetNbinsX();
  
  Double_t IsoCutEff[4][n];
  
  for (int i_eff=1; i_eff < n; i_eff++){
    
    IsoCutEff[0][i_eff-1] = (*histos)[sampleConfigs[1].sampleName][hisname]->GetBinContent(i_eff);
    IsoCutEff[1][i_eff-1] = (*histos)[sampleConfigs[0].sampleName][hisname]->GetBinContent(i_eff);
    
    IsoCutEff[2][i_eff-1] = (*histos)[sampleConfigs[1].sampleName][hisname]->GetBinError(i_eff);
    IsoCutEff[3][i_eff-1] = (*histos)[sampleConfigs[0].sampleName][hisname]->GetBinError(i_eff);
    
  }// for(bins);
  
  TGraphErrors *gEffIsoCut;
  
  gEffIsoCut = new TGraphErrors(n-1,IsoCutEff[0], IsoCutEff[1], IsoCutEff[2], IsoCutEff[3]);
  
  gEffIsoCut->SetNameTitle(sampleConfigs[1].sampleName + "_Vs_" + sampleConfigs[0].sampleName + "_" + hisname,
			   sampleConfigs[1].sampleName + " Vs " + sampleConfigs[0].sampleName + ": " + hisname);

  gEffIsoCut->GetXaxis()->SetTitle(sampleConfigs[1].sampleName);
  gEffIsoCut->GetYaxis()->SetTitle(sampleConfigs[0].sampleName);
  return gEffIsoCut;
  
}

