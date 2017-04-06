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
#include "TGraphAsymmErrors.h"
#include "THStack.h"
#include "TLegend.h"
#include "TRegexp.h"
#include "TTree.h"
#include "TBranch.h"
#include "TLorentzVector.h"

#include "../src/MuonPogTree.h"
#include "Utils.h"
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
// 1. SampleConfig : configuration class containing sample information
// 2. TagAndProbeConfig : configuration class containing TnP cuts information
// 3. Plotter : class containing the plot definition and defining the plot filling 
//              for a given sample <= CB modify this to add new variables
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
    std::vector<int> runs;
        
    SampleConfig() {};
    
#ifndef __MAKECINT__ // CB CINT doesn't like boost :'-(    
    SampleConfig(boost::property_tree::ptree::value_type & vt); 
#endif

    ~SampleConfig() {};

  private:
    std::vector<int> toArray(const std::string & entries); 
    
  };

  class TagAndProbeConfig {

  public :
    
    // config parameters (public for direct access)
    
    Float_t     gen_DrCut;
    Float_t     gen_DrCutVeto;
    
    std::string muon_trackType; // applies to both, tag and probe     
    TString muon_cuts;
  
    Float_t probe_minPt;      
   
    TagAndProbeConfig() {};
    
#ifndef __MAKECINT__ // CB CINT doesn't like boost :'-(    
    TagAndProbeConfig(boost::property_tree::ptree::value_type & vt); 
#endif

    ~TagAndProbeConfig() {};
    
  private:
    std::vector<TString> toArray(const std::string & entries); 
  
  };


  class Plotter {

  public :

    enum HistoType {GLB=0, TRK, LOOSE, TIGHT, MEDIUM, SOFT, HIGHpT};
      
    Plotter(muon_pog::TagAndProbeConfig tnpConfig, muon_pog::SampleConfig & sampleConfig) :
      m_tnpConfig(tnpConfig) , m_sampleConfig(sampleConfig) {};
    ~Plotter() {};
    
    void book(TFile *outFile);
    void fill(const std::vector<muon_pog::Muon> & muons, const muon_pog::HLT & hlt, const muon_pog::Event & ev, float weight);
    void fillGen(const std::vector<muon_pog::GenParticle> & genpars, const muon_pog::Event & ev);

    std::map<TString, std::map<TString, std::map<TString, TH1F *> > >        m_plots;
    std::map<TString, std::map<TString, std::map<TString, TH2F *> > >        m_2Dplots;
    std::map<TString, std::map<TString, std::map<TString, TProfile *> > >    m_prof;
    std::map<TString, std::map<TString, std::map<TString, TEfficiency *> > > m_effs;

    std::map<TString, TH2F *> m_2Dyields;

    TagAndProbeConfig m_tnpConfig;
    SampleConfig      m_sampleConfig;
        
  };
  
}

// Helper classes defintion *****
// 1. parseConfig : parse the full cfg file
// 1. comparisonPlot : make a plot overlayng data and MC for a given plot
// ******************************

namespace muon_pog {
  void parseConfig(const std::string configFile, TagAndProbeConfig & tpConfig,
		   std::vector<SampleConfig> & sampleConfigs);
  
  void comparisonPlots(std::vector<Plotter> & plotters,
		       TFile *outFile, TString &  outputDir);

  void copyPhp(const TString &  outputDir);
  Bool_t IsMediumHIP(const muon_pog::Muon & muon);
  Bool_t IsSoftHIP  (const muon_pog::Muon & muon);
  // void setTProfY(TProfile &prof1, TProfile &prof2);
}



// The main program******** *****
// 1. Get configuration file and produces configurations
// 2. Create Plotters and loop on the event to fill them
// 3. Writes results in cnfigurable outuput file
// ******************************

int main(int argc, char* argv[]){
  using namespace muon_pog;


  if (argc != 4) 
    {
      std::cout << "Usage : "
		<< argv[0] << " PATH_TO_CONFIG_FILE PATH_TO_OUTPUT_DIR OUTPUT_FILE_NAME\n";
      exit(100);
    }

  std::string configFile(argv[1]);
  
  std::cout << "[" << argv[0] << "] Using config file " << configFile << std::endl;

  // Output directory
  TString dirName    = argv[2];
  TString outputName = argv[3];
  system("mkdir -p " + dirName);
  TFile* outputFile = TFile::Open(dirName + "/histos_" + outputName  + ".root","RECREATE"); // CB find a better name for output file  

  // Set it to kTRUE if you do not run interactively
  gROOT->SetBatch(kTRUE); 

  // Initialize Root application
  TRint* app = new TRint("CMS Root Application", &argc, argv);

  setTDRStyle();
 
  TagAndProbeConfig tnpConfig;
  std::vector<SampleConfig> sampleConfigs;
  
  parseConfig(configFile,tnpConfig,sampleConfigs);

  std::vector<Plotter> plotters;

  for (auto sampleConfig : sampleConfigs)
    {

      Plotter plotter(tnpConfig, sampleConfig);
      plotter.book(outputFile);
      
      plotters.push_back(plotter);
    }
 
  for (auto plotter : plotters)
    {

      TString fileName = plotter.m_sampleConfig.fileName;
      std::cout << "[" << argv[0] << "] Processing file "
		<< fileName.Data() << std::endl;  
  
      // Initialize pointers to summary and full event structure

      muon_pog::Event*   ev   = new muon_pog::Event();

      TTree* tree;
      TBranch* evBranch;

      // Open file, get tree, set branches

      TFile* inputFile = TFile::Open(fileName,"READONLY");
      tree = (TTree*)inputFile->Get("MUONPOGTREE");
      if (!tree) inputFile->GetObject("MuonPogTree/MUONPOGTREE",tree);

      evBranch = tree->GetBranch("event");
      evBranch->SetAddress(&ev);

      // Watch number of entries
      int nEntries = plotter.m_sampleConfig.nEvents > 0 ? plotter.m_sampleConfig.nEvents : tree->GetEntriesFast();
      std::cout << "[" << argv[0] << "] Number of entries = " << nEntries << std::endl;

      int nFilteredEvents = 0;

      for (Long64_t iEvent=0; iEvent<nEntries; ++iEvent) {
	  if (tree->LoadTree(iEvent)<0) break;

	  if (iEvent % 25000 == 0 )
	    std::cout << "[" << argv[0] << "] processing event : " << iEvent << "\r" << std::flush;
	    
	  evBranch->GetEntry(iEvent);
	  float weight = ev->genInfos.size() > 0 ?
	    ev->genInfos[0].genWeight/fabs(ev->genInfos[0].genWeight) : 1.;
	 
	  // CB to be fixed when kown what to do for MC !!!

	  //	  if(plotter.m_sampleConfig.applyReweighting==true)
	  //  weight *= ev->nVtx < 60 ? PUweight[ev->nVtx] : 0;
	  
	  plotter.fill(ev->muons, ev->hlt, (*ev), weight);

	  plotter.fillGen(ev->genParticles, (*ev));
	}
      
      delete ev;
      inputFile->Close();
      std::cout << std::endl;
	   
    }
  
  muon_pog::comparisonPlots(plotters,outputFile,dirName);
  muon_pog::copyPhp(dirName);
  
  outputFile->Write();
  
  if (!gROOT->IsBatch()) app->Run();

  return 0;

}


muon_pog::TagAndProbeConfig::TagAndProbeConfig(boost::property_tree::ptree::value_type & vt)
{

  try
    {

      gen_DrCut      = vt.second.get<Float_t>("gen_DrCut");
      gen_DrCutVeto  = vt.second.get<Float_t>("gen_DrCutVeto");
     
      muon_trackType = vt.second.get<std::string>("muon_trackType");

      probe_minPt  = vt.second.get<Float_t>("probe_minPt");

      muon_cuts = TString(vt.second.get<std::string>("muon_cuts").c_str());

    }

  catch (boost::property_tree::ptree_bad_data bd)
    {
      std::cout << "[TagAndProbeConfig] Can't get data : has error : "
		<< bd.what() << std::endl;
      throw std::runtime_error("Bad INI variables");
    }

}

muon_pog::SampleConfig::SampleConfig(boost::property_tree::ptree::value_type & vt)
{
 
 try
    {
      fileName     = TString(vt.second.get<std::string>("fileName").c_str());
      sampleName   = TString(vt.first.c_str());
      cSection = vt.second.get<Float_t>("cSection");
      nEvents = vt.second.get<Float_t>("nEvents"); //CB do we really need this? can't we take nEvents from the file itself?
      applyReweighting = vt.second.get<Bool_t>("applyReweighting");
      runs = toArray(vt.second.get<std::string>("runs"));
    }
  
  catch (boost::property_tree::ptree_bad_data bd)
    {
      std::cout << "[TagAndProbeConfig] Can't get data : has error : "
		<< bd.what() << std::endl;
      throw std::runtime_error("Bad INI variables");
    }
  
}

std::vector<int> muon_pog::SampleConfig::toArray(const std::string& entries)
{
  
  std::vector<int> result;
  std::stringstream sentries(entries);
  std::string item;
  while(std::getline(sentries, item, ','))
    result.push_back(atoi(item.c_str()));
  return result;

}


std::vector<TString> muon_pog::TagAndProbeConfig::toArray(const std::string& entries)
{
  
  std::vector<TString> result;
  std::stringstream sentries(entries);
  std::string item;
  while(std::getline(sentries, item, ','))
    result.push_back(TString(item));
  return result;

}


void muon_pog::Plotter::book(TFile *outFile)
{

  // const double etaBins[18] = {-2.4,-2.1,-1.6,-1.2,-1.05,-0.9,-0.6,-0.3,-0.2,
  //  			      0.2, 0.3, 0.6, 0.9, 1.05, 1.2, 1.6, 2.1, 2.4};
  const double etaBins[9] = {-2.4,-1.8,-1.2,-0.6,0.0,
   			      0.6, 1.2, 1.8, 2.4};

  //const double pTBins[24] = {5,6,7,8,9,10,12,14,16,18,20,22,24,26,28,30,35,40,45,50,60,70,80,100};
  //const double pTBins[13] = {5,7,9,12,16,20,25,30,35,40,50,70,100};
  const double pTBins[7] = {5,15,25,35,50,70,100};

  const double VtxBins[10] = {0,10,20,24,28,32,36,40,50,60};
	  
  TString sampleTag = m_sampleConfig.sampleName;
  
  outFile->cd("/");
  outFile->mkdir(sampleTag);
  outFile->cd(sampleTag);
  
  outFile->mkdir(sampleTag + "/efficiencies");      
  outFile->mkdir(sampleTag + "/kinematical_variables");
  outFile->mkdir(sampleTag + "/Yields");
 
  TString IDName[8] = {"GLB",  "TRK", "LOOSE", "MEDIUM", "TIGHT", "SOFT", "HIGHpT", "TRKNoGLB"};
  TString region[4] = {"Full", "Barrel", "Endcap", "Overlap"};
  TString Origin[4] = {"bGEN", "WGEN", "SIGN", "FAKE"};

  TH1::SetDefaultSumw2(kTRUE);
  
  // -- Yields
  for (unsigned nor=0; nor<5; nor++){
    TString hNameOrigin;
    if (nor<4) hNameOrigin = Origin[nor];
    else hNameOrigin = "TOTAL";
    outFile->cd(sampleTag + "/Yields/");
    m_2Dyields[hNameOrigin] = new TH2F("Yields_" + hNameOrigin ,"Yields " + sampleTag + " " + hNameOrigin,10,0,10,4,0,4);
    m_2Dyields[hNameOrigin]->SetOption("COLTEXT"); 
    for (unsigned nid=0; nid<8; nid++) m_2Dyields[hNameOrigin]->GetXaxis()->SetBinLabel(nid+1,IDName[nid]);
    m_2Dyields[hNameOrigin]->GetXaxis()->SetBinLabel(9, "Total Pass Cut");
    m_2Dyields[hNameOrigin]->GetXaxis()->SetBinLabel(10,"Total");
    for (unsigned nre=0; nre<4; nre++) m_2Dyields[hNameOrigin]->GetYaxis()->SetBinLabel(nre+1,region[nre]);

    m_2Dyields["GenMuon"] = new TH2F("Yields_NGenMuons" ,"Yields GenMuons " + sampleTag,4,0,4,4,0,4);
    m_2Dyields["GenMuon"]->SetOption("COLTEXT"); 
    m_2Dyields["GenMuon"]->GetXaxis()->SetBinLabel(1, "Total Events");
    m_2Dyields["GenMuon"]->GetXaxis()->SetBinLabel(2, "Gen Muon W");
    m_2Dyields["GenMuon"]->GetXaxis()->SetBinLabel(3, "Gen Muon Others");
    m_2Dyields["GenMuon"]->GetXaxis()->SetBinLabel(4, "Gen Muon Total");
    for (unsigned nre=0; nre<4; nre++) m_2Dyields["GenMuon"]->GetYaxis()->SetBinLabel(nre+1,region[nre]);

  }
  
  for (unsigned nid=0; nid<8; nid++){
    outFile->mkdir(sampleTag + "/kinematical_variables/" + IDName[nid]);
    outFile->mkdir(sampleTag + "/efficiencies/"          + IDName[nid]);      
    for (unsigned nre=0; nre<4; nre++){
      outFile->mkdir(sampleTag + "/kinematical_variables/" + IDName[nid] + "/" + region[nre]);
      outFile->mkdir(sampleTag + "/efficiencies/"          + IDName[nid] + "/" + region[nre]);      
      for (unsigned nor=0; nor<4; nor++){
	outFile->cd(sampleTag+"/kinematical_variables/" + IDName[nid] + "/" + region[nre]);
	m_plots[IDName[nid]][region[nre]][Origin[nor] + "Vtx"] = new TH1F ("Vtx_" + Origin[nor] + "_" + IDName[nid] + "_" + region[nre],  IDName[nid] + " Vtx " + region[nre] + "; N. Vtx; # entries", 60,0.,60.);
	m_plots[IDName[nid]][region[nre]][Origin[nor] + "Pt"]  = new TH1F ("Pt_"  + Origin[nor] + "_" + IDName[nid] + "_" + region[nre],  IDName[nid] + " Pt "  + region[nre] + "; p_{T} [GeV]; # entries", 100,0.,100);
	m_plots[IDName[nid]][region[nre]][Origin[nor] + "Eta"] = new TH1F ("Eta_" + Origin[nor] + "_" + IDName[nid] + "_" + region[nre],  IDName[nid] + " Eta " + region[nre] + "; #eta; # entries", 48,-2.4, 2.4);
	m_plots[IDName[nid]][region[nre]][Origin[nor] + "VHF"] = new TH1F ("VHF_" + Origin[nor] + "_" + IDName[nid] + "_" + region[nre],  IDName[nid] + " Number of Valid Hit Fraction " + region[nre] + "; #N. Valid Hit Fraction; # entries", 55, 0.0, 1.1);
	
	outFile->cd(sampleTag+"/efficiencies/" + IDName[nid] + "/" + region[nre]);      
	m_effs[IDName[nid]][region[nre]][Origin[nor] + "EffVsPt"]    = new TEfficiency("EffVsPt_"  + Origin[nor] + "_" + IDName[nid] + "_" + region[nre], IDName[nid] + " Eff. Vs Pt "   + region[nre] + ";p_{T} [GeV]; Eff.",  6, pTBins);
	m_effs[IDName[nid]][region[nre]][Origin[nor] + "EffVsEta"]   = new TEfficiency("EffVsEta_" + Origin[nor] + "_" + IDName[nid] + "_" + region[nre], IDName[nid] + " Eff. Vs Eta "  + region[nre] + ";#eta; Eff.", 8, etaBins);
	m_effs[IDName[nid]][region[nre]][Origin[nor] + "EffVsVtx"]  = new TEfficiency("EffVsVtx_" + Origin[nor] + "_" + IDName[nid] + "_" + region[nre], IDName[nid] + " Eff. Vs Vtx "  + region[nre] + ";Vtx; Eff.", 9, VtxBins);
      } // for(nor)
    } // for(nre)
  } // for(nid)
  
} // void::book 

void muon_pog::Plotter::fillGen(const std::vector<muon_pog::GenParticle> & genpars, const muon_pog::Event & ev){
  
  float etaBarrel = 1.2;  
  float pTeta = 20.0;
  TString region[4] = {"Full","Barrel","Endcap","Overlap"};
  
  m_2Dyields["GenMuon"]->Fill(0.0, 0.0);
  
  //std::cout << " ++++ EVENT ++++ " << std::endl;
  for (auto & genpar : genpars){

    bool genWmuon     = false;
    bool genOthermuon = false;

    if(abs(genpar.pdgId) == 13 &&
       genpar.pt          > m_tnpConfig.probe_minPt &&
       fabs(genpar.eta)   < 2.4 ){

      bool IsMuonpTCorr = false;
      for (auto motherId : genpar.mothers){
	if(abs(motherId) == 13) IsMuonpTCorr =true;
      }      
      if (IsMuonpTCorr) continue;

      // IsLastCopy to avoid muon double counting
      // if (!genpar.flags[13]) continue;
      // -- Debug
      // std::cout << " ----------------- " << std::endl;
      // std::cout << "pT = " << genpar.pt << " - eta = " << genpar.eta << std::endl;
      // for (auto motherId : genpar.mothers) std::cout << "Mother = " << motherId << " - " << genpar.pt << std::endl;
      // std::cout << "Flags = " << std::endl;
      // int idf = 0;
      // for (auto muonflag : genpar.flags){
      // 	std::cout << idf << ": " <<  muonflag << " - " ;
      // 	idf++;
      // }
      // std::cout << "Is Last Copy = " << genpar.flags[13] << std::endl;

      genWmuon = hasMother(genpar, 24);
      if(!genWmuon) genOthermuon = true;
    }
    if(genWmuon || genOthermuon){
      
      //if(genOthermuon) for (auto motherId : genpar.mothers) std::cout << "Mother = " << motherId << " - " << genpar.pt << std::endl;
      for (unsigned int nre = 0; nre<4; nre++){
	
	bool etaRegion = false;
	if (nre == 0) etaRegion = true;
	if (nre == 1 && std::abs(genpar.eta) <  etaBarrel) etaRegion = true;
	if (nre == 2 && std::abs(genpar.eta) >= etaBarrel) etaRegion = true;
	if (nre == 3 && (std::abs(genpar.eta) > 0.9 && std::abs(genpar.eta) < 1.3)) etaRegion = true;
	
	if (!etaRegion)     continue;
	
	if (genWmuon)     m_2Dyields["GenMuon"]->Fill(1,nre);
	if (genOthermuon) m_2Dyields["GenMuon"]->Fill(2,nre);
	m_2Dyields["GenMuon"]->Fill(3,nre);


      }// for(nre)
    }//if(GenMuon)
    
  }// for(genpar)  
}

void muon_pog::Plotter::fill(const std::vector<muon_pog::Muon> & muons,
			     const muon_pog::HLT & hlt, const muon_pog::Event & ev, float weight)
{

  float etaBarrel = 1.2;  
  float pTeta = 20.0;
  
  for (auto & muon : muons){

    // Total Number of Muons
    m_2Dyields["SIGN"]->Fill(9, 0);

    // if (muon.isCloneFilter) std::cout << "-- Clone Muon " << std::endl; 
    // if (muon.isBadFilter)   std::cout << "-- Bad Muon " << std::endl;    
    // if (muon.isBadFilter || muon.isCloneFilter)   std::cout << "-- Bad/Clone Muon " << std::endl;    
    // else continue;
 

    // General Probe Muons	      
    if(muon.pt > m_tnpConfig.probe_minPt && 
       fabs(muon.eta) < 2.4){

      // Number of Muons Passing (pT,eta) cuts
      m_2Dyields["SIGN"]->Fill(8, 0);

      // if(muon.isLoose){
	Bool_t dXYZ  = std::abs(muon.dxyBest) < 0.2 && std::abs(muon.dzBest)  < 0.5;
	Bool_t MuIso = muon.isoPflow03 < 0.4; 
	
	if(m_tnpConfig.muon_cuts == "Iso"   && !(MuIso))         continue;
	if(m_tnpConfig.muon_cuts == "IsoIP" && !(MuIso && dXYZ)) continue;
	
	bool IsbGEN = false;
	bool IsWGEN = false;
	bool IsSIGN = false;
	bool IsFAKE = false;

	// Check first if muon comes from W then from b.
	IsWGEN = muon_pog::hasGenMatch(muon, ev.genParticles, m_tnpConfig.gen_DrCut, 24);
	if(!IsWGEN) IsbGEN = muon_pog::hasbGenMatch(muon,ev.genParticles, 0.2, 3.);
	IsSIGN = IsbGEN || IsWGEN;
	if (!IsSIGN) IsFAKE = muon_pog::hasNoGenMatch(muon,ev.genParticles,m_tnpConfig.gen_DrCutVeto, 3.);

	TString IDName[8] = {"GLB",  "TRK", "LOOSE", "MEDIUM", "TIGHT", "SOFT", "HIGHpT","TRKNoGLB"};
	bool IsMuonID[8];
	for (int bmid = 0; bmid<8; bmid++) IsMuonID[bmid] = false;
	// {"GLB", "TRK", "LOOSE", "MEDIUM", "TIGHT", "SOFT", "HIGHpT", "TRKNoGLB"};
	IsMuonID[0] = muon.isGlobal ; 
	IsMuonID[1] = muon.isTracker ; 
	IsMuonID[2] = muon.isLoose ; 
	IsMuonID[3] = IsMediumHIP(muon) ; 
	IsMuonID[4] = muon.isTight ; 
	IsMuonID[5] = IsSoftHIP(muon) ; 
	IsMuonID[6] = muon.isHighPt ; 
	IsMuonID[7] = muon.isTracker && !(muon.isGlobal) ; 

	TString region[4] = {"Full","Barrel","Endcap","Overlap"};


	for (unsigned int mid = 0; mid<8; mid++){

	  for (unsigned int nre = 0; nre<4; nre++){
	    	    
	    TLorentzVector muTk(muon_pog::muonTk(muon, m_tnpConfig.muon_trackType));
	    bool etaRegion = false;
	    if (nre == 0) etaRegion = true;
	    if (nre == 1 && std::abs(muTk.Eta()) <  etaBarrel) etaRegion = true;
	    if (nre == 2 && std::abs(muTk.Eta()) >= etaBarrel) etaRegion = true;
	    if (nre == 3 && (std::abs(muTk.Eta()) > 0.9 && std::abs(muTk.Eta()) < 1.3)) etaRegion = true;
	    
	    // Debug
	    // std::cout << "---------------------------" << std::endl;
	    // std::cout << "Vtx= " << ev.nVtx << " pT= " << muTk.Pt() << " eta= " << muTk.Eta() << std::endl;
	    // std::cout << "ID: " << IDName[mid] << std::endl;
	    // std::cout << "Region: " << region[nre] << std::endl;
	    // std::cout << "IsID: " << IsMuonID[mid] << std::endl;

	    if (!IsMuonID[mid]) continue;
	    if (!etaRegion)     continue;

	    if(IsSIGN){	    
	      m_2Dyields["SIGN"]->Fill(mid,nre);
	      m_plots[IDName[mid]][region[nre]]["SIGNVtx"]->Fill(ev.nVtx);
	      m_plots[IDName[mid]][region[nre]]["SIGNPt"] ->Fill(muTk.Pt());
	      if(muTk.Pt() > pTeta) m_plots[IDName[mid]][region[nre]]["SIGNEta"]->Fill(muTk.Eta());
	      
	      if(IsbGEN){
		m_2Dyields["bGEN"]->Fill(mid,nre);
		m_plots[IDName[mid]][region[nre]]["bGENVtx"]->Fill(ev.nVtx);
		m_plots[IDName[mid]][region[nre]]["bGENPt"] ->Fill(muTk.Pt());
		if(muTk.Pt() > pTeta) m_plots[IDName[mid]][region[nre]]["bGENEta"]->Fill(muTk.Eta());	      
	      }
	      else{
		m_2Dyields["WGEN"]->Fill(mid,nre);
		m_plots[IDName[mid]][region[nre]]["WGENVtx"]->Fill(ev.nVtx);
		m_plots[IDName[mid]][region[nre]]["WGENPt"] ->Fill(muTk.Pt());
		if(muTk.Pt() > pTeta) m_plots[IDName[mid]][region[nre]]["WGENEta"]->Fill(muTk.Eta());
	      }
	    } // if(IsSIGN)
	    else{
	      m_2Dyields["FAKE"]->Fill(mid,nre);
	      m_plots[IDName[mid]][region[nre]]["FAKEVtx"]->Fill(ev.nVtx);
	      m_plots[IDName[mid]][region[nre]]["FAKEPt"] ->Fill(muTk.Pt());
	      if(muTk.Pt() > pTeta) m_plots[IDName[mid]][region[nre]]["FAKEEta"]->Fill(muTk.Eta());
	    } // (IsFAKE)
	    
	    m_2Dyields["TOTAL"]->Fill(mid,nre);

	    m_effs[IDName[mid]][region[nre]]["bGENEffVsPt"]->Fill(IsbGEN, muTk.Pt());
	    if(muTk.Pt() > pTeta) m_effs[IDName[mid]][region[nre]]["bGENEffVsEta"]->Fill(IsbGEN, muTk.Eta());
	    m_effs[IDName[mid]][region[nre]]["bGENEffVsVtx"]->Fill(IsbGEN, ev.nVtx);
	    
	    m_effs[IDName[mid]][region[nre]]["WGENEffVsPt"]->Fill(IsWGEN, muTk.Pt());
	    if(muTk.Pt() > pTeta) m_effs[IDName[mid]][region[nre]]["WGENEffVsEta"]->Fill(IsWGEN, muTk.Eta());
	    m_effs[IDName[mid]][region[nre]]["WGENEffVsVtx"]->Fill(IsWGEN, ev.nVtx);
	    
	    m_effs[IDName[mid]][region[nre]]["SIGNEffVsPt"]->Fill(IsSIGN, muTk.Pt());
	    if(muTk.Pt() > pTeta) m_effs[IDName[mid]][region[nre]]["SIGNEffVsEta"]->Fill(IsSIGN, muTk.Eta());
	    m_effs[IDName[mid]][region[nre]]["SIGNEffVsVtx"]->Fill(IsSIGN, ev.nVtx);
	    
	    m_effs[IDName[mid]][region[nre]]["FAKEEffVsPt"]->Fill(IsFAKE, muTk.Pt());
	    if(muTk.Pt() > pTeta) m_effs[IDName[mid]][region[nre]]["FAKEEffVsEta"]->Fill(IsFAKE, muTk.Eta());
	    m_effs[IDName[mid]][region[nre]]["FAKEEffVsVtx"]->Fill(IsFAKE, ev.nVtx);
	    
	  } // for(etaRegion)
	} // for(muonID)
	
    } // if(muon pT eta)
  }// for(Muons)
  

} // void::fill


//Functions

void muon_pog::parseConfig(const std::string configFile, muon_pog::TagAndProbeConfig & tpConfig,
			   std::vector<muon_pog::SampleConfig> & sampleConfigs)
{

  boost::property_tree::ptree pt;
  
  try
    {
      boost::property_tree::ini_parser::read_ini(configFile, pt);
    }
  catch (boost::property_tree::ini_parser::ini_parser_error iniParseErr)
    {
      std::cout << "[TagAndProbeConfig] Can't open : " << iniParseErr.filename()
		<< "\n\tin line : " << iniParseErr.line()
		<< "\n\thas error :" << iniParseErr.message()
		<< std::endl;
      throw std::runtime_error("Bad INI parsing");
    }

  for( auto vt : pt )
    {
      if (vt.first.find("TagAndProbe") != std::string::npos)
	tpConfig = muon_pog::TagAndProbeConfig(vt);
      else
	sampleConfigs.push_back(muon_pog::SampleConfig(vt));
    }
}

void muon_pog::comparisonPlots(std::vector<muon_pog::Plotter> & plotters,
 			       TFile *outFile, TString &  outputDir)
{

  outFile->cd("/");

}

void muon_pog::copyPhp(const TString &  outputDir)
{

  system("cp index.php " + outputDir);
  
  boost::filesystem::directory_iterator dirIt(outputDir.Data());
  boost::filesystem::directory_iterator dirEnd;
  for (;dirIt != dirEnd; ++ dirIt)
    {
      if (boost::filesystem::is_directory(dirIt->status()))
	copyPhp(TString(dirIt->path().string()));
    }

}

Bool_t muon_pog::IsMediumHIP(const muon_pog::Muon & muon){
  Bool_t step0 = muon.isLoose && 
    muon.trkValidHitFrac > 0.49;
  Bool_t step1 = muon.isGlobal && 
    muon.glbNormChi2 < 3.      && 
    muon.trkStaChi2 < 12.      && 
    muon.trkKink < 20          && 
    muon.muSegmComp > 0.303;
  Bool_t step2 = muon.muSegmComp > 0.451;
  
  Bool_t medium = step0  && (step1 || step2);    
  return medium;
}

Bool_t muon_pog::IsSoftHIP(const muon_pog::Muon & muon){
  Bool_t soft = (muon.isTrkMuOST)         && 
    (muon.trkTrackerLayersWithMeas > 5)   &&
    (muon.trkPixelLayersWithMeas   > 0)   &&
    (std::abs(muon.dxyInner)       < 0.3) &&
    (std::abs(muon.dzInner)        < 20.0);

  return soft;
}
