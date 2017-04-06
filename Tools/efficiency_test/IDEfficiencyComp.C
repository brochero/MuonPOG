void SetUpEff(TEfficiency *Eff, TString titleaxis = "", float MarkerSize = 0.02, int hColor = 1){

  Eff->SetTitle(titleaxis);
  Eff->SetMarkerStyle(20);
  Eff->SetMarkerSize(MarkerSize);
  Eff->SetMarkerColor(hColor);
  Eff->SetLineColor(hColor);

}


void IDEfficiencyCompE(TString plot        = "EffVsPt",     
		       TString region      = "Full",   
		       TString ID          = "TRK",   
		       TString muon        = "WGEN",   
		       TString input_title = "Histogram",   
		       TString x_label     = "",       
		       TString y_label     = "", 
		       float y_min         = -1,
		       float y_max         = -1,
		       bool    LogScale    = false){
  
  int hcol[7];
  hcol[0] = (1);
  hcol[1] = (2);
  hcol[2] = (kGreen+2);
  hcol[3] = (4);
  hcol[4] = (kOrange+2);
  hcol[5] = (kViolet+2);
  hcol[6] = (kAzure+7);
  hcol[7] = (kYellow+1);

  TString DrawHistoOption = "HIST PE"; // H=Histogram, E=Error

  TString fname   = "histos_FullHIPStudiesOverlapingNewBinKinVarSumw2";
  TString dirname = "MuonResults/" + fname;  

  TString inputfile0 = dirname + ".root" ;  
 
  TFile *hfile0 = NULL;
  hfile0 = TFile::Open(inputfile0);
  cout << "file loaded: " << inputfile0 << endl;

  TString dirfigname_pdf = "figures_" + fname + "/";
  // make a dir if it does not exist!!
  gSystem->mkdir(dirfigname_pdf, kTRUE);

  TEfficiency *hEff_0, *hEff_1, *hEff_2, *hEff_3, *hEff_4, *hEff_5, *hEff_6, *hEff_7;

  cout << "Loading histograms... " << endl;

  hEff_0  = (TEfficiency*)hfile0->Get("TTbarHIPNoSIMU-NoMITI/efficiencies/" + ID + "/" + region + "/" + plot + "_" + muon + "_" + ID + "_" + region);
  hEff_1  = (TEfficiency*)hfile0->Get("TTbarHIPNoSIMU-MITI/efficiencies/" + ID + "/" + region + "/" + plot + "_" + muon + "_" + ID + "_" + region);
  hEff_2  = (TEfficiency*)hfile0->Get("TTbarHIPSIMU-NoMITI_0p6/efficiencies/" + ID + "/" + region + "/" + plot + "_" + muon + "_" + ID + "_" + region);
  hEff_3  = (TEfficiency*)hfile0->Get("TTbarHIPSIMU-MITI_0p6/efficiencies/" + ID + "/" + region + "/" + plot + "_" + muon + "_" + ID + "_" + region);
  hEff_4  = (TEfficiency*)hfile0->Get("TTbarHIPSIMU-NoMITI_0p8/efficiencies/" + ID + "/" + region + "/" + plot + "_" + muon + "_" + ID + "_" + region);
  hEff_5  = (TEfficiency*)hfile0->Get("TTbarHIPSIMU-MITI_0p8/efficiencies/" + ID + "/" + region + "/" + plot + "_" + muon + "_" + ID + "_" + region);
  hEff_6  = (TEfficiency*)hfile0->Get("TTbarHIPSIMU-NoMITI_1p2/efficiencies/" + ID + "/" + region + "/" + plot + "_" + muon + "_" + ID + "_" + region);
  hEff_7  = (TEfficiency*)hfile0->Get("TTbarHIPSIMU-MITI_1p2/efficiencies/" + ID + "/" + region + "/" + plot + "_" + muon + "_" + ID + "_" + region);


  TCanvas *cPlots;//histos
  cPlots = new TCanvas("cPlots" ,"Plots");
  cPlots->Divide(2,1);

  TPad *cPad[2];
  cPad[0] = (TPad*)cPlots->GetPad(1);
  cPad[0]->SetPad(0.05, 0.01, 0.92, 0.99);
  if(LogScale) cPad[0]->SetLogy();

  cPad[1] = (TPad*)cPlots->GetPad(2);
  cPad[1]->SetPad(0.84, 0.01, 0.99, 0.99);
  

  cPad[0]->cd();

  TString titleaxis = ";" + x_label + ";" + y_label ;

  float MarkerSize = 0.7;

  SetUpEff(hEff_0, titleaxis, MarkerSize, hcol[0]);
  SetUpEff(hEff_1, titleaxis, MarkerSize, hcol[1]);
  SetUpEff(hEff_2, titleaxis, MarkerSize, hcol[2]);
  SetUpEff(hEff_3, titleaxis, MarkerSize, hcol[3]);
  SetUpEff(hEff_4, titleaxis, MarkerSize, hcol[4]);
  SetUpEff(hEff_5, titleaxis, MarkerSize, hcol[5]);
  SetUpEff(hEff_6, titleaxis, MarkerSize, hcol[6]);
  SetUpEff(hEff_7, titleaxis, MarkerSize, hcol[7]);



  hEff_0->Draw(); 
  //hEff_1->Draw(); 
  //hEff_1->Draw("SAME"); 
  //hEff_2->Draw("SAME"); 
  //hEff_3->Draw("SAME"); 
  hEff_4->Draw("SAME"); 
  hEff_5->Draw("SAME"); 
  //hEff_6->Draw("SAME"); 
  //hEff_7->Draw("SAME"); 



  gPad->Update();

  if(y_max > 0.0) hEff_0->GetPaintedGraph()->SetMaximum(y_max);
  if(y_min > 0.0) hEff_0->GetPaintedGraph()->SetMinimum(y_min);


  //if(y_max > 0.0) hEff_1->GetPaintedGraph()->SetMaximum(y_max);
  //if(y_min > 0.0) hEff_1->GetPaintedGraph()->SetMinimum(y_min);

  /***********************
      CMS Legend
  **********************/
  TLatex *titlePr;
  titlePr  = new TLatex(-20.,50.,"CMS Preliminary");
  titlePr->SetNDC();
  titlePr->SetTextAlign(12);
  titlePr->SetX(0.15);
  titlePr->SetY(0.93);
  titlePr->SetTextColor(2);
  titlePr->SetTextFont(42);
  titlePr->SetTextSize(0.05);
  titlePr->SetTextSizePixels(24);
  titlePr->Draw("SAME");
  
  TLatex *title;
  title  = new TLatex(-20.,50.,input_title);
  title->SetNDC();
  title->SetTextAlign(12);
  title->SetX(0.15);
  title->SetY(0.83);
  title->SetTextFont(42);
  title->SetTextSize(0.05);
  title->SetTextSizePixels(24);
  title->Draw("SAME");


  cPad[1]->cd();

  TLegend *leg1;
  //leg1 = new TLegend(0.70,.74,0.89,0.89);
  leg1 = new TLegend(0.01,0.1,0.9,0.9);
  leg1->SetLineColor(1);
  leg1->SetTextSize(0.065);

  leg1->AddEntry((TObject*)0,"HIP No SIM","");
  leg1->AddEntry(hEff_0, "No MITI","PL");
  //leg1->AddEntry(hEff_1, "MITI","PL");
  //leg1->AddEntry((TObject*)0,"HIP-SIM 0.6e34","");
  //leg1->AddEntry(hEff_2, "No MITI","PL");
  //leg1->AddEntry(hEff_3, "MITI","PL");
  leg1->AddEntry((TObject*)0,"HIP-SIM 0.8e34","");
  leg1->AddEntry(hEff_4, "No MITI","PL");
  leg1->AddEntry(hEff_5, "MITI","PL");
  // leg1->AddEntry((TObject*)0,"HIP-SIM 1.2e34","");
  // leg1->AddEntry(hEff_6, "No MITI","PL");
  // leg1->AddEntry(hEff_7, "MITI","PL");

  TString Scenario = "NOMITI_Plus_0p8";

  leg1->Draw();

  TString hname= hEff_0->GetName(); 
  if(!LogScale) cPlots->SaveAs(dirfigname_pdf + hname + Scenario + ".pdf");
  else cPlots->SaveAs(dirfigname_pdf + hname + "_" + Scenario + "log.pdf");

  cPlots->Clear();

}

void IDEfficiencyCompA(TString hnID = "", TString hnMu = ""){
  TString LongName = "";
  if(hnMu == "SIGN") LongName = "Prompt ";
  if(hnMu == "FAKE") LongName = "Fake ";

  IDEfficiencyCompE("EffVsPt", "Full", hnID, hnMu, LongName + hnID + " muons ", "pT[GeV]", "Purity", -1, -1, false);
  IDEfficiencyCompE("EffVsPt", "Barrel", hnID, hnMu, LongName + hnID + " muons (Barrel |#eta|<1.2)", "pT[GeV]", "Purity", -1, -1, false);
  IDEfficiencyCompE("EffVsPt", "Endcap", hnID, hnMu, LongName + hnID + " muons (Endcap |#eta|>1.2)", "pT[GeV]", "Purity", -1, -1, false);
  IDEfficiencyCompE("EffVsPt", "Overlap", hnID, hnMu, LongName + hnID + " muons (Overlap 0.9<|#eta|<1.3)", "pT[GeV]", "Purity", -1, -1, false);

  IDEfficiencyCompE("EffVsVtx", "Full", hnID, hnMu, LongName + hnID + " muons ", "Vertices", "Purity", -1, -1, false);
  IDEfficiencyCompE("EffVsVtx", "Barrel", hnID, hnMu, LongName + hnID + " muons (Barrel |#eta|<1.2)", "Verties", "Purity", -1, -1, false);
  IDEfficiencyCompE("EffVsVtx", "Endcap", hnID, hnMu, LongName + hnID + " muons (Endcap |#eta|>1.2)", "Vertices", "Purity", -1, -1, false);
  IDEfficiencyCompE("EffVsVtx", "Overlap", hnID, hnMu, LongName + hnID + " muons (Overlap 0.9<|#eta|<1.3)", "Vertices", "Purity", -1, -1, false);

  if(hnMu == "SIGN") IDEfficiencyCompE("EffVsEta", "Full", hnID, hnMu, LongName + hnID + " muons ", "#eta", "Purity", 0.4, 0.9, false);
  else if(hnMu == "FAKE" && hnID == "TRK") IDEfficiencyCompE("EffVsEta", "Full", hnID, hnMu, LongName + hnID + " muons ", "#eta", "Purity", 0.0, 0.4, false);
  else if(hnMu == "FAKE" && hnID == "GLB") IDEfficiencyCompE("EffVsEta", "Full", hnID, hnMu, LongName + hnID + " muons ", "#eta", "Purity", 0.0, 0.15, false);
  else if(hnMu == "FAKE" && hnID == "SOFT") IDEfficiencyCompE("EffVsEta", "Full", hnID, hnMu, LongName + hnID + " muons ", "#eta", "Purity", 0.0, 0.15, false);
  //else if(hnMu == "FAKE" && hnID == "LOOSE") 
  else {IDEfficiencyCompE("EffVsEta", "Full", hnID, hnMu, LongName + hnID + " muons ", "#eta", "Purity", 0.0, 0.1, false);}

}


void IDEfficiencyComp(){

  IDEfficiencyCompA("TRK","SIGN");
  IDEfficiencyCompA("TRK","FAKE");

  IDEfficiencyCompA("GLB","SIGN");
  IDEfficiencyCompA("GLB","FAKE");

  IDEfficiencyCompA("LOOSE","SIGN");
  IDEfficiencyCompA("LOOSE","FAKE");

  IDEfficiencyCompA("MEDIUM","SIGN");
  IDEfficiencyCompA("MEDIUM","FAKE");

  IDEfficiencyCompA("TIGHT","SIGN");
  IDEfficiencyCompA("TIGHT","FAKE");

  IDEfficiencyCompA("SOFT","SIGN");
  IDEfficiencyCompA("SOFT","FAKE");

  IDEfficiencyCompA("HIGHpT","SIGN");
  IDEfficiencyCompA("HIGHpT","FAKE");

}
