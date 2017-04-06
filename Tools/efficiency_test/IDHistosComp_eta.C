
void SetUpRatio(TH1F *Ratio, TString x_label = ""){

  Ratio->SetMarkerStyle(20);
  Ratio->SetMarkerSize(0.5);
  Ratio->SetMarkerColor(1);
  Ratio->SetLineWidth(2);
  Ratio->SetLineStyle(1);
  Ratio->SetMaximum(2.0);
  Ratio->SetMinimum(0);
  Ratio->SetTitle("");

  Ratio->GetYaxis()->SetTitle("");
  Ratio->GetYaxis()->CenterTitle();
  Ratio->GetYaxis()->SetTitleOffset(0.45);
  Ratio->GetYaxis()->SetTitleSize(0.16);
  Ratio->GetYaxis()->SetLabelSize(0.1);
  Ratio->GetYaxis()->SetNdivisions(506);
  //Ratio->GetXaxis()->SetNdivisions(509); 
  Ratio->GetXaxis()->SetTitleOffset(0.6);
  Ratio->GetXaxis()->SetLabelSize(0.1);
  Ratio->GetXaxis()->SetTitleSize(0.16);
  Ratio->GetXaxis()->SetTitle(x_label);
 
  Ratio->SetMinimum(0.6);
  Ratio->SetMaximum(1.4);

}

void SetUpHisto(TH1 *histo, TString titleaxis = "", float MarkerSize = 0.02, int hColor = 1){

  histo->SetTitle(titleaxis);
  histo->SetMarkerStyle(20);
  histo->SetMarkerSize(MarkerSize);
  histo->SetMarkerColor(hColor);
  histo->SetLineColor(hColor);
  histo->GetXaxis()->SetLabelSize(0.);

}

void IDHistosCompF(TString plot        = "Eta",     
		   TString region      = "Full",   
		   TString ID          = "TRK",   
		   TString muon        = "FAKE",   
		   TString input_title = "Histogram",   
		   TString x_label     = "",       
		   TString y_label     = "", 
		   float y_min         = -1,
		   float y_max         = -1,
		   bool    LogScale    = false,
		   bool    Norm        = false,
		   bool    Rebinning   = false){
  
  int hcol[8];
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

  TH1F *h0, *h1, *h2, *h3, *h4, *h5, *h6, *h7;
  TH1 *h_0, *h_1, *h_2, *h_3, *h_4, *h_5, *h_6, *h_7;
  TH1F *Ratio_01, *Ratio_02, *Ratio_03, *Ratio_04,*Ratio_05, *Ratio_06, *Ratio_07;

  cout << "Loading histograms... " << endl;

  h0  = (TH1F*)hfile0->Get("TTbarHIPNoSIMU-NoMITI/kinematical_variables/" + ID + "/" + region + "/" + plot + "_" + muon + "_" + ID + "_" + region);
  h1  = (TH1F*)hfile0->Get("TTbarHIPNoSIMU-MITI/kinematical_variables/" + ID + "/" + region + "/" + plot + "_" + muon + "_" + ID + "_" + region);
  h2  = (TH1F*)hfile0->Get("TTbarHIPSIMU-NoMITI_0p6/kinematical_variables/" + ID + "/" + region + "/" + plot + "_" + muon + "_" + ID + "_" + region);
  h3  = (TH1F*)hfile0->Get("TTbarHIPSIMU-MITI_0p6/kinematical_variables/" + ID + "/" + region + "/" + plot + "_" + muon + "_" + ID + "_" + region);
  h4  = (TH1F*)hfile0->Get("TTbarHIPSIMU-NoMITI_0p8/kinematical_variables/" + ID + "/" + region + "/" + plot + "_" + muon + "_" + ID + "_" + region);
  h5  = (TH1F*)hfile0->Get("TTbarHIPSIMU-MITI_0p8/kinematical_variables/" + ID + "/" + region + "/" + plot + "_" + muon + "_" + ID + "_" + region);
  h6  = (TH1F*)hfile0->Get("TTbarHIPSIMU-NoMITI_1p2/kinematical_variables/" + ID + "/" + region + "/" + plot + "_" + muon + "_" + ID + "_" + region);
  h7  = (TH1F*)hfile0->Get("TTbarHIPSIMU-MITI_1p2/kinematical_variables/" + ID + "/" + region + "/" + plot + "_" + muon + "_" + ID + "_" + region);
  
 
  h_0 = (TH1*) h0;
  h_1 = (TH1*) h1;
  h_2 = (TH1*) h2;
  h_3 = (TH1*) h3;
  h_4 = (TH1*) h4;
  h_5 = (TH1*) h5;
  h_6 = (TH1*) h6;
  h_7 = (TH1*) h7;

  if(Rebinning){
    cout << "Rebinning histogram............" << endl;
    // int NewBins=17;
    // Double_t xbins[18] = {-2.4,-2.1,-1.6,-1.2,-1.0,-0.9,-0.6,-0.3,-0.2,                                                                                                                
    // 			 0.2, 0.3, 0.6, 0.9, 1.0, 1.2, 1.6, 2.1, 2.4};
    // int NewBins=24;
    // Double_t xbins[25] = {-2.4,-2.2,-2.0,-1.8,-1.6,-1.4,-1.2,-1.0,-0.8,-0.6,-0.4,-0.2,
    // 			  0.0,
    // 			   0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4};

    // int NewBins=12;
    // Double_t xbins[13] = {-2.4,-2.0,-1.6,-1.2,-0.8,-0.4,
    // 			  0.0,
    // 			   0.4, 0.8, 1.2, 1.6, 2.0, 2.4};

    int NewBins=8;
    Double_t xbins[9] = {-2.4,-1.8,-1.2,-0.6,
    			  0.0,
    			   0.6, 1.2, 1.8, 2.4};

    h_0 = h0->Rebin(NewBins,"h_0", xbins); 
    h_0->SetName(h0->GetName());
    h_1 = h1->Rebin(NewBins,"h_1", xbins); 
    h_2 = h2->Rebin(NewBins,"h_2", xbins); 
    h_3 = h3->Rebin(NewBins,"h_3", xbins); 
    h_4 = h4->Rebin(NewBins,"h_4", xbins); 
    h_5 = h5->Rebin(NewBins,"h_5", xbins); 
    h_6 = h6->Rebin(NewBins,"h_6", xbins); 
    h_7 = h7->Rebin(NewBins,"h_7", xbins); 
  }

  if(Norm){
    if(h_0->Integral() != 0.0) h_0->Scale(1.0/(h_0->Integral()));
    if(h_1->Integral() != 0.0) h_1->Scale(1.0/(h_1->Integral()));
    if(h_2->Integral() != 0.0) h_2->Scale(1.0/(h_2->Integral()));
    if(h_3->Integral() != 0.0) h_3->Scale(1.0/(h_3->Integral()));
    if(h_4->Integral() != 0.0) h_4->Scale(1.0/(h_4->Integral()));
    if(h_5->Integral() != 0.0) h_5->Scale(1.0/(h_5->Integral()));
    if(h_6->Integral() != 0.0) h_6->Scale(1.0/(h_6->Integral()));
    if(h_7->Integral() != 0.0) h_7->Scale(1.0/(h_7->Integral()));
  }

  TCanvas *cPlots;//histos
  cPlots = new TCanvas("cPlots" ,"Plots");
  cPlots->Divide(2,1);

  TLegend *leg1;
  leg1 = new TLegend(0.01,0.1,0.9,0.9);
  leg1->SetLineColor(1);
  leg1->SetTextSize(0.065);

  TPad *cPad[4];
  cPad[0] = (TPad*)cPlots->GetPad(1);
  cPad[0]->SetPad(0.05, 0.01, 0.92, 0.99);
  cPad[0]->Divide(1,2);
  cPad[2] = (TPad*)cPad[0]->GetPad(1);
  cPad[2]->SetPad(0.04, 0.23, 0.92, 0.99);
  cPad[2]->SetLeftMargin(0.12);
  //cPad[2]->SetGridx();
  //cPad[2]->SetGridy();
  cPad[3] = (TPad*)cPad[0]->GetPad(2);
  cPad[3]->SetPad(0.04, 0.05, 0.92, 0.3);
  cPad[3]->SetBottomMargin(0.25);
  cPad[3]->SetTopMargin(0.02);
  cPad[3]->SetLeftMargin(0.12);
  cPad[3]->SetGridx();
  cPad[3]->SetGridy();

  if(LogScale) cPad[2]->SetLogy();

  cPad[1] = (TPad*)cPlots->GetPad(2);
  cPad[1]->SetPad(0.84, 0.01, 0.99, 0.99);
  

  cPad[2]->cd();

  TString titleaxis = ";" + x_label + ";" + y_label ;

  float MarkerSize = 0.7;

  SetUpHisto(h_0, titleaxis, MarkerSize, hcol[0]);
  SetUpHisto(h_1, titleaxis, MarkerSize, hcol[1]);
  SetUpHisto(h_2, titleaxis, MarkerSize, hcol[2]);
  SetUpHisto(h_3, titleaxis, MarkerSize, hcol[3]);
  SetUpHisto(h_4, titleaxis, MarkerSize, hcol[4]);
  SetUpHisto(h_5, titleaxis, MarkerSize, hcol[5]);
  SetUpHisto(h_6, titleaxis, MarkerSize, hcol[6]);
  SetUpHisto(h_7, titleaxis, MarkerSize, hcol[7]);


  Ratio_01 = (TH1F*)h_1->Clone();
  Ratio_01->SetLineColor(hcol[1]);
  Ratio_02 = (TH1F*)h_2->Clone();
  Ratio_02->SetLineColor(hcol[2]);
  Ratio_03 = (TH1F*)h_3->Clone();
  Ratio_03->SetLineColor(hcol[3]);
  Ratio_04 = (TH1F*)h_4->Clone();
  Ratio_04->SetLineColor(hcol[4]);
  Ratio_05 = (TH1F*)h_5->Clone();
  Ratio_05->SetLineColor(hcol[5]);
  Ratio_06 = (TH1F*)h_6->Clone();
  Ratio_06->SetLineColor(hcol[6]);
  Ratio_07 = (TH1F*)h_7->Clone();
  Ratio_07->SetLineColor(hcol[7]);


  Ratio_01->Divide(h_0);
  SetUpRatio(Ratio_01, x_label);
  Ratio_02->Divide(h_0);
  SetUpRatio(Ratio_02, x_label);
  Ratio_03->Divide(h_0);
  SetUpRatio(Ratio_03, x_label);
  Ratio_04->Divide(h_0);
  SetUpRatio(Ratio_04, x_label);
  Ratio_05->Divide(h_0);
  SetUpRatio(Ratio_05, x_label);
  Ratio_06->Divide(h_0);
  SetUpRatio(Ratio_06, x_label);
  Ratio_07->Divide(h_0);
  SetUpRatio(Ratio_07, x_label);

  h_0->Draw(); 
  //h_2->Draw(); 
  //h_1->Draw("SAME"); 
  h_2->Draw("SAME"); 
  h_3->Draw("SAME"); 
  //h_4->Draw("SAME"); 
  //h_5->Draw("SAME"); 
  //h_6->Draw("SAME"); 
  //h_7->Draw("SAME"); 

  cPad[3]->cd();
  //Ratio_01->Draw();
  Ratio_02->Draw();
  Ratio_03->Draw("SAME");
  //Ratio_04->Draw();
  //Ratio_05->Draw("SAME");
  //Ratio_06->Draw();
  //Ratio_07->Draw("SAME");

  cPad[2]->cd();

  TString Scenario = "NOMITI_Plus_0p6";

  leg1->AddEntry((TObject*)0,"HIP No SIM","");
  leg1->AddEntry(h_0, "No MITI","PL");
  //leg1->AddEntry(h_1, "MITI","PL");
  leg1->AddEntry((TObject*)0,"HIP-SIM 0.6e34","");
  leg1->AddEntry(h_2, "No MITI","PL");
  leg1->AddEntry(h_3, "MITI","PL");
  //leg1->AddEntry((TObject*)0,"HIP-SIM 0.8e34","");
  //leg1->AddEntry(h_4, "No MITI","PL");
  //leg1->AddEntry(h_5, "MITI","PL");
  //leg1->AddEntry((TObject*)0,"HIP-SIM 1.2e34","");
  //leg1->AddEntry(h_6, "No MITI","PL");
  //leg1->AddEntry(h_7, "MITI","PL");



  gPad->Update();
  if(y_max > 0.0) h_0->SetMaximum(y_max);
  else{h_0->SetMaximum( h_0->GetMaximum()*1.3);}
  if(y_min > 0.0) h_0->SetMinimum(y_min);
  // if(y_max > 0.0) h_2->SetMaximum(y_max);
  // if(y_min > 0.0) h_2->SetMinimum(y_min);

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

  leg1->Draw();

  TString hname= h_0->GetName(); 

  if(LogScale)  Scenario = Scenario + "Log";
  if(Norm)      Scenario = Scenario + "Norm";
  if(Rebinning) Scenario = Scenario + "Rebin";
  
  cPlots->SaveAs(dirfigname_pdf + hname + "_" + Scenario + ".pdf");
  cPlots->Clear();
}

void IDHistosCompR(TString hnID = "", TString hnMu = ""){
  TString LongName = "";
  if(hnMu == "SIGN") LongName = "Prompt ";
  if(hnMu == "FAKE") LongName = "Fake ";

  IDHistosCompF("Eta","Full", hnID, hnMu, LongName + hnID + " muons","#eta","N.Muons",-1,-1,false,false,false);
  IDHistosCompF("Eta","Full", hnID, hnMu, LongName + hnID + " muons","#eta","N.Muons",-1,-1,false,false,true);
  IDHistosCompF("Eta","Full", hnID, hnMu, LongName + hnID + " muons","#eta","N.Muons",-1,-1,false,true, true);
  IDHistosCompF("Eta","Full", hnID, hnMu, LongName + hnID + " muons","#eta","N.Muons",-1,-1,false,true, false);


}

void IDHistosComp_eta(){
  IDHistosCompR ("TRK","SIGN");
  IDHistosCompR ("TRK","FAKE");

  IDHistosCompR ("GLB","SIGN");
  IDHistosCompR ("GLB","FAKE");

  IDHistosCompR ("LOOSE","SIGN");
  IDHistosCompR ("LOOSE","FAKE");

  IDHistosCompR ("MEDIUM","SIGN");
  IDHistosCompR ("MEDIUM","FAKE");

  IDHistosCompR ("TIGHT","SIGN");
  IDHistosCompR ("TIGHT","FAKE");

  IDHistosCompR ("SOFT","SIGN");
  IDHistosCompR ("SOFT","FAKE");

  IDHistosCompR ("HIGHpT","SIGN");
  IDHistosCompR ("HIGHpT","FAKE");

}

