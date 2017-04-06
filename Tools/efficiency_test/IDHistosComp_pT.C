
void IDHistosCompF(TString plot        = "Pt",     
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
  TH1F *Ratio_01, *Ratio_02;

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
    int NewBins=42;
    Double_t xbins[43] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,28,30,32,34,36,38,40,44,48,52,56,60,65,70,80,100};
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

  float MarkerSize = 0.4;
  h_0->SetTitle(titleaxis);
  h_0->SetMarkerStyle(20);
  h_0->SetMarkerSize(MarkerSize);
  h_0->SetMarkerColor(hcol[0]);
  h_0->SetLineColor(hcol[0]);
  h_0->GetXaxis()->SetLabelSize(0.);

  h_1->SetMarkerStyle(20);
  h_1->SetTitle(titleaxis);
  h_1->SetMarkerColor(hcol[1]);
  h_1->SetMarkerSize(MarkerSize);
  h_1->SetLineColor(hcol[1]);
  h_1->GetXaxis()->SetLabelSize(0.);

  h_2->SetTitle(titleaxis);
  h_2->SetMarkerStyle(21);
  h_2->SetMarkerColor(hcol[2]);
  h_2->SetMarkerSize(MarkerSize);
  h_2->SetLineColor(hcol[2]);
  h_2->GetXaxis()->SetLabelSize(0.);

  h_3->SetTitle(titleaxis);
  h_3->SetMarkerStyle(21);
  h_3->SetMarkerColor(hcol[3]);
  h_3->SetMarkerSize(MarkerSize);
  h_3->SetLineColor(hcol[3]);
  h_3->GetXaxis()->SetLabelSize(0.);

  h_4->SetTitle(titleaxis);
  h_4->SetMarkerStyle(22);
  h_4->SetMarkerColor(hcol[4]);
  h_4->SetMarkerSize(MarkerSize);
  h_4->SetLineColor(hcol[4]);
  h_4->GetXaxis()->SetLabelSize(0.);

  h_5->SetTitle(titleaxis);
  h_5->SetMarkerStyle(22);
  h_5->SetMarkerColor(hcol[5]);
  h_5->SetMarkerSize(MarkerSize);
  h_5->SetLineColor(hcol[5]);
  h_5->GetXaxis()->SetLabelSize(0.);

  h_6->SetTitle(titleaxis);
  h_6->SetMarkerStyle(23);
  h_6->SetMarkerColor(hcol[6]);
  h_6->SetMarkerSize(MarkerSize);
  h_6->SetLineColor(hcol[6]);
  h_6->GetXaxis()->SetLabelSize(0.);


  h_7->SetTitle(titleaxis);
  h_7->SetMarkerStyle(23);
  h_7->SetMarkerColor(hcol[7]);
  h_7->SetMarkerSize(MarkerSize);
  h_7->SetLineColor(hcol[7]);
  h_7->GetXaxis()->SetLabelSize(0.);


  h_0->Draw(); 
  //h_2->Draw(); 
  //h_1->Draw("SAME"); 
  //h_2->Draw("SAME"); 
  //h_3->Draw("SAME"); 
  //h_4->Draw("SAME"); 
  //h_5->Draw("SAME"); 
  h_6->Draw("SAME"); 
  h_7->Draw("SAME"); 

  Ratio_01 = (TH1F*)h_6->Clone();
  Ratio_01->SetLineColor(hcol[6]);
  Ratio_02 = (TH1F*)h_7->Clone();
  Ratio_02->SetLineColor(hcol[7]);

  Ratio_01->Divide(h_0);
  Ratio_02->Divide(h_0);

  TString Scenario = "NOMITI_Plus_1p2";



  gPad->Update();
  if(y_max > 0.0) h_0->SetMaximum(y_max);
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

  cPad[3]->cd();

  Ratio_01->SetMarkerStyle(20);
  Ratio_01->SetMarkerSize(0.5);
  Ratio_01->SetMarkerColor(1);
  Ratio_01->SetLineWidth(2);
  Ratio_01->SetLineStyle(1);
  Ratio_01->SetMaximum(2.0);
  Ratio_01->SetMinimum(0);
  Ratio_01->SetTitle("");

  Ratio_01->GetYaxis()->SetTitle("");
  Ratio_01->GetYaxis()->CenterTitle();
  Ratio_01->GetYaxis()->SetTitleOffset(0.45);
  Ratio_01->GetYaxis()->SetTitleSize(0.16);
  Ratio_01->GetYaxis()->SetLabelSize(0.1);
  Ratio_01->GetYaxis()->SetNdivisions(506);
  //Ratio_01->GetXaxis()->SetNdivisions(509); 
  Ratio_01->GetXaxis()->SetTitleOffset(0.6);
  Ratio_01->GetXaxis()->SetLabelSize(0.1);
  Ratio_01->GetXaxis()->SetTitleSize(0.16);
  Ratio_01->GetXaxis()->SetTitle(x_label);
 
  Ratio_01->SetMinimum(0.5);
  Ratio_01->SetMaximum(2.0);

  Ratio_01->Draw();
  Ratio_02->Draw("SAME");

  cPad[1]->cd();

  TLegend *leg1;
  //leg1 = new TLegend(0.70,.74,0.89,0.89);
  leg1 = new TLegend(0.01,0.1,0.9,0.9);
  leg1->SetLineColor(1);
  leg1->SetTextSize(0.065);

  leg1->AddEntry((TObject*)0,"HIP No SIM","");
  leg1->AddEntry(h_0, "No MITI","PL");
  //leg1->AddEntry(h_1, "MITI","PL");
  //leg1->AddEntry((TObject*)0,"HIP-SIM 0.6e34","");
  //leg1->AddEntry(h_2, "No MITI","PL");
  //leg1->AddEntry(h_3, "MITI","PL");
  //leg1->AddEntry((TObject*)0,"HIP-SIM 0.8e34","");
  //leg1->AddEntry(h_4, "No MITI","PL");
  //leg1->AddEntry(h_5, "MITI","PL");
  leg1->AddEntry((TObject*)0,"HIP-SIM 1.2e34","");
  leg1->AddEntry(h_6, "No MITI","PL");
  leg1->AddEntry(h_7, "MITI","PL");

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

  IDHistosCompF("Pt","Full", hnID, hnMu, LongName + hnID + " muons","pT[GeV]","N.Muons",-1,-1,false,false,false);
  IDHistosCompF("Pt","Full", hnID, hnMu, LongName + hnID + " muons","pT[GeV]","N.Muons",-1,-1,false,false,true);
  IDHistosCompF("Pt","Full", hnID, hnMu, LongName + hnID + " muons","pT[GeV]","N.Muons",-1,-1,false,true, true);
  IDHistosCompF("Pt","Full", hnID, hnMu, LongName + hnID + " muons","pT[GeV]","N.Muons",-1,-1,false,true, false);
  IDHistosCompF("Pt","Full", hnID, hnMu, LongName + hnID + " muons","pT[GeV]","N.Muons",-1,-1,true,false,false);
  IDHistosCompF("Pt","Full", hnID, hnMu, LongName + hnID + " muons","pT[GeV]","N.Muons",-1,-1,true,false,true);
  IDHistosCompF("Pt","Full", hnID, hnMu, LongName + hnID + " muons","pT[GeV]","N.Muons",-1,-1,true,true, true);
  IDHistosCompF("Pt","Full", hnID, hnMu, LongName + hnID + " muons","pT[GeV]","N.Muons",-1,-1,true,true, false);

  IDHistosCompF("Pt","Barrel", hnID, hnMu, LongName + hnID + " muons (Barrel |#eta|<1.2)","pT[GeV]","N.Muons",-1,-1,false,false,false);
  IDHistosCompF("Pt","Barrel", hnID, hnMu, LongName + hnID + " muons (Barrel |#eta|<1.2)","pT[GeV]","N.Muons",-1,-1,false,false,true);
  IDHistosCompF("Pt","Barrel", hnID, hnMu, LongName + hnID + " muons (Barrel |#eta|<1.2)","pT[GeV]","N.Muons",-1,-1,false,true, true);
  IDHistosCompF("Pt","Barrel", hnID, hnMu, LongName + hnID + " muons (Barrel |#eta|<1.2)","pT[GeV]","N.Muons",-1,-1,false,true, false);
  IDHistosCompF("Pt","Barrel", hnID, hnMu, LongName + hnID + " muons (Barrel |#eta|<1.2)","pT[GeV]","N.Muons",-1,-1,true,false,false);
  IDHistosCompF("Pt","Barrel", hnID, hnMu, LongName + hnID + " muons (Barrel |#eta|<1.2)","pT[GeV]","N.Muons",-1,-1,true,false,true);
  IDHistosCompF("Pt","Barrel", hnID, hnMu, LongName + hnID + " muons (Barrel |#eta|<1.2)","pT[GeV]","N.Muons",-1,-1,true,true, true);
  IDHistosCompF("Pt","Barrel", hnID, hnMu, LongName + hnID + " muons (Barrel |#eta|<1.2)","pT[GeV]","N.Muons",-1,-1,true,true, false);

  IDHistosCompF("Pt","Endcap", hnID, hnMu, LongName + hnID + " muons (Endcap |#eta|>1.2)","pT[GeV]","N.Muons",-1,-1,false,false,false);
  IDHistosCompF("Pt","Endcap", hnID, hnMu, LongName + hnID + " muons (Endcap |#eta|>1.2)","pT[GeV]","N.Muons",-1,-1,false,false,true);
  IDHistosCompF("Pt","Endcap", hnID, hnMu, LongName + hnID + " muons (Endcap |#eta|>1.2)","pT[GeV]","N.Muons",-1,-1,false,true, true);
  IDHistosCompF("Pt","Endcap", hnID, hnMu, LongName + hnID + " muons (Endcap |#eta|>1.2)","pT[GeV]","N.Muons",-1,-1,false,true, false);
  IDHistosCompF("Pt","Endcap", hnID, hnMu, LongName + hnID + " muons (Endcap |#eta|>1.2)","pT[GeV]","N.Muons",-1,-1,true,false,false);
  IDHistosCompF("Pt","Endcap", hnID, hnMu, LongName + hnID + " muons (Endcap |#eta|>1.2)","pT[GeV]","N.Muons",-1,-1,true,false,true);
  IDHistosCompF("Pt","Endcap", hnID, hnMu, LongName + hnID + " muons (Endcap |#eta|>1.2)","pT[GeV]","N.Muons",-1,-1,true,true, true);
  IDHistosCompF("Pt","Endcap", hnID, hnMu, LongName + hnID + " muons (Endcap |#eta|>1.2)","pT[GeV]","N.Muons",-1,-1,true,true, false);

  IDHistosCompF("Pt","Overlap", hnID, hnMu, LongName + hnID + " muons (Overlap 0.9<|#eta|<1.3)","pT[GeV]","N.Muons",-1,-1,false,false,false);
  IDHistosCompF("Pt","Overlap", hnID, hnMu, LongName + hnID + " muons (Overlap 0.9<|#eta|<1.3)","pT[GeV]","N.Muons",-1,-1,false,false,true);
  IDHistosCompF("Pt","Overlap", hnID, hnMu, LongName + hnID + " muons (Overlap 0.9<|#eta|<1.3)","pT[GeV]","N.Muons",-1,-1,false,true, true);
  IDHistosCompF("Pt","Overlap", hnID, hnMu, LongName + hnID + " muons (Overlap 0.9<|#eta|<1.3)","pT[GeV]","N.Muons",-1,-1,false,true, false);
  IDHistosCompF("Pt","Overlap", hnID, hnMu, LongName + hnID + " muons (Overlap 0.9<|#eta|<1.3)","pT[GeV]","N.Muons",-1,-1,true,false,false);
  IDHistosCompF("Pt","Overlap", hnID, hnMu, LongName + hnID + " muons (Overlap 0.9<|#eta|<1.3)","pT[GeV]","N.Muons",-1,-1,true,false,true);
  IDHistosCompF("Pt","Overlap", hnID, hnMu, LongName + hnID + " muons (Overlap 0.9<|#eta|<1.3)","pT[GeV]","N.Muons",-1,-1,true,true, true);
  IDHistosCompF("Pt","Overlap", hnID, hnMu, LongName + hnID + " muons (Overlap 0.9<|#eta|<1.3)","pT[GeV]","N.Muons",-1,-1,true,true, false);




}

void IDHistosComp_pT(){
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

