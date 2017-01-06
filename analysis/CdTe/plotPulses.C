#include <iostream>
#include <fstream> 
#include <sstream>
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TF1.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TVirtualFitter.h"
#include "TPaveText.h"
#include "TProfile.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TColor.h"
#include "TGraphErrors.h"
#include "TLatex.h"
#include <math.h> 
//#include "SiliconPadUtils.h"




void plotPulses() {

  TFile *file_100GeV = TFile::Open( "/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Timing/Nov2016CERN/ntuples_v6/analysis_5568.root" , "READ");
  TFile *file_200GeV = TFile::Open( "/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Timing/Nov2016CERN/ntuples_v6/analysis_5570.root" , "READ");
  TFile *file_50GeV = TFile::Open( "/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Timing/Nov2016CERN/ntuples_v6/analysis_5571-5576.root" , "READ");
  TTree *tree_100GeV = (TTree*)(file_100GeV->Get("t1065"));
  TTree *tree_200GeV = (TTree*)(file_200GeV->Get("t1065"));
  TTree *tree_50GeV = (TTree*)(file_50GeV->Get("t1065"));

  TH1F *pulse_100GeV = new TH1F("pulse_100GeV", ";Time [ns];Amplitude [V]", 1024, -0.5*0.2, 1023.5*0.2);
  tree_100GeV->Draw("t0*0.2>>pulse_100GeV", "-1*raw[1]*(1.0/4096.)*3.16228*(1.0/63.0957)*(event==187 && amp[0]>0.2116 && amp[0]<0.2117)");
  TH1F *pulse_200GeV = new TH1F("pulse_200GeV", ";Time [ns];Amplitude [V]", 1024, -0.5*0.2, 1023.5*0.2);
  tree_200GeV->Draw("t0*0.2>>pulse_200GeV", "-1*raw[1]*(1.0/4096.)*3.16228*(1.0/63.0957)*(event==496 && amp[0]>0.4025 && amp[0]<0.4026)");
  TH1F *pulse_50GeV = new TH1F("pulse_50GeV", ";Time [ns];Amplitude [V]", 1024, -0.5*0.2, 1023.5*0.2);
  tree_50GeV->Draw("t0*0.2>>pulse_50GeV", "-1*raw[1]*(1.0/4096.)*3.16228*(1.0/63.0957)*(event==340 && amp[0]>0.0810 && amp[0]<0.0811)");
 

  TCanvas *cv = new TCanvas("cv","cv", 800, 600);
  pulse_200GeV->Draw("hist");
  pulse_200GeV->SetStats(false);
  pulse_200GeV->GetYaxis()->SetTitleOffset(1.5);
  pulse_200GeV->GetXaxis()->SetTitleOffset(1.1);
  pulse_200GeV->GetXaxis()->SetTitleSize(0.05);
  pulse_200GeV->GetXaxis()->SetLabelSize(0.04);
  pulse_200GeV->GetYaxis()->SetTitleSize(0.05);
  pulse_200GeV->GetYaxis()->SetLabelSize(0.04);
  pulse_200GeV->SetLineWidth(2);
  cv->SetLeftMargin(0.15);
  cv->SetBottomMargin(0.12);

  pulse_100GeV->SetLineWidth(2);
  pulse_100GeV->SetLineColor(kRed);

  pulse_100GeV->Draw("histsame");

  pulse_50GeV->SetLineColor(kBlack);
  pulse_50GeV->Draw("histsame");
  pulse_50GeV->SetLineWidth(2);

  TLegend *legend = new TLegend( 0.6,0.6,0.8,0.8);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->SetTextSize(0.05);
  legend->AddEntry(pulse_200GeV ,"200 GeV");
  legend->AddEntry(pulse_100GeV ,"100 GeV");
  legend->AddEntry(pulse_50GeV ,"50 GeV");
  legend->Draw();
 
  cv->SaveAs("TypicalPulses.C");
  cv->SaveAs("TypicalPulses.png");
  cv->SaveAs("TypicalPulses.pdf");

  pulse_200GeV->GetXaxis()->SetRangeUser(30,60);
  cv->SaveAs("TypicalPulses_Zoom.C");
  cv->SaveAs("TypicalPulses_Zoom.png");
  cv->SaveAs("TypicalPulses_Zoom.pdf");
  



}
