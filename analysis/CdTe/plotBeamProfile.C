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




void plotBeamProfile( int option = 100 ) {

  TFile *file = 0;
  if (option == 100) {
    file = TFile::Open( "/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Timing/Nov2016CERN/ntuples_v6/analysis_5568.root" , "READ");
  }
  if (option == 50) {
    file = TFile::Open( "/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Timing/Nov2016CERN/ntuples_v6/analysis_5571-5576.root" , "READ");
  }
  if (option == 200) {
    file = TFile::Open( "/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Timing/Nov2016CERN/ntuples_v6/analysis_5570.root" , "READ");
  }

  TTree *tree = (TTree*)(file->Get("t1065"));

  TH2F *beam = new TH2F("beam", " ; Vertical Beam Position [mm] ; Horizontal Beam Position [mm]; Number of Events", 35,-20,15,35,-10,25);
  TH2F *beamSensor = new TH2F("beamSensor", " ; Vertical Beam Position [mm] ; Horizontal Beam Position [mm]; Amplitude Weighted Events", 35,-20,15,35,-10,25);
  TH2F *beamXVsAmp = new TH2F("beamXVsAmp", " ; Vertical Beam Position [mm] ; Signal Amplitude [V]; Number of Events", 30,-5,25,4000,0,1.0);
  TH2F *beamYVsAmp = new TH2F("beamYVsAmp", " ; Horizontal Beam Position [mm] ; Signal Amplitude [V]; Number of Events", 30,-15,15,4000,0,1.0);
  tree->Draw("TDCx:TDCy>>beam","","colz");
  tree->Draw("TDCx:TDCy>>beamSensor","amp[1]","colz");

  // For 100 GeV run
  if (option == 100) {
    tree->Draw("amp[1]*3.16228*(1.0/63.0957):TDCx>>beamXVsAmp","TDCy>-8.0 && TDCy < 2.0","colz");
    tree->Draw("amp[1]*3.16228*(1.0/63.0957):TDCy>>beamYVsAmp","TDCx>2.5 && TDCx < 13.5","colz");
  }

  //For 200 GeV run
  if (option == 200) {
     tree->Draw("amp[1]*3.16228*(1.0/63.0957):TDCx>>beamXVsAmp","TDCy>-8.0 && TDCy < 2.0","colz");
     tree->Draw("amp[1]*3.16228*(1.0/63.0957):TDCy>>beamYVsAmp","TDCx>2.5 && TDCx < 13.5","colz");
  }

  //For 50 GeV run
  if (option == 50) {
    tree->Draw("amp[1]*3.16228*(1.0/63.0957):TDCx>>beamXVsAmp","TDCy>-8.0 && TDCy < 2.0","colz");
    tree->Draw("amp[1]*3.16228*(1.0/63.0957):TDCy>>beamYVsAmp","TDCx>2.5 && TDCx < 13.5","colz");
  }

  TCanvas *cv = new TCanvas("cv","cv", 800, 800);
  beam->Draw("colz");
  beam->SetStats(false);
  beam->GetYaxis()->SetTitleOffset(1.5);
  beam->GetXaxis()->SetTitleOffset(1.2);

  cv->SetRightMargin(0.15);
  cv->SetLeftMargin(0.12);
  cv->SetBottomMargin(0.12);

  cv->SaveAs("BeamProfile.png");
  cv->SaveAs("BeamProfile.pdf");

  cv = new TCanvas("cv","cv", 800, 800);
  beamSensor->Draw("colz");
  beamSensor->SetStats(false);
  beamSensor->GetYaxis()->SetTitleOffset(1.5);
  beamSensor->GetXaxis()->SetTitleOffset(1.2);

  cv->SetRightMargin(0.15);
  cv->SetLeftMargin(0.12);
  cv->SetBottomMargin(0.12);

  cv->SaveAs("BeamSensorProfile.png");
  cv->SaveAs("BeamSensorProfile.pdf");

  cv = new TCanvas("cv","cv", 800, 800);
  TProfile *beamXProfile = beamXVsAmp->ProfileX();
  beamXProfile->Draw();
  beamXProfile->SetLineWidth(2);
  beamXProfile->SetStats(false);
  beamXProfile->GetYaxis()->SetTitleOffset(2.5);
  beamXProfile->GetYaxis()->SetTitle("Mean Signal Amplitude [V]");
  beamXProfile->GetXaxis()->SetTitleOffset(1.2);

  cv->SetLeftMargin(0.2);
  cv->SetBottomMargin(0.12);

  cv->SaveAs("SensorXProfile.png");
  cv->SaveAs("SensorXProfile.pdf");

  cv = new TCanvas("cv","cv", 800, 800);
  TProfile *beamYProfile = beamYVsAmp->ProfileX();
  beamYProfile->Draw();
  beamYProfile->SetLineWidth(2);
  beamYProfile->SetStats(false);
  beamYProfile->GetYaxis()->SetTitleOffset(2.5);
  beamYProfile->GetYaxis()->SetTitle("Mean Signal Amplitude [V]");
  beamYProfile->GetXaxis()->SetTitleOffset(1.2);

  cv->SetLeftMargin(0.2);
  cv->SetBottomMargin(0.12);

  cv->SaveAs("SensorYProfile.png");
  cv->SaveAs("SensorYProfile.pdf");



}
