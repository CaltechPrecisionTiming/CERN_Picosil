//#ifdef __MAKECINT__
//#pragma link C++ class vector<vector<float> >+;
//#endif

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
#include "SiliconPadUtils.h"




void MakeAmplitudeVsBeamEnergyGraph() {

  //use beam energy for xaxis
  const int nPoints = 4;
  float x[nPoints] = { 4.0, 8.0, 16.0, 32.0 };
  float xerr[nPoints] = { 0.027*4.0, 0.023*8.0, 0.045*16.0, 0.05*32.0 };
  float y_charge[nPoints] = {  54.8, 106, 217, 350 }; 
  float yerr_charge[nPoints] = { 19, 46, 46, 92 };
  float y_MIP[nPoints] = {  0.0, 0.0, 0.0, 0.0 }; 
  float yerr_MIP[nPoints] = { 0.0, 0.0, 0.0, 0.0 };

  double chargePerMIP = 6.5;
  for (int i=0; i<nPoints; ++i) {
    y_MIP[i] = y_charge[i] / chargePerMIP;
    yerr_MIP[i] =  yerr_charge[i] / chargePerMIP;
  }

  TGraphErrors *graphChargeVsEnergy = new TGraphErrors(nPoints,x,y_charge,xerr,yerr_charge);
  graphChargeVsEnergy->SetLineWidth(3);
  TGraphErrors *graphMIPVsEnergy = new TGraphErrors(nPoints,x,y_MIP,xerr,yerr_MIP);
  graphMIPVsEnergy->SetLineWidth(3);

  TCanvas *c = 0;
  TVirtualFitter *fitter = 0;

  c = new TCanvas("c","c",800,600);
  graphChargeVsEnergy->Draw("AP");
  graphChargeVsEnergy->SetTitle("");
  graphChargeVsEnergy->GetXaxis()->SetTitle("Electron Beam Energy [GeV/c^{2}]");
  graphChargeVsEnergy->GetXaxis()->SetTitleSize(0.045);
  graphChargeVsEnergy->GetXaxis()->SetLabelSize(0.045);
  graphChargeVsEnergy->GetXaxis()->SetTitleOffset(1.0);
  graphChargeVsEnergy->GetYaxis()->SetTitle("Integrated Charge [pC]");
  graphChargeVsEnergy->GetYaxis()->SetTitleOffset(1.0);
  graphChargeVsEnergy->GetYaxis()->SetTitleSize(0.05);
  graphChargeVsEnergy->GetYaxis()->SetLabelSize(0.045);
  graphChargeVsEnergy->GetXaxis()->SetRangeUser(0,40);
  graphChargeVsEnergy->GetYaxis()->SetRangeUser(0,500);

  graphChargeVsEnergy->Fit("pol1","","");
  fitter = TVirtualFitter::GetFitter();
  
  c->SaveAs( "ChargeVsEnergyAt6X0.gif" );
  c->SaveAs( "ChargeVsEnergyAt6X0.pdf" );


  c = new TCanvas("c","c",800,600);
  c->SetLeftMargin(0.11);
  graphMIPVsEnergy->Draw("AP");
  graphMIPVsEnergy->SetTitle("");
  graphMIPVsEnergy->GetXaxis()->SetTitle("Electron Beam Energy [GeV/c^{2}]");
  graphMIPVsEnergy->GetXaxis()->SetTitleSize(0.045);
  graphMIPVsEnergy->GetXaxis()->SetLabelSize(0.045);
  graphMIPVsEnergy->GetXaxis()->SetTitleOffset(1.0);
  graphMIPVsEnergy->GetYaxis()->SetTitle("Integrated Charge [ Q_{MIP} ]");
  graphMIPVsEnergy->GetYaxis()->SetTitleOffset(0.9);
  graphMIPVsEnergy->GetYaxis()->SetTitleSize(0.05);
  graphMIPVsEnergy->GetYaxis()->SetLabelSize(0.045);
  graphMIPVsEnergy->GetXaxis()->SetRangeUser(0,40);
  graphMIPVsEnergy->GetYaxis()->SetRangeUser(0,80);

  graphMIPVsEnergy->Fit("pol1","","");
  fitter = TVirtualFitter::GetFitter();
  
  c->SaveAs( "MIPVsEnergyAt6X0.gif" );
  c->SaveAs( "MIPVsEnergyAt6X0.pdf" );


}

void MakeAmplitudeVsAbsorberThicknessGraph() {

  //use beam energy for xaxis
  const int nPoints = 5;
  float x[nPoints] = { 0.0, 1.0, 2.0, 4.0, 6.0 };
  float xerr[nPoints] = { 0.1, 0.25, 0.25, 0.25, 0.25 };
  float y_charge[nPoints] = {  6.5, 33, 57.0, 109, 106 }; 
  float yerr_charge[nPoints] = { 1.4, 24, 22, 60, 46 };
  float y_MIP[nPoints] = {  0.0, 0.0, 0.0, 0.0, 0.0 }; 
  float yerr_MIP[nPoints] = { 0.0, 0.0, 0.0, 0.0, 0.0 };

  double chargePerMIP = 6.5;
  for (int i=0; i<nPoints; ++i) {
    y_MIP[i] = y_charge[i] / chargePerMIP;
    yerr_MIP[i] =  yerr_charge[i] / chargePerMIP;
  }

  TGraphErrors *graphChargeVsAbsorber = new TGraphErrors(nPoints,x,y_charge,xerr,yerr_charge);
  graphChargeVsAbsorber->SetLineWidth(3);
  TGraphErrors *graphMIPVsAbsorber = new TGraphErrors(nPoints,x,y_MIP,xerr,yerr_MIP);
  graphMIPVsAbsorber->SetLineWidth(3);

  TCanvas *c = 0;
  TVirtualFitter *fitter = 0;

  c = new TCanvas("c","c",800,600);
  c->SetBottomMargin(0.12);
  graphChargeVsAbsorber->Draw("AP");
  graphChargeVsAbsorber->SetTitle("");
  graphChargeVsAbsorber->GetXaxis()->SetTitle("Tungsten Absorber Thickness [X_{0}]");
  graphChargeVsAbsorber->GetXaxis()->SetTitleSize(0.045);
  graphChargeVsAbsorber->GetXaxis()->SetLabelSize(0.045);
  graphChargeVsAbsorber->GetXaxis()->SetTitleOffset(1.2);
  graphChargeVsAbsorber->GetYaxis()->SetTitle("Integrated Charge [pC]");
  graphChargeVsAbsorber->GetYaxis()->SetTitleOffset(1.0);
  graphChargeVsAbsorber->GetYaxis()->SetTitleSize(0.05);
  graphChargeVsAbsorber->GetYaxis()->SetLabelSize(0.045);
  graphChargeVsAbsorber->GetXaxis()->SetRangeUser(0.0,8.0);
  graphChargeVsAbsorber->GetYaxis()->SetRangeUser(0,200);

  // graphChargeVsAbsorber->Fit("pol2","","");
  // fitter = TVirtualFitter::GetFitter();
  
  c->SaveAs( "ChargeVsAbsorberAt8GeV.gif" );
  c->SaveAs( "ChargeVsAbsorberAt8GeV.pdf" );


  c = new TCanvas("c","c",800,600);
  c->SetBottomMargin(0.12);
  c->SetLeftMargin(0.11);

  graphMIPVsAbsorber->Draw("AP");
  graphMIPVsAbsorber->SetTitle("");
  graphMIPVsAbsorber->GetXaxis()->SetTitle("Tungsten Absorber Thickness [X_{0}]");
  graphMIPVsAbsorber->GetXaxis()->SetTitleSize(0.045);
  graphMIPVsAbsorber->GetXaxis()->SetLabelSize(0.045);
  graphMIPVsAbsorber->GetXaxis()->SetTitleOffset(1.2);
  graphMIPVsAbsorber->GetYaxis()->SetTitle("Integrated Charge [ Q_{MIP} ]");
  graphMIPVsAbsorber->GetYaxis()->SetTitleOffset(0.9);
  graphMIPVsAbsorber->GetYaxis()->SetTitleSize(0.05);
  graphMIPVsAbsorber->GetYaxis()->SetLabelSize(0.045);
  graphMIPVsAbsorber->GetXaxis()->SetRangeUser(0.0, 8.0);
  graphMIPVsAbsorber->GetYaxis()->SetRangeUser(0, 30);

  // graphMIPVsAbsorber->Fit("pol1","","");
  // fitter = TVirtualFitter::GetFitter();
  
  c->SaveAs( "MIPVsAbsorberAt8GeV.gif" );
  c->SaveAs( "MIPVsAbsorberAt8GeV.pdf" );


}



void MakeAmplitudeVsBiasVoltageGraph() {

  //With the 600V point
  // const int nPoints = 5;
  // float x[nPoints] = { 200, 300, 400, 500, 600 };
  // float xerr[nPoints] = { 0.01*200, 0.01*300, 0.01*400, 0.01*500, 0.01*600 };
  // float y_charge[nPoints] = {  58.4, 68.0, 73.8, 77.1, 78.6 }; 
  // float yerr_charge[nPoints] = { 1.0, 1.4, 1.3, 0.5, 1.4 };
  // float y_MIP[nPoints] = {  0.0, 0.0, 0.0, 0.0, 0.0 }; 
  // float yerr_MIP[nPoints] = { 0.0, 0.0, 0.0, 0.0, 0.0 };

  //Without the 600V point
  const int nPoints = 5;
  float x[nPoints] = { 200, 300, 400, 500 };
  float xerr[nPoints] = { 0.01*200, 0.01*300, 0.01*400, 0.01*500 };
  float y_charge[nPoints] = {  58.4, 68.0, 73.8, 77.1 }; 
  float yerr_charge[nPoints] = { 1.0, 1.4, 1.3, 0.5 };
  float y_MIP[nPoints] = {  0.0, 0.0, 0.0, 0.0 }; 
  float yerr_MIP[nPoints] = { 0.0, 0.0, 0.0, 0.0 };


  double chargePerMIP = 1.3;
  for (int i=0; i<nPoints; ++i) {
    y_MIP[i] = y_charge[i] / chargePerMIP;
    yerr_MIP[i] =  yerr_charge[i] / chargePerMIP;
  }

  TGraphErrors *graphChargeVsBiasVoltage = new TGraphErrors(nPoints,x,y_charge,xerr,yerr_charge);
  graphChargeVsBiasVoltage->SetLineWidth(3);
  TGraphErrors *graphMIPVsBiasVoltage = new TGraphErrors(nPoints,x,y_MIP,xerr,yerr_MIP);
  graphMIPVsBiasVoltage->SetLineWidth(3);

  TCanvas *c = 0;
  TVirtualFitter *fitter = 0;

  c = new TCanvas("c","c",800,600);
  graphChargeVsBiasVoltage->Draw("AP");
  graphChargeVsBiasVoltage->SetTitle("");
  graphChargeVsBiasVoltage->GetXaxis()->SetTitle("Bias Voltage [V]");
  graphChargeVsBiasVoltage->GetXaxis()->SetTitleSize(0.045);
  graphChargeVsBiasVoltage->GetXaxis()->SetLabelSize(0.045);
  graphChargeVsBiasVoltage->GetXaxis()->SetTitleOffset(1.0);
  graphChargeVsBiasVoltage->GetYaxis()->SetTitle("Integrated Charge [pC]");
  graphChargeVsBiasVoltage->GetYaxis()->SetTitleOffset(1.02);
  graphChargeVsBiasVoltage->GetYaxis()->SetTitleSize(0.045);
  graphChargeVsBiasVoltage->GetYaxis()->SetLabelSize(0.045);
  graphChargeVsBiasVoltage->GetXaxis()->SetRangeUser(0.0,700);
  graphChargeVsBiasVoltage->GetYaxis()->SetRangeUser(0,90);

  // graphChargeVsBiasVoltage->Fit("pol2","","");
  // fitter = TVirtualFitter::GetFitter();
  
  c->SaveAs( "ChargeVsBiasVoltageAt6X0_16GeV.gif" );
  c->SaveAs( "ChargeVsBiasVoltageAt6X0_16GeV.pdf" );


  c = new TCanvas("c","c",800,600);
  c->SetLeftMargin(0.11);
  graphMIPVsBiasVoltage->Draw("AP");
  graphMIPVsBiasVoltage->SetTitle("");
  graphMIPVsBiasVoltage->GetXaxis()->SetTitle("Bias Voltage [V]");
  graphMIPVsBiasVoltage->GetXaxis()->SetTitleSize(0.045);
  graphMIPVsBiasVoltage->GetXaxis()->SetLabelSize(0.045);
  graphMIPVsBiasVoltage->GetXaxis()->SetTitleOffset(1.0);
  graphMIPVsBiasVoltage->GetYaxis()->SetTitle("Integrated Charge [ Q_{MIP} ]");
  graphMIPVsBiasVoltage->GetYaxis()->SetTitleOffset(1.02);
  graphMIPVsBiasVoltage->GetYaxis()->SetTitleSize(0.045);
  graphMIPVsBiasVoltage->GetYaxis()->SetLabelSize(0.045);
  graphMIPVsBiasVoltage->GetXaxis()->SetRangeUser(0.0, 700);
  graphMIPVsBiasVoltage->GetYaxis()->SetRangeUser(0, 70);

  // graphMIPVsBiasVoltage->Fit("pol1","","");
  // fitter = TVirtualFitter::GetFitter();
  
  c->SaveAs( "MIPVsBiasVoltageAt6X0_16GeV.gif" );
  c->SaveAs( "MIPVsBiasVoltageAt6X0_16GeV.pdf" );


}



void makeEMShowerEnergyPlots() {

  //*************************************
  // Charge Vs Energy
  //*************************************

  //Noise Runs
  //MakeChargePlot("root://eoscms:///eos/cms/store/group/phys_susy/razor/Timing/t1065-dec-2015/reco/v1/t1065_dec2015_run120-0.dat-full.root", "NoiseNoBeam", -99, sqrt(10.0), -10, 10, -4,4 , true);
  //MakeChargePlot("root://eoscms:///eos/cms/store/group/phys_susy/razor/Timing/t1065-dec-2015/reco/v1/t1065_dec2015_run121-0.dat-full.root", "NoiseNoBeamNoAttenuator", -99, 1, -10, 10, -4,4 , true);


  //Protons
  //MakeChargePlot("root://eoscms:///eos/cms/store/group/phys_susy/razor/Timing/t1065-dec-2015/reco/v1/t1065_dec2015_run125.root", "Proton", -99, 1, -5, 50, 0,10, true );


  //Electrons
  //MakeChargePlot("root://eoscms:///eos/cms/store/group/phys_susy/razor/Timing/t1065-dec-2015/reco/v1/t1065_dec2015_run123-124.root", "Electron_0X0", -99, 1, -5, 50, 0,10, true );
  
  //MakeChargePlot("root://eoscms:///eos/cms/store/group/phys_susy/razor/Timing/t1065-dec-2015/reco/v1/t1065_dec2015_run112-500V_123.root", "Electron_6X0_4GeV", 0.02, 10.0, 0, 200, 20, 80 );
  //MakeChargePlot("root://eoscms:///eos/cms/store/group/phys_susy/razor/Timing/t1065-dec-2015/reco/v1/t1065_dec2015_run109.root","Electron_6X0_8GeV", 0.05, 10.0, 0, 300, 50, 175);
  //MakeChargePlot("root://eoscms://eos/cms/store/group/phys_susy/razor/Timing/t1065-dec-2015/reco/v1/t1065_dec2015_run110.root","Electron_6X0_16GeV", 0.1, 10.0, 0, 500, 150, 300);
  MakeChargePlot("root://eoscms://eos/cms/store/group/phys_susy/razor/Timing/t1065-dec-2015/reco/v1/t1065_dec2015_run113-114.root","Electron_6X0_32GeV", 0.1, 10.0, 0, 700, 250, 450);
 
  //MakeAmplitudeVsBeamEnergyGraph();
  

  //*************************************
  // Charge Vs Absorber Thickness
  //*************************************
  //MakeChargePlot("root://eoscms:///eos/cms/store/group/phys_susy/razor/Timing/t1065-dec-2015/reco/v1/t1065_dec2015_run122.root", "Electron_1X0_8GeV", 0.05, 1.0, 0, 100, 15, 60 );
  //MakeChargePlot("root://eoscms:///eos/cms/store/group/phys_susy/razor/Timing/t1065-dec-2015/reco/v1/t1065_dec2015_run119.root", "Electron_1X0_8GeV_10dbAttenuator", 0.02, sqrt(10.0), 0, 100, 15, 60 );
  //MakeChargePlot("root://eoscms:///eos/cms/store/group/phys_susy/razor/Timing/t1065-dec-2015/reco/v1/t1065_dec2015_run118.root", "Electron_2X0_8GeV", 0.03, sqrt(10.0), 0, 200, 25, 90 );
  //MakeChargePlot("root://eoscms:///eos/cms/store/group/phys_susy/razor/Timing/t1065-dec-2015/reco/v1/t1065_dec2015_run117.root", "Electron_4X0_8GeV", 0.05, sqrt(10.0), 0, 300, 40, 175 );
  //MakeChargePlot("root://eoscms:///eos/cms/store/group/phys_susy/razor/Timing/t1065-dec-2015/reco/v1/t1065_dec2015_run116.root", "Electron_6X0_8GeV", 0.05, sqrt(10.0), 0, 300, 30, 150 );

  //MakeAmplitudeVsAbsorberThicknessGraph();


  //*************************************
  // Bias Voltage Scan
  //*************************************
  //MakeChargePlot("root://eoscms:///eos/cms/store/group/phys_susy/razor/Timing/t1065-dec-2015/reco/v1/t1065_dec2015_run111_600V.root", "Electron_6X0_8GeV_600V", 0.1, 10.0, 0, 300, 30, 150 );
  //MakeChargePlot("root://eoscms:///eos/cms/store/group/phys_susy/razor/Timing/t1065-dec-2015/reco/v1/t1065_dec2015_run110.root", "Electron_6X0_8GeV_500V", 0.1, 10.0, 0, 300, 30, 150 );
  //MakeChargePlot("root://eoscms:///eos/cms/store/group/phys_susy/razor/Timing/t1065-dec-2015/reco/v1/t1065_dec2015_run111_400V.root", "Electron_6X0_8GeV_400V", 0.1, 10.0, 0, 300, 30, 140 );
  //MakeChargePlot("root://eoscms:///eos/cms/store/group/phys_susy/razor/Timing/t1065-dec-2015/reco/v1/t1065_dec2015_run111_300V.root", "Electron_6X0_8GeV_300V", 0.1, 10.0, 0, 300, 30, 130 );
  //MakeChargePlot("root://eoscms:///eos/cms/store/group/phys_susy/razor/Timing/t1065-dec-2015/reco/v1/t1065_dec2015_run111_200V.root", "Electron_6X0_8GeV_200V", 0.1, 10.0, 0, 300, 20, 100 );

  //MakeAmplitudeVsBiasVoltageGraph();

}
