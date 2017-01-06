
#include <iostream>
#include <fstream> 
#include <sstream>
#include <map>
#include <utility>
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
#include <time.h>
#include <algorithm>
#include <functional>
#include "TRandom3.h"

void Fitter(TH1F *hist) {
  //Helper function for fitting Gaussian, 2RMS around mean
  double xmin = hist->GetMean() - 2.0*hist->GetRMS();
  double xmax = hist->GetMean() + 2.0*hist->GetRMS();
  hist->Fit("gaus","QMLES","",xmin,xmax); // Q suppresses fit results
  gStyle->SetOptFit(1);
}

void makeOutputFile( string filename, float photekAmpCut, float photekChargeCut, float centerAmpCut, float centerChargeCut, float TDCxmin, float TDCxmax, float TDCymin, float TDCymax, float centerAtn, float ring1Atn, float ring2Atn, float photekAtn ) {

  TFile *inputfile = TFile::Open(filename.c_str(),"READ");
  TTree *tree = (TTree*)inputfile->Get("t1065");

  // get the variables from the ntuple
  float amp[36];
  float integral[36];
  float gauspeak[36];
  float linearTime45[36];
  float TDCx[1];
  float TDCy[1];

  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("gauspeak",1);
  tree->SetBranchStatus("linearTime45",1);
  tree->SetBranchStatus("amp",1);
  tree->SetBranchStatus("int",1);
  tree->SetBranchStatus("TDCx",1);
  tree->SetBranchStatus("TDCy",1);
  
  tree->SetBranchAddress("gauspeak",gauspeak);
  tree->SetBranchAddress("linearTime45",linearTime45);
  tree->SetBranchAddress("amp",amp);
  tree->SetBranchAddress("int",integral);
  tree->SetBranchAddress("TDCx",TDCx);
  tree->SetBranchAddress("TDCy",TDCy);

  //Create histograms
  float width = 0.8;
  int bins = 80;
  string axes = Form(";#Deltat (ns);Entries/(%.3f ns)",(width*2/bins));
  string axesHGC = Form(";#Deltat_{HGC} (ns);Entries/(%.3f ns)",(width*2/bins));

  TH1F *histDeltaTPicoSil[19];
  TH1F *histDeltaTPicoSilAt0[19];
  TH1F *histCharges[19]; // collects charge values for picosil pixels in every event in which they pass the cuts.

  for(int i=0; i<19; i++) {
    histDeltaTPicoSil[i] = new TH1F(Form("histDeltaTPicoSil_%d",i),";#Deltat (ns);Entries", 100, -50, 50); //DeltaT of PicoSil pixels
    histDeltaTPicoSilAt0[i] = new TH1F(Form("histDeltaTPicoSilAt0_%d",i), axes.c_str(), bins, -width, width); 
    histCharges[i] = new TH1F( Form("histCharges_%d",i),";Charge (pC);Entries/(3.0 pC)", 50, 0, 150);
  }
  TH1F *histDeltaTCenter = new TH1F("histDeltaTCenter",";#Deltat (ns);Entries", 50, -20, 20); //DeltaT of center picosil pixel
  TH1F *histDeltaTPicoSilAt0TotalCharge = new TH1F("histDeltaTPicoSilAt0TotalCharge", axesHGC.c_str(), bins, -width, width); //All pixels combined
  TH1F *histDeltaTPicoSilAt0EventCharge = new TH1F("histDeltaTPicoSilAt0EventCharge", axesHGC.c_str(), bins, -width, width);
  TH1F *histDeltaTPicoSilAt0LandauCharge = new TH1F("histDeltaTPicoSilAt0LandauCharge", axesHGC.c_str(), bins, -width, width);// uses charge MPV
  TH1F *histDeltaTCenterAt0 = new TH1F("histDeltaTCenterAt0", axes.c_str(), bins, -width, width); //shifted to be centered at zero

  float totalPicoSilCharge[19] = {0.};


  //read all entries and fill the histogram
  Long64_t nentries = tree->GetEntries();

  //Loop through every event in .root file
  std::cout<<"Number of events in sample: "<<nentries<<std::endl;  
  for (Long64_t iEntry=0; iEntry<nentries; iEntry++) {
    if (iEntry %1000 == 0) cout << "Processing Event " << iEntry << "\n";
    tree->GetEntry(iEntry); 
  
    float photekTimeGauss0 = gauspeak[0]; // 0-->Channel Group 0
    //UNUSED: float photekAmp0 = photekAtn * amp[0];
    //UNUSED: float photekCharge0 = photekAtn * integral[0];

    float photekTimeGauss1 = gauspeak[8]; // 1-->Channel Group 1
    //UNUSED: float photekAmp1 = photekAtn * amp[8];
    //UNUSED: float photekCharge1 = photekAtn * integral[8];

    float photekTimeGauss2 = gauspeak[16]; // 2-->Channel Group 2
    float photekAmp2 = photekAtn * amp[16];
    float photekCharge2 = photekAtn * integral[16];

    float photekTimeGauss3 = gauspeak[24]; // 3-->Channel Group 3
    //UNUSED: float photekAmp3 = photekAtn * amp[24];
    //UNUSED: float photekCharge3 = photekAtn * integral[24];

    float centerAmp = centerAtn * amp[17]; //accounts for 10dB attenuator
    float centerCharge = centerAtn * integral[17];
    float centerTime = linearTime45[17];
    float centerTDCx = TDCx[0];
    float centerTDCy = TDCy[0];

    // APPLY EVENT CUTS:
    //require photek (for electron selection)
    if( !(photekAmp2 > photekAmpCut && photekCharge2 > photekChargeCut) ) continue;
    //require signal in the central pixel
    if( !(centerAmp > centerAmpCut && centerCharge > centerChargeCut) ) continue;
    //require wire chamber
    if( !(centerTDCx > TDCxmin && centerTDCx < TDCxmax) ) continue;
    if( !(centerTDCy > TDCymin && centerTDCy < TDCymax) ) continue;

    //Calculates the Delta T's if the event passes the cuts:
    // (Will be used to fill the histogram at every event)
    float DeltaTPicoSil[19];
    std::fill(DeltaTPicoSil, DeltaTPicoSil+19, -99);
    // Center pixel
    DeltaTPicoSil[0] = photekTimeGauss2 - centerTime;
    totalPicoSilCharge[0] += centerCharge;
    histCharges[0]->Fill(centerCharge);
    // Ring 1: Cuts on additional pixels to determine whether they get incorporated.
    for ( int j = 1; j <= 6; j++){
      if ( ring1Atn*amp[j+17] > 0.01 && ring1Atn*integral[j+17] > 1 ) { //j+17 --> want elts 18-23 for 1st ring.
        DeltaTPicoSil[j] = photekTimeGauss2 - linearTime45[j+17]; 
        totalPicoSilCharge[j] += ring1Atn * integral[j+17];
        histCharges[j]->Fill( ring1Atn * integral[j+17] );
      }
    }
    // Ring 2:
    for ( int j = 7; j <= 13; j++){  // Group 3 pixels
      if ( ring2Atn*amp[j+18] > 0.01 && ring2Atn*integral[j+18] > 1 ) { //j+18 --> want elts 25-31 for 2nd ring.
        DeltaTPicoSil[j] = photekTimeGauss3 - linearTime45[j+18];
        totalPicoSilCharge[j] += ring2Atn * integral[j+18];
        histCharges[j]->Fill( ring2Atn * integral[j+18] );
      }
    }
    for ( int j = 14; j <= 16; j++){ // Group 0 pixels
      if ( ring2Atn*amp[j-9] > 0.01 && ring2Atn*integral[j-9] > 1 ) { //j-9 --> want elts 5-7 for 2nd ring.
        DeltaTPicoSil[j] = photekTimeGauss0 - linearTime45[j-9];
        totalPicoSilCharge[j] += ring2Atn * integral[j-9];
        histCharges[j]->Fill( ring2Atn * integral[j-9] );
      }
    }
    for ( int j = 17; j <= 18; j++){ // Group 1 pixels
      if ( ring2Atn*amp[j-4] > 0.01 && ring2Atn*integral[j-4] > 1 ) { //j-4 --> want elts 13-14 for 2nd ring.
        DeltaTPicoSil[j] = photekTimeGauss1 - linearTime45[j-4];
        totalPicoSilCharge[j] += ring2Atn * integral[j-4];
        histCharges[j]->Fill( ring2Atn * integral[j-4] );
      }
    }



    histDeltaTCenter->Fill(DeltaTPicoSil[0]);
    for(int k=0; k<=18; k++) { if (DeltaTPicoSil[k] != -99.) histDeltaTPicoSil[k]->Fill(DeltaTPicoSil[k]); }
  }




  // Do Gaussian fit of delta T distributions from (mean-2RMS) to (mean+2RMS)
  Fitter(histDeltaTCenter);


  TF1 *flandau[19];
  double MPVlandau[19]; // Seeing if we can get better results by using Landau mean as weighting charge value.
  for(int i=1; i<=18; i++) {
    double mean = histCharges[i]->GetMean();
    double rms = histCharges[i]->GetRMS();
    double xmin = mean-2.0*rms;
    double xmax = mean+2.0*rms;
    flandau[i] = new TF1( Form("flandau_%d",i), "landau", xmin, xmax); // 1-D landau func
    histCharges[i]->Fit( Form("flandau_%d",i), "QMLES","", xmin, xmax);
    gStyle->SetOptFit(1);
    MPVlandau[i] = flandau[i]->GetMaximumX(); // In order to find MPV.
  }

  double mean = histCharges[0]->GetMean();
  double rms = histCharges[0]->GetRMS();
  double xmin = mean-2.0*rms;
  double xmax = mean+2.0*rms;
  flandau[0] = new TF1("flandau_0","gaus", xmin, xmax); // storing a Gaus fit for the central pixel with the other landau fits
  histCharges[0]->Fit("flandau_0","QMLES","", xmin, xmax);
  gStyle->SetOptFit(1);
  MPVlandau[0] = flandau[0]->GetParameter(1);
  //cout<<"\nhistCharges[0]\nGauss Mean: "<<flandau[0]->GetParameter(1)<<
  //                      "\nUncertainy: "<<flandau[0]->GetParError(1)<<
  //"\n"<<endl;


  double meanPicoSil[19];
  for (int i = 0; i < 19; i++) meanPicoSil[i] = histDeltaTPicoSil[i]->GetMean();




  //Loop through again. This time, subtracting the respective means. Process is almost the same.
  for (Long64_t iEntry=0; iEntry<nentries; iEntry++) {
    if (iEntry %1000 == 0) cout << "Processing Event " << iEntry << "\n";
    tree->GetEntry(iEntry);    
  
    float photekTimeGauss0 = gauspeak[0];
    //UNUSED: float photekAmp0 = photekAtn * amp[0]; //accounts for attenuator
    //UNUSED: float photekCharge0 = photekAtn * integral[0];

    float photekTimeGauss1 = gauspeak[8]; 
    //UNUSED: float photekAmp1 = photekAtn * amp[8];
    //UNUSED: float photekCharge1 = photekAtn * integral[8];

    float photekTimeGauss2 = gauspeak[16];
    float photekAmp2 = photekAtn * amp[16];
    float photekCharge2 = photekAtn * integral[16];

    float photekTimeGauss3 = gauspeak[24];
    //UNUSED: float photekAmp3 = photekAtn * amp[24];
    //UNUSED: float photekCharge3 = photekAtn * integral[24];

    float centerAmp = centerAtn * amp[17]; //accounts for attenuator
    float centerCharge = centerAtn * integral[17];
    float centerTime = linearTime45[17];
    float centerTDCx = TDCx[0];
    float centerTDCy = TDCy[0];


    // APPLY EVENT CUTS:
    if( !(photekAmp2 > photekAmpCut && photekCharge2 > photekChargeCut) ) continue;
    if( !(centerAmp > centerAmpCut && centerCharge > centerChargeCut) ) continue;
    if( !(centerTDCx > TDCxmin && centerTDCx < TDCxmax) ) continue;
    if( !(centerTDCy > TDCymin && centerTDCy < TDCymax) ) continue;



    float DeltaTPicoSil[19];
    std::fill(DeltaTPicoSil, DeltaTPicoSil+19, -99);
    // Center pixel:
    DeltaTPicoSil[0] = photekTimeGauss2 - centerTime - meanPicoSil[0];
    histDeltaTPicoSilAt0[0]->Fill(DeltaTPicoSil[0]);
    // Ring 1: Cuts on additional pixels to determine whether they get incorporated.
    for ( int j = 1; j <= 6; j++){
      if ( ring1Atn*amp[j+17] > 0.01 && ring1Atn*integral[j+17] > 1 ) { //j+17--> want elts 18-23 for first ring.
        DeltaTPicoSil[j] = photekTimeGauss2 - linearTime45[j+17] - meanPicoSil[j];
        histDeltaTPicoSilAt0[j]->Fill(DeltaTPicoSil[j]);
      }
    }
    // Ring 2:
    for ( int j = 7; j <= 13; j++){  // Group 3 pixels
      if ( ring2Atn*amp[j+18] > 0.01 && ring2Atn*integral[j+18] > 1 ) { //j+18 --> want elts 25-31 for 2nd ring.
        DeltaTPicoSil[j] = photekTimeGauss3 - linearTime45[j+18] - meanPicoSil[j];
        histDeltaTPicoSilAt0[j]->Fill(DeltaTPicoSil[j]);
      }
    }
    for ( int j = 14; j <= 16; j++){ // Group 0 pixels
      if ( ring2Atn*amp[j-9] > 0.01 && ring2Atn*integral[j-9] > 1 ) { //j-9 --> want elts 5-7 for 2nd ring.
        DeltaTPicoSil[j] = photekTimeGauss0 - linearTime45[j-9] - meanPicoSil[j];
        histDeltaTPicoSilAt0[j]->Fill(DeltaTPicoSil[j]);
      }
    }
    for ( int j = 17; j <= 18; j++){ // Group 1 pixels
      if ( ring2Atn*amp[j-4] > 0.01 && ring2Atn*integral[j-4] > 1 ) { //j-4 --> want elts 13-14 for 2nd ring.
        DeltaTPicoSil[j] = photekTimeGauss1 - linearTime45[j-4] - meanPicoSil[j];
        histDeltaTPicoSilAt0[j]->Fill(DeltaTPicoSil[j]);
      }
    }




    float ChargeTotal = totalPicoSilCharge[0];
    float WeightTotal = totalPicoSilCharge[0] * DeltaTPicoSil[0];
    float ChargeEvent = centerCharge;
    float WeightEvent = centerCharge * DeltaTPicoSil[0];
    float ChargeLandau = MPVlandau[0];
    float WeightLandau = MPVlandau[0] * DeltaTPicoSil[0];
    for (int j = 1; j <= 6; j++) { //Ring 1 -- all Group 2
      if (DeltaTPicoSil[j] != -99.) {
        ChargeTotal += totalPicoSilCharge[j];
        WeightTotal += totalPicoSilCharge[j]   * DeltaTPicoSil[j];
        ChargeEvent += ring1Atn*integral[j+17]; 
        WeightEvent += ring1Atn*integral[j+17] * DeltaTPicoSil[j]; 
        ChargeLandau += MPVlandau[j];
        WeightLandau += MPVlandau[j]           * DeltaTPicoSil[j];
      }
    }
    for (int j = 7; j <= 13; j++) { //Ring 2, Group 3 pixels
      if (DeltaTPicoSil[j] != -99.) {
        ChargeTotal += totalPicoSilCharge[j];
        WeightTotal += totalPicoSilCharge[j]   * DeltaTPicoSil[j];
        ChargeEvent += ring2Atn*integral[j+18];
        WeightEvent += ring2Atn*integral[j+18] * DeltaTPicoSil[j];
        ChargeLandau += MPVlandau[j];
        WeightLandau += MPVlandau[j]           * DeltaTPicoSil[j];
      }
    }
    for (int j = 14; j <= 16; j++) {
      if (DeltaTPicoSil[j] != -99.) { //Ring 2, Group 0 pixels
        ChargeTotal += totalPicoSilCharge[j];
        WeightTotal += totalPicoSilCharge[j]   * DeltaTPicoSil[j];
        ChargeEvent += ring2Atn*integral[j-9];
        WeightEvent += ring2Atn*integral[j-9] * DeltaTPicoSil[j];
        ChargeLandau += MPVlandau[j];
        WeightLandau += MPVlandau[j]           * DeltaTPicoSil[j];
      }
    }
    for (int j = 17; j <= 18; j++) {
      if (DeltaTPicoSil[j] != -99.) { //Ring 2, Group 1 pixels
        ChargeTotal += totalPicoSilCharge[j];
        WeightTotal += totalPicoSilCharge[j]   * DeltaTPicoSil[j];
        ChargeEvent += ring2Atn*integral[j-4];
        WeightEvent += ring2Atn*integral[j-4] * DeltaTPicoSil[j];
        ChargeLandau += MPVlandau[j];
        WeightLandau += MPVlandau[j]           * DeltaTPicoSil[j];
      }
    }




    histDeltaTPicoSilAt0TotalCharge->Fill(WeightTotal/ChargeTotal); // This incorporates the center, ring 1, AND ring 2 pixels. Maybe add an intermediate histogram that just incorporates ring 1?
    histDeltaTPicoSilAt0EventCharge->Fill(WeightEvent/ChargeEvent);
    histDeltaTPicoSilAt0LandauCharge->Fill(WeightLandau/ChargeLandau);
    histDeltaTCenterAt0->Fill(DeltaTPicoSil[0]);
  }

  for (int i = 0; i < 19; i++) std::cout<<"mean "<<i<<":   "<<meanPicoSil[i]<<"\tsubtracted mean: "<<histDeltaTPicoSilAt0[i]->GetMean()<<std::endl;

  // Add Gaussian fit
  Fitter(histDeltaTCenterAt0);
  Fitter(histDeltaTPicoSilAt0TotalCharge);
  Fitter(histDeltaTPicoSilAt0EventCharge);
  Fitter(histDeltaTPicoSilAt0LandauCharge);
  for (int i = 0; i < 19; i++) Fitter(histDeltaTPicoSilAt0[i]);




  // Creates output root file
  TFile *file = TFile::Open(("output"+filename).c_str(), "RECREATE");
  file->cd();
  file->WriteTObject(histDeltaTCenterAt0,"histDeltaTCenter", "WriteDelete");
  file->WriteTObject(histDeltaTPicoSilAt0EventCharge,"histDeltaTPicoSilEventCharge", "WriteDelete");
  file->WriteTObject(histDeltaTPicoSilAt0TotalCharge,"histDeltaTPicoSilTotalCharge", "WriteDelete");
  file->WriteTObject(histDeltaTPicoSilAt0LandauCharge,"histDeltaTPicoSilLandauCharge", "WriteDelete");
  for(int i=0; i<=18; i++) file->WriteTObject(histDeltaTPicoSilAt0[i], Form("histDeltaTPicoSil[%d]",i),"WriteDelete");
  for(int i=0; i<=18; i++) file->WriteTObject(histCharges[i],Form("histCharges[%d]",i),"WriteDelete");
  // Above are in separate loops to be organized in the TBrowser


  gStyle->SetOptFit(1);
  gStyle->SetOptStat(1);
  TH1F *histPhotekAmpCut = new TH1F("histPhotekAmpCut",";Amplitude (mA);Entries/(0.04 mA)", 75, 0, 3);
  TH1F *histPhotekChargeCut = new TH1F("histPhotekChargeCut",";Charge (pC);Entries/(0.4 pC)", 75, 0, 30);
  TH1F *histCenterAmpCut = new TH1F("histCenterAmpCut",";Amplitude (mA);Entries/(0.03 mA)", 100, 0, 3);
  TH1F *histCenterChargeCut = new TH1F("histCenterChargeCut",";Charge (pC);Entries/(1.25 pC)", 80, 0, 150);

  tree->Draw( Form("%f*amp[16]>>histPhotekAmpCut",photekAtn), Form("%f*amp[16]>%f && TDCx > %f && TDCx < %f && TDCy > %f && TDCy < %f",photekAtn,photekAmpCut,TDCxmin,TDCxmax,TDCymin,TDCymax) );
  tree->Draw( Form("%f*int[16]>>histPhotekChargeCut",photekAtn), Form("%f*int[16]>%f && TDCx > %f && TDCx < %f && TDCy > %f && TDCy < %f",photekAtn,photekChargeCut,TDCxmin,TDCxmax,TDCymin,TDCymax));
  tree->Draw( Form("%f*amp[17]>>histCenterAmpCut",centerAtn), Form("%f*amp[17]>%f && TDCx > %f && TDCx < %f && TDCy > %f && TDCy < %f",centerAtn,centerAmpCut,TDCxmin,TDCxmax,TDCymin,TDCymax));
  tree->Draw( Form("%f*int[17]>>histCenterChargeCut",centerAtn), Form("%f*int[17]>%f && TDCx > %f && TDCx < %f && TDCy > %f && TDCy < %f",centerAtn,centerChargeCut,TDCxmin,TDCxmax,TDCymin,TDCymax));

  TH1F *histPhotekAmp = new TH1F("histPhotekAmp",";Amplitude (mA);Entries/(0.04 mA)", 75, 0, 3);
  TH1F *histPhotekCharge = new TH1F("histPhotekCharge",";Charge (pC);Entries/(0.4 pC)", 75, 0, 30);
  TH1F *histCenterAmp = new TH1F("histCenterAmp",";Amplitude (mA);Entries/(0.03 mA)", 100, 0, 3);
  TH1F *histCenterCharge = new TH1F("histCenterCharge",";Charge (pC);Entries/(1.25 pC)", 80, 0, 150);

  tree->Draw( Form("%f*amp[16]>>histPhotekAmp",photekAtn) );
  tree->Draw( Form("%f*int[16]>>histPhotekCharge",photekAtn) );
  tree->Draw( Form("%f*amp[17]>>histCenterAmp",centerAtn) );
  tree->Draw( Form("%f*int[17]>>histCenterCharge",centerAtn) );

  file->WriteTObject(histPhotekAmp, "Photek Amp", "WriteDelete");
  file->WriteTObject(histPhotekAmpCut, "Cut on Photek Amp with X,Y cuts", "WriteDelete");
  file->WriteTObject(histPhotekCharge, "Photek Charge", "WriteDelete");
  file->WriteTObject(histPhotekChargeCut, "Cut on Photek Charge with X,Y cuts", "WriteDelete");
  file->WriteTObject(histCenterAmp, "Center Amp", "WriteDelete");
  file->WriteTObject(histCenterAmpCut, "Cut on Center Amp with X,Y cuts", "WriteDelete");
  file->WriteTObject(histCenterCharge, "Center Charge", "WriteDelete");
  file->WriteTObject(histCenterChargeCut, "Cut on Center Charge with X,Y cuts", "WriteDelete");


  TH1F *histPicoSilAmp = new TH1F("histPicoSilAmp",";Amplitude (mA);Entries/(0.04 mA)", 75, -0.5, 3.5);
  TH1F *histPicoSilCharge = new TH1F("histPicoSilCharge",";Charge (pC);Entries/(0.4 pC)", 75, -50, 250);
  tree->Draw("sqrt(10)*amp[17]+sqrt(10)*amp[18]+sqrt(10)*amp[19]+sqrt(10)*amp[20]+sqrt(10)*amp[21]+sqrt(10)*amp[22]+sqrt(10)*amp[23]+amp[25]+amp[26]+amp[27]+amp[28]+amp[29]+amp[30]+amp[31]+amp[5]+amp[6]+amp[7]+amp[13]+amp[14]>>histPicoSilAmp");
  tree->Draw("sqrt(10)*int[17]+sqrt(10)*int[18]+sqrt(10)*int[19]+sqrt(10)*int[20]+sqrt(10)*int[21]+sqrt(10)*int[22]+sqrt(10)*int[23]+int[25]+int[26]+int[27]+int[28]+int[29]+int[30]+int[31]+int[5]+int[6]+int[7]+int[13]+int[14]>>histPicoSilCharge");
  file->WriteTObject(histPicoSilAmp, "Center + Ring 1 + Ring 2 Amps--no cuts", "WriteDelete");
  file->WriteTObject(histPicoSilCharge, "Center + Ring 1 + Ring 2 Charges--no cuts", "WriteDelete");

 
  file->Close();
  delete file;

}




void PlotDeltaTPDF(TCanvas *c, TLatex *tex, TH1F *hist, string outfile) {
  hist->Draw();
  gStyle->SetOptFit(0); //Hides the parameter box
  gStyle->SetOptStat(0);
  double mean = hist->GetMean();
  double rms = hist->GetRMS();
  TF1 *gausfit = new TF1("gausfit","gaus", mean - 2.0*rms, mean + 2.0*rms);//1-D gaus function defined around hist peak
  hist->Fit("gausfit","QMLES","", mean - 2.0*rms, mean + 2.0*rms);// Fit the hist; Q-quiet, L-log likelihood method, E-Minos errors technique, M-improve fit results
  //hist->GetXaxis()->SetTitle("#Deltat (ns)");
  if(1000*gausfit->GetParError(2)>2) tex->DrawLatex(0.59, 0.83, Form("#sigma = %.0f #pm %.0f ps", 1000*gausfit->GetParameter(2), 1000*gausfit->GetParError(2)));
  else tex->DrawLatex(0.59, 0.83, Form("#sigma = %.1f #pm %.1f ps", 1000*gausfit->GetParameter(2), 1000*gausfit->GetParError(2)));
  c->SaveAs(outfile.c_str()); //outfile should end in .pdf

  //if(outfile == "deltaTPicoSilLandauCharge.pdf") c->SaveAs("deltaTPicoSilLandauCharge.C");
}


void makePlots( string filename ) {

  gStyle->SetOptFit(1);
  gStyle->SetOptStat(1);

  TFile *_file = TFile::Open( ("output"+filename).c_str() ); //Should be .root

  //Create variables containing hists:
  TH1F *histDeltaTCenter = (TH1F*)_file->Get("histDeltaTCenter"); //Picosil center pixel
  TH1F *histDeltaTPicoSilTotalCharge = (TH1F*)_file->Get("histDeltaTPicoSilTotalCharge");//Picosil, total charge
  TH1F *histDeltaTPicoSilEventCharge = (TH1F*)_file->Get("histDeltaTPicoSilEventCharge");//Picosil, event charge
  TH1F *histDeltaTPicoSilLandauCharge= (TH1F*)_file->Get("histDeltaTPicoSilLandauCharge");//Picosil, landau charge
  TH1F *histPhotekAmp = (TH1F*)_file->Get("Photek Amp");
  TH1F *histPhotekAmpCut = (TH1F*)_file->Get("Cut on Photek Amp with X,Y cuts");
  TH1F *histPhotekCharge = (TH1F*)_file->Get("Photek Charge");
  TH1F *histPhotekChargeCut = (TH1F*)_file->Get("Cut on Photek Charge with X,Y cuts");
  TH1F *histDeltaTPicoSil[19];
  for(int i=0; i<19; i++) histDeltaTPicoSil[i] = (TH1F*)_file->Get( Form("histDeltaTPicoSil[%d]",i) ); //Already wrote center pixel


  TCanvas *c = new TCanvas ("c","c",800, 600); 
  TLatex *tex = new TLatex();
  tex->SetNDC(); // Sets coords such that (0,0) is bottom left & (1,1) is top right.
  tex->SetTextSize(0.060);
  tex->SetTextFont(42);
  tex->SetTextColor(kBlack);

  c->cd();

  PlotDeltaTPDF(c, tex, histDeltaTCenter, "deltaTCenter.pdf");
  PlotDeltaTPDF(c, tex, histDeltaTPicoSilEventCharge, "deltaTPicoSilEventCharge.pdf");
  PlotDeltaTPDF(c, tex, histDeltaTPicoSilTotalCharge, "deltaTPicoSilTotalCharge.pdf");
  PlotDeltaTPDF(c, tex, histDeltaTPicoSilLandauCharge, "deltaTPicoSilLandauCharge.pdf");
  for(int i=0; i<19; i++) PlotDeltaTPDF(c, tex, histDeltaTPicoSil[i], Form("deltaTPicoSilPixel%d.pdf",i+1) );

  c->SetLogy();
  histPhotekAmp->Draw();
  c->SaveAs( "PhotekAmp.pdf" );
  histPhotekAmpCut->Draw();
  c->SaveAs( "PhotekAmpCut.pdf" );
  histPhotekCharge->Draw();
  c->SaveAs( "PhotekCharge.pdf" );
  histPhotekChargeCut->Draw();
  c->SaveAs( "PhotekChargeCut.pdf" );

  c->Close();
}


void PicosilStudy() {
  gStyle->SetTitleOffset(0.8,"x");
  gStyle->SetTitleOffset(0.85,"y");
  gStyle->SetTitleSize(0.055,"x");
  gStyle->SetTitleSize(0.055,"y");
  gStyle->SetLabelSize(0.045,"x");
  gStyle->SetLabelSize(0.045,"y");
  gStyle->SetOptFit(1);
  gStyle->SetOptStat(1);

  string infile = "analysis_5488.root";
  float photekAmpCut = sqrt(10)*0.01; //THESE ARE THE CUT VALUES AFTER ADJUSTING FOR ATTENUATORS
  float photekChargeCut = 0;  // sqrt(10)*0.1;
  float centerAmpCut = sqrt(10)*0.01;
  float centerChargeCut = 0;  //sqrt(10)*1.0;
  float TDCxmin = 4;  // Can write a function to calculate optimal region.
  float TDCxmax = 14;  // Leaving as manual inputs for now.
  float TDCymin = -3;
  float TDCymax = 9;
  float centerAtn = sqrt(10); // 10dB attenuator
  float ring1Atn = sqrt(10);
  float ring2Atn = 1; // No attenuator
  float photekAtn= sqrt(10);

  makeOutputFile( infile.c_str(), photekAmpCut, photekChargeCut, centerAmpCut, centerChargeCut, TDCxmin, TDCxmax, TDCymin, TDCymax, centerAtn, ring1Atn, ring2Atn, photekAtn );//Outputs a ROOT file with histograms
  makePlots( infile.c_str() ); // Outputs PDF files with histograms
}
