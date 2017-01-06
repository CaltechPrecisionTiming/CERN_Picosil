
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

// This code generates the combined Delta T histogram for the picosil (center pixel and first ring of pixels) and the Photonis MCP for the re-cabled configuration **with respect to the Photek**. The SiPad is not used. 
// An update to the code includes use of the first ring pixels in the picosil.
// Another update gives the Delta T histogram between the picosil and the Photonis MCP.
// Another update gives the SKIROC emulation.
// Author: Daniel Gawerc


// Histogram names start with histDeltaT. Then they say the devices they incorporate: either MCP, PicoSil center pixel, or PicoSil with all pixels. The PicoSil with all the Pixels will either specify equal, total charge, or event charge; this signifies all pixels weighted equally (1/7), or by charge (either weighting differently event-by-event or the same weights using total charge in the each pixel throughout the run). Then the histogram will specify equal, total charge, or event charge, with the same meanings as above, except for weighting the PicoSil delta T against the MCP delta T. Finally, a last option No Shift indicates the component histograms hadn't been shifted to zero by having their means subtracted, which should result in a histogram with a high timing resolution.
// Having At0 indicates the mean has been subtracted from the original hist.

void Fitter(TH1F *hist) {
  //Helper function for fitting Gaussian
  double xmin = hist->GetMean() - 2.0*hist->GetRMS();
  double xmax = hist->GetMean() + 2.0*hist->GetRMS();
  hist->Fit("gaus","QMLES","",xmin,xmax); // Q suppresses fit results
  gStyle->SetOptFit(0);
}

void SKIROCPlotPDF(TCanvas *c, TLatex *tex, TH1F *hist, string outfile, int *pixelsUsed) {

  // Draw hist; Remove Stat and Fit parameter boxes:
  hist->Draw();
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  // Add text box:
  TPaveText *txtbox = new TPaveText(0.7, 0.48, 0.9, 0.78, "NDC");
  txtbox->SetBorderSize(1);
  txtbox->SetFillColor(10);
  txtbox->AddText("Pixels Passing Cuts");
  txtbox->AddLine(0., 0.875, 1., 0.875);
  for (int i = 0; i<7; i++) txtbox->AddText( Form("%d Pixel(s):  %d Events", i+1, pixelsUsed[i]) );
  // Add Sigma Parameter text:
  double mean = hist->GetMean();
  double rms = hist->GetRMS();
  TF1 *gausfit = new TF1("gausfit","gaus", mean - 2.0*rms, mean + 2.0*rms);//1-D gaus function defined around hist peak
  hist->Fit("gausfit","QMLES","", mean - 2.0*rms, mean + 2.0*rms);// Fit the hist; Q-quiet, L-log likelihood method, E-Minos errors technique, M-improve fit results
  hist->GetXaxis()->SetTitle("#Deltat (ns)");
  txtbox->Draw();
  if(1000*gausfit->GetParError(2)>2) tex->DrawLatex(0.59, 0.83, Form("#sigma = %.0f #pm %.0f ps", 1000*gausfit->GetParameter(2), 1000*gausfit->GetParError(2)));
  else tex->DrawLatex(0.59, 0.83, Form("#sigma = %.1f #pm %.1f ps", 1000*gausfit->GetParameter(2), 1000*gausfit->GetParError(2)));
  c->SaveAs( (outfile + ".pdf").c_str() ); 
  c->SaveAs( (outfile +   ".C").c_str() );
}


void DoMultiDeviceStudy( string filename, float photekAmpCut, float photekChargeCut, float centerAmpCut, float centerChargeCut, float MCPAmpCut) {

  srand (time(NULL));
  int seed = rand();

  TFile *inputfile = TFile::Open(filename.c_str(),"READ");
  TTree *tree = (TTree*)inputfile->Get("pulse");

  // get the variables from the ntuple
  float amp[36];
  float integral[36];
  float gauspeak[36];
  float linearTime45[36];

  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("gauspeak",1);
  tree->SetBranchStatus("linearTime45",1);
  tree->SetBranchStatus("amp",1);
  tree->SetBranchStatus("int",1);
  
  tree->SetBranchAddress("gauspeak",gauspeak);
  tree->SetBranchAddress("linearTime45",linearTime45);
  tree->SetBranchAddress("amp",amp);
  tree->SetBranchAddress("int",integral);

  //Create histograms
  float width = 0.15;
  float smearWidth = 0.3;
  int bins = 75;
  int smearBins = 75;
  float pixelSmear = 0.050; // in ns
  float MCPSmear = 0.045;

  TH1F *histDeltaT_Center_MCP_Equal = new TH1F("histDeltaT_Center_MCP_Equal",";#Deltat (ns);Entries/(0.004 ns)", bins, -width, width); // Weights MCP and PicoSil center pixel 50-50
  TH1F *histDeltaT_PicoSilEventCharge_MCP_Equal = new TH1F("histDeltaT_PicoSilEventCharge_MCP_Equal",";#Deltat (ns);Entries/(0.004 ns)", bins, -width, width);// Picosil delta T found by weighting with event pixel charge, but then overall delta T fund by weighting PicoSil and MCP equally.
  TH1F *histDeltaT_PicoSilTotalCharge_MCP_Equal = new TH1F("histDeltaT_PicoSilTotalCharge_MCP_Equal",";#Deltat (ns);Entries/(0.004 ns)", bins, -width, width);// Picosil delta T found by weighting with total pixel charge, but then overall delta T fund by weighting PicoSil and MCP equally.
  TH1F *histDeltaT_PicoSilLandauCharge_MCP_Equal = new TH1F("histDeltaT_PicoSilLandauCharge_MCP_Equal",";#Deltat (ns);Entries/(0.004 ns)", bins, -width, width);
  TH1F *histDeltaT_PicoSilLandauCharge_MCP_Equal_PicoSilSmear = new TH1F("histDeltaT_PicoSilLandauCharge_MCP_Equal_PicoSilSmear",";#Deltat (ns);Entries/(0.008 ns)", smearBins, -smearWidth, smearWidth);//Smear
  TH1F *histDeltaT_PicoSilLandauCharge_MCP_Equal_BothSmear = new TH1F("histDeltaT_PicoSilLandauCharge_MCP_Equal_BothSmear",";#Deltat (ns);Entries/(0.008 ns)", smearBins, -smearWidth, smearWidth);//Smear
  TH1F *histDeltaT_PicoSilEqual_MCP_Equal = new TH1F("histDeltaT_PicoSilEqual_MCP_Equal",";#Deltat (ns);Entries/(0.004 ns)", bins, -width, width);//Picosil delta T found by weighting pixels equally, and then overall delta T fund by weighting PicoSil and MCP equally. Thus MCP is weighted by 1/2 and each pixel is weighted 1/14.
  TH1F *histDeltaT_PicoSilEqual_MCP_Equal_BothSmear = new TH1F("histDeltaT_PicoSilEqual_MCP_Equal_BothSmear",";#Deltat (ns);Entries/(0.008 ns)", smearBins, -smearWidth, smearWidth);
  TH1F *histDeltaT_PicoSil_MCP_Equal = new TH1F("histDeltaT_PicoSil_MCP_Equal",";#Deltat (ns);Entries/(0.004 ns)", bins, -width, width);// delta T found by placing equal emphasis on the pixels as on the MCP. Thus everything is weighted 1/8.
  TH1F *histDeltaT_PicoSil_MCP_EventCharge = new TH1F("histDeltaT_PicoSil_MCP_EventCharge",";#Deltat (ns);Entries/(0.004 ns)", bins, -width, width);// delta T found by weighting with event charge in each device.
  TH1F *histDeltaT_PicoSil_MCP_TotalCharge = new TH1F("histDeltaT_PicoSil_MCP_TotalCharge",";#Deltat (ns);Entries/(0.004 ns)", bins, -width, width);// delta T found by constant weighting with charge in device over entire run.
  TH1F *histDeltaT_Center_MCP_EventCharge_NoShift = new TH1F("histDeltaT_Center_MCP_EventCharge_NoShift",";#Deltat (ns);Entries/(0.04 ns)", 100, 2, 6); //Weights each event based on charge. Poor resolution.
  TH1F *histDeltaTPicoSil[7];
  TH1F *histDeltaTPicoSilAt0[7];
  TH1F *histDeltaTPicoSilSmear[7];
  TH1F *histDeltaTPicoSilSmearAt0[7];
  for(int i=0; i<7; i++) {
    histDeltaTPicoSil[i] = new TH1F(Form("histDeltaTPicoSil_%d",i),";#Deltat (ns);Entries/(0.06 ns)", 50, 3, 6); //DeltaT of PicoSil pixels
    histDeltaTPicoSilAt0[i] = new TH1F(Form("histDeltaTPicoSilAt0_%d",i),";#Deltat (ns);Entries/(0.008 ns)", bins/2, -width, width); 
    histDeltaTPicoSilSmear[i] = new TH1F(Form("histDeltaTPicoSilSmear_%d",i),";#Deltat (ns);Entries/(0.06 ns)", 50, 3, 6);
    histDeltaTPicoSilSmearAt0[i] = new TH1F(Form("histDeltaTPicoSilSmearAt0_%d",i),";#Deltat (ns);Entries/(0.008 ns)", smearBins, -smearWidth, smearWidth); 
  }
  TH1F *histDeltaTCenter = new TH1F("histDeltaTCenter",";#Deltat (ns);Entries/(0.02 ns)", 50, 4, 5); //DeltaT of center picosil pixel
  TH1F *histDeltaTMCP = new TH1F("histDeltaTMCP",";#Deltat (ns);Entries/(0.02 ns)", 50, 2, 3); //DeltaT of MCP
  TH1F *histDeltaTMCPSmear = new TH1F("histDeltaTMCPSmear",";#Deltat (ns);Entries/(0.02 ns)", 50, 2, 3); //DeltaT of MCP  
  TH1F *histDeltaTPicoSilAt0TotalCharge = new TH1F("histDeltaTPicoSilAt0TotalCharge",";#Deltat_{HGC} (ns);Entries/(0.004 ns)", bins, -width, width); //All pixels combined
  TH1F *histDeltaTPicoSilAt0EventCharge = new TH1F("histDeltaTPicoSilAt0EventCharge",";#Deltat_{HGC} (ns);Entries/(0.004 ns)", bins, -width, width);
  TH1F *histDeltaTPicoSilAt0LandauCharge = new TH1F("histDeltaTPicoSilAt0LandauCharge",";#Deltat_{HGC} (ns);Entries/(0.004 ns)", bins, -width, width);// uses charge MPV
  TH1F *histDeltaTPicoSilAt0LandauChargeSmear = new TH1F("histDeltaTPicoSilAt0LandauChargeSmear", ";#Deltat_{HGC} (ns);Entries/(0.008 ns)", smearBins, -smearWidth, smearWidth);//SKIROC
  TH1F *histDeltaTPicoSilAt0EqualSmear = new TH1F("histDeltaTPicoSilAt0EqualSmear", ";#Deltat_{HGC} (ns);Entries/(0.008 ns)", smearBins, -smearWidth, smearWidth);//SKIROC
  TH1F *histDeltaTPicoSilAt0EqualSmear_nEventsCombine[7];//Central pixel smear is [0], central + second highest nEvents is [1], etc...
  for(int i=0; i<7; i++) histDeltaTPicoSilAt0EqualSmear_nEventsCombine[i] = 
      new TH1F(Form("histDeltaTPicoSilAt0EqualSmear_nEventsCombine_%dPixels",i+1), ";#Deltat (ns);Entries/(0.008 ns)", smearBins, -smearWidth, smearWidth);//SKIROC
  TH1F *histDeltaTCenterAt0 = new TH1F("histDeltaTCenterAt0",";#Deltat (ns);Entries/(0.004 ns)", bins, -width, width); //shifted to be centered at zero
  TH1F *histDeltaTMCPAt0 = new TH1F("histDeltaTMCPAt0",";#Deltat (ns);Entries/(0.004 ns)", bins, -width, width); //shifted to be centered at zero
  TH1F *histDeltaTMCPAt0Smear = new TH1F("histDeltaTMCPAt0Smear",";#Deltat (ns);Entries/(0.008 ns)", smearBins, -smearWidth, smearWidth); //shifted to be centered at zero

  TH1F *histDeltaT_PicoSil_vs_MCP_EventCharge = new TH1F("histDeltaT_PicoSil_vs_MCP_EventCharge",";#Deltat (ns);Entries/(0.005 ns)", 80, -0.2, 0.2);
  TH1F *histDeltaT_PicoSil_vs_MCP_TotalCharge = new TH1F("histDeltaT_PicoSil_vs_MCP_TotalCharge",";#Deltat (ns);Entries/(0.005 ns)", 80, -0.2, 0.2);
  TH1F *histDeltaT_PicoSil_vs_MCP[7];
  for(int i=0; i<7; i++) histDeltaT_PicoSil_vs_MCP[i] = new TH1F(Form("histDeltaT_PicoSil_vs_MCP_%d",i),";#Deltat (ns);Entries/(0.02 ns)", 100, -3, -1); // DeltaT between PicoSil and MCP instead of Photek.


  TH1F *histCharges[7]; // collects charge values for picosil pixels in every event in which they pass the cuts.
  for(int i=0; i<7; i++) histCharges[i] = new TH1F( Form("histCharges_%d",i),";Charge (pC);Entries/(1.6 pC)", 50, 0, 80);
  //TH1F *histAllCharges = new TH1F("histAllCharges",";Charge (pC);Entries/(1.6 pC)",50,0,80); // histogram filled with charge values from every pixel *after* cuts. Basically just adding all histCharges[i].
  //Above hist removed because cutting at 20pC is basically just histCharges[0] -->The center pixel charge dist 




  float totalPicoSilCharge[7] = {0.};
  float totalCenterCharge = 0; // Should be same value as totalPicoSilCharge[0]
  float totalMCPCharge = 0;


  //read all entries and fill the histogram
  Long64_t nentries = tree->GetEntries();

  //Loop through every event in .root file
  TRandom3 *rando = new TRandom3(seed);
  std::cout<<"Number of events in Sample: "<<nentries<<std::endl;  
  for (Long64_t iEntry=0; iEntry<nentries; iEntry++) {
    if (iEntry %1000 == 0) cout << "Processing Event " << iEntry << "\n";
    tree->GetEntry(iEntry);    
  
    float photekTimeGauss0 = gauspeak[0];
    float photekAmp0 = sqrt(10)*amp[0]; //accounts for 10dB attenuator
    float photekCharge0 = sqrt(10)*integral[0];

    float photekTimeGauss1 = gauspeak[9];
    float photekAmp1 = sqrt(10)*amp[9];
    float photekCharge1 = sqrt(10)*integral[9];

    float centerAmp = 2*amp[1]; //accounts for 6dB attenuator
    float centerCharge = 2*integral[1];
    float centerTime = linearTime45[1];

    float MCPAmp = amp[11];
    float MCPCharge = integral[11];
    float MCPTime = linearTime45[11];


    // APPLY EVENT CUTS:
    //require photek (for electron selection)
    if( !(photekAmp0 > photekAmpCut && photekCharge0 > photekChargeCut) ) continue;
    //require signal in the central pixel
    if( !(centerAmp > centerAmpCut && centerCharge > centerChargeCut) ) continue;
    //require MCP minimum amplitude. Possibly will include cut for charge at some point...
    if( !( MCPAmp > MCPAmpCut) ) continue;

    
    double linearTime45Smear[7];
    for (int j = 0; j < 7; j++)  linearTime45Smear[j] = rando->Gaus(linearTime45[j+1], pixelSmear); //Samples from smear
    double MCPTimeSmear = rando->Gaus(MCPTime, MCPSmear);


    //Calculates the Delta T's if the event passes the cuts:
    // (Will be used to fill the histogram at every event)
    float DeltaTCenter = photekTimeGauss0 - centerTime; 
    float DeltaTMCP = photekTimeGauss1 - MCPTime;
    float DeltaTMCPSmear = photekTimeGauss1 - MCPTimeSmear;
    float DeltaTPicoSil[7]; 
    float DeltaTPicoSilSmear[7];
    float DeltaTPicoSil_vs_MCP[7];
    std::fill(DeltaTPicoSil, DeltaTPicoSil+7, -99);
    std::fill(DeltaTPicoSilSmear, DeltaTPicoSilSmear+7, -99);
    std::fill(DeltaTPicoSil_vs_MCP, DeltaTPicoSil_vs_MCP+7, -99);
    for ( int j = 1; j <= 7; j++){
      if ( (amp[j] > 0.01 && integral[j] > 1) || j == 1 ) {
        DeltaTPicoSil[j-1] = photekTimeGauss0 - linearTime45[j];
        DeltaTPicoSilSmear[j-1] = photekTimeGauss0 - linearTime45Smear[j-1];
        DeltaTPicoSil_vs_MCP[j-1] = (linearTime45[j] - photekTimeGauss0) - (MCPTime - photekTimeGauss1); //subtracting photek accts for MCP and HGC in diff groups
	if ( j == 1 ) {
	  totalPicoSilCharge[j-1] += centerCharge;
	  histCharges[j-1]->Fill(centerCharge);
	}
	else {
	  totalPicoSilCharge[j-1] += integral[j];
	  histCharges[j-1]->Fill(integral[j]);
	}
      }
    }

    float DeltaT_Center_MCP_EventCharge_NoShift = (DeltaTCenter*centerCharge + DeltaTMCP*MCPCharge) / (centerCharge+MCPCharge);
    //The above weights each individual event being added by the charge of each event. The resolution will be large, as the mean should be subtracted from each of the events prior to weighting them by charge.

    totalCenterCharge += centerCharge;
    totalMCPCharge += MCPCharge;

    histDeltaT_Center_MCP_EventCharge_NoShift->Fill(DeltaT_Center_MCP_EventCharge_NoShift);
    histDeltaTCenter->Fill(DeltaTCenter);
    histDeltaTMCP->Fill(DeltaTMCP);
    histDeltaTMCPSmear->Fill(DeltaTMCPSmear);
    for(int k=0; k<7; k++) {
      if (DeltaTPicoSil[k] != -99.) {
        histDeltaTPicoSil[k]->Fill(DeltaTPicoSil[k]);
        histDeltaTPicoSilSmear[k]->Fill(DeltaTPicoSilSmear[k]);
        histDeltaT_PicoSil_vs_MCP[k]->Fill(DeltaTPicoSil_vs_MCP[k]);
      }
    }
  }


  float ringChargeTotal = 0;
  for(int i=1; i<=6; i++) ringChargeTotal += totalPicoSilCharge[i];

  // Do Gaussian fit of delta T distributions from (mean-2RMS) to (mean+2RMS)
  Fitter(histDeltaT_Center_MCP_EventCharge_NoShift); 
  Fitter(histDeltaTCenter);
  Fitter(histDeltaTMCP);
  Fitter(histDeltaTMCPSmear);

  double meanCenter = histDeltaTCenter->GetMean();
  double meanMCP = histDeltaTMCP->GetMean();
  double meanMCPSmear = histDeltaTMCPSmear->GetMean();

  TF1 *flandau[7];
  double MPVlandau[7]; // Seeing if we can get better results by using Landau mean as weighting charge value.
  for(int i=1; i<=6; i++) {
  double mean = histCharges[i]->GetMean();
  double rms = histCharges[i]->GetRMS();
  double xmin = mean-2.0*rms;
  double xmax = mean+2.0*rms;
  flandau[i] = new TF1( Form("flandau_%d",i), "landau", xmin, xmax); // 1-D landau func
  histCharges[i]->Fit( Form("flandau_%d",i), "QMLES","", xmin, xmax);
  gStyle->SetOptFit(0);
  MPVlandau[i] = flandau[i]->GetMaximumX(); // In order to find MPV.
  }
  double mean = histCharges[0]->GetMean();
  double rms = histCharges[0]->GetRMS();
  double xmin = mean-2.0*rms;
  double xmax = mean+2.0*rms;
  flandau[0] = new TF1("flandau_0","gaus", xmin, xmax); // storing a Gaus fit for the central pixel with the other landau fits
  histCharges[0]->Fit("flandau_0","QMLES","", xmin, xmax);
  gStyle->SetOptFit(0);
  MPVlandau[0] = flandau[0]->GetParameter(1);
  cout<<"\nhistCharges[0]\nGauss Mean: "<<flandau[0]->GetParameter(1)<<
                        "\nUncertainy: "<<flandau[0]->GetParError(1)<<
  "\n"<<endl;


  double meanPicoSil[7];
  double meanPicoSilSmear[7];
  int DeltaTPicoSilSmear_Events[7];
  double meanPicoSil_vs_MCP[7];
  for (int i = 0; i < 7; i++){
    meanPicoSil[i] = histDeltaTPicoSil[i]->GetMean();
    meanPicoSilSmear[i] = histDeltaTPicoSilSmear[i]->GetMean();
    DeltaTPicoSilSmear_Events[i] = histDeltaTPicoSilSmear[i]->GetEntries();
    meanPicoSil_vs_MCP[i] = histDeltaT_PicoSil_vs_MCP[i]->GetMean();
  }
  int DeltaTPicoSilSmear_Events_Sorted[7];
  int DeltaTPicoSilSmear_Events_SortedIndices[7];
  std::copy(DeltaTPicoSilSmear_Events, DeltaTPicoSilSmear_Events+7, DeltaTPicoSilSmear_Events_Sorted);
  std::sort(DeltaTPicoSilSmear_Events_Sorted, DeltaTPicoSilSmear_Events_Sorted+7, std::greater<int>() );
  for (int i = 0; i < 7; i++){
    for (int j = 0; j < 7; j++) {
      if( DeltaTPicoSilSmear_Events[j] == DeltaTPicoSilSmear_Events_Sorted[i] ) {
        DeltaTPicoSilSmear_Events_SortedIndices[i] = j;
        break;
      }
    }
  }

  int pixelsUsed[7][7] = {0}; // [nPixels added] [nPixels actually used]




  TRandom3 *rando2 = new TRandom3(seed);
  //Loop through again. This time, subtracting the respective means. Process is almost the same.
  for (Long64_t iEntry=0; iEntry<nentries; iEntry++) {
    if (iEntry %1000 == 0) cout << "Processing Event " << iEntry << "\n";
    tree->GetEntry(iEntry);    
  
    float photekTimeGauss0 = gauspeak[0];
    float photekAmp0 = sqrt(10)*amp[0]; //accounts for 10dB attenuator
    float photekCharge0 = sqrt(10)*integral[0];

    float photekTimeGauss1 = gauspeak[9];
    float photekAmp1 = sqrt(10)*amp[9];
    float photekCharge1 = sqrt(10)*integral[9];

    float centerAmp = 2*amp[1]; //accounts for 6dB attenuator
    float centerCharge = 2*integral[1];
    float centerTime = linearTime45[1];

    float MCPAmp = amp[11];
    float MCPCharge = integral[11];
    float MCPTime = linearTime45[11];

    if( !(photekAmp0 > photekAmpCut && photekCharge0 > photekChargeCut) ) continue;
    if( !(centerAmp > centerAmpCut && centerCharge > centerChargeCut) ) continue;
    if( !( MCPAmp > MCPAmpCut) ) continue;

    double linearTime45Smear[7];
    for (int j = 0; j < 7; j++)  linearTime45Smear[j] = rando2->Gaus(linearTime45[j+1], pixelSmear);
    double MCPTimeSmear = rando2->Gaus(MCPTime, MCPSmear);

    float DeltaTPicoSil[7]; 
    float DeltaTPicoSilSmear[7];
    float DeltaTPicoSil_vs_MCP[7];
    std::fill(DeltaTPicoSil, DeltaTPicoSil+7, -99);
    std::fill(DeltaTPicoSilSmear, DeltaTPicoSilSmear+7, -99);
    std::fill(DeltaTPicoSil_vs_MCP, DeltaTPicoSil_vs_MCP+7, -99);
    for ( int j = 1; j <= 7; j++){
      if ( (amp[j] > 0.01 && integral[j] > 1) || j == 1 ) {
        DeltaTPicoSil[j-1] = photekTimeGauss0 - linearTime45[j] - meanPicoSil[j-1];
        histDeltaTPicoSilAt0[j-1]->Fill(DeltaTPicoSil[j-1]);
        DeltaTPicoSilSmear[j-1] = photekTimeGauss0 - linearTime45Smear[j-1] - meanPicoSilSmear[j-1];
        histDeltaTPicoSilSmearAt0[j-1]->Fill(DeltaTPicoSilSmear[j-1]);
        DeltaTPicoSil_vs_MCP[j-1] = (linearTime45[j] - photekTimeGauss0) - (MCPTime - photekTimeGauss1) - meanPicoSil_vs_MCP[j-1];
      }
    }


    float ringWeightTotal = 0;
    float ringWeightEvent = 0;
    float ringChargeEvent = 0;
    for(int ii = 2; ii <= 7; ii++){
      if(DeltaTPicoSil[ii-1] != -99.) {
        ringChargeEvent += integral[ii];
        ringWeightTotal += totalPicoSilCharge[ii-1] * DeltaTPicoSil[ii-1];
        ringWeightEvent += integral[ii] * DeltaTPicoSil[ii-1];
      }
    }


    float PicoSil_vs_MCP_WeightTotal = 0;
    float PicoSil_vs_MCP_WeightEvent = 0;
    float PicoSil_vs_MCP_ChargeTotal = 0;
    float PicoSil_vs_MCP_ChargeEvent = 0;
    for (int jj = 1; jj <= 7; jj++) {
      if (DeltaTPicoSil_vs_MCP[jj-1] != -99.) {
        PicoSil_vs_MCP_WeightTotal += totalPicoSilCharge[jj-1] * DeltaTPicoSil_vs_MCP[jj-1];
        PicoSil_vs_MCP_ChargeTotal += totalPicoSilCharge[jj-1];
      	if ( jj != 1 ) {
          PicoSil_vs_MCP_ChargeEvent += integral[jj];
          PicoSil_vs_MCP_WeightEvent += integral[jj] * DeltaTPicoSil_vs_MCP[jj-1];
        }
      	else {
          PicoSil_vs_MCP_ChargeEvent += centerCharge;
          PicoSil_vs_MCP_WeightEvent += centerCharge * DeltaTPicoSil_vs_MCP[jj-1];
        }
      }
    }

    
    //Here are the corrections that center the histograms at 0.
    float DeltaTCenter = photekTimeGauss0 - centerTime - meanCenter; 
    float DeltaTMCP = photekTimeGauss1 - MCPTime - meanMCP;
    float DeltaTMCPSmear = photekTimeGauss1 - MCPTimeSmear - meanMCPSmear;

    float DeltaT_Center_MCP_Equal = (DeltaTCenter + DeltaTMCP)/2; //MCP and center pixel evenly weighted
    float DeltaT_PicoSil_MCP_EventCharge = (DeltaTCenter*centerCharge + DeltaTMCP*MCPCharge + ringWeightEvent) / (centerCharge+MCPCharge+ringChargeEvent);
    float DeltaT_PicoSil_MCP_TotalCharge = (DeltaTCenter*totalCenterCharge + DeltaTMCP*totalMCPCharge + ringWeightTotal) / (totalCenterCharge+totalMCPCharge+ringChargeTotal);

    float DeltaT_PicoSilTotalCharge_MCP_Equal = 0.5*DeltaTMCP; //Add PicoSil on next line.
    float temp_pstotalweight = 0;
    float temp_pstotalcharge = 0;
    for (int j = 0; j <= 6; j++) {
      if (DeltaTPicoSil[j] != -99.) {
        temp_pstotalcharge += totalPicoSilCharge[j];
        temp_pstotalweight += DeltaTPicoSil[j]*totalPicoSilCharge[j];
      }
    }
    DeltaT_PicoSilTotalCharge_MCP_Equal += 0.5*temp_pstotalweight/temp_pstotalcharge;

    float DeltaT_PicoSilEventCharge_MCP_Equal = 0.5*DeltaTMCP;
    float temp_pseventweight = 0;
    float temp_pseventcharge = 0;
    if (DeltaTPicoSil[0] != -99.) {
      temp_pseventweight += DeltaTPicoSil[0]*centerCharge;
      temp_pseventcharge += centerCharge;
    }
    for (int j = 1; j <= 6; j++) {
      if (DeltaTPicoSil[j] != -99) {
        temp_pseventweight += DeltaTPicoSil[j]*integral[j+1];
        temp_pseventcharge += integral[j+1];
      }
    }
    DeltaT_PicoSilEventCharge_MCP_Equal += 0.5*temp_pseventweight/temp_pseventcharge;

    float DeltaT_PicoSilLandauCharge_MCP_Equal = 0.5*DeltaTMCP;
    float temp_pslandauweight = 0;
    float temp_pslandaucharge = 0;
    for (int j = 0; j <= 6; j++) {
      if (DeltaTPicoSil[j] != -99.) {
        temp_pslandauweight += DeltaTPicoSil[j]*MPVlandau[j];
        temp_pslandaucharge += MPVlandau[j];
      }
    }
    DeltaT_PicoSilLandauCharge_MCP_Equal += 0.5*temp_pslandauweight/temp_pslandaucharge;

    float DeltaT_PicoSilLandauCharge_MCP_Equal_PicoSilSmear = 0.5*DeltaTMCP; // MCP, photek not smeared
    float temp_pslandauweight_smear = 0;
    float temp_pslandaucharge_smear = 0;
    for (int j = 0; j <= 6; j++) {
      if (DeltaTPicoSilSmear[j] != -99.) {
        temp_pslandauweight_smear += DeltaTPicoSilSmear[j]*MPVlandau[j];
        temp_pslandaucharge_smear += MPVlandau[j];
      }
    }
    DeltaT_PicoSilLandauCharge_MCP_Equal_PicoSilSmear += 0.5*temp_pslandauweight_smear/temp_pslandaucharge_smear;


    float DeltaT_PicoSilLandauCharge_MCP_Equal_BothSmear = 0.5*temp_pslandauweight_smear/temp_pslandaucharge_smear + 0.5*DeltaTMCPSmear; // Photek not smeared


    float DeltaT_PicoSil_MCP_Equal = DeltaTMCP;
    int inc = 0; //Divide by inc at the end, because some PicoSil pixels may not have passed cuts.
    for (int j = 0; j <= 6; j++) {
      if (DeltaTPicoSil[j] != -99.) {
        DeltaT_PicoSil_MCP_Equal += DeltaTPicoSil[j];
        inc += 1;
      }
    }
    DeltaT_PicoSil_MCP_Equal /= inc;

    float DeltaT_PicoSilEqual_MCP_Equal = 0.5*DeltaTMCP;
    inc = 0;
    float temp_psdeltat = 0;
    for (int j = 0; j <= 6; j++) {
      if (DeltaTPicoSil[j] != -99.) {
        inc += 1;
        temp_psdeltat += DeltaTPicoSil[j];
      }
    }
    temp_psdeltat /= inc;
    DeltaT_PicoSilEqual_MCP_Equal += 0.5*temp_psdeltat;

    float DeltaTPicoSilAt0EqualSmear = 0;
    inc = 0;
    for (int j = 0; j <= 6; j++) {
      if (DeltaTPicoSilSmear[j] != -99.) {
        inc += 1;
        DeltaTPicoSilAt0EqualSmear += DeltaTPicoSilSmear[j];
      }
    }
    DeltaTPicoSilAt0EqualSmear /= inc;


    float DeltaT_PicoSil_vs_MCP_TotalCharge = PicoSil_vs_MCP_WeightTotal / PicoSil_vs_MCP_ChargeTotal;
    float DeltaT_PicoSil_vs_MCP_EventCharge = PicoSil_vs_MCP_WeightEvent / PicoSil_vs_MCP_ChargeEvent;

    // For adding in each pixel in descending order by number of events.
    for(int j = 0; j < 7; j++) {
      double fill = 0;
      inc = 0;
      for(int k = 0; k <= j; k++) {
        if( DeltaTPicoSilSmear[ DeltaTPicoSilSmear_Events_SortedIndices[k] ] != -99. ) {
          fill += DeltaTPicoSilSmear[ DeltaTPicoSilSmear_Events_SortedIndices[k] ];
          inc += 1;
        }
      }
      fill /= inc;
      pixelsUsed[j][inc - 1] += 1;
      histDeltaTPicoSilAt0EqualSmear_nEventsCombine[j]->Fill( fill );
    }


    histDeltaT_Center_MCP_Equal->Fill(DeltaT_Center_MCP_Equal);
    histDeltaT_PicoSil_MCP_EventCharge->Fill(DeltaT_PicoSil_MCP_EventCharge); 
    histDeltaT_PicoSil_MCP_TotalCharge->Fill(DeltaT_PicoSil_MCP_TotalCharge); 
    histDeltaT_PicoSilTotalCharge_MCP_Equal->Fill(DeltaT_PicoSilTotalCharge_MCP_Equal);
    histDeltaT_PicoSilEventCharge_MCP_Equal->Fill(DeltaT_PicoSilEventCharge_MCP_Equal);
    histDeltaT_PicoSilLandauCharge_MCP_Equal->Fill(DeltaT_PicoSilLandauCharge_MCP_Equal);
    histDeltaT_PicoSilLandauCharge_MCP_Equal_PicoSilSmear->Fill(DeltaT_PicoSilLandauCharge_MCP_Equal_PicoSilSmear);
    histDeltaT_PicoSilLandauCharge_MCP_Equal_BothSmear->Fill(DeltaT_PicoSilLandauCharge_MCP_Equal_BothSmear);
    histDeltaT_PicoSil_MCP_Equal->Fill(DeltaT_PicoSil_MCP_Equal);
    histDeltaT_PicoSilEqual_MCP_Equal->Fill(DeltaT_PicoSilEqual_MCP_Equal);
    histDeltaT_PicoSilEqual_MCP_Equal_BothSmear->Fill(0.5*DeltaTPicoSilAt0EqualSmear + 0.5*DeltaTMCPSmear);

    histDeltaTPicoSilAt0TotalCharge->Fill(temp_pstotalweight/temp_pstotalcharge);
    histDeltaTPicoSilAt0EventCharge->Fill(temp_pseventweight/temp_pseventcharge);
    histDeltaTPicoSilAt0LandauCharge->Fill(temp_pslandauweight/temp_pslandaucharge);
    histDeltaTPicoSilAt0LandauChargeSmear->Fill(temp_pslandauweight_smear/temp_pslandaucharge_smear);
    histDeltaTPicoSilAt0EqualSmear->Fill(DeltaTPicoSilAt0EqualSmear);
    histDeltaTCenterAt0->Fill(DeltaTCenter);
    histDeltaTMCPAt0->Fill(DeltaTMCP);
    histDeltaTMCPAt0Smear->Fill(DeltaTMCPSmear);

    histDeltaT_PicoSil_vs_MCP_TotalCharge->Fill(DeltaT_PicoSil_vs_MCP_TotalCharge);
    histDeltaT_PicoSil_vs_MCP_EventCharge->Fill(DeltaT_PicoSil_vs_MCP_EventCharge);
  }

  //for(int i=0;i<7;i++) { histAllCharges->Add(histCharges[i]); }

  // Add Gaussian fit
  Fitter(histDeltaT_Center_MCP_Equal);
  Fitter(histDeltaT_PicoSil_MCP_Equal);
  Fitter(histDeltaT_PicoSilEqual_MCP_Equal);
  Fitter(histDeltaT_PicoSilEqual_MCP_Equal_BothSmear);
  Fitter(histDeltaT_PicoSil_MCP_EventCharge);
  Fitter(histDeltaT_PicoSil_MCP_TotalCharge);
  Fitter(histDeltaT_PicoSilTotalCharge_MCP_Equal);
  Fitter(histDeltaT_PicoSilEventCharge_MCP_Equal);
  Fitter(histDeltaT_PicoSilLandauCharge_MCP_Equal);
  Fitter(histDeltaT_PicoSilLandauCharge_MCP_Equal_PicoSilSmear);
  Fitter(histDeltaT_PicoSilLandauCharge_MCP_Equal_BothSmear);
  Fitter(histDeltaTCenterAt0);
  Fitter(histDeltaTMCPAt0);
  Fitter(histDeltaTMCPAt0Smear);
  Fitter(histDeltaTPicoSilAt0TotalCharge);
  Fitter(histDeltaTPicoSilAt0EventCharge);
  Fitter(histDeltaTPicoSilAt0LandauCharge);
  Fitter(histDeltaTPicoSilAt0LandauChargeSmear);
  Fitter(histDeltaTPicoSilAt0EqualSmear);
  Fitter(histDeltaT_PicoSil_vs_MCP_TotalCharge);
  Fitter(histDeltaT_PicoSil_vs_MCP_EventCharge);
  //Fitter(histAllCharges);

  for (int i = 0; i < 7; i++) {
    Fitter(histDeltaTPicoSilSmearAt0[i]);
    Fitter(histDeltaTPicoSilAt0[i]);
    Fitter(histDeltaTPicoSilAt0EqualSmear_nEventsCombine[i]);
    Fitter(histDeltaT_PicoSil_vs_MCP[i]);
  }


  TCanvas *c = new TCanvas ("c","c",800, 600); 
  TLatex *tex = new TLatex();
  tex->SetNDC(); // Sets coords such that (0,0) is bottom left & (1,1) is top right.
  tex->SetTextSize(0.060);
  tex->SetTextFont(42);
  tex->SetTextColor(kBlack);

  c->cd();

  string pixels_added[7];
  pixels_added[0] = Form("#Deltat of Equal-Weight HGC Pixels in Order of # Events. Smeared Pixels: %d", DeltaTPicoSilSmear_Events_SortedIndices[0]);
  for (int i = 1; i<7; i++) pixels_added[i] = pixels_added[i-1] + Form(",%d",DeltaTPicoSilSmear_Events_SortedIndices[i]);
  for (int i = 0; i<7; i++){
    //histDeltaTPicoSilAt0EqualSmear_nEventsCombine[i]->SetTitle( pixels_added[i].c_str() ); <-- COMMENT OUT TITLE
    int tempArray[7] = {0};
    for (int j = 0; j<7; j++) tempArray[j] = pixelsUsed[i][j];
    SKIROCPlotPDF(c, tex, histDeltaTPicoSilAt0EqualSmear_nEventsCombine[i], Form("SKIROC_%d_Pixels",i+1), tempArray );
  }
  c->Close();

  // Creates output root file
  TFile *file = TFile::Open(("output"+filename).c_str(), "RECREATE");
  file->cd();
  file->WriteTObject(histDeltaTCenterAt0,"histDeltaTCenter", "WriteDelete");
  file->WriteTObject(histDeltaTPicoSilAt0EventCharge,"histDeltaTPicoSilEventCharge", "WriteDelete");
  file->WriteTObject(histDeltaTPicoSilAt0TotalCharge,"histDeltaTPicoSilTotalCharge", "WriteDelete");
  file->WriteTObject(histDeltaTPicoSilAt0LandauCharge,"histDeltaTPicoSilLandauCharge", "WriteDelete");
  file->WriteTObject(histDeltaTPicoSilAt0LandauChargeSmear,"histDeltaTPicoSilLandauChargeSmear", "WriteDelete");
  file->WriteTObject(histDeltaTPicoSilAt0EqualSmear,"histDeltaTPicoSilEqualSmear", "WriteDelete");
  file->WriteTObject(histDeltaTMCPAt0,"histDeltaTMCP", "WriteDelete");
  file->WriteTObject(histDeltaTMCPAt0Smear,"histDeltaTMCPSmear", "WriteDelete");
  file->WriteTObject(histDeltaT_Center_MCP_EventCharge_NoShift,"histDeltaT_Center_MCP_EventCharge_NoShift", "WriteDelete");
  file->WriteTObject(histDeltaT_Center_MCP_Equal,"histDeltaT_Center_MCP_Equal", "WriteDelete");
  file->WriteTObject(histDeltaT_PicoSil_MCP_Equal,"histDeltaT_PicoSil_MCP_Equal", "WriteDelete");
  file->WriteTObject(histDeltaT_PicoSilEqual_MCP_Equal,"histDeltaT_PicoSilEqual_MCP_Equal", "WriteDelete");
  file->WriteTObject(histDeltaT_PicoSilEqual_MCP_Equal_BothSmear,"histDeltaT_PicoSilEqual_MCP_Equal_BothSmear", "WriteDelete");
  file->WriteTObject(histDeltaT_PicoSilEventCharge_MCP_Equal,"histDeltaT_PicoSilEventCharge_MCP_Equal", "WriteDelete");
  file->WriteTObject(histDeltaT_PicoSilTotalCharge_MCP_Equal,"histDeltaT_PicoSilTotalCharge_MCP_Equal", "WriteDelete");
  file->WriteTObject(histDeltaT_PicoSilLandauCharge_MCP_Equal,"histDeltaT_PicoSilLandauCharge_MCP_Equal", "WriteDelete");
  file->WriteTObject(histDeltaT_PicoSilLandauCharge_MCP_Equal_PicoSilSmear,"histDeltaT_PicoSilLandauCharge_MCP_Equal_PicoSilSmear", "WriteDelete");
  file->WriteTObject(histDeltaT_PicoSilLandauCharge_MCP_Equal_BothSmear,"histDeltaT_PicoSilLandauCharge_MCP_Equal_BothSmear", "WriteDelete");
  file->WriteTObject(histDeltaT_PicoSil_MCP_EventCharge,"histDeltaT_PicoSil_MCP_EventCharge", "WriteDelete");
  file->WriteTObject(histDeltaT_PicoSil_MCP_TotalCharge,"histDeltaT_PicoSil_MCP_TotalCharge", "WriteDelete");
  file->WriteTObject(histDeltaT_PicoSil_vs_MCP_TotalCharge,"histDeltaT_PicoSil_vs_MCP_TotalCharge", "WriteDelete");
  file->WriteTObject(histDeltaT_PicoSil_vs_MCP_EventCharge,"histDeltaT_PicoSil_vs_MCP_EventCharge", "WriteDelete");
  //file->WriteTObject(histAllCharges,"histAllCharges","WriteDelete");
  for(int i=0; i<=6; i++) file->WriteTObject(histDeltaTPicoSilAt0[i], Form("histDeltaTPicoSil[%d]",i),"WriteDelete");
  for(int i=0; i<=6; i++) file->WriteTObject(histCharges[i],Form("histCharges[%d]",i),"WriteDelete");
  for(int i=0; i<=6; i++) file->WriteTObject(histDeltaTPicoSilSmearAt0[i],Form("histDeltaTPicoSilSmear[%d]",i),"WriteDelete");
  for(int i=0; i<=6; i++) file->WriteTObject(histDeltaTPicoSilAt0EqualSmear_nEventsCombine[i], pixels_added[i].c_str(), "WriteDelete");
  for(int i=0; i<=6; i++) file->WriteTObject(histDeltaT_PicoSil_vs_MCP[i],Form("histDeltaT_PicoSil_vs_MCP[%d]",i),"WriteDelete");
  // Above are in separate loops to be organized in the TBrowser


  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  TH1F *histPhotekAmpCut = new TH1F("histPhotekAmpCut",";Amplitude (mA);Entries/(0.04 mA)", 75, 0, 3);
  TH1F *histPhotekChargeCut = new TH1F("histPhotekChargeCut",";Charge (pC);Entries/(0.4 pC)", 75, 0, 30);
  TH1F *histCenterAmpCut = new TH1F("histCenterAmpCut",";Amplitude (mA);Entries/(0.015 mA)", 100, 0, 1.5);
  TH1F *histCenterChargeCut = new TH1F("histCenterChargeCut",";Charge (pC);Entries/(0.8 pC)", 75, 0, 60);
  TH1F *histMCPAmpCut = new TH1F("histMCPAmpCut",";Amplitude (mA);Entries/(0.008 mA)", 100, 0, 0.8);

  tree->Draw("sqrt(10)*amp[0]>>histPhotekAmpCut", Form("sqrt(10)*amp[0]>%f",photekAmpCut) );
  tree->Draw("sqrt(10)*int[0]>>histPhotekChargeCut", Form("sqrt(10)*int[0]>%f",photekChargeCut));
  tree->Draw("2*amp[1]>>histCenterAmpCut", Form("2*amp[1]>%f",centerAmpCut));
  tree->Draw("2*int[1]>>histCenterChargeCut", Form("2*int[1]>%f",centerChargeCut));
  tree->Draw("amp[11]>>histMCPAmpCut", Form("amp[11]>%f",MCPAmpCut));

  TH1F *histPhotekAmp = new TH1F("histPhotekAmp",";Amplitude (mA);Entries/(0.04 mA)", 75, 0, 3);
  TH1F *histPhotekCharge = new TH1F("histPhotekCharge",";Charge (pC);Entries/(0.4 pC)", 75, 0, 30);
  TH1F *histCenterAmp = new TH1F("histCenterAmp",";Amplitude (mA);Entries/(0.015 mA)", 100, 0, 1.5);
  TH1F *histCenterCharge = new TH1F("histCenterCharge",";Charge (pC);Entries/(0.8 pC)", 75, 0, 60);
  TH1F *histMCPAmp = new TH1F("histMCPAmp",";Amplitude (mA);Entries/(0.008 mA)", 100, 0, 0.8);

  tree->Draw("sqrt(10)*amp[0]>>histPhotekAmp", "", " " );
  tree->Draw("sqrt(10)*int[0]>>histPhotekCharge");
  tree->Draw("2*amp[1]>>histCenterAmp");
  tree->Draw("2*int[1]>>histCenterCharge");
  tree->Draw("amp[11]>>histMCPAmp");

  file->WriteTObject(histPhotekAmp, "Photek Amp", "WriteDelete");
  file->WriteTObject(histPhotekAmpCut, "Cut on Photek Amp", "WriteDelete");
  file->WriteTObject(histPhotekCharge, "Photek Charge", "WriteDelete");
  file->WriteTObject(histPhotekChargeCut, "Cut on Photek Charge", "WriteDelete");
  file->WriteTObject(histCenterAmp, "Center Amp", "WriteDelete");
  file->WriteTObject(histCenterAmpCut, "Cut on Center Amp", "WriteDelete");
  file->WriteTObject(histCenterCharge, "Center Charge", "WriteDelete");
  file->WriteTObject(histCenterChargeCut, "Cut on Center Charge", "WriteDelete");
  file->WriteTObject(histMCPAmp, "MCP Amp", "WriteDelete");
  file->WriteTObject(histMCPAmpCut, "Cut on MCP Amp", "WriteDelete");

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
  // Just use original X axis titles:
  //hist->GetXaxis()->SetTitle("#Deltat (ns)");
  if(1000*gausfit->GetParError(2)>2) tex->DrawLatex(0.59, 0.83, Form("#sigma = %.0f #pm %.0f ps", 1000*gausfit->GetParameter(2), 1000*gausfit->GetParError(2)));
  else tex->DrawLatex(0.59, 0.83, Form("#sigma = %.1f #pm %.1f ps", 1000*gausfit->GetParameter(2), 1000*gausfit->GetParError(2)));
  c->SaveAs(outfile.c_str()); //outfile should end in .pdf
}


void makeTimeResolution( string filename, float photekAmpCut, float photekChargeCut, float centerAmpCut, float centerChargeCut, float MCPAmpCut ) {

  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);

  DoMultiDeviceStudy( filename.c_str(), photekAmpCut, photekChargeCut, centerAmpCut, centerChargeCut, MCPAmpCut );

  TFile *_file = TFile::Open( ("output"+filename).c_str() ); //Should be .root

  //Create variables containing hists:
  TH1F *histDeltaTCenter = (TH1F*)_file->Get("histDeltaTCenter"); //Picosil center pixel
  TH1F *histDeltaTPicoSilTotalCharge = (TH1F*)_file->Get("histDeltaTPicoSilTotalCharge");//Picosil, total charge
  TH1F *histDeltaTPicoSilEventCharge = (TH1F*)_file->Get("histDeltaTPicoSilEventCharge");//Picosil, event charge
  TH1F *histDeltaTPicoSilLandauCharge= (TH1F*)_file->Get("histDeltaTPicoSilLandauCharge");//Picosil, landau charge
  TH1F *histDeltaTPicoSilLandauChargeSmear = (TH1F*)_file->Get("histDeltaTPicoSilLandauChargeSmear");//Picosil, landau charge
  TH1F *histDeltaTPicoSilEqualSmear = (TH1F*)_file->Get("histDeltaTPicoSilEqualSmear");//Picosil, landau charge
  TH1F *histDeltaTMCP = (TH1F*)_file->Get("histDeltaTMCP"); // MCP
  TH1F *histDeltaTMCPSmear = (TH1F*)_file->Get("histDeltaTMCPSmear"); // MCP
  TH1F *histDeltaT_Center_MCP_Equal = (TH1F*)_file->Get("histDeltaT_Center_MCP_Equal"); // each device weighted equally
  TH1F *histDeltaT_PicoSil_MCP_Equal = (TH1F*)_file->Get("histDeltaT_PicoSil_MCP_Equal");
  TH1F *histDeltaT_PicoSilEqual_MCP_Equal = (TH1F*)_file->Get("histDeltaT_PicoSilEqual_MCP_Equal");
  TH1F *histDeltaT_PicoSilEqual_MCP_Equal_BothSmear = (TH1F*)_file->Get("histDeltaT_PicoSilEqual_MCP_Equal_BothSmear");
  TH1F *histDeltaT_PicoSilEventCharge_MCP_Equal = (TH1F*)_file->Get("histDeltaT_PicoSilEventCharge_MCP_Equal");
  TH1F *histDeltaT_PicoSilTotalCharge_MCP_Equal = (TH1F*)_file->Get("histDeltaT_PicoSilTotalCharge_MCP_Equal");
  TH1F *histDeltaT_PicoSilLandauCharge_MCP_Equal= (TH1F*)_file->Get("histDeltaT_PicoSilLandauCharge_MCP_Equal");
  TH1F *histDeltaT_PicoSilLandauCharge_MCP_Equal_PicoSilSmear = (TH1F*)_file->Get("histDeltaT_PicoSilLandauCharge_MCP_Equal_PicoSilSmear");
  TH1F *histDeltaT_PicoSilLandauCharge_MCP_Equal_BothSmear = (TH1F*)_file->Get("histDeltaT_PicoSilLandauCharge_MCP_Equal_BothSmear");
  TH1F *histDeltaT_PicoSil_MCP_EventCharge = (TH1F*)_file->Get("histDeltaT_PicoSil_MCP_EventCharge"); //Combination of device Delta T's after shifting distributions around 0 and then weighting event by event.
  TH1F *histDeltaT_PicoSil_MCP_TotalCharge = (TH1F*)_file->Get("histDeltaT_PicoSil_MCP_TotalCharge"); //Combination after shifting around 0 and weighting with total charge.
  TH1F *histDeltaT_PicoSil_vs_MCP_TotalCharge = (TH1F*)_file->Get("histDeltaT_PicoSil_vs_MCP_TotalCharge");
  TH1F *histDeltaT_PicoSil_vs_MCP_EventCharge = (TH1F*)_file->Get("histDeltaT_PicoSil_vs_MCP_EventCharge");
  TH1F *histPhotekAmp = (TH1F*)_file->Get("Photek Amp");
  TH1F *histPhotekAmpCut = (TH1F*)_file->Get("Cut on Photek Amp");
  TH1F *histPhotekCharge = (TH1F*)_file->Get("Photek Charge");
  TH1F *histPhotekChargeCut = (TH1F*)_file->Get("Cut on Photek Charge");
  TH1F *histDeltaTPicoSil[6];
  for(int i=1; i<=6; i++) histDeltaTPicoSil[i-1] = (TH1F*)_file->Get( Form("histDeltaTPicoSil[%d]",i) ); //Already wrote center pixel


  TCanvas *c = new TCanvas ("c","c",800, 600); 
  TLatex *tex = new TLatex();
  tex->SetNDC(); // Sets coords such that (0,0) is bottom left & (1,1) is top right.
  tex->SetTextSize(0.060);
  tex->SetTextFont(42);
  tex->SetTextColor(kBlack);

  c->cd();

  /* REMOVING TITLES:
  histDeltaTCenter->SetTitle("HGC Center Pixel: TOF");
  histDeltaTPicoSilEventCharge->SetTitle("HGC: TOF w/ Event Charge Weighting");
  histDeltaTPicoSilTotalCharge->SetTitle("HGC: TOF w/ Total Charge Weighting");
  histDeltaTPicoSilLandauCharge->SetTitle("HGC: TOF w/ Landau MPV Charge Weighting");
  histDeltaTPicoSilLandauChargeSmear->SetTitle("#splitline{SKIROC Emulation: HGC TOF w/ Landau}{MPV Charge Weighting}");
  histDeltaTPicoSilEqualSmear->SetTitle("SKIROC Emulation: HGC TOF w/ Equal Weighting");
  histDeltaTMCP->SetTitle("MCP: TOF");
  histDeltaTMCPSmear->SetTitle("SKIROC Emulation: MCP TOF");
  histDeltaT_Center_MCP_Equal->SetTitle("1/2 HGC Center Pixel, 1/2 MCP: TOF");
  histDeltaT_PicoSil_MCP_Equal->SetTitle("HGC 7 Pixels, MCP: TOF w/ Weighting 1/8");
  histDeltaT_PicoSilEqual_MCP_Equal->SetTitle("1/14 HGC 7 Pixels, 1/2 MCP: TOF");
  histDeltaT_PicoSilEqual_MCP_Equal_BothSmear->SetTitle("1/14 HGC Smeared 7 Pixels, 1/2 MCP Smeared: TOF");
  histDeltaT_PicoSil_MCP_EventCharge->SetTitle("HGC, MCP: TOF w/ Event Charge Weighting");
  histDeltaT_PicoSil_MCP_TotalCharge->SetTitle("HGC, MCP: TOF w/ Total Charge Weighting");
  histDeltaT_PicoSilEventCharge_MCP_Equal->SetTitle("1/2 HGC w/ Event Charge Weighting, 1/2 MCP: TOF");
  histDeltaT_PicoSilTotalCharge_MCP_Equal->SetTitle("1/2 HGC w/ Total Charge Weighting, 1/2 MCP: TOF");
  histDeltaT_PicoSilLandauCharge_MCP_Equal->SetTitle("1/2 HGC w/ Landau MPV Charge Weighting, 1/2 MCP: TOF");
  histDeltaT_PicoSilLandauCharge_MCP_Equal_PicoSilSmear->SetTitle("#splitline{SKIROC Emulation: 1/2 Smeared HGC w/ Landau}{MPV Charge Weighting, 1/2 MCP: TOF}");
  histDeltaT_PicoSilLandauCharge_MCP_Equal_BothSmear->SetTitle("#splitline{SKIROC Emulation: 1/2 Smeared HGC w/ Landau}{MPV Charge Weighting, 1/2 Smeared MCP: TOF}");
  histDeltaT_PicoSil_vs_MCP_TotalCharge->SetTitle("#Deltat b/t HGC and Photonis -- Total Charge Weighted");
  histDeltaT_PicoSil_vs_MCP_EventCharge->SetTitle("#Deltat b/t HGC and Photonis -- Event Charge Weighted");
  for(int i=0; i<6; i++) histDeltaTPicoSil[i]->SetTitle( Form("HGC Pixel %d: TOF",i+1) ); //pixel 0 is center pixel
  */


  PlotDeltaTPDF(c, tex, histDeltaTCenter, "deltaTCenter.pdf");
  PlotDeltaTPDF(c, tex, histDeltaTPicoSilEventCharge, "deltaTPicoSilEventCharge.pdf");
  PlotDeltaTPDF(c, tex, histDeltaTPicoSilTotalCharge, "deltaTPicoSilTotalCharge.pdf");
  PlotDeltaTPDF(c, tex, histDeltaTPicoSilLandauCharge, "deltaTPicoSilLandauCharge.pdf");
  PlotDeltaTPDF(c, tex, histDeltaTPicoSilLandauChargeSmear, "deltaTPicoSilLandauChargeSmear.pdf");
  PlotDeltaTPDF(c, tex, histDeltaTPicoSilEqualSmear, "deltaTPicoSilEqualSmear.pdf");
  PlotDeltaTPDF(c, tex, histDeltaTMCP, "deltaTMCP.pdf");
  PlotDeltaTPDF(c, tex, histDeltaTMCPSmear, "deltaTMCPSmear.pdf");
  PlotDeltaTPDF(c, tex, histDeltaT_Center_MCP_Equal, "deltaT_Center_MCP_Equal.pdf");
  PlotDeltaTPDF(c, tex, histDeltaT_PicoSil_MCP_Equal, "deltaT_PicoSil_MCP_Equal.pdf");
  PlotDeltaTPDF(c, tex, histDeltaT_PicoSilEqual_MCP_Equal, "deltaT_PicoSilEqual_MCP_Equal.pdf");
  PlotDeltaTPDF(c, tex, histDeltaT_PicoSilEqual_MCP_Equal_BothSmear, "deltaT_PicoSilEqual_MCP_Equal_BothSmear.pdf");
  PlotDeltaTPDF(c, tex, histDeltaT_PicoSil_MCP_EventCharge, "deltaT_PicoSil_MCP_EventCharge.pdf");
  PlotDeltaTPDF(c, tex, histDeltaT_PicoSil_MCP_TotalCharge, "deltaT_PicoSil_MCP_TotalCharge.pdf");
  PlotDeltaTPDF(c, tex, histDeltaT_PicoSilEventCharge_MCP_Equal, "deltaT_PicoSilEventCharge_MCP_Equal.pdf");
  PlotDeltaTPDF(c, tex, histDeltaT_PicoSilTotalCharge_MCP_Equal, "deltaT_PicoSilTotalCharge_MCP_Equal.pdf");
  PlotDeltaTPDF(c, tex, histDeltaT_PicoSilLandauCharge_MCP_Equal, "deltaT_PicoSilLandauCharge_MCP_Equal.pdf");
  PlotDeltaTPDF(c, tex, histDeltaT_PicoSilLandauCharge_MCP_Equal_PicoSilSmear, "deltaT_PicoSilLandauCharge_MCP_Equal_PicoSilSmear.pdf");
  PlotDeltaTPDF(c, tex, histDeltaT_PicoSilLandauCharge_MCP_Equal_BothSmear, "deltaT_PicoSilLandauCharge_MCP_Equal_BothSmear.pdf");
  PlotDeltaTPDF(c, tex, histDeltaT_PicoSil_vs_MCP_TotalCharge, "deltaT_PicoSil_vs_MCP_TotalCharge.pdf");
  PlotDeltaTPDF(c, tex, histDeltaT_PicoSil_vs_MCP_EventCharge, "deltaT_PicoSil_vs_MCP_EventCharge.pdf");
  for(int i=0; i<6; i++) PlotDeltaTPDF(c, tex, histDeltaTPicoSil[i], Form("deltaTPicoSilPixel%d.pdf",i+1) );

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


void MultiDeviceStudy_PicosilMCP() {
  gStyle->SetTitleOffset(0.8,"x");
  gStyle->SetTitleOffset(0.85,"y");
  gStyle->SetTitleSize(0.055,"x");
  gStyle->SetTitleSize(0.055,"y");
  gStyle->SetLabelSize(0.045,"x");
  gStyle->SetLabelSize(0.045,"y");

  string infile = "104-116except111-114.root";
  float photekAmpCut = sqrt(10)*0.1; //THESE ARE THE CUT VALUES AFTER ADJUSTING FOR ATTENUATORS
  float photekChargeCut = sqrt(10)*2;
  float centerAmpCut = 2*0.15;
  float centerChargeCut = 2*11;
  float MCPAmpCut = 0.08;
  makeTimeResolution(infile.c_str(), photekAmpCut, photekChargeCut, centerAmpCut, centerChargeCut, MCPAmpCut); // Outputs PDFs with histograms
  //Un-comment following lines to make all output files at once:
  /*cout<<"\n\n 65-83:"<<endl;
  makeTimeResolution("65-83.root",                sqrt(10)*0.1,   sqrt(10)*2,    2*0.15, 2*10,  0.05);

  cout<<"\n\n 84-93:"<<endl;
  makeTimeResolution("84-93.root",                sqrt(10)*0.1,   sqrt(10)*2,    2*0.03, 2*2.5, 0.05);

  cout<<"\n\n 94-103:"<<endl;
  makeTimeResolution("94-103.root",               sqrt(10)*0.1,   sqrt(10)*2.5,  2*0.1,  2*7,   0.05);

  cout<<"\n\n 104-110,115-116:"<<endl;
  makeTimeResolution("104-116except111-114.root", sqrt(10)*0.1,   sqrt(10)*2,    2*0.15, 2*11,  0.08);

  cout<<"\n\n 117-122:"<<endl;
  makeTimeResolution("117-122.root",              sqrt(10)*0.09,  sqrt(10)*2,    2*0.05, 2*3,   0.055);

  cout<<"\n\n 129-138:"<<endl;
  makeTimeResolution("129-138.root",              sqrt(10)*0.1,   sqrt(10)*2,    2*0.1,  2*8,   0.075);

  cout<<"\n\n 144-155:"<<endl;
  makeTimeResolution("144-155.root",              sqrt(10)*0.03,  sqrt(10)*0.8,  2*0.07, 2*6,   0.025);

  cout<<"\n\n 167-171:"<<endl;
  makeTimeResolution("167-171.root",              sqrt(10)*0.015, sqrt(10)*0.4,  2*0.01, 2*2.5, 0.01);

  cout<<"\n\n 178-185:"<<endl;
  makeTimeResolution("178-185.root",              sqrt(10)*0.015, sqrt(10)*0.3,  2*0.01, 2*1,   0.03);

  cout<<"\n\n 186-200:"<<endl;
  makeTimeResolution("186-200.root",              sqrt(10)*0.03,  sqrt(10)*0.75, 2*0.02, 2*1,   0.05);*/

}
