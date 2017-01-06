
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

// This code generates the combined Delta T histogram for the picosil center pixel, the SiPad, and the MCP for the re-cabled configuration. 
// In order to use in early runs (prior to recabling), the array numbers and reference photek signals must be changed.

// Author: Daniel Gawerc


void Fitter(TH1F *hist) {
  //Helper function for fitting Gaussian
  double xmin = hist->GetMean() - 2.0*hist->GetRMS();
  double xmax = hist->GetMean() + 2.0*hist->GetRMS();
  hist->Fit("gaus","QMLES","",xmin,xmax); // Q suppresses fit results
  gStyle->SetOptFit(1);
  gStyle->SetOptStat(0);
}


void DoMultiDeviceStudy( string filename, float photekAmpCut, float photekChargeCut, float centerAmpCut, float centerChargeCut, float MCPAmpCut, float SiPadAmpCut ) {


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
  float SiPadWidth = 1.5;
  int bins = 75;

  TH1F *histDeltaTCenter = new TH1F("histDeltaTCenter",";#Deltat (ns);Entries/(0.012 ns)", 500, 0, 6); //DeltaT of center HGC pixel
  TH1F *histDeltaTCenterAt0 = new TH1F("histDeltaTCenterAt0",";#Deltat (ns);Entries/(0.004 ns)", bins, -width, width); 
  TH1F *histDeltaTMCP = new TH1F("histDeltaTMCP",";#Deltat (ns);Entries/(0.01 ns)", 500, 0, 5); //DeltaT of MCP
  TH1F *histDeltaTMCPAt0 = new TH1F("histDeltaTMCPAt0",";#Deltat (ns);Entries/(0.004 ns)", bins, -width, width); 
  TH1F *histDeltaTSiPad = new TH1F("histDeltaTSiPad",";#Deltat (ns);Entries/(0.06 ns)", 500, -15, 15); //DeltaT of SiPad
  TH1F *histDeltaTSiPadAt0 = new TH1F("histDeltaTSiPadAt0",";#Deltat (ns);Entries/(0.04 ns)", bins, -SiPadWidth, SiPadWidth);
  TH1F *histDeltaT_Center_MCP_SiPad_Equal = new TH1F("histDeltaT_Center_MCP_SiPad_Equal",";#Deltat (ns);Entries/(0.012 ns)", bins, -3*width, 3*width); // Weights HGC center pixel, MCP, and SiPad equally
  TH1F *histDeltaT_Center_MCP_SiPad_EventCharge = new TH1F("histDeltaT_Center_MCP_SiPad_EventCharge",";#Deltat (ns);Entries/(0.012 ns)", bins, -3*width, 3*width); //Event charge weighting
TH1F *histDeltaT_Center_MCP_SiPad_TotalCharge = new TH1F("histDeltaT_Center_MCP_SiPad_TotalCharge",";#Deltat (ns);Entries/(0.012 ns)", bins, -3*width, 3*width); //Total Charge Weighting

  

  float totalCenterCharge = 0;
  float totalSiPadCharge = 0;
  float totalMCPCharge = 0;

  //read all entries and fill the histogram
  Long64_t nentries = tree->GetEntries();

  //Loop through every event in .root file
  std::cout<<"Number of events in sample: "<<nentries<<std::endl;  
  for (Long64_t iEntry=0; iEntry<nentries; iEntry++) {
    if (iEntry %1000 == 0) cout << "Processing Event " << iEntry << "\n";
    tree->GetEntry(iEntry);    
  
    float photekTimeGauss0 = gauspeak[0];
    float photekAmp0 = sqrt(10)*amp[0]; //Acct for attenuator
    float photekCharge0 = sqrt(10)*integral[0];

    float photekTimeGauss1 = gauspeak[9];
    float photekAmp1 = sqrt(10)*amp[9];
    float photekCharge1 = sqrt(10)*integral[9];

    float centerAmp = 2*amp[1];
    float centerCharge = 2*integral[1];
    float centerTime = linearTime45[1];

    float SiPadAmp = sqrt(10)*amp[10];
    float SiPadCharge = sqrt(10)*integral[10];
    float SiPadTime = linearTime45[10];

    float MCPAmp = amp[11];
    float MCPCharge = integral[11];
    float MCPTime = linearTime45[11];


    // APPLY EVENT CUTS:
    //require photek (for electron selection)
    if( !(photekAmp0 > photekAmpCut && photekCharge0 > photekChargeCut) ) continue;
    
    //require signal in the central pixel
    if( !(centerAmp > centerAmpCut && centerCharge > centerChargeCut) ) continue;

    //require Si Pad and MCP minimum amplitude. Possibly will include cut for charge at some point...
    if( !(SiPadAmp > MCPAmpCut && SiPadAmp > SiPadAmpCut) ) continue;

    //Calculates the Delta T's if the event passes the cuts:
    float DeltaTCenter = photekTimeGauss0 - centerTime; // Will be used to fill the histogram at every event
    float DeltaTSiPad = photekTimeGauss1 - SiPadTime;
    float DeltaTMCP = photekTimeGauss1 - MCPTime;


    totalCenterCharge += centerCharge;
    totalSiPadCharge += SiPadCharge;
    totalMCPCharge += MCPCharge;

    histDeltaTCenter->Fill(DeltaTCenter);
    histDeltaTSiPad->Fill(DeltaTSiPad);
    histDeltaTMCP->Fill(DeltaTMCP);
    
  }


  double meanCenter = histDeltaTCenter->GetMean();
  double meanMCP = histDeltaTMCP->GetMean();
  double meanSiPad = histDeltaTSiPad->GetMean();



  //Loop through again. This time, subtracting the respective means. Process is almost the same.
  for (Long64_t iEntry=0; iEntry<nentries; iEntry++) {
    if (iEntry %1000 == 0) cout << "Processing Event " << iEntry << "\n";
    tree->GetEntry(iEntry);    
  
    float photekTimeGauss0 = gauspeak[0];
    float photekAmp0 = sqrt(10)*amp[0]; //Acct for attenuator
    float photekCharge0 = sqrt(10)*integral[0];

    float photekTimeGauss1 = gauspeak[9];
    float photekAmp1 = sqrt(10)*amp[9];
    float photekCharge1 = sqrt(10)*integral[9];

    float centerAmp = 2*amp[1];
    float centerCharge = 2*integral[1];
    float centerTime = linearTime45[1];

    float SiPadAmp = sqrt(10)*amp[10];
    float SiPadCharge = sqrt(10)*integral[10];
    float SiPadTime = linearTime45[10];

    float MCPAmp = amp[11];
    float MCPCharge = integral[11];
    float MCPTime = linearTime45[11];

    if( !(photekAmp0 > photekAmpCut && photekCharge0 > photekChargeCut) ) continue;
    if( !(centerAmp > centerAmpCut && centerCharge > centerChargeCut) ) continue;
    if( !(SiPadAmp > MCPAmpCut && SiPadAmp > SiPadAmpCut) ) continue;

    //shift histograms to 0.
    float DeltaTCenter = photekTimeGauss0 - centerTime - meanCenter; 
    float DeltaTSiPad = photekTimeGauss1 - SiPadTime - meanSiPad;
    float DeltaTMCP = photekTimeGauss1 - MCPTime - meanMCP;
    float DeltaT_Center_MCP_SiPad_Equal = (DeltaTCenter + DeltaTMCP + DeltaTSiPad)/3;
    float DeltaT_Center_MCP_SiPad_EventCharge = (DeltaTCenter*centerCharge + DeltaTSiPad*SiPadCharge + DeltaTMCP*MCPCharge) / (centerCharge+SiPadCharge+MCPCharge);
    float DeltaT_Center_MCP_SiPad_TotalCharge = (DeltaTCenter*totalCenterCharge + DeltaTSiPad*totalSiPadCharge + DeltaTMCP*totalMCPCharge) / (totalCenterCharge+totalSiPadCharge+totalMCPCharge);


    histDeltaTCenterAt0->Fill(DeltaTCenter);
    histDeltaTSiPadAt0->Fill(DeltaTSiPad);
    histDeltaTMCPAt0->Fill(DeltaTMCP);
    histDeltaT_Center_MCP_SiPad_Equal->Fill(DeltaT_Center_MCP_SiPad_Equal);
    histDeltaT_Center_MCP_SiPad_EventCharge->Fill(DeltaT_Center_MCP_SiPad_EventCharge);
    histDeltaT_Center_MCP_SiPad_TotalCharge->Fill(DeltaT_Center_MCP_SiPad_TotalCharge);
  }

  // Do Gaussian fit of delta T distributions from (mean-2RMS) to (mean+2RMS)
  TCanvas *c1 = new TCanvas ("c1","c1",800, 600);
  Fitter(histDeltaTCenterAt0);
  Fitter(histDeltaTMCPAt0);
  Fitter(histDeltaTSiPadAt0);
  Fitter(histDeltaT_Center_MCP_SiPad_Equal);
  Fitter(histDeltaT_Center_MCP_SiPad_EventCharge);
  Fitter(histDeltaT_Center_MCP_SiPad_TotalCharge);



  // Creates output root file
  TFile *file = TFile::Open(("output"+filename).c_str(), "RECREATE");
  file->cd();
  file->WriteTObject(histDeltaTCenterAt0,"histDeltaTCenter", "WriteDelete");
  file->WriteTObject(histDeltaTMCPAt0,"histDeltaTMCP", "WriteDelete");
  file->WriteTObject(histDeltaTSiPadAt0,"histDeltaTSiPad", "WriteDelete");
  file->WriteTObject(histDeltaT_Center_MCP_SiPad_Equal, "histDeltaT_Center_MCP_SiPad_Equal", "WriteDelete");
  file->WriteTObject(histDeltaT_Center_MCP_SiPad_EventCharge, "histDeltaT_Center_MCP_SiPad_EventCharge", "WriteDelete");
  file->WriteTObject(histDeltaT_Center_MCP_SiPad_TotalCharge, "histDeltaT_Center_MCP_SiPad_TotalCharge", "WriteDelete");


  c1->Close();
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
  if(1000*gausfit->GetParError(2)>2) tex->DrawLatex(0.59, 0.83, Form("#sigma = %.0f #pm %.0f ps", 1000*gausfit->GetParameter(2), 1000*gausfit->GetParError(2)));
  else tex->DrawLatex(0.59, 0.83, Form("#sigma = %.1f #pm %.1f ps", 1000*gausfit->GetParameter(2), 1000*gausfit->GetParError(2)));
  c->SaveAs(outfile.c_str()); //outfile should end in .pdf
}






void makeTimeResolution( string filename, float photekAmpCut, float photekChargeCut, float centerAmpCut, float centerChargeCut, float MCPAmpCut, float SiPadAmpCut ) {

  DoMultiDeviceStudy( filename.c_str(), photekAmpCut, photekChargeCut, centerAmpCut, centerChargeCut, MCPAmpCut, SiPadAmpCut );

  TFile *_file = TFile::Open( ("output"+filename).c_str() ); //Should be .root file

  //Create variables containing hists:
  TH1F *histDeltaTCenter = (TH1F*)_file->Get("histDeltaTCenter"); //Picosil center pixel
  TH1F *histDeltaTMCP = (TH1F*)_file->Get("histDeltaTMCP"); // MCP
  TH1F *histDeltaTSiPad = (TH1F*)_file->Get("histDeltaTSiPad"); 
  TH1F *histDeltaT_Center_MCP_SiPad_Equal = (TH1F*)_file->Get("histDeltaT_Center_MCP_SiPad_Equal");
  TH1F *histDeltaT_Center_MCP_SiPad_EventCharge = (TH1F*)_file->Get("histDeltaT_Center_MCP_SiPad_EventCharge");
  TH1F *histDeltaT_Center_MCP_SiPad_TotalCharge = (TH1F*)_file->Get("histDeltaT_Center_MCP_SiPad_TotalCharge");


  TCanvas *c = new TCanvas ("c","c",800, 600); 
  TLatex *tex = new TLatex();
  tex->SetNDC(); // Sets coords such that (0,0) is bottom left & (1,1) is top right.
  tex->SetTextSize(0.060);
  tex->SetTextFont(42);
  tex->SetTextColor(kBlack);

  c->cd();

 
  PlotDeltaTPDF(c, tex, histDeltaTCenter, "deltaTCenter.pdf");
  PlotDeltaTPDF(c, tex, histDeltaTMCP, "deltaTMCP.pdf");
  PlotDeltaTPDF(c, tex, histDeltaTSiPad, "deltaTSiPad.pdf");
  PlotDeltaTPDF(c, tex, histDeltaT_Center_MCP_SiPad_Equal, "deltaT_Center_MCP_SiPad_Equal.pdf");
  PlotDeltaTPDF(c, tex, histDeltaT_Center_MCP_SiPad_EventCharge, "deltaT_Center_MCP_SiPad_EventCharge.pdf");
  PlotDeltaTPDF(c, tex, histDeltaT_Center_MCP_SiPad_TotalCharge, "deltaT_Center_MCP_SiPad_TotalCharge.pdf");
  
}


void MultiDeviceStudy() {
  gStyle->SetTitleOffset(0.8,"x");
  gStyle->SetTitleOffset(0.85,"y");
  gStyle->SetTitleSize(0.055,"x");
  gStyle->SetTitleSize(0.055,"y");
  gStyle->SetLabelSize(0.045,"x");
  gStyle->SetLabelSize(0.045,"y");


  string infile = "65-83.root";
  float photekAmpCut = sqrt(10)*0.1; 
  float photekChargeCut = sqrt(10)*2;
  float centerAmpCut = 2*0.15;
  float centerChargeCut = 2*10;
  float MCPAmpCut = 0.05;
  float SiPadAmpCut = 0.0001;
  makeTimeResolution(infile.c_str(), photekAmpCut, photekChargeCut, centerAmpCut, centerChargeCut, MCPAmpCut, SiPadAmpCut); // Outputs PDFs with histograms
}
