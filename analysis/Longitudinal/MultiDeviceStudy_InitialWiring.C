
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
  gStyle->SetOptFit(1);
}


void DoMultiDeviceStudy( string filename, float photekAmpCut, float photekChargeCut, float centerAmpCut, float centerChargeCut, float MCPAmpCut, float SiPadAmpCut) {

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
  float width = 1.5;
  float smearWidth = 0.75;
  int bins = 100;
  int smearBins = 75;
  float pixelSmear = 0.050; // in ns
  float MCPSmear = 0.045;
  float SiPadSmear = 0.050;

  TH1F *histDeltaTCenter = new TH1F("histDeltaTCenter","; Time [ns];Number of Events", 50, 4, 5); //DeltaT of center HGC pixel
  TH1F *histDeltaTCenterSmear = new TH1F("histDeltaTCenterSmear","; Time [ns];Number of Events", 50, 4, 5); 
  TH1F *histDeltaTCenterAt0 = new TH1F("histDeltaTCenterAt0","; Time [ns];Number of Events", bins, -width, width); 
  TH1F *histDeltaTCenterAt0Smear = new TH1F("histDeltaTCenterAt0Smear","; Time [ns];Number of Events", smearBins, -smearWidth, smearWidth);
  TH1F *histDeltaTMCP = new TH1F("histDeltaTMCP","; Time [ns];Number of Events", 50, 2, 3); //DeltaT of MCP
  TH1F *histDeltaTMCPSmear = new TH1F("histDeltaTMCPSmear","; Time [ns];Number of Events", 50, 2, 3); 
  TH1F *histDeltaTMCPAt0 = new TH1F("histDeltaTMCPAt0","; Time [ns];Number of Events", bins, -width, width); 
  TH1F *histDeltaTMCPAt0Smear = new TH1F("histDeltaTMCPAt0Smear","; Time [ns];Number of Events", smearBins, -smearWidth, smearWidth); 
  TH1F *histDeltaTSiPad = new TH1F("histDeltaTSiPad","; Time [ns];Number of Events", 500, -15, 15); //DeltaT of SiPad
  TH1F *histDeltaTSiPadSmear = new TH1F("histDeltaTSiPadSmear","; Time [ns];Number of Events", 500, -15, 15);
  TH1F *histDeltaTSiPadAt0 = new TH1F("histDeltaTSiPadAt0","; Time [ns];Number of Events", bins, -width, width);
  TH1F *histDeltaTSiPadAt0Smear = new TH1F("histDeltaTSiPadAt0Smear","; Time [ns];Number of Events", smearBins, -smearWidth, smearWidth);
  TH1F *histDeltaT_Center_MCP_SiPad_Equal = new TH1F("histDeltaT_Center_MCP_SiPad_Equal","; Time [ns];Number of Events", bins, -width, width); // Weights HGC center pixel, MCP, and SiPad equally
  TH1F *histDeltaT_Center_MCP_SiPad_EqualSmear = new TH1F("histDeltaT_Center_MCP_SiPad_EqualSmear","; Time [ns];Number of Events", bins, -width, width); // Weights HGC center pixel, MCP, and SiPad equally

  float totalCenterCharge = 0; 
  float totalMCPCharge = 0;
  float totalSiPadCharge = 0;


  //read all entries and fill the histogram
  Long64_t nentries = tree->GetEntries();

  //Loop through every event in .root file
  TRandom3 *rando = new TRandom3(seed);
  std::cout<<"Number of Events in Sample: "<<nentries<<std::endl;  
  for (Long64_t iEntry=0; iEntry<nentries; iEntry++) {
    if (iEntry %1000 == 0) cout << "Processing Event " << iEntry << "\n";
    tree->GetEntry(iEntry);    
  
    float photekTimeGauss0 = gauspeak[0];
    float photekAmp0 = amp[0]; //no attenuator
    float photekCharge0 = integral[0];

    float photekTimeGauss3 = gauspeak[27];
    float photekAmp3 = amp[27];
    float photekCharge3 = integral[27];

    float centerAmp = amp[34]; //no attenuator
    float centerCharge = integral[34];
    float centerTime = linearTime45[34];

    float SiPadAmp = amp[1]; //Runs 1-33 SiPad has no attenuator. After 33, it does.
    float SiPadCharge = integral[1];
    float SiPadTime = linearTime45[1];

    float MCPAmp = amp[2];
    float MCPCharge = integral[2];
    float MCPTime = linearTime45[2];


    // APPLY EVENT CUTS:
    //require photek (for electron selection)
    if( !(photekAmp0 > photekAmpCut && photekCharge0 > photekChargeCut) ) continue;
    //require signal in the central pixel
    if( !(centerAmp > centerAmpCut && centerCharge > centerChargeCut) ) continue;
    //require MCP minimum amplitude. Possibly will include cut for charge at some point...
    if( !( MCPAmp > MCPAmpCut && SiPadAmp > SiPadAmpCut) ) continue;

    
    double centerTimeSmear = rando->Gaus(centerTime, pixelSmear); //Samples from smear
    double MCPTimeSmear = rando->Gaus(MCPTime, MCPSmear);
    double SiPadTimeSmear = rando->Gaus(SiPadTime, SiPadSmear);


    //Calculates the Delta T's
    // (Will be used to fill the histogram at every event)
    float DeltaTCenter = photekTimeGauss3 - centerTime;
    float DeltaTCenterSmear = photekTimeGauss3 - centerTimeSmear;
    float DeltaTMCP = photekTimeGauss0 - MCPTime;
    float DeltaTMCPSmear = photekTimeGauss0 - MCPTimeSmear;
    float DeltaTSiPad = photekTimeGauss0 - SiPadTime;
    float DeltaTSiPadSmear = photekTimeGauss0 - SiPadTimeSmear;
    
    totalCenterCharge += centerCharge;
    totalMCPCharge += MCPCharge;
    totalSiPadCharge += SiPadCharge; // These aren't currently being used for any combinations

    histDeltaTCenter->Fill(DeltaTCenter);
    histDeltaTCenterSmear->Fill(DeltaTCenterSmear);
    histDeltaTMCP->Fill(DeltaTMCP);
    histDeltaTMCPSmear->Fill(DeltaTMCPSmear);
    histDeltaTSiPad->Fill(DeltaTSiPad);
    histDeltaTSiPadSmear->Fill(DeltaTSiPadSmear);
  }



  double meanCenter = histDeltaTCenter->GetMean();
  double meanCenterSmear = histDeltaTCenterSmear->GetMean();
  double meanMCP = histDeltaTMCP->GetMean();
  double meanMCPSmear = histDeltaTMCPSmear->GetMean();
  double meanSiPad = histDeltaTSiPad->GetMean();
  double meanSiPadSmear = histDeltaTSiPadSmear->GetMean();




  TRandom3 *rando2 = new TRandom3(seed);
  //Loop through again. This time, subtracting the respective means. Process is almost the same.
  for (Long64_t iEntry=0; iEntry<nentries; iEntry++) {
    if (iEntry %1000 == 0) cout << "Processing Event " << iEntry << "\n";
    tree->GetEntry(iEntry);    
  
    float photekTimeGauss0 = gauspeak[0];
    float photekAmp0 = amp[0]; // no attenuator
    float photekCharge0 = integral[0];

    float photekTimeGauss3 = gauspeak[27];
    float photekAmp3 = amp[27];
    float photekCharge3 = integral[27];

    float centerAmp = amp[34]; // no attenuator
    float centerCharge = integral[34];
    float centerTime = linearTime45[34];

    float SiPadAmp = amp[1]; //Runs 1-33 SiPad has no attenuator. After 33, it does.
    float SiPadCharge = integral[1];
    float SiPadTime = linearTime45[1];

    float MCPAmp = amp[2];
    float MCPCharge = integral[2];
    float MCPTime = linearTime45[2];

    if( !(photekAmp0 > photekAmpCut && photekCharge0 > photekChargeCut) ) continue;
    if( !(centerAmp > centerAmpCut && centerCharge > centerChargeCut) ) continue;
    if( !( MCPAmp > MCPAmpCut && SiPadAmp > SiPadAmpCut) ) continue;


    double centerTimeSmear = rando2->Gaus(centerTime, pixelSmear); //Samples from smear
    double MCPTimeSmear = rando2->Gaus(MCPTime, MCPSmear);
    double SiPadTimeSmear = rando2->Gaus(SiPadTime, SiPadSmear);

    // Shift histograms to 0
    float DeltaTCenter = photekTimeGauss3 - centerTime - meanCenter;
    float DeltaTCenterSmear = photekTimeGauss3 - centerTimeSmear - meanCenterSmear;
    float DeltaTMCP = photekTimeGauss0 - MCPTime - meanMCP;
    float DeltaTMCPSmear = photekTimeGauss0 - MCPTimeSmear - meanMCPSmear;
    float DeltaTSiPad = photekTimeGauss0 - SiPadTime - meanSiPad;
    float DeltaTSiPadSmear = photekTimeGauss0 - SiPadTimeSmear - meanSiPadSmear;

    float DeltaT_Center_MCP_SiPad_Equal = (DeltaTCenter + DeltaTMCP + DeltaTSiPad)/3;
    float DeltaT_Center_MCP_SiPad_EqualSmear = (DeltaTCenterSmear + DeltaTMCPSmear + DeltaTSiPadSmear)/3;


    histDeltaTCenterAt0->Fill(DeltaTCenter);
    histDeltaTCenterAt0Smear->Fill(DeltaTCenterSmear);
    histDeltaTMCPAt0->Fill(DeltaTMCP);
    histDeltaTMCPAt0Smear->Fill(DeltaTMCPSmear);
    histDeltaTSiPadAt0->Fill(DeltaTSiPad);
    histDeltaTSiPadAt0Smear->Fill(DeltaTSiPadSmear);
    histDeltaT_Center_MCP_SiPad_Equal->Fill(DeltaT_Center_MCP_SiPad_Equal);
    histDeltaT_Center_MCP_SiPad_EqualSmear->Fill(DeltaT_Center_MCP_SiPad_EqualSmear);
  }

  // Do Gaussian fit of delta T distributions from (mean-2RMS) to (mean+2RMS)
  Fitter(histDeltaTCenterAt0);
  Fitter(histDeltaTCenterAt0Smear);
  Fitter(histDeltaTMCPAt0);
  Fitter(histDeltaTMCPAt0Smear);
  Fitter(histDeltaTSiPadAt0);
  Fitter(histDeltaTSiPadAt0Smear);
  Fitter(histDeltaT_Center_MCP_SiPad_Equal);
  Fitter(histDeltaT_Center_MCP_SiPad_EqualSmear);


  // Creates output root file
  TFile *file = TFile::Open(("output"+filename).c_str(), "RECREATE");
  file->cd();
  file->WriteTObject(histDeltaTCenterAt0,"histDeltaTCenter", "WriteDelete");
  file->WriteTObject(histDeltaTCenterAt0Smear, "histDeltaTCenterSmear", "WriteDelete");
  file->WriteTObject(histDeltaTMCPAt0,"histDeltaTMCP", "WriteDelete");
  file->WriteTObject(histDeltaTMCPAt0Smear, "histDeltaTMCPSmear", "WriteDelete");
  file->WriteTObject(histDeltaTSiPad,"SiPad", "WriteDelete"); //delete after running -- use to find range
  file->WriteTObject(histDeltaTSiPadAt0,"histDeltaTSiPad", "WriteDelete");
  file->WriteTObject(histDeltaTSiPadAt0Smear, "histDeltaTSiPadSmear", "WriteDelete");
  file->WriteTObject(histDeltaT_Center_MCP_SiPad_Equal, "histDeltaT_Center_MCP_SiPad_Equal", "WriteDelete");
  file->WriteTObject(histDeltaT_Center_MCP_SiPad_EqualSmear, "histDeltaT_Center_MCP_SiPad_EqualSmear", "WriteDelete");


  TH1F *histPhotekAmpCut = new TH1F("histPhotekAmpCut","; Amp;Number of Events", 400, 0, 2.5);
  TH1F *histPhotekChargeCut = new TH1F("histPhotekChargeCut","; Charge;Number of Events", 400, 0, 30);
  TH1F *histCenterAmpCut = new TH1F("histCenterAmpCut","; Amp;Number of Events", 200, 0, 1.5);
  TH1F *histCenterChargeCut = new TH1F("histCenterChargeCut","; Charge;Number of Events", 400, 0, 60);
  TH1F *histMCPAmpCut = new TH1F("histMCPAmpCut","; Amp;Number of Events", 100, 0, 0.75);
  TH1F *histSiPadAmpCut = new TH1F("histSiPadAmpCut","; Amp;Number of Events", 100, 0, 0.75);

  tree->Draw("amp[0]>>histPhotekAmpCut", Form("amp[0]>%f",photekAmpCut) );
  tree->Draw("int[0]>>histPhotekChargeCut", Form("int[0]>%f",photekChargeCut));
  tree->Draw("amp[34]>>histCenterAmpCut", Form("amp[34]>%f",centerAmpCut));
  tree->Draw("int[34]>>histCenterChargeCut", Form("int[34]>%f",centerChargeCut));
  tree->Draw("amp[2]>>histMCPAmpCut", Form("amp[2]>%f",MCPAmpCut));
  tree->Draw("amp[1]>>histSiPadAmpCut", Form("amp[1]>%f",SiPadAmpCut) );

  TH1F *histPhotekAmp = new TH1F("histPhotekAmp","; Amp;Number of Events", 400, 0, 2.5);
  TH1F *histPhotekCharge = new TH1F("histPhotekCharge","; Charge;Number of Events", 400, 0, 30);
  TH1F *histCenterAmp = new TH1F("histCenterAmp","; Amp;Number of Events", 200, 0, 1.5);
  TH1F *histCenterCharge = new TH1F("histCenterCharge","; Charge;Number of Events", 400, 0, 60);
  TH1F *histMCPAmp = new TH1F("histMCPAmp","; Amp;Number of Events", 100, 0, 0.75);
  TH1F *histSiPadAmp = new TH1F("histSiPadAmp","; Amp;Number of Events", 100, 0, 0.75);

  tree->Draw("amp[0]>>histPhotekAmp");
  tree->Draw("int[0]>>histPhotekCharge");
  tree->Draw("amp[34]>>histCenterAmp");
  tree->Draw("int[34]>>histCenterCharge");
  tree->Draw("amp[2]>>histMCPAmp");
  tree->Draw("amp[1]>>histSiPadAmp");

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
  file->WriteTObject(histSiPadAmp, "SiPad Amp", "WriteDelete");
  file->WriteTObject(histSiPadAmpCut, "Cut on SiPad Amp", "WriteDelete");

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
  hist->GetXaxis()->SetTitle("Time Resolution [ns]");
  if(1000*gausfit->GetParError(2)>2) tex->DrawLatex(0.6, 0.8, Form("#sigma = %.0f #pm %.0f ps", 1000*gausfit->GetParameter(2), 1000*gausfit->GetParError(2)));
  else tex->DrawLatex(0.6, 0.8, Form("#sigma = %.1f #pm %.1f ps", 1000*gausfit->GetParameter(2), 1000*gausfit->GetParError(2)));
  c->SaveAs(outfile.c_str()); //outfile should end in .pdf
}


void makeTimeResolution( string filename, float photekAmpCut, float photekChargeCut, float centerAmpCut, float centerChargeCut, float MCPAmpCut, float SiPadAmpCut ) {

  DoMultiDeviceStudy( filename.c_str(), photekAmpCut, photekChargeCut, centerAmpCut, centerChargeCut, MCPAmpCut, SiPadAmpCut );

  TFile *_file = TFile::Open( ("output"+filename).c_str() ); //Should be .root

  //Create variables containing hists:
  TH1F *histDeltaTCenter = (TH1F*)_file->Get("histDeltaTCenter"); //Picosil center pixel
  TH1F *histDeltaTCenterSmear = (TH1F*)_file->Get("histDeltaTCenterSmear");
  TH1F *histDeltaTMCP = (TH1F*)_file->Get("histDeltaTMCP"); // MCP
  TH1F *histDeltaTMCPSmear = (TH1F*)_file->Get("histDeltaTMCPSmear");
  TH1F *histDeltaTSiPad = (TH1F*)_file->Get("histDeltaTSiPad"); 
  TH1F *histDeltaTSiPadSmear = (TH1F*)_file->Get("histDeltaTSiPadSmear"); 
  TH1F *histDeltaT_Center_MCP_SiPad_Equal = (TH1F*)_file->Get("histDeltaT_Center_MCP_SiPad_Equal");
  TH1F *histDeltaT_Center_MCP_SiPad_EqualSmear = (TH1F*)_file->Get("histDeltaT_Center_MCP_SiPad_EqualSmear");


  TCanvas *c = new TCanvas ("c","c",800, 600); 
  TLatex *tex = new TLatex();
  tex->SetNDC(); // Sets coords such that (0,0) is bottom left & (1,1) is top right.
  tex->SetTextSize(0.060);
  tex->SetTextFont(42);
  tex->SetTextColor(kBlack);

  c->cd();

  histDeltaTCenter->SetTitle("HGC Center Pixel: TOF");
  histDeltaTCenterSmear->SetTitle("SKIROC Emulation: HGC Center Pixel TOF");
  histDeltaTMCP->SetTitle("MCP: TOF");
  histDeltaTMCPSmear->SetTitle("SKIROC Emulation: MCP TOF");
  histDeltaTSiPad->SetTitle("SiPad: TOF");
  histDeltaTSiPadSmear->SetTitle("SKIROC Emulation: SiPad TOF");
  histDeltaT_Center_MCP_SiPad_Equal->SetTitle("#frac{1}{3} HGC Center Pixel, #frac{1}{3} MCP, #frac{1}{3} SiPad: TOF");
  histDeltaT_Center_MCP_SiPad_EqualSmear->SetTitle("#splitline{SKIROC Emulation: #frac{1}{3} Smeared HGC Center Pixel,}{#frac{1}{3} MCP, #frac{1}{3} SiPad: TOF}");


  PlotDeltaTPDF(c, tex, histDeltaTCenter, "deltaTCenter.pdf");
  PlotDeltaTPDF(c, tex, histDeltaTCenterSmear, "deltaTCenterSmear.pdf");
  PlotDeltaTPDF(c, tex, histDeltaTMCP, "deltaTMCP.pdf");
  PlotDeltaTPDF(c, tex, histDeltaTMCPSmear, "deltaTMCPSmear.pdf");
  PlotDeltaTPDF(c, tex, histDeltaTSiPad, "deltaTSiPad.pdf");
  PlotDeltaTPDF(c, tex, histDeltaTSiPadSmear, "deltaTSiPadSmear.pdf");
  PlotDeltaTPDF(c, tex, histDeltaT_Center_MCP_SiPad_Equal, "histDeltaT_Center_MCP_SiPad_Equal.pdf");
  PlotDeltaTPDF(c, tex, histDeltaT_Center_MCP_SiPad_EqualSmear, "histDeltaT_Center_MCP_SiPad_EqualSmear.pdf");

}


void MultiDeviceStudy_InitialWiring() {

  string infile = "t1065-jun-2016-9.dat-full.root";
  float photekAmpCut = 0.01; //no attenuators
  float photekChargeCut = 0.2;
  float centerAmpCut = 0.001;
  float centerChargeCut = 0.05;
  float MCPAmpCut = 0.01;
  float SiPadAmpCut = 0.0001;
  makeTimeResolution(infile.c_str(), photekAmpCut, photekChargeCut, centerAmpCut, centerChargeCut, MCPAmpCut, SiPadAmpCut); // Outputs PDFs with histograms

}
