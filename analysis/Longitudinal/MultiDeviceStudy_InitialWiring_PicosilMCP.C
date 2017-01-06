
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

// This code generates the combined Delta T histogram for the picosil center pixel and the MCP for the initial configuration. Runs 1-52 have this configuration.

void DoMultiDeviceStudy_InitialWiring( string filename, string outputFilename ) {

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

  //Create histogram
  TH1F *histDeltaT = new TH1F("histDeltaT","; Time [ns];Number of Events", 50, 2, 3); // Weights 50-50
  TH1F *histDeltaTWeighted = new TH1F("histDeltaTWeighted","; Time [ns];Number of Events", 50, 3,5); //Weights each event based on charge
  TH1F *histDeltaTCenter = new TH1F("histDeltaTCenter","; Time [ns];Number of Events", 50, 4, 5); //DeltaT of center picosil pixel
  TH1F *histDeltaTMCP = new TH1F("histDeltaTMCP","; Time [ns];Number of Events", 50, 2, 3); //DeltaT of MCP
  

  float totalCenterCharge = 0;
  float totalMCPCharge = 0;

  //read all entries and fill the histogram
  Long64_t nentries = tree->GetEntries();

  //Loop through every event in .root file
  std::cout<<"Number of events in Sample: "<<nentries<<std::endl;  
  for (Long64_t iEntry=0; iEntry<nentries; iEntry++) {
    if (iEntry %1000 == 0) cout << "Processing Event " << iEntry << "\n";
    tree->GetEntry(iEntry);    
  
    float photekTimeGauss0 = gauspeak[0];
    float photekAmp0 = amp[0];
    float photekCharge0 = integral[0];

    float photekTimeGauss3 = gauspeak[27];
    float photekAmp3 = amp[27];
    float photekCharge3 = integral[27];

    float centerAmp = amp[34];
    float centerCharge = integral[34];
    float centerTime = linearTime45[34];

    float MCPAmp = amp[2];
    float MCPCharge = integral[2];
    float MCPTime = linearTime45[2];


    // APPLY EVENT CUTS:
    //require photek (for electron selection)
    if( !(photekAmp0 > 0.05 && photekCharge0 > 2) ) continue;
    
    //require signal in the central pixel
    if( !(centerAmp > 0.03 && centerCharge > 2) ) continue;

    //require Si Pad and MCP minimum amplitude. Possibly will include cut for charge at some point...
    if( !( MCPAmp > 0.03) ) continue;

    //Calculates the Delta T's if the event passes the cuts:
    float DeltaTCenter = photekTimeGauss3 - centerTime; // Will be used to fill the histogram at every event
    float DeltaTMCP = photekTimeGauss0 - MCPTime;

    float DeltaT = (DeltaTCenter + DeltaTMCP)/3; //Evenly weighted
    float DeltaTWeighted = (DeltaTCenter*centerCharge + DeltaTMCP*MCPCharge) / (centerCharge+MCPCharge);
    //The above weights each individual event being added by the charge of each event. 

    totalCenterCharge += centerCharge;
    totalMCPCharge += MCPCharge;

    histDeltaT->Fill(DeltaT);
    histDeltaTWeighted->Fill(DeltaTWeighted);
    histDeltaTCenter->Fill(DeltaTCenter);
    histDeltaTMCP->Fill(DeltaTMCP);
    
    
  }

  // Do Gaussian fit of delta T distributions from (mean-2RMS) to (mean+2RMS)
  double mean = histDeltaT->GetMean();
  double rms = histDeltaT->GetRMS();
  double xmin = mean-2.0*rms;
  double xmax = mean+2.0*rms;
  histDeltaT->Fit("gaus","MLES","",xmin,xmax);

  mean = histDeltaTWeighted->GetMean();
  rms = histDeltaTWeighted->GetRMS();
  xmin = mean-2.0*rms;
  xmax = mean+2.0*rms;
  histDeltaTWeighted->Fit("gaus","MLES","",xmin,xmax);

  double meanCenter = histDeltaTCenter->GetMean();
  rms = histDeltaTCenter->GetRMS();
  xmin = meanCenter-2.0*rms;
  xmax = meanCenter+2.0*rms;
  histDeltaTCenter->Fit("gaus","MLES","",xmin,xmax);

  double meanMCP = histDeltaTMCP->GetMean();
  rms = histDeltaTMCP->GetRMS();
  xmin = meanMCP-2.0*rms;
  xmax = meanMCP+2.0*rms;
  histDeltaTMCP->Fit("gaus","MLES","",xmin,xmax);


  TH1F *histDeltaTWeightsCorrected = new TH1F("histDeltaTWeightsCorrected","; Time [ns];Number of Events", 50, -0.3, 0.3);//Corrects for histograms not initially centered at 0.
  TH1F *histDeltaTWeightsCorrected_totalcharge = new TH1F("histDeltaTWeightsCorrected_totalcharge","; Time [ns];Number of Events", 50, -0.3, 0.3);

  TH1F *histDeltaTCenterAt0 = new TH1F("histDeltaTCenterAt0","; Time [ns];Number of Events", 30, -0.3, 0.3); //shifted to be centered at zero
  TH1F *histDeltaTMCPAt0 = new TH1F("histDeltaTMCPAt0","; Time [ns];Number of Events", 30, -0.3, 0.3); //shifted to be centered at zero


  //Loop through again. This time, subtracting the respective means. Process is almost the same.
  for (Long64_t iEntry=0; iEntry<nentries; iEntry++) {
    if (iEntry %1000 == 0) cout << "Processing Event " << iEntry << "\n";
    tree->GetEntry(iEntry);    
  
    float photekTimeGauss0 = gauspeak[0];
    float photekAmp0 = amp[0];
    float photekCharge0 = integral[0];

    float photekTimeGauss3 = gauspeak[27];
    float photekAmp3 = amp[27];
    float photekCharge3 = integral[27];

    float centerAmp = amp[34];
    float centerCharge = integral[34];
    float centerTime = linearTime45[34];

    float MCPAmp = amp[2];
    float MCPCharge = integral[2];
    float MCPTime = linearTime45[2];

    if( !(photekAmp0 > 0.05 && photekCharge0 > 2) ) continue;
    if( !(centerAmp > 0.03 && centerCharge > 2) ) continue;
    if( !( MCPAmp > 0.03) ) continue;

    //Here are the corrections that center the histograms at 0.
    float DeltaTCenter = photekTimeGauss3 - centerTime - meanCenter; 
    float DeltaTMCP = photekTimeGauss0 - MCPTime - meanMCP;
    float DeltaTWeighted = (DeltaTCenter*centerCharge + DeltaTMCP*MCPCharge) / (centerCharge+MCPCharge);
    float DeltaTWeightsCorrected_totalcharge = (DeltaTCenter*totalCenterCharge + DeltaTMCP*totalMCPCharge) / (totalCenterCharge+totalMCPCharge);

    histDeltaTWeightsCorrected->Fill(DeltaTWeighted); // Uses charge in each device over each individual event as weighting.
    histDeltaTWeightsCorrected_totalcharge->Fill(DeltaTWeightsCorrected_totalcharge); // Uses total charge in each device over whole run as weighting.

    histDeltaTCenterAt0->Fill(DeltaTCenter);
    histDeltaTMCPAt0->Fill(DeltaTMCP);
  }

  mean = histDeltaTWeightsCorrected->GetMean();
  rms = histDeltaTWeightsCorrected->GetRMS();
  xmin = mean-2.0*rms;
  xmax = mean+2.0*rms;
  histDeltaTWeightsCorrected->Fit("gaus","MLES","",xmin,xmax);

  mean = histDeltaTWeightsCorrected_totalcharge->GetMean();
  rms = histDeltaTWeightsCorrected_totalcharge->GetRMS();
  xmin = mean-2.0*rms;
  xmax = mean+2.0*rms;
  histDeltaTWeightsCorrected_totalcharge->Fit("gaus","MLES","",xmin,xmax);

  mean = histDeltaTCenterAt0->GetMean();
  rms = histDeltaTCenterAt0->GetRMS();
  xmin = mean-2.0*rms;
  xmax = mean+2.0*rms;
  histDeltaTCenterAt0->Fit("gaus","MLES","",xmin,xmax);

  mean = histDeltaTMCPAt0->GetMean();
  rms = histDeltaTMCPAt0->GetRMS();
  xmin = mean-2.0*rms;
  xmax = mean+2.0*rms;
  histDeltaTMCPAt0->Fit("gaus","MLES","",xmin,xmax);


  // Creates output root file
  TFile *file = TFile::Open(outputFilename.c_str(), "RECREATE");
  file->cd();
  file->WriteTObject(histDeltaT,"histDeltaT", "WriteDelete");
  file->WriteTObject(histDeltaTWeighted,"histDeltaTWeighted", "WriteDelete");
  file->WriteTObject(histDeltaTCenterAt0,"histDeltaTCenter", "WriteDelete"); //Uses hist centered near 0
  file->WriteTObject(histDeltaTMCPAt0,"histDeltaTMCP0", "WriteDelete"); //Uses hist centered near 0
  file->WriteTObject(histDeltaTWeightsCorrected,"histDeltaTWeightsCorrected", "WriteDelete");
  file->WriteTObject(histDeltaTWeightsCorrected_totalcharge,"histDeltaTWeightsCorrected_totalcharge", "WriteDelete");


  file->Close();
  delete file;

}




void makeTimeResolution( string filename ) {

  TFile *_file = TFile::Open( filename.c_str() ); //Should be .root file

  //Create variables containing hists:
  TH1F* histDeltaT = (TH1F*)_file->Get("histDeltaT"); // each device weighted equally
  // Not doing histDeltaTWeighted, since it gives bad results.
  TH1F* histDeltaTCenter = (TH1F*)_file->Get("histDeltaTCenter"); //Picosil center pixel
  TH1F* histDeltaTMCP = (TH1F*)_file->Get("histDeltaTMCP"); // MCP
  TH1F* histDeltaTWeightsCorrected = (TH1F*)_file->Get("histDeltaTWeightsCorrected"); //Combination of device Delta T's after shifting distributions around 0 and then weighting event by event.
  TH1F* histDeltaTWeightsCorrected_totalcharge = (TH1F*)_file->Get("histDeltaTWeightsCorrected_totalcharge"); //Combination after shifting around 0 and weighting with total charge.

  TCanvas *c = new TCanvas ("c","c",800, 600); 
  TLatex *tex = new TLatex();
  tex->SetNDC(); // Sets coords such that (0,0) is bottom left & (1,1) is top right.
  tex->SetTextSize(0.060);
  tex->SetTextFont(42);
  tex->SetTextColor(kBlack);

  c->cd();
  gStyle->SetOptStat(0); // Hides the parameter box


  histDeltaT->Draw();
  double mean = histDeltaT->GetMean();
  double rms = histDeltaT->GetRMS();
  TF1* fgausEqual = new TF1("fgausEqual","gaus", mean - 2.0*rms, mean + 2.0*rms); // 1-D gaus func defined around hist peak
  histDeltaT->Fit("fgausEqual","QMLE","", mean - 2.0*rms, mean + 2.0*rms); // Fit the hist; Q-quiet, L-log likelihood method, E-Minos errors technique, M-improve fit results
  histDeltaT->GetXaxis()->SetTitle("Time Resolution [ns]");
  histDeltaT->SetTitle("");
  tex->DrawLatex(0.6, 0.80, Form("#sigma = %.0f #pm %.2f ps", 1000*fgausEqual->GetParameter(2), 1000*fgausEqual->GetParError(2)));
  c->SaveAs("deltaTEqual.pdf");

  histDeltaTCenter->Draw();
  mean = histDeltaTCenter->GetMean();
  rms = histDeltaTCenter->GetRMS();
  TF1* fgausCenter = new TF1("fgausCenter","gaus", mean - 2.0*rms, mean + 2.0*rms); 
  histDeltaTCenter->Fit("fgausCenter","QMLE","", mean - 2.0*rms, mean + 2.0*rms);
  histDeltaTCenter->GetXaxis()->SetTitle("Time Resolution [ns]");
  histDeltaTCenter->SetTitle("");
  tex->DrawLatex(0.6, 0.80, Form("#sigma = %.0f #pm %.2f ps", 1000*fgausCenter->GetParameter(2), 1000*fgausCenter->GetParError(2)));
  c->SaveAs("deltaTCenter.pdf");

  histDeltaTMCP->Draw();
  mean = histDeltaTMCP->GetMean();
  rms = histDeltaTMCP->GetRMS();
  TF1* fgausMCP = new TF1("fgausMCP","gaus", mean - 2.0*rms, mean + 2.0*rms); 
  histDeltaTMCP->Fit("fgausMCP","QMLE","", mean - 2.0*rms, mean + 2.0*rms);
  histDeltaTMCP->GetXaxis()->SetTitle("Time Resolution [ns]");
  histDeltaTMCP->SetTitle("");
  tex->DrawLatex(0.6, 0.80, Form("#sigma = %.0f #pm %.2f ps", 1000*fgausMCP->GetParameter(2), 1000*fgausMCP->GetParError(2)));
  c->SaveAs("deltaTMCP.pdf");
  
  histDeltaTWeightsCorrected->Draw();
  mean = histDeltaTWeightsCorrected->GetMean();
  rms = histDeltaTWeightsCorrected->GetRMS();
  TF1* fgausWC = new TF1("fgausWC","gaus", mean - 2.0*rms, mean + 2.0*rms); 
  histDeltaTWeightsCorrected->Fit("fgausWC","QMLE","", mean - 2.0*rms, mean + 2.0*rms);
  histDeltaTWeightsCorrected->GetXaxis()->SetTitle("Time Resolution [ns]");
  histDeltaTWeightsCorrected->SetTitle("");
  tex->DrawLatex(0.6, 0.80, Form("#sigma = %.0f #pm %.2f ps", 1000*fgausWC->GetParameter(2), 1000*fgausWC->GetParError(2)));
  c->SaveAs("deltaTWeightsCorrected.pdf");

  histDeltaTWeightsCorrected_totalcharge->Draw();
  mean = histDeltaTWeightsCorrected_totalcharge->GetMean();
  rms = histDeltaTWeightsCorrected_totalcharge->GetRMS();
  TF1* fgausWCtc = new TF1("fgausWCtc","gaus", mean - 2.0*rms, mean + 2.0*rms); 
  histDeltaTWeightsCorrected_totalcharge->Fit("fgausWCtc","QMLE","", mean - 2.0*rms, mean + 2.0*rms);
  histDeltaTWeightsCorrected_totalcharge->GetXaxis()->SetTitle("Time Resolution [ns]");
  histDeltaTWeightsCorrected_totalcharge->SetTitle("");
  tex->DrawLatex(0.6, 0.80, Form("#sigma = %.0f #pm %.2f ps", 1000*fgausWCtc->GetParameter(2), 1000*fgausWCtc->GetParError(2)));
  c->SaveAs("deltaTWeightsCorrected_totalcharge.pdf");
  
}


void MultiDeviceStudy_InitialWiring_PicosilMCP() {
  DoMultiDeviceStudy_InitialWiring("t1065-jun-2016-46.dat-full.root","output46_PicosilMCP.root");
  makeTimeResolution( "output46_PicosilMCP.root" ); // Outputs PDFs with histograms
}
