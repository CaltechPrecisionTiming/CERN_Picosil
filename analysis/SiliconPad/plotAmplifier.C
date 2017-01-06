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



void MakeAmplificationGraph() {

  //use beam energy for xaxis
  const int nPoints = 24;
  float x[nPoints] = { 11.7, 13.5, 16.7, 20.7, 24.6, 33.6, 47.4, 59, 73, 87, 
		       101, 116, 131, 144, 163, 197, 235, 306, 391, 466, 
		       544, 669, 831, 937 };		       
  float xerr[nPoints] = { 0.8, 0.8, 1, 1, 1, 1, 1.2, 1.2, 1.4, 1.4, 
			  2, 2, 2, 2, 2, 2, 2, 2, 4, 5,
			  5, 5, 7, 7};

  //correcting the extrapolated points by measured / extrapolated ratio
  float y_extrapCorrected[nPoints] = { 27.9/1.37 , 26.8/1.37, 31.6/1.37, 33.9/1.37, 37.3/1.37, 40.2/1.37, 43.5/1.37, 43.1/1.37, 44.1/1.37, 45.5/1.37, 
				       34, 35.6, 37, 40, 42, 45, 48,53, 58, 61, 
				       63, 66, 68, 66 };  
  //uncorrected gain values
  float y_extrap[nPoints] = { 27.9 , 26.8, 31.6, 33.9, 37.3, 40.2, 43.5, 43.1, 44.1, 45.5, 
		       34, 35.6, 37, 40, 42, 45, 48,53, 58, 61, 
		       63, 66, 68, 66 };  

  float yerr[nPoints];
  for (int i=0; i<nPoints; ++i) {   
    y_extrapCorrected[i] =  y_extrapCorrected[i] * 10;
    y_extrap[i] =  y_extrap[i] * 10;
    yerr[i] =  y_extrapCorrected[i] * 0.1;
  }

  //Do interpolation
  const int nBins = 1000;
  float xFine[nBins];
  float yFine[nBins];
  for(int i=0; i< nBins; ++i) {
    xFine[i] = i;
    yFine[i] = GetAmplificationFactor(xFine[i]);
  }

  //Do interpolation
  float xInput[nBins];
  float xOutput[nBins];
  for(int i=0; i< nBins; ++i) {
    xOutput[i] = i*10;
    xInput[i] = i*10 / GetAmplificationFactor(xOutput[i]);
  }



  TGraphErrors *graphAmplificationVsOutputAmplitude = new TGraphErrors(nPoints,x,y_extrapCorrected,xerr,yerr);
  graphAmplificationVsOutputAmplitude->SetLineWidth(3);
  TGraph *graphAmplificationVsOutputAmplitudeInterpolated = new TGraphErrors(nBins,xFine,yFine);
  graphAmplificationVsOutputAmplitudeInterpolated->SetLineWidth(3);
  graphAmplificationVsOutputAmplitudeInterpolated->SetLineColor(kRed);

  TCanvas *c = 0;
  TVirtualFitter *fitter = 0;

  c = new TCanvas("c","c",800,600);
  graphAmplificationVsOutputAmplitude->Draw("AP");
  graphAmplificationVsOutputAmplitude->SetTitle("");
  graphAmplificationVsOutputAmplitude->GetXaxis()->SetTitle("Output Amplitude [mV]");
  graphAmplificationVsOutputAmplitude->GetXaxis()->SetTitleSize(0.045);
  graphAmplificationVsOutputAmplitude->GetXaxis()->SetLabelSize(0.045);
  graphAmplificationVsOutputAmplitude->GetXaxis()->SetTitleOffset(1.0);
  graphAmplificationVsOutputAmplitude->GetYaxis()->SetTitle("Amplification Factor");
  graphAmplificationVsOutputAmplitude->GetYaxis()->SetTitleOffset(1.02);
  graphAmplificationVsOutputAmplitude->GetYaxis()->SetTitleSize(0.045);
  graphAmplificationVsOutputAmplitude->GetYaxis()->SetLabelSize(0.045);
  graphAmplificationVsOutputAmplitude->GetXaxis()->SetRangeUser(0,1000);
  graphAmplificationVsOutputAmplitude->GetYaxis()->SetRangeUser(0,1000);

  graphAmplificationVsOutputAmplitudeInterpolated->Draw("Lsame");
  // graphAmplificationVsOutputAmplitude->Fit("pol1","","");
  // fitter = TVirtualFitter::GetFitter();
  
  c->SaveAs( "AmplificationVsOutputAmplitude.gif" );
  c->SaveAs( "AmplificationVsOutputAmplitude.pdf" );


  c = new TCanvas("c","c",800,600);
  TGraph *graphInputVsOutput = new TGraphErrors(nBins,xInput,xOutput);
  graphInputVsOutput->Draw("APL");
  graphInputVsOutput->SetTitle("");
  graphInputVsOutput->GetXaxis()->SetTitle("Input [mV]");
  graphInputVsOutput->GetXaxis()->SetTitleSize(0.045);
  graphInputVsOutput->GetXaxis()->SetLabelSize(0.045);
  graphInputVsOutput->GetXaxis()->SetTitleOffset(1.0);
  graphInputVsOutput->GetYaxis()->SetTitle("Output [mV]");
  graphInputVsOutput->GetYaxis()->SetTitleOffset(1.02);
  graphInputVsOutput->GetYaxis()->SetTitleSize(0.045);
  graphInputVsOutput->GetYaxis()->SetLabelSize(0.045);
  //graphInputVsOutput->GetXaxis()->SetRangeUser(0,1000);
  graphInputVsOutput->GetYaxis()->SetRangeUser(0,10000);
  graphInputVsOutput->SetMarkerStyle(20);
  graphInputVsOutput->SetMarkerSize(0.2);

  c->SaveAs( "AmplifierInputOutput.gif" );
  c->SaveAs( "AmplifierInputOutput.pdf" );


}


void plotAmplifier() {

  MakeAmplificationGraph();

}
