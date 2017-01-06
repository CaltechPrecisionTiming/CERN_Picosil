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

static const int nPoints = 24;
static float outputAmplitude[nPoints] = { 11.7, 13.5, 16.7, 20.7, 24.6, 33.6, 47.4, 59, 73, 87, 
					  101, 116, 131, 144, 163, 197, 235, 306, 391, 466, 
					  544, 669, 831, 937 };		       

//correcting the extrapolated points by measured / extrapolated ratio
static float amplificationFactor[nPoints] = { 27.9/1.37*10 , 26.8/1.37*10, 31.6/1.37*10, 33.9/1.37*10, 37.3/1.37*10, 40.2/1.37*10, 43.5/1.37*10, 43.1/1.37*10, 44.1/1.37*10, 45.5/1.37*10, 
					      34*10, 35.6*10, 37*10, 40*10, 42*10, 45*10, 48*10, 53*10, 58*10, 61*10, 
					      63*10, 66*10, 68*10, 66*10 };  

double GetAmplificationFactor ( double measuredAmplitude ) {
  
  int index_firstBinAboveInput = -1;
  for (int i=0; i < nPoints; ++i) {
    index_firstBinAboveInput = i;
    if (measuredAmplitude < outputAmplitude[i]) break;
  }
  
  double answer = 0; 

  if (measuredAmplitude > outputAmplitude[21]) answer =amplificationFactor[21];
  else if (index_firstBinAboveInput == 0) answer = amplificationFactor[0];
  else {
    
    //cout << "index_firstBinAboveInput = " << index_firstBinAboveInput << " : "
    //	 << amplificationFactor[index_firstBinAboveInput-1] << " " << outputAmplitude[index_firstBinAboveInput]
    //	 << "\n";
    double x = measuredAmplitude - outputAmplitude[index_firstBinAboveInput-1];
    double y = amplificationFactor[index_firstBinAboveInput-1] + x * (amplificationFactor[index_firstBinAboveInput] - amplificationFactor[index_firstBinAboveInput-1]) / (outputAmplitude[index_firstBinAboveInput] - outputAmplitude[index_firstBinAboveInput-1]);
    //cout << "x = " << x << " , y = " << y << "\n";
    answer = y;
  }

  //cout << measuredAmplitude << " " << answer << "\n";

  return answer;
  
}




void MakeChargePlot(string filename, string plotname, double ampCutOnPhotek, double attenuationFactor, 
		    double xmin, double xmax, 
		    double fitmin, double fitmax, bool forMIP = false) {
  // Get the tree

  TFile *inputfile = TFile::Open(filename.c_str(),"READ");
  TTree *tree = (TTree*)inputfile->Get("pulse");

  // get the variables from the ntuple
  float amp[36];
  float integral[36];
  float gauspeak[36];

  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("gauspeak",1);
  tree->SetBranchStatus("amp",1);
  tree->SetBranchStatus("int",1);
  tree->SetBranchAddress("gauspeak",gauspeak);
  tree->SetBranchAddress("amp",amp);
  tree->SetBranchAddress("int",integral);

  //create histograms
  TH1F *histIntCharge;
  histIntCharge = new TH1F("histIntCharge","; Integrated Charge [fC];Number of Events", 50, xmin,xmax);
  TH1F *histIntMIP;
  histIntMIP = new TH1F("histIntMIP","; Integrated Charge / Charge for MIP ; Number of Events", 50, xmin/6.5,xmax/6.5);
  
  //read all entries and fill the histograms
  Long64_t nentries = tree->GetEntries();

  std::cout<<"Number of events in Sample: "<<nentries<<std::endl;  
  for (Long64_t iEntry=0;iEntry<nentries;iEntry++) {
    if (iEntry %1000 == 0) cout << "Processing Event " << iEntry << "\n";
    tree->GetEntry(iEntry);    
  
    float photekTimeGauss = gauspeak[18];
    float siliconTimeGauss = gauspeak[21];
    float photekAmp = amp[18];
    float siliconAmp = amp[21];
    float photekIntegral = integral[18];
    float siliconIntegral = integral[21];
       
    //use photek amplitude cut for electron ID
    //cout << "test: " << photekAmp << " " << siliconIntegral << "\n";
    if( !(photekAmp > ampCutOnPhotek)) continue;
    
    //Fill histogram
    double amplificationFactor = GetAmplificationFactor( 1000 * amp[21] * (attenuationFactor/10) );
    if (forMIP) {

      //don't fill overflow bins
      if (1000* siliconIntegral * attenuationFactor / amplificationFactor / 1.37 > xmax) continue;

      //for MIPs only, we use the amplification factor without the 1.37 correction factor
      histIntCharge->Fill(1000* siliconIntegral * attenuationFactor / amplificationFactor / 1.37 ); 
      histIntMIP->Fill(1000* siliconIntegral * attenuationFactor / amplificationFactor / 6.5 / 1.37); 
    } else {

      //don't fill overflow bins
      if (1000* siliconIntegral * attenuationFactor / amplificationFactor > xmax) continue;

      histIntCharge->Fill(1000* siliconIntegral * attenuationFactor / amplificationFactor );
      histIntMIP->Fill(1000* siliconIntegral * attenuationFactor / amplificationFactor / 6.5);
      //cout << 1000* amp[21] << " : " << amplificationFactor << " : " << siliconIntegral * attenuationFactor / amplificationFactor << "\n";
    }

  }


  TCanvas * c = 0;


  //Energy plot
  c = new TCanvas("c","c",600,600);  
  c->SetRightMargin(0.05);
  c->SetLeftMargin(0.17);
  histIntCharge->SetAxisRange(xmin,xmax,"X");
  histIntCharge->SetTitle("");
  histIntCharge->GetXaxis()->SetTitle("Integrated Charge [fC]");
  histIntCharge->GetXaxis()->SetTitleSize(0.045);
  histIntCharge->GetXaxis()->SetLabelSize(0.045);
  histIntCharge->GetYaxis()->SetTitle("Number of Events");
  histIntCharge->GetYaxis()->SetTitleOffset(1.3);
  histIntCharge->GetYaxis()->SetTitleSize(0.05);
  histIntCharge->GetYaxis()->SetLabelSize(0.045);
  histIntCharge->GetYaxis()->SetLabelOffset(0.015);
  histIntCharge->GetYaxis()->SetTitleOffset(1.7);
  histIntCharge->SetMaximum(1.2*histIntCharge->GetMaximum());
  histIntCharge->Draw();
  histIntCharge->SetStats(0);
  histIntCharge->Fit("gaus","","",fitmin,fitmax);
  TVirtualFitter * fitter = TVirtualFitter::GetFitter();
  
  TLatex *tex = new TLatex();
  tex->SetNDC();
  tex->SetTextSize(0.050);
  tex->SetTextFont(42);
  tex->SetTextColor(kBlack);
  /* tex->DrawLatex(0.45, 0.85, Form("Mean = %.1f #pm %.1f %s",fitter->GetParameter(1),TMath::Max(0.01,fitter->GetParError(1)),"fC")); */
  /* tex->DrawLatex(0.45, 0.80, Form("#sigma = %.1f #pm %.1f %s",fitter->GetParameter(2),TMath::Max(0.01,fitter->GetParError(2)),"fC")); */
  tex->DrawLatex(0.45, 0.85, Form("Mean = %.1f %s",fitter->GetParameter(1),"fC"));
  tex->DrawLatex(0.45, 0.80, Form("#sigma = %.1f %s",fitter->GetParameter(2),"fC"));
  //tex->DrawLatex(0.15, 0.92, Form("Attenuation Factor = %.3f",attenuationFactor));
  
  c->SaveAs( Form("%s_charge.gif", plotname.c_str()) );
  c->SaveAs( Form("%s_charge.pdf", plotname.c_str()) );
 
  if (!forMIP) {
    //Energy plot
    c = new TCanvas("c","c",600,600);  
    c->SetRightMargin(0.05);
    c->SetLeftMargin(0.15);
    c->SetBottomMargin(0.13);
    histIntMIP->SetAxisRange(xmin/6.5,xmax/6.5,"X");
    histIntMIP->SetTitle("");
    histIntMIP->GetXaxis()->SetTitle("Integrated Charge [ Q_{MIP} ]");
    histIntMIP->GetXaxis()->SetTitleSize(0.045);
    histIntMIP->GetXaxis()->SetLabelSize(0.045);
    histIntMIP->GetXaxis()->SetTitleOffset(1.15);
    histIntMIP->GetYaxis()->SetTitle("Number of Events");
    histIntMIP->GetYaxis()->SetTitleSize(0.05);
    histIntMIP->GetYaxis()->SetLabelSize(0.045);
    histIntMIP->GetYaxis()->SetLabelOffset(0.015);
    histIntMIP->GetYaxis()->SetTitleOffset(1.5);
    histIntMIP->SetMaximum(1.2*histIntMIP->GetMaximum());
    histIntMIP->Draw();
    histIntMIP->SetStats(0);
    histIntMIP->Fit("gaus","","",fitmin/6.5,fitmax/6.5);
    fitter = TVirtualFitter::GetFitter();
  
    tex = new TLatex();
    tex->SetNDC();
    tex->SetTextSize(0.040);
    tex->SetTextFont(42);
    tex->SetTextColor(kBlack);
    tex->DrawLatex(0.45, 0.85, Form("Mean = %.0f %s",fitter->GetParameter(1),"MIPs"));
    tex->DrawLatex(0.45, 0.80, Form("#sigma = %.0f %s",fitter->GetParameter(2),"MIPs"));
    //tex->DrawLatex(0.15, 0.92, Form("Attenuation Factor = %.3f",attenuationFactor));
  
    c->SaveAs( Form("%s_chargeMIP.gif", plotname.c_str()) );
    c->SaveAs( Form("%s_chargeMIP.pdf", plotname.c_str()) );
  }

}




