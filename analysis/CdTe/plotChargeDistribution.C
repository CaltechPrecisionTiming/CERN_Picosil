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



void MakeAmplitudeVsBeamEnergyGraph() {

  //use beam energy for xaxis
  const int nPoints_H2 = 3;
  float x_H2[nPoints_H2] = { 50.0 , 100.0, 200.0 };
  float xerr_H2[nPoints_H2] = { 0, 0, 0 };
  float y_charge_H2[nPoints_H2] = {  3.6 , 6.8, 11.5 }; 
  float yerr_charge_H2[nPoints_H2] = { 0.09, 0.05, 0.12 };
  float yResolution_charge_H2[nPoints_H2] = { 0.6, 1.2, 2.2 };


  TGraphErrors *graphChargeVsEnergyAt6X0_Resolution = new TGraphErrors(nPoints_H2,x_H2,y_charge_H2,xerr_H2,yResolution_charge_H2);
  TGraphErrors *graphChargeVsEnergyAt6X0 = new TGraphErrors(nPoints_H2,x_H2,y_charge_H2,xerr_H2,yerr_charge_H2);
  graphChargeVsEnergyAt6X0_Resolution->SetMarkerStyle(20);
  graphChargeVsEnergyAt6X0_Resolution->SetMarkerSize(1.25);
  graphChargeVsEnergyAt6X0_Resolution->SetLineColor(kGreen+1);
  graphChargeVsEnergyAt6X0_Resolution->SetLineWidth(2);
  graphChargeVsEnergyAt6X0->SetLineColor(kBlue);
  graphChargeVsEnergyAt6X0->SetLineWidth(2);


  TCanvas *c = 0;
  TVirtualFitter *fitter = 0;
  TLatex *tex = 0;

  c = new TCanvas("c","c",800,800);
  graphChargeVsEnergyAt6X0_Resolution->Draw("AP");
  graphChargeVsEnergyAt6X0_Resolution->SetTitle("");
  graphChargeVsEnergyAt6X0_Resolution->GetXaxis()->SetTitle("Electron Beam Energy [GeV]");
  graphChargeVsEnergyAt6X0_Resolution->GetXaxis()->SetTitleSize(0.045);
  graphChargeVsEnergyAt6X0_Resolution->GetXaxis()->SetLabelSize(0.045);
  graphChargeVsEnergyAt6X0_Resolution->GetXaxis()->SetTitleOffset(1.05);
  graphChargeVsEnergyAt6X0_Resolution->GetYaxis()->SetTitle("Integrated Charge [pC]");
  graphChargeVsEnergyAt6X0_Resolution->GetYaxis()->SetTitleOffset(1.2);
  graphChargeVsEnergyAt6X0_Resolution->GetYaxis()->SetTitleSize(0.05);
  graphChargeVsEnergyAt6X0_Resolution->GetYaxis()->SetLabelSize(0.045);
  //graphChargeVsEnergyAt6X0_Resolution->GetXaxis()->SetRangeUser(0,40);
  // graphChargeVsEnergyAt6X0_Resolution->SetMarkerStyle(20);
  // graphChargeVsEnergyAt6X0_Resolution->SetMarkerSize(1.5);
  graphChargeVsEnergyAt6X0_Resolution->GetYaxis()->SetRangeUser(0,15);

  // graphChargeVsEnergyAt6X0_Resolution->Fit("pol1","","");
  // fitter = TVirtualFitter::GetFitter();
  c->SetLeftMargin(0.15);
  c->SetBottomMargin(0.12);

  tex = new TLatex();
  tex->SetNDC();
  tex->SetTextSize(0.050);
  tex->SetTextFont(42);
  tex->SetTextColor(kBlack);
  tex->DrawLatex(0.30, 0.93, "Absorber : 6 X_{0} W-Pb");


  c->SaveAs( "ChargeVsEnergyAt6X0.C" );
  c->SaveAs( "ChargeVsEnergyAt6X0.gif" );
  c->SaveAs( "ChargeVsEnergyAt6X0.pdf" );


  //use beam energy for xaxis
  const int nPoints_T9 = 5;
  float x_T9[nPoints_T9] = {0, 2, 3.5, 5, 7 };
  float xerr_T9[nPoints_T9] = {0, 0, 0, 0, 0 };

  //Computed from charge integral
  //float y_charge_T9[nPoints_T9] = {  0.16, 0.26, 0.26, 0.40 }; 
  //float yerr_charge_T9[nPoints_T9] = { 0.008, 0.015, 0.03, 0.026 };
  //float yResolution_charge_T9[nPoints_T9] = { 0.08, 0.15, 0.17, 0.20 };

  //Computed from 
  float y_charge_T9[nPoints_T9] = {-100,  0.076, 0.104, 0.161, 0.243 }; 
  float yerr_charge_T9[nPoints_T9] = { 0, 0.006, 0.003, 0.005, 0.008 };
  float yResolution_charge_T9[nPoints_T9] = {0, 0.062, 0.070, 0.094, 0.118 };

  TGraphErrors *graphChargeVsEnergyAt2X0_Resolution = new TGraphErrors(nPoints_T9,x_T9,y_charge_T9,xerr_T9,yResolution_charge_T9);
  TGraphErrors *graphChargeVsEnergyAt2X0 = new TGraphErrors(nPoints_T9,x_T9,y_charge_T9,xerr_T9,yerr_charge_T9);
  graphChargeVsEnergyAt2X0_Resolution->SetMarkerStyle(20);
  graphChargeVsEnergyAt2X0_Resolution->SetMarkerSize(1.25);
  graphChargeVsEnergyAt2X0_Resolution->SetLineWidth(3);
  graphChargeVsEnergyAt2X0_Resolution->SetLineColor(kGreen+1);
  graphChargeVsEnergyAt2X0->SetLineColor(kBlue);
  graphChargeVsEnergyAt2X0->SetLineWidth(2);

  c = new TCanvas("c","c",800,800);
  graphChargeVsEnergyAt2X0_Resolution->Draw("AP");
  graphChargeVsEnergyAt2X0_Resolution->SetTitle("");
  graphChargeVsEnergyAt2X0_Resolution->GetXaxis()->SetTitle("Electron Beam Energy [GeV]");
  graphChargeVsEnergyAt2X0_Resolution->GetXaxis()->SetTitleSize(0.045);
  graphChargeVsEnergyAt2X0_Resolution->GetXaxis()->SetLabelSize(0.045);
  graphChargeVsEnergyAt2X0_Resolution->GetXaxis()->SetTitleOffset(1.05);
  graphChargeVsEnergyAt2X0_Resolution->GetYaxis()->SetTitle("Integrated Charge [pC]");
  graphChargeVsEnergyAt2X0_Resolution->GetYaxis()->SetTitleOffset(1.5);
  graphChargeVsEnergyAt2X0_Resolution->GetYaxis()->SetTitleSize(0.05);
  graphChargeVsEnergyAt2X0_Resolution->GetYaxis()->SetLabelSize(0.045);
  //graphChargeVsEnergyAt2X0_Resolution->GetXaxis()->SetRangeUser(0,40);
  graphChargeVsEnergyAt2X0_Resolution->GetYaxis()->SetRangeUser(0,0.4);

  graphChargeVsEnergyAt2X0->Draw("PE1same");
  // graphChargeVsEnergyAt2X0->Fit("pol1","","");
  // fitter = TVirtualFitter::GetFitter();

  c->SetLeftMargin(0.17);
  c->SetBottomMargin(0.12);

  tex = new TLatex();
  tex->SetNDC();
  tex->SetTextSize(0.050);
  tex->SetTextFont(42);
  tex->SetTextColor(kBlack);
  tex->DrawLatex(0.35, 0.93, "Absorber : 2 X_{0} Lead");

  c->SaveAs( "ChargeVsEnergyAt2X0.C" );
  c->SaveAs( "ChargeVsEnergyAt2X0.gif" );
  c->SaveAs( "ChargeVsEnergyAt2X0.pdf" );


}



void makeChargeDistributionH2(string filename, string plotname, string plotTitle,
			      double ampCutOnPhotek, 
			      double beamXMin, double beamXMax, double beamYMin, double beamYMax,
			      int nbins, double xmin, double xmax, double fitmin, double fitmax) {
  

  TFile *inputfile = TFile::Open(filename.c_str(),"READ");
  
  TTree *tree = (TTree*)inputfile->Get("t1065");

  // get the variables from the ntuple
  float amp[36];
  float integral[36];
  float gauspeak[36];
  float linearTime30[36];
  float beamX;
  float beamY;

  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("gauspeak",1);
  tree->SetBranchStatus("amp",1);
  tree->SetBranchStatus("intfull",1);
  tree->SetBranchStatus("linearTime30",1);
  tree->SetBranchStatus("TDCx",1);
  tree->SetBranchStatus("TDCy",1);
  tree->SetBranchAddress("gauspeak",gauspeak);
  tree->SetBranchAddress("amp",amp);
  tree->SetBranchAddress("intfull",integral);
  tree->SetBranchAddress("linearTime30",linearTime30);
  tree->SetBranchAddress("TDCx",&beamX);
  tree->SetBranchAddress("TDCy",&beamY);

  //create histograms
  TH1F *histIntCharge;
  histIntCharge = new TH1F("histIntCharge","; Integrated Charge [pC];Number of Events", nbins, xmin, xmax);
  TH1F *histAmp;
  histAmp = new TH1F("histAmp","; Amplitude [mV];Number of Events", 100, 0, 1.0);

  
  //read all entries and fill the histograms
  Long64_t nentries = tree->GetEntries();

  std::cout<<"Number of events in Sample: "<<nentries<<std::endl;  
  for (Long64_t iEntry=0;iEntry<nentries;iEntry++) {
    if (iEntry %1000 == 0) 
      cout << "Processing Event " << iEntry << "\n";
    tree->GetEntry(iEntry);    
    // cout << "here1\n";
    float photekTimeGauss = gauspeak[0];
    float CdTeTime = linearTime30[1];
    float photekAmp = amp[0]*(10.0);
    float CdTeAmp = amp[1]*3.16228*(1.0/63.0957);
    float photekCharge = integral[0]*(10.0);
    float CdTeCharge = integral[1]*3.16228*(1.0/63.0957);
    // cout << "here2\n";
       
    //use photek amplitude cut for electron ID
    //cout << "test: " << photekAmp << " " << siliconIntegral << "\n";
    if( !(photekAmp > 0.5)) continue;
    // if(!(beamX > 4.5 && beamX < 12.5)) continue;
    // if(!(beamY > -6.5 && beamY < 1.5)) continue;
    if(!(beamX > beamXMin && beamX < beamXMax)) continue;
    if(!(beamY > beamYMin && beamY < beamYMax)) continue;
    // cout << "here3\n";

    //don't fill overflow bins
    //if (1000* siliconIntegral * attenuationFactor / amplificationFactor > xmax) continue;
    
    histIntCharge->Fill( CdTeCharge );
    histAmp->Fill( amp[1] );

    //cout << CdTeCharge << " " << beamX << " " << beamY << " " << CdTeAmp << "\n";

    //cout << 1000* amp[21] << " : " << amplificationFactor << " : " << siliconIntegral * attenuationFactor / amplificationFactor << "\n";
 
  }


  TCanvas * c = 0;


  //Energy plot
  c = new TCanvas("c","c",600,600);  
  c->SetRightMargin(0.05);
  c->SetLeftMargin(0.17);
  histIntCharge->SetAxisRange(xmin,xmax,"X");
  histIntCharge->SetTitle("");
  histIntCharge->GetXaxis()->SetTitle("Integrated Charge [pC]");
  histIntCharge->GetXaxis()->SetTitleSize(0.045);
  histIntCharge->GetXaxis()->SetLabelSize(0.045);
  histIntCharge->GetYaxis()->SetTitle("Number of Events");
  histIntCharge->GetYaxis()->SetTitleOffset(1.2);
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
  tex->DrawLatex(0.45, 0.85, Form("Mean = %.1f %s",fitter->GetParameter(1),"pC"));
  tex->DrawLatex(0.45, 0.80, Form("#sigma = %.1f %s",fitter->GetParameter(2),"pC"));

  tex->DrawLatex(0.13, 0.93, Form("%s", plotTitle.c_str()));
    
  c->SaveAs( Form("%s_charge.C", plotname.c_str()) );
  c->SaveAs( Form("%s_charge.gif", plotname.c_str()) );
  c->SaveAs( Form("%s_charge.pdf", plotname.c_str()) );
 
  // c = new TCanvas("c","c",600,600);  
  // histAmp->Draw();

}

void makeChargeDistributionT9(string filename, string plotname, string plotTitle,
			      double ampCutOnMCP, double ampCutOnLYSO, double ampCutOnTrigger,
			    double beamXMin, double beamXMax, double beamYMin, double beamYMax,
			    int nbins, double xmin, double xmax, double fitmin, double fitmax) {


  TFile *inputfile = TFile::Open(filename.c_str(),"READ");
  
  TTree *tree = (TTree*)inputfile->Get("t1065");

  // get the variables from the ntuple
  float amp[36];
  float integral[36];
  float gauspeak[36];
  float linearTime30[36];
  float beamX;
  float beamY;

  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("gauspeak",1);
  tree->SetBranchStatus("amp",1);
  tree->SetBranchStatus("intfull",1);
  tree->SetBranchStatus("linearTime30",1);
  tree->SetBranchStatus("TDCx",1);
  tree->SetBranchStatus("TDCy",1);
  tree->SetBranchAddress("gauspeak",gauspeak);
  tree->SetBranchAddress("amp",amp);
  tree->SetBranchAddress("intfull",integral);
  tree->SetBranchAddress("linearTime30",linearTime30);
  tree->SetBranchAddress("TDCx",&beamX);
  tree->SetBranchAddress("TDCy",&beamY);

  //create histograms
  TH1F *histIntCharge;
  histIntCharge = new TH1F("histIntCharge","; Integrated Charge [pC];Number of Events", nbins, xmin, xmax);

  
  //read all entries and fill the histograms
  Long64_t nentries = tree->GetEntries();

  std::cout<<"Number of events in Sample: "<<nentries<<std::endl;  
  for (Long64_t iEntry=0;iEntry<nentries;iEntry++) {
    if (iEntry %1000 == 0) 
      cout << "Processing Event " << iEntry << "\n";
    tree->GetEntry(iEntry);    
    // cout << "here1\n";
    float MCPTimeGauss = gauspeak[0];
    float CdTeTime = linearTime30[1];
    float MCPAmp = amp[0];
    float LYSOAmp = amp[2];
    float TriggerAmp = amp[3];
    float CherenkovAmp = amp[7];
    float CdTeAmp = amp[1]*(1.0/63.0957);
    float MCPCharge = integral[0];
    float CdTeCharge = integral[1]*(1.0/63.0957);
    // cout << "here2\n";
       
    //use MCP amplitude cut for electron ID
    //cout << "test: " << MCPAmp << " " << siliconIntegral << "\n";
    if( !(MCPAmp > ampCutOnMCP)) continue;
    if( !(TriggerAmp > ampCutOnTrigger)) continue;
    if( !(LYSOAmp > ampCutOnLYSO)) continue;
    if(!(beamX > beamXMin && beamX < beamXMax)) continue;
    if(!(beamY > beamYMin && beamY < beamYMax)) continue;
     // cout << "here3\n";

    //don't fill overflow bins
    //if (1000* siliconIntegral * attenuationFactor / amplificationFactor > xmax) continue;
    
    histIntCharge->Fill( CdTeCharge );

    //cout << CdTeCharge << " " << beamX << " " << beamY << " " << CdTeAmp << "\n";

    //cout << 1000* amp[21] << " : " << amplificationFactor << " : " << siliconIntegral * attenuationFactor / amplificationFactor << "\n";
 
  }


  TCanvas * c = 0;


  //Energy plot
  c = new TCanvas("c","c",600,600);  
  c->SetRightMargin(0.05);
  c->SetLeftMargin(0.17);
  histIntCharge->SetAxisRange(xmin,xmax,"X");
  histIntCharge->SetTitle("");
  histIntCharge->GetXaxis()->SetTitle("Integrated Charge [pC]");
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
  tex->DrawLatex(0.45, 0.85, Form("Mean = %.2f %s",fitter->GetParameter(1),"pC"));
  tex->DrawLatex(0.45, 0.80, Form("#sigma = %.2f %s",fitter->GetParameter(2),"pC"));

  tex->DrawLatex(0.18, 0.93, Form("%s", plotTitle.c_str()));

  c->SaveAs( Form("%s_charge.gif", plotname.c_str()) );
  c->SaveAs( Form("%s_charge.pdf", plotname.c_str()) );
 

}



void plotChargeDistribution(double energy = -1) {

  if (energy == 0) {
    makeChargeDistributionH2( "/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Timing/Nov2016CERN/ntuples_v6/analysis_5568.root", 
			    "100GeV", "100 GeV Electrons, 6 X_{0} W-Pb Absorber", 0.5,
			    2.5,13.5,-8.0,2.0,
			    50, 0, 12, 4, 8.5
			    );
  }

  if (energy == 100) {
    makeChargeDistributionH2( "/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Timing/Nov2016CERN/ntuples_v6/analysis_5568.root", 
			    "100GeV", "100 GeV Electrons, 6 X_{0} W-Pb Absorber", 0.5,
			    // 2.5,13.5,-8.0,2.0,
			    6.5,9.5,-3.5,-0.5,
			    50, 0, 12, 5, 11.0
			    );
  }

  if (energy == 200) {
    makeChargeDistributionH2( "/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Timing/Nov2016CERN/ntuples_v6/analysis_5570.root", 
			    "200GeV", "200 GeV Electrons, 6 X_{0} W-Pb Absorber", 0.5,
			    // 2.5,13.5,-8,2,
			    6.5,9.5,-3.5,-0.5,
			    50, 0, 25, 8.0, 20.0
			    );
  }
  
  if (energy == 50) {
    makeChargeDistributionH2( "/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Timing/Nov2016CERN/ntuples_v6/analysis_5571-5576.root", 
			    "50GeV", "50 GeV Electrons, 6 X_{0} Tungsten Absorber" , 0.5,
			    // 2.5,13.5,-8.0,2.0,
			    6.5,9.5,-3.5,-0.5,
			    25, 0,7,2.5,6.0
			    );
  }

  if (energy == 2) {
    makeChargeDistributionT9( "/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Timing/Nov2016CERN/ntuples_v6/analysis_4627.root", 
			      "2GeV", "2 GeV Electrons, 2 X_{0} Lead Absorber", 0.025, 0.85, 0.15,
			      -20,20, -10,10,
			      25, 0,0.7,0.06,0.28
			      );
  }

  if (energy == 3.5) {
    makeChargeDistributionT9( "/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Timing/Nov2016CERN/ntuples_v6/analysis_4626.root", 
			      "3p5GeV", "3.5 GeV Electrons, 2 X_{0} Lead Absorber", 0.025, 0.85, 0.15,
			      -20,20, -10,10,
			      20, 0,1.0,0.05,0.5
			      );
  }

  if (energy == 5) {
    makeChargeDistributionT9( "/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Timing/Nov2016CERN/ntuples_v6/analysis_4629.root", 
			      "5GeV", "5 GeV Electrons, 2 X_{0} Lead Absorber", 0.025, 0.80, 0.15,
			      -20,20, -10,10,
			      20, 0,1.0,0.05,0.55
			      );
  }

  if (energy == 7) {
    makeChargeDistributionT9( "/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Timing/Nov2016CERN/ntuples_v6/analysis_4632.root", 
			      "7GeV", "7 GeV Electrons, 2 X_{0} Lead Absorber", 0.025, 0.87, 0.15,
			      -20,20, -10,10,
			      20, 0,2.0,0.15,0.75
			      );
  }


  if (energy == -1) {
    MakeAmplitudeVsBeamEnergyGraph();
  }
 

}
