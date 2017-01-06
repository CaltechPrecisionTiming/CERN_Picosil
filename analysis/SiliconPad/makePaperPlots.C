#include "TFile.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TF1.h"
#include "TGraphErrors.h"


// make the plots as a function of beam energy
void makeTimeResolutionVsBeamEnergy(){

  TFile *_file0 = TFile::Open("doc/papers/SiliconPad/data/run112_lin30_DeltaT.root");
  TFile *_file1 = TFile::Open("doc/papers/SiliconPad/data/run109_lin30_DeltaT.root");
  TFile *_file2 = TFile::Open("doc/papers/SiliconPad/data/run110_lin30_DeltaT.root");
  TFile *_file3 = TFile::Open("doc/papers/SiliconPad/data/run114_lin30_DeltaT.root");

  TH1F* deltaT0 = (TH1F*)_file0->Get("DeltaT");
  TH1F* deltaT1 = (TH1F*)_file1->Get("DeltaT");
  TH1F* deltaT2 = (TH1F*)_file2->Get("DeltaT");
  TH1F* deltaT3 = (TH1F*)_file3->Get("DeltaT");
  deltaT0->SetStats(false);
  deltaT1->SetStats(false);
  deltaT2->SetStats(false);
  deltaT3->SetStats(false);

  TCanvas *c = 0;
  TLatex *tex = 0;

  

  c = new TCanvas ("c","c",800, 800);
  c->SetRightMargin(0.05);
  c->SetLeftMargin(0.15);
  deltaT0->Draw();
  deltaT0->GetXaxis()->SetRangeUser(-0.5,0.5);
  TF1* fgaus0 = new TF1("fgaus0","gaus", deltaT0->GetMean() - 1.5*deltaT0->GetRMS(), deltaT0->GetMean() + 1.5*deltaT0->GetRMS());
  deltaT0->Fit("fgaus0","Q","", deltaT0->GetMean() - 1.5*deltaT0->GetRMS(), deltaT0->GetMean() + 1.5*deltaT0->GetRMS());
  float res0 = fgaus0->GetParameter(2);
  float error0 = fgaus0->GetParError(2);
  deltaT0->GetYaxis()->SetTitle("Number of Events");
  deltaT0->GetYaxis()->SetTitleSize(0.045);
  deltaT0->GetYaxis()->SetTitleOffset(1.7);
  deltaT0->GetYaxis()->SetLabelSize(0.045);
  deltaT0->GetYaxis()->SetLabelOffset(0.015);
  deltaT0->GetYaxis()->SetTitle("Number of Events");
  deltaT0->GetXaxis()->SetTitle("#Delta t [ns]");
  deltaT0->GetXaxis()->SetTitleSize(0.045);
  deltaT0->GetXaxis()->SetLabelSize(0.045);
  deltaT0->SetTitle("");
  tex = new TLatex();
  tex->SetNDC();
  tex->SetTextSize(0.050);
  tex->SetTextFont(42);
  tex->SetTextColor(kBlack);
  tex->DrawLatex(0.62, 0.80, Form("#sigma = %.0f #pm %.1f ps", 1000*fgaus0->GetParameter(2), 1000*fgaus0->GetParError(2)));
  tex->DrawLatex(0.15, 0.93, "4 GeV Electron Beam, 6 X_{0} Absorber");
  c->SaveAs("deltaT_4GeV_6X0.pdf");

  c = new TCanvas ("c","c",800, 800);
  c->SetRightMargin(0.05);
  c->SetLeftMargin(0.15);
  deltaT1->Draw();
  deltaT1->GetXaxis()->SetRangeUser(-0.2,0.2);
  TF1* fgaus1 = new TF1("fgaus1","gaus", deltaT1->GetMean() - 1.5*deltaT1->GetRMS(), deltaT1->GetMean() + 1.5*deltaT1->GetRMS());
  deltaT1->Fit("fgaus1","Q","", deltaT1->GetMean() - 1.5*deltaT1->GetRMS(), deltaT1->GetMean() + 1.5*deltaT1->GetRMS());
  float res1 = fgaus1->GetParameter(2);
  float error1 = fgaus1->GetParError(2);
  deltaT1->GetYaxis()->SetTitle("Number of Events");
  deltaT1->GetYaxis()->SetTitleSize(0.045);
  deltaT1->GetYaxis()->SetTitleOffset(1.7);
  deltaT1->GetYaxis()->SetLabelSize(0.045);
  deltaT1->GetYaxis()->SetLabelOffset(0.015);
  deltaT1->GetYaxis()->SetTitle("Number of Events");
  deltaT1->GetXaxis()->SetTitle("#Delta t [ns]");
  deltaT1->GetXaxis()->SetTitleSize(0.045);
  deltaT1->GetXaxis()->SetLabelSize(0.045);
  deltaT1->SetTitle("");
  tex = new TLatex();
  tex->SetNDC();
  tex->SetTextSize(0.050);
  tex->SetTextFont(42);
  tex->SetTextColor(kBlack);
  tex->DrawLatex(0.62, 0.80, Form("#sigma = %.0f #pm %.1f ps", 1000*fgaus1->GetParameter(2), 1000*fgaus1->GetParError(2)));
  tex->DrawLatex(0.15, 0.93, "8 GeV Electron Beam, 6 X_{0} Absorber");
  c->SaveAs("deltaT_8GeV_6X0.pdf");


  c = new TCanvas ("c","c",800, 800);
  c->SetRightMargin(0.05);
  c->SetLeftMargin(0.15);
  deltaT2->Draw();
  deltaT2->GetXaxis()->SetRangeUser(-0.2,0.2);
  TF1* fgaus2 = new TF1("fgaus2","gaus", deltaT2->GetMean() - 1.5*deltaT2->GetRMS(), deltaT2->GetMean() + 1.5*deltaT2->GetRMS());
  deltaT2->Fit("fgaus2","Q","", deltaT2->GetMean() - 1.5*deltaT2->GetRMS(), deltaT2->GetMean() + 1.5*deltaT2->GetRMS());
  float res2 = fgaus2->GetParameter(2);
  float error2 = fgaus2->GetParError(2);
  deltaT2->GetYaxis()->SetTitle("Number of Events");
  deltaT2->GetYaxis()->SetTitleSize(0.045);
  deltaT2->GetYaxis()->SetTitleOffset(1.7);
  deltaT2->GetYaxis()->SetLabelSize(0.045);
  deltaT2->GetYaxis()->SetLabelOffset(0.015);
  deltaT2->GetYaxis()->SetTitle("Number of Events");
  deltaT2->GetXaxis()->SetTitle("#Delta t [ns]");
  deltaT2->GetXaxis()->SetTitleSize(0.045);
  deltaT2->GetXaxis()->SetLabelSize(0.045);
  deltaT2->SetTitle("");
  tex = new TLatex();
  tex->SetNDC();
  tex->SetTextSize(0.050);
  tex->SetTextFont(42);
  tex->SetTextColor(kBlack);
  tex->DrawLatex(0.62, 0.80, Form("#sigma = %.0f #pm %.1f ps", 1000*fgaus2->GetParameter(2), 1000*fgaus2->GetParError(2)));
  tex->DrawLatex(0.15, 0.93, "16 GeV Electron Beam, 6 X_{0} Absorber");
  c->SaveAs("deltaT_16GeV_6X0.pdf");


  c = new TCanvas ("c","c",800, 800);
  c->SetRightMargin(0.05);
  c->SetLeftMargin(0.15);
  deltaT3->Draw();
  deltaT3->GetXaxis()->SetRangeUser(3.7,4.1);
  TF1* fgaus3 = new TF1("fgaus3","gaus", deltaT3->GetMean() - 1.5*deltaT3->GetRMS(), deltaT3->GetMean() + 1.5*deltaT3->GetRMS());
  deltaT3->Fit("fgaus3","Q","", deltaT3->GetMean() - 1.5*deltaT3->GetRMS(), deltaT3->GetMean() + 1.5*deltaT3->GetRMS());
  float res3 = fgaus3->GetParameter(2);
  float error3 = fgaus3->GetParError(2);
  deltaT3->GetYaxis()->SetTitle("Number of Events");
  deltaT3->GetYaxis()->SetTitleSize(0.045);
  deltaT3->GetYaxis()->SetTitleOffset(1.7);
  deltaT3->GetYaxis()->SetLabelSize(0.045);
  deltaT3->GetYaxis()->SetLabelOffset(0.015);
  deltaT3->GetYaxis()->SetTitle("Number of Events");
  deltaT3->GetXaxis()->SetTitle("#Delta t [ns]");
  deltaT3->GetXaxis()->SetTitleSize(0.045);
  deltaT3->GetXaxis()->SetLabelSize(0.045);
  deltaT3->SetTitle("");
  tex = new TLatex();
  tex->SetNDC();
  tex->SetTextSize(0.050);
  tex->SetTextFont(42);
  tex->SetTextColor(kBlack);
  tex->DrawLatex(0.62, 0.80, Form("#sigma = %.0f #pm %.1f ps", 1000*fgaus3->GetParameter(2), 1000*fgaus3->GetParError(2)));
  tex->DrawLatex(0.15, 0.93, "32 GeV Electron Beam, 6 X_{0} Absorber");
  c->SaveAs("deltaT_32GeV_6X0.pdf");


  return;



  float res[4]    = {1000*res0, 1000*res1, 1000*res2, 1000*res3}; 
  float charge[4] = {4., 8, 16, 32};
  float errorX[4] = {0.}; 
  float errorY[4] = {1000*error0, 1000*error1, 1000*error2, 1000*error3};

  
  c = new TCanvas ("c","c",800, 600);
  c->SetGridx();
  c->SetGridy();
  TGraphErrors* gr = new TGraphErrors( 4, charge, res, errorX, errorY );
  gr -> SetTitle("");

  gr -> SetMarkerStyle(8);
  gr -> SetLineColor(kBlue);
  gr -> SetMarkerColor(kBlue);
  gr -> SetMarkerSize(0.9);
  gr -> SetLineWidth(2);

  gr -> Draw("APL");

  gr->GetXaxis()->SetTitle("Beam Energy [GeV]");
  gr->GetYaxis()->SetTitle("Time Resolution [ps]");
  gr->GetXaxis()->SetTitleSize(0.045);
  gr->GetXaxis()->SetLabelSize(0.045);
  gr->GetXaxis()->SetTitleOffset(1.0);
  gr->GetYaxis()->SetTitleOffset(0.95);
  gr->GetYaxis()->SetTitleSize(0.05);
  gr->GetYaxis()->SetLabelSize(0.045);
  gr->GetXaxis()->SetRangeUser(0.0, 700);
  gr->GetYaxis()->SetRangeUser(0., 80.);
 
  c->SaveAs("SigmaT_vs_BeamEnergy_lin30Stamp.pdf");

}


// make the plots as a function of X0
void makeTimeResolutionVsAbsorber(){

  TFile *_file0 = TFile::Open("doc/papers/SiliconPad/data/run119_lin30_DeltaT.root");
  TFile *_file1 = TFile::Open("doc/papers/SiliconPad/data/run118_lin30_DeltaT.root");
  TFile *_file2 = TFile::Open("doc/papers/SiliconPad/data/run117_lin30_DeltaT.root");
  TFile *_file3 = TFile::Open("doc/papers/SiliconPad/data/run116_lin30_DeltaT.root");

  TH1F* deltaT0 = (TH1F*)_file0->Get("DeltaT");
  TH1F* deltaT1 = (TH1F*)_file1->Get("DeltaT");
  TH1F* deltaT2 = (TH1F*)_file2->Get("DeltaT");
  TH1F* deltaT3 = (TH1F*)_file3->Get("DeltaT");

  TCanvas *c = new TCanvas ("c","c",800, 600);
  TLatex *tex = new TLatex();
  tex->SetNDC();
  tex->SetTextSize(0.060);
  tex->SetTextFont(42);
  tex->SetTextColor(kBlack);

  c->cd();
  //gStyle->SetOptStat(0);
 

  deltaT0->Draw();
  TF1* fgaus0 = new TF1("fgaus0","gaus", deltaT0->GetMean() - 1.5*deltaT0->GetRMS(), deltaT0->GetMean() + 1.5*deltaT0->GetRMS());
  deltaT0->Fit("fgaus0","Q","", deltaT0->GetMean() - 1.5*deltaT0->GetRMS(), deltaT0->GetMean() + 1.5*deltaT0->GetRMS());
  float res0 = fgaus0->GetParameter(2);
  float error0 = fgaus0->GetParError(2);
  deltaT0->GetXaxis()->SetTitle("Time Resolution [ps]");
  deltaT0->SetTitle("");
  tex->DrawLatex(0.6, 0.80, Form("#sigma = %.0f #pm %.2f ps", 1000*fgaus0->GetParameter(2), 1000*fgaus0->GetParError(2)));
  c->SaveAs("deltaT_1X0.pdf");

  deltaT1->Draw();
  TF1* fgaus1 = new TF1("fgaus1","gaus", deltaT1->GetMean() - 1.5*deltaT1->GetRMS(), deltaT1->GetMean() + 1.5*deltaT1->GetRMS());
  deltaT1->Fit("fgaus1","Q","", deltaT1->GetMean() - 1.5*deltaT1->GetRMS(), deltaT1->GetMean() + 1.5*deltaT1->GetRMS());
  float res1 = fgaus1->GetParameter(2);
  float error1 = fgaus1->GetParError(2);
  deltaT1->GetXaxis()->SetTitle("Time Resolution [ps]");
  deltaT1->SetTitle("");
  tex->DrawLatex(0.6, 0.80, Form("#sigma = %.0f #pm %.2f ps", 1000*fgaus1->GetParameter(2), 1000*fgaus1->GetParError(2)));
  c->SaveAs("deltaT_2X0.pdf");

  deltaT2->Draw();
  TF1* fgaus2 = new TF1("fgaus2","gaus", deltaT2->GetMean() - 1.5*deltaT2->GetRMS(), deltaT2->GetMean() + 1.5*deltaT2->GetRMS());
  deltaT2->Fit("fgaus2","Q","", deltaT2->GetMean() - 1.5*deltaT2->GetRMS(), deltaT2->GetMean() + 1.5*deltaT2->GetRMS());
  float res2 = fgaus2->GetParameter(2);
  float error2 = fgaus2->GetParError(2);
  deltaT2->GetXaxis()->SetTitle("Time Resolution [ps]");
  deltaT2->SetTitle("");
  tex->DrawLatex(0.6, 0.80, Form("#sigma = %.0f #pm %.2f ps", 1000*fgaus2->GetParameter(2), 1000*fgaus2->GetParError(2)));
  c->SaveAs("deltaT_4X0.pdf");

  deltaT3->Draw();
  TF1* fgaus3 = new TF1("fgaus3","gaus", deltaT3->GetMean() - 1.5*deltaT3->GetRMS(), deltaT3->GetMean() + 1.5*deltaT3->GetRMS());
  deltaT3->Fit("fgaus3","Q","", deltaT3->GetMean() - 1.5*deltaT3->GetRMS(), deltaT3->GetMean() + 1.5*deltaT3->GetRMS());
  float res3 = fgaus3->GetParameter(2);
  float error3 = fgaus3->GetParError(2);
  deltaT3->GetXaxis()->SetTitle("Time Resolution [ps]");
  deltaT3->SetTitle("");
  tex->DrawLatex(0.6, 0.80, Form("#sigma = %.0f #pm %.2f ps", 1000*fgaus3->GetParameter(2), 1000*fgaus3->GetParError(2)));
  c->SaveAs("deltaT_6X0.pdf");


  c = new TCanvas ("c","c",800, 600);
  c->SetGridx();
  c->SetGridy();
  c->SetBottomMargin(0.12);

  float res[4]    = {1000*res0, 1000*res1, 1000*res2, 1000*res3}; 
  float charge[4] = {1., 2, 4, 6};
  float errorX[4] = {0.}; 
  float errorY[4] = {1000*error0, 1000*error1, 1000*error2, 1000*error3};

  TGraphErrors* gr = new TGraphErrors( 4, charge, res, errorX, errorY );
  gr -> SetTitle("");

  gr -> SetMarkerStyle(8);
  gr -> SetLineColor(kBlue);
  gr -> SetMarkerColor(kBlue);
  gr -> SetMarkerSize(0.9);
  gr -> SetLineWidth(2);

  gr -> Draw("APL");

  gr->GetXaxis()->SetTitle("Tungsten Absober Thickness [X_{0}]");
  gr->GetYaxis()->SetTitle("Time Resolution [ps]");
  gr->GetXaxis()->SetTitleSize(0.045);
  gr->GetXaxis()->SetLabelSize(0.045);
  gr->GetXaxis()->SetTitleOffset(1.2);
  gr->GetYaxis()->SetTitleOffset(0.95);
  gr->GetYaxis()->SetTitleSize(0.05);
  gr->GetYaxis()->SetLabelSize(0.045);
  gr->GetXaxis()->SetRangeUser(0.0, 700);
  gr->GetYaxis()->SetRangeUser(0., 150.);
 
  c->SaveAs("SigmaT_vs_X0_lin30Stamp.pdf");
}


// make the plots as a function of Silicon Voltage
void makeTimeResolutionVsVoltage(){

  TFile *_file0 = TFile::Open("doc/papers/SiliconPad/data/run111-234_lin30_DeltaT.root");
  TFile *_file1 = TFile::Open("doc/papers/SiliconPad/data/run111-1_lin30_DeltaT.root");
  TFile *_file2 = TFile::Open("doc/papers/SiliconPad/data/run111-0_lin30_DeltaT.root");
  TFile *_file3 = TFile::Open("doc/papers/SiliconPad/data/run110_lin30_DeltaT.root");

  TH1F* deltaT0 = (TH1F*)_file0->Get("DeltaT");
  TH1F* deltaT1 = (TH1F*)_file1->Get("DeltaT");
  TH1F* deltaT2 = (TH1F*)_file2->Get("DeltaT");
  TH1F* deltaT3 = (TH1F*)_file3->Get("DeltaT");

  TCanvas *c = new TCanvas ("c","c",800, 600);
  TLatex *tex = new TLatex();
  tex->SetNDC();
  tex->SetTextSize(0.060);
  tex->SetTextFont(42);
  tex->SetTextColor(kBlack);

  c->cd();
  //gStyle->SetOptStat(0);

  deltaT0->Draw();
  TF1* fgaus0 = new TF1("fgaus0","gaus", deltaT0->GetMean() - 1.0*deltaT0->GetRMS(), deltaT0->GetMean() + 1.0*deltaT0->GetRMS());
  deltaT0->Fit("fgaus0","Q","", deltaT0->GetMean() - 1.0*deltaT0->GetRMS(), deltaT0->GetMean() + 1.0*deltaT0->GetRMS());
  float res0 = fgaus0->GetParameter(2);
  float error0 = fgaus0->GetParError(2);
  deltaT0->GetXaxis()->SetTitle("Time Resolution [ps]");
  deltaT0->SetTitle("");
  tex->DrawLatex(0.6, 0.80, Form("#sigma = %.0f #pm %.2f ps", 1000*fgaus0->GetParameter(2), 1000*fgaus0->GetParError(2)));
  c->SaveAs("deltaT_200V.pdf");

  deltaT1->Draw();
  TF1* fgaus1 = new TF1("fgaus1","gaus", deltaT1->GetMean() - 1.0*deltaT1->GetRMS(), deltaT1->GetMean() + 1.0*deltaT1->GetRMS());
  deltaT1->Fit("fgaus1","Q","", deltaT1->GetMean() - 1.0*deltaT1->GetRMS(), deltaT1->GetMean() + 1.0*deltaT1->GetRMS());
  float res1 = fgaus1->GetParameter(2);
  float error1 = fgaus1->GetParError(2);
  deltaT1->GetXaxis()->SetTitle("Time Resolution [ps]");
  deltaT1->SetTitle("");
  tex->DrawLatex(0.6, 0.80, Form("#sigma = %.0f #pm %.2f ps", 1000*fgaus1->GetParameter(2), 1000*fgaus1->GetParError(2)));
  c->SaveAs("deltaT_300V.pdf");

  deltaT2->Draw();
  TF1* fgaus2 = new TF1("fgaus2","gaus", deltaT2->GetMean() - 1.0*deltaT2->GetRMS(), deltaT2->GetMean() + 1.0*deltaT2->GetRMS());
  deltaT2->Fit("fgaus2","Q","", deltaT2->GetMean() - 1.0*deltaT2->GetRMS(), deltaT2->GetMean() + 1.0*deltaT2->GetRMS());
  float res2 = fgaus2->GetParameter(2);
  float error2 = fgaus2->GetParError(2);
  deltaT2->GetXaxis()->SetTitle("Time Resolution [ps]");
  deltaT2->SetTitle("");
  tex->DrawLatex(0.6, 0.80, Form("#sigma = %.0f #pm %.2f ps", 1000*fgaus2->GetParameter(2), 1000*fgaus2->GetParError(2)));
  c->SaveAs("deltaT_400V.pdf");

  deltaT3->Draw();
  TF1* fgaus3 = new TF1("fgaus3","gaus", deltaT3->GetMean() - 1.0*deltaT3->GetRMS(), deltaT3->GetMean() + 1.0*deltaT3->GetRMS());
  deltaT3->Fit("fgaus3","Q","", deltaT3->GetMean() - 1.0*deltaT3->GetRMS(), deltaT3->GetMean() + 1.0*deltaT3->GetRMS());
  float res3 = fgaus3->GetParameter(2);
  float error3 = fgaus3->GetParError(2);
  deltaT3->GetXaxis()->SetTitle("Time Resolution [ps]");
  deltaT3->SetTitle("");
  tex->DrawLatex(0.6, 0.80, Form("#sigma = %.0f #pm %.2f ps", 1000*fgaus3->GetParameter(2), 1000*fgaus3->GetParError(2)));
  c->SaveAs("deltaT_500V.pdf");


  float res[4]    = {1000*res0, 1000*res1, 1000*res2, 1000*res3}; 
  float charge[4] = {200., 300., 400., 500.};
  float errorX[4] = {0.}; 
  float errorY[4] = {1000*error0, 1000*error1, 1000*error2, 1000*error3};

  TGraphErrors* gr = new TGraphErrors( 4, charge, res, errorX, errorY );
  gr -> SetTitle("");

  gr -> SetMarkerStyle(8);
  gr -> SetLineColor(kBlue);
  gr -> SetMarkerColor(kBlue);
  gr -> SetMarkerSize(0.9);
  gr -> SetLineWidth(2);

  gr -> Draw("APL");

  gr->GetXaxis()->SetTitle("Bias Voltage [V] ");
  gr->GetYaxis()->SetTitle("Time Resolution [ps]");
  gr->GetXaxis()->SetTitleSize(0.045);
  gr->GetXaxis()->SetLabelSize(0.045);
  gr->GetXaxis()->SetTitleOffset(1.0);
  gr->GetYaxis()->SetTitleOffset(0.95);
  gr->GetYaxis()->SetTitleSize(0.05);
  gr->GetYaxis()->SetLabelSize(0.045);
  gr->GetXaxis()->SetRangeUser(0.0, 700);
  gr->GetYaxis()->SetRangeUser(0., 80.);
 
  c->SaveAs("SigmaT_vs_DV_lin30Stamp.pdf");
}

void makePaperPlots() {

  //makeTimeResolutionVsBeamEnergy();
  makeTimeResolutionVsAbsorber();
  //makeTimeResolutionVsVoltage();

}
