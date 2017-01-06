{
//=========Macro generated from canvas: CdTeTiming2/CdTeTiming2
//=========  (Thu Nov 24 18:10:49 2016) by ROOT version5.34/07
   TCanvas *CdTeTiming2 = new TCanvas("CdTeTiming2", "CdTeTiming2",254,29,500,500);
   gStyle->SetOptFit(1);
   gStyle->SetOptStat(0);
   CdTeTiming2->Range(-1.5,1.249999,13.5,88.75);
   CdTeTiming2->SetFillColor(0);
   CdTeTiming2->SetBorderMode(0);
   CdTeTiming2->SetBorderSize(2);
   CdTeTiming2->SetFrameBorderMode(0);
   CdTeTiming2->SetFrameBorderMode(0);
   CdTeTiming2->SetLeftMargin(0.15);
   CdTeTiming2->SetBottomMargin(0.12);


   TH1F *timingvsr = new TH1F("timingvsr","timingvsr",12,0,12);
   timingvsr->SetBinContent(1,20.18473);
   timingvsr->SetBinContent(2,24.91624);
   timingvsr->SetBinContent(3,23.7178);
   timingvsr->SetBinContent(4,24.40884);
   timingvsr->SetBinContent(5,26.2355);
   timingvsr->SetBinContent(6,28.0548);
   timingvsr->SetBinContent(7,31.54168);
   timingvsr->SetBinContent(8,42.36994);
   timingvsr->SetBinContent(9,88.27947);
   timingvsr->SetBinContent(10,356.4376);
   timingvsr->SetBinContent(11,70.81522);
   timingvsr->SetBinContent(12,113.9672);
   timingvsr->SetBinError(1,1.05596);
   timingvsr->SetBinError(2,0.8719087);
   timingvsr->SetBinError(3,0.691177);
   timingvsr->SetBinError(4,0.598555);
   timingvsr->SetBinError(5,0.6738172);
   timingvsr->SetBinError(6,0.829842);
   timingvsr->SetBinError(7,1.167989);
   timingvsr->SetBinError(8,3.701801);
   timingvsr->SetBinError(9,45.97428);
   timingvsr->SetBinError(10,272.1961);
   timingvsr->SetBinError(11,51.89134);
   timingvsr->SetBinError(12,225.4484);
   timingvsr->SetMinimum(0);
   timingvsr->SetMaximum(80);
   timingvsr->SetEntries(12);
   timingvsr->SetTitle("");
   timingvsr->GetXaxis()->SetTitle("Distance to wire-bond [mm]");
   timingvsr->GetXaxis()->SetTitleSize(0.05);
   timingvsr->GetXaxis()->SetLabelSize(0.04);
   timingvsr->GetYaxis()->SetTitle("Time Resolution [ps]");
   timingvsr->GetYaxis()->SetTitleSize(0.05);
   timingvsr->GetYaxis()->SetTitleOffset(1.25);
   timingvsr->GetYaxis()->SetLabelSize(0.04);
   timingvsr->GetYaxis()->SetTitle("Time Resolution [ps]");
   timingvsr->SetLineColor(kBlack);
   timingvsr->SetLineWidth(2);
 
   timingvsr->Draw("");
   
   TPaveText *pt = new TPaveText(0.01,0.945,0.04008064,0.995,"blNDC");
   pt->SetName("title");
   pt->SetBorderSize(1);
   pt->SetFillColor(0);
   TText *text = pt->AddText("timingvsr");
   pt->Draw();
   CdTeTiming2->Modified();
   CdTeTiming2->cd();
   CdTeTiming2->SetSelected(CdTeTiming2);
}
