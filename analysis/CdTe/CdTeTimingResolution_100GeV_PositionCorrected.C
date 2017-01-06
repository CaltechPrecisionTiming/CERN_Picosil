{
//=========Macro generated from canvas: CdTeTiming/CdTeTiming
//=========  (Fri Nov 25 16:40:54 2016) by ROOT version5.34/07
  TCanvas *CdTeTiming = new TCanvas("CdTeTiming", "CdTeTiming",800, 800);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   CdTeTiming->Range(-4.449655,-31.85495,-3.949793,284.547);
   CdTeTiming->SetFillColor(0);
   CdTeTiming->SetBorderMode(0);
   CdTeTiming->SetBorderSize(2);
   CdTeTiming->SetFrameBorderMode(0);
   CdTeTiming->SetFrameBorderMode(0);
   CdTeTiming->SetLeftMargin(0.15);
   CdTeTiming->SetBottomMargin(0.12);


   TH1F *h2000__1 = new TH1F("h2000__1","(time_first-(time_ref2)/1.0)*0.2-0.0096*TDCx+0.011*TDCy-0.0379*(int(time_first+0.5)-time_first)+0.000106*max_first {(time_first-(time_ref1+time_ref2)/2.0)*0.2>-5 && (time_first-(time_ref1+time_ref2)/2.0)*0.2<-3.15 && TDCx>10 && TDCx<12 && TDCy>-6.0 && TDCy<-4.0 && max_first>100 && max_first<1600}",30,-4.4,-4);
   h2000__1->SetBinContent(0,1);
   h2000__1->SetBinContent(1,1);
   h2000__1->SetBinContent(6,2);
   h2000__1->SetBinContent(7,2);
   h2000__1->SetBinContent(8,3);
   h2000__1->SetBinContent(9,5);
   h2000__1->SetBinContent(10,17);
   h2000__1->SetBinContent(11,34);
   h2000__1->SetBinContent(12,35);
   h2000__1->SetBinContent(13,81);
   h2000__1->SetBinContent(14,134);
   h2000__1->SetBinContent(15,225);
   h2000__1->SetBinContent(16,241);
   h2000__1->SetBinContent(17,217);
   h2000__1->SetBinContent(18,115);
   h2000__1->SetBinContent(19,65);
   h2000__1->SetBinContent(20,34);
   h2000__1->SetBinContent(21,16);
   h2000__1->SetBinContent(22,8);
   h2000__1->SetBinContent(23,3);
   h2000__1->SetBinContent(24,2);
   h2000__1->SetEntries(1241);
   h2000__1->SetDirectory(0);
   
   h2000__1->GetXaxis()->SetTitle("#Delta t [ns]");
   h2000__1->GetXaxis()->SetRange(3,150);
   h2000__1->GetXaxis()->SetLabelFont(42);
   h2000__1->GetXaxis()->SetTitleSize(0.05);
   h2000__1->GetXaxis()->SetTitleFont(42);
   h2000__1->GetYaxis()->SetTitle("Number of Events");
   h2000__1->GetYaxis()->SetLabelFont(42);
   h2000__1->GetYaxis()->SetTitleSize(0.05);
   h2000__1->GetYaxis()->SetTitleFont(42);
   h2000__1->GetZaxis()->SetLabelFont(42);
   h2000__1->GetZaxis()->SetLabelSize(0.035);
   h2000__1->GetZaxis()->SetTitleSize(0.035);
   h2000__1->GetZaxis()->SetTitleFont(42);
   h2000__1->GetYaxis()->SetTitleOffset(1.3);
   h2000__1->SetLineColor(kBlack);
   h2000__1->SetLineWidth(2);
   h2000__1->Draw("");

   h2000__1->Fit("gaus","","",-4.25,-4.15);
   TVirtualFitter * fitter = TVirtualFitter::GetFitter();
 
   TLatex *tex = new TLatex();
   tex->SetNDC();
   tex->SetTextSize(0.050);
   tex->SetTextFont(42);
   tex->SetTextColor(kBlack);
   tex->DrawLatex(0.20, 0.80, Form("#sigma = %.0f %s",1000*fitter->GetParameter(2),"ps"));   
   tex->DrawLatex(0.06, 0.93, "100 GeV Electrons, 6 X_{0}, Position Corrected");

   CdTeTiming->SaveAs("CdTeTimingResolution_100GeV_PositionCorrected.pdf");

   CdTeTiming->Modified();
   CdTeTiming->cd();
   CdTeTiming->SetSelected(CdTeTiming);
}
