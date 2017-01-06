{
//=========Macro generated from canvas: Canvas_1/Canvas_1
//=========  (Thu Nov 24 17:16:53 2016) by ROOT version5.34/36
   TCanvas *Canvas_1 = new TCanvas("Canvas_1", "Canvas_1",572,176,748,579);
   Canvas_1->Range(1.375,-6.25,17.625,56.25);
   Canvas_1->SetFillColor(0);
   Canvas_1->SetBorderMode(0);
   Canvas_1->SetBorderSize(2);
   Canvas_1->SetFrameBorderMode(0);
   Canvas_1->SetFrameBorderMode(0);
   
   TH1F *timingvsx = new TH1F("timingvsx","",12,3,16);
   timingvsx->SetBinContent(1,102.8503);
   timingvsx->SetBinContent(2,27.44033);
   timingvsx->SetBinContent(3,26.71203);
   timingvsx->SetBinContent(4,24.99561);
   timingvsx->SetBinContent(5,25.28536);
   timingvsx->SetBinContent(6,24.81045);
   timingvsx->SetBinContent(7,22.61036);
   timingvsx->SetBinContent(8,22.30497);
   timingvsx->SetBinContent(9,24.55057);
   timingvsx->SetBinContent(10,27.56612);
   timingvsx->SetBinContent(11,108.7242);
   timingvsx->SetBinContent(12,352.7209);
   timingvsx->SetBinError(1,117.9207);
   timingvsx->SetBinError(2,6.207026);
   timingvsx->SetBinError(3,2.088589);
   timingvsx->SetBinError(4,1.258595);
   timingvsx->SetBinError(5,1.319373);
   timingvsx->SetBinError(6,1.096787);
   timingvsx->SetBinError(7,0.9656131);
   timingvsx->SetBinError(8,0.975975);
   timingvsx->SetBinError(9,1.378177);
   timingvsx->SetBinError(10,3.018471);
   timingvsx->SetBinError(11,121.1112);
   timingvsx->SetBinError(12,352.8225);
   timingvsx->SetMinimum(0);
   timingvsx->SetMaximum(50);
   timingvsx->SetEntries(12);
   timingvsx->SetStats(0);
   timingvsx->SetLineWidth(2);
   timingvsx->SetMarkerStyle(20);
   timingvsx->GetXaxis()->SetTitle("Vertical Beam Position [mm]");
   timingvsx->GetYaxis()->SetTitle("Time Resolution [ps]");
   timingvsx->Draw("");
   Canvas_1->Modified();
   Canvas_1->cd();
   Canvas_1->SetSelected(Canvas_1);
}
