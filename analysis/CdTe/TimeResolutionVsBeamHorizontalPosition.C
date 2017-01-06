{
//=========Macro generated from canvas: Canvas_1/Canvas_1
//=========  (Thu Nov 24 17:17:55 2016) by ROOT version5.34/36
   TCanvas *Canvas_1 = new TCanvas("Canvas_1", "Canvas_1",572,176,748,579);
   Canvas_1->Range(-9.625,-6.25,6.625,56.25);
   Canvas_1->SetFillColor(0);
   Canvas_1->SetBorderMode(0);
   Canvas_1->SetBorderSize(2);
   Canvas_1->SetFrameBorderMode(0);
   Canvas_1->SetFrameBorderMode(0);
   
   TH1F *timingvsy = new TH1F("timingvsy","",12,-8,5);
   timingvsy->SetBinContent(1,76.78795);
   timingvsy->SetBinContent(2,26.19271);
   timingvsy->SetBinContent(3,24.10236);
   timingvsy->SetBinContent(4,21.07939);
   timingvsy->SetBinContent(5,21.43319);
   timingvsy->SetBinContent(6,22.79189);
   timingvsy->SetBinContent(7,24.39891);
   timingvsy->SetBinContent(8,24.35709);
   timingvsy->SetBinContent(9,28.47322);
   timingvsy->SetBinContent(10,30.06272);
   timingvsy->SetBinContent(11,77.07889);
   timingvsy->SetBinContent(12,103.981);
   timingvsy->SetBinError(1,60.21739);
   timingvsy->SetBinError(2,3.119071);
   timingvsy->SetBinError(3,1.605079);
   timingvsy->SetBinError(4,1.219044);
   timingvsy->SetBinError(5,1.617441);
   timingvsy->SetBinError(6,1.500618);
   timingvsy->SetBinError(7,1.322545);
   timingvsy->SetBinError(8,1.587265);
   timingvsy->SetBinError(9,2.93994);
   timingvsy->SetBinError(10,4.746157);
   timingvsy->SetBinError(11,93.20974);
   timingvsy->SetBinError(12,139.1813);
   timingvsy->SetMinimum(0);
   timingvsy->SetMaximum(50);
   timingvsy->SetEntries(12);
   timingvsy->SetStats(0);
   timingvsy->SetLineWidth(2);
   timingvsy->SetMarkerStyle(20);
   timingvsy->GetXaxis()->SetTitle("Horizontal Beam Position [mm]");
   timingvsy->GetYaxis()->SetTitle("Time Resolution [ns]");
   timingvsy->Draw("");
   Canvas_1->Modified();
   Canvas_1->cd();
   Canvas_1->SetSelected(Canvas_1);
}
