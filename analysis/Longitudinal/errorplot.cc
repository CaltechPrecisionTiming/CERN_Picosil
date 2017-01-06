{
   TCanvas *c = new TCanvas("c","c",200,10,600,400);
   double x[] = {8, 16, 32};
   double y[] = {10.08, 20.9, 42.7};
   double ex[] = {0, 0, 0};
   double ey[] = {0.05, 0.10, 0.4};
   TGraphErrors* ge = new TGraphErrors(3, x, y, ex, ey);
   ge->Draw("ap");

   Xaxis = ge->GetXaxis();
   Xaxis->SetTitle("Electron Beam Energy (GeV/c^{2})");
   Xaxis->SetTitleSize(0.055);
   Xaxis->SetLabelSize(0.045);
   Xaxis->SetTitleOffset(0.8);

   Yaxis = ge->GetYaxis();
   Yaxis->SetTitle("Integrated Charge (pC)"); 
   Yaxis->SetTitleSize(0.055);
   Yaxis->SetLabelSize(0.045);
   Yaxis->SetTitleOffset(0.85);

   ge->SetTitle("");
   ge->SetLineWidth(3);

   TF1 *linfit = new TF1("linfit", "[0]*x + [1]",0,33);
   ge->Fit("linfit","QMES",0,33);
   TLatex *tex = new TLatex();
   tex->SetNDC();
   tex->SetTextSize(0.080);
   tex->SetTextFont(42);
   tex->SetTextColor(kBlack);
   if(linfit->GetParameter(1)>0) 
   tex->DrawLatex(0.15, 0.8, Form("y = %.2fx + %.2f", linfit->GetParameter(0), linfit->GetParameter(1)));
   else tex->DrawLatex(0.15, 0.8, Form("y = %.2fx - %.2f", linfit->GetParameter(0), -linfit->GetParameter(1) ));

   c->SaveAs("chargemean.pdf");
/////////////////////////////////////////////////////////////////////////////
//  Values from PicoSil combined with Landau Charge weighted equally with MCP
   
   double res[] = {21.6,14.2,11.3};
   double energy[] = {8,16,32};
   double Eres[] = {0.3,0.2,0.3};
   double Eenergy[] = {0,0,0};

   TGraphErrors* resplot = new TGraphErrors(3, energy, res, Eenergy, Eres);
   resplot->Draw("ap");
   
   Xaxis = resplot->GetXaxis();
   Xaxis->SetTitle("Electron Beam Energy (GeV/c^{2})");
   Xaxis->SetTitleSize(0.055);
   Xaxis->SetLabelSize(0.045);
   Xaxis->SetTitleOffset(0.8);

   Yaxis = resplot->GetYaxis();
   Yaxis->SetTitle("2-Layer Time Resolution (ps)");
   Yaxis->SetTitleSize(0.055);
   Yaxis->SetLabelSize(0.045);
   Yaxis->SetTitleOffset(0.85);

   resplot->SetTitle("");
   resplot->SetLineWidth(3);

   TF1 *fit = new TF1("fit", "[0]/sqrt(x) + [1]",0,33);
   resplot->Fit("fit","QMES",0,33);
   TLatex *tex = new TLatex();
   tex->SetNDC();
   tex->SetTextSize(0.080);
   tex->SetTextFont(42);
   tex->SetTextColor(kBlack);
   if(fit->GetParameter(1)>0)
   tex->DrawLatex(0.45, 0.8, Form("y = #frac{%.2f}{x} + %.2f", fit->GetParameter(0), fit->GetParameter(1)));
   else tex->DrawLatex(0.45, 0.8, Form("y = #frac{%.2f}{x} - %.2f", fit->GetParameter(0), -fit->GetParameter(1)));

   c->SaveAs("resplot.pdf");
//////////////////////////////////////////////////////////////////////////
//  Smeared values from equally-combined PicoSil weighted equally with MCP

   double res[] = {41.9,37,33.2};
   double energy[] = {8,16,32};
   double Eres[] = {0.8,0.5,0.5};
   double Eenergy[] = {0,0,0};

   TGraphErrors* smear = new TGraphErrors(3, energy, res, Eenergy, Eres);
   smear->Draw("ap");

   Xaxis = smear->GetXaxis();
   Xaxis->SetTitle("Electron Beam Energy (GeV/c^{2})");
   Xaxis->SetTitleSize(0.055);
   Xaxis->SetLabelSize(0.045);
   Xaxis->SetTitleOffset(0.8);

   Yaxis = smear->GetYaxis();
   Yaxis->SetTitle("2-Layer Smear Time Resolution (ps)");
   Yaxis->SetTitleSize(0.055);
   Yaxis->SetLabelSize(0.045);
   Yaxis->SetTitleOffset(0.85);

   smear->SetTitle("");
   smear->SetLineWidth(3);

   TF1 *fit = new TF1("fit", "[0]/sqrt(x) + [1]",0,33);
   smear->Fit("fit","QMES",0,33);
   TLatex *tex = new TLatex();
   tex->SetNDC();
   tex->SetTextSize(0.080);
   tex->SetTextFont(42);
   tex->SetTextColor(kBlack);
   if(fit->GetParameter(1)>0)
   tex->DrawLatex(0.45, 0.8, Form("y = #frac{%.2f}{x} + %.2f", fit->GetParameter(0), fit->GetParameter(1)));
   else tex->DrawLatex(0.45, 0.8, Form("y = #frac{%.2f}{x} - %.2f", fit->GetParameter(0), -fit->GetParameter(1)));

   c->SaveAs("resplotsmear.pdf");
}
