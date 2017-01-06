{
  TCanvas *c[5];
  for (int i = 0; i<5; i++){
  c[i] = new TCanvas(Form("c_%d",i),Form("canvas_%d",i),1000,600);
  c[i]->Divide(1,2,0,0);
  }
  TRandom3 *r = new TRandom3(0);
  TH1F *histChi[5];
  int nEvents[5] = {10, 100, 1000, 10000, 100000};
  gStyle->SetOptFit(1);
  gStyle->SetOptStat(1000000001);
  gStyle->SetStatW(0.13);
  gStyle->SetStatH(0.4);

  for ( int j = 0; j < 5; j++ ){
    histChi[j] = new TH1F( Form("histChi%d",j), "", 50, -5,5);
    for ( int i = 0; i < nEvents[j]; i++ )histChi[j]->Fill(r->Gaus(0,1));
    std::cout<<"\n\nSAMPLING FROM A GAUSSIAN NORMAL DISTRIBUTION: "<< nEvents[j] <<" EVENTS"<<std::endl;
    TH1F *histLikelihood = (TH1F*)histChi[j]->Clone(Form("histLikelihood_%d",j));
    c[j]->cd(1);
    std::cout<<"\nCHI-SQUARED FIT RESULTS: \n"<<std::endl;
    histChi[j]->Draw();
    histChi[j]->Fit("gaus","","");
    c[j]->cd(2);
    std::cout<<"\n\nLOG-LIKELIHOOD FIT RESULTS:\n"<<std::endl;
    histLikelihood->Draw();
    histLikelihood->Fit("gaus","MLES","");
    c[j]->SaveAs( Form("ChiSqd_vs_Likelihood_%d.pdf",j) );
  }
}
