
#include <iostream>
#include <fstream> 
#include <sstream>
#include <map>
#include <utility>
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TF1.h"
#include "TF2.h"
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
#include "TRandom.h"
#include "TRandom3.h"



//*************************************************************************************************
//Normalize Hist
//*************************************************************************************************
TH1F* NormalizeHist(TH1F *originalHist) {
  TH1F* hist = (TH1F*)originalHist->Clone((string(originalHist->GetName())+"_normalized").c_str());
  Double_t norm = 0;
  hist->SetTitle("");
  for (UInt_t b=0; int(b)<hist->GetXaxis()->GetNbins()+2; ++b) {
    norm += hist->GetBinContent(b);
  }
  for (UInt_t b=0; int(b)<hist->GetXaxis()->GetNbins()+2; ++b) {
    hist->SetBinContent(b,hist->GetBinContent(b) / norm);
    hist->SetBinError(b,hist->GetBinError(b) / norm);
  }

  return hist;
}



struct Pixel
{
  int   index;
  float charge;
  float time;
  float amp;
  float time_smear;
};

void DoMultiChannelStudy( string filename , string outputFilename) {
  // TFile *inputfile = TFile::Open("t1065-jun-2016-65To80.root","READ"); //lead 6X0, 3cm away, 32GeV
  // TFile *inputfile = TFile::Open("t1065-jun-2016-84To89.root","READ"); //lead 6X0, 0cm away, 32GeV
  // TFile *inputfile = TFile::Open("t1065-jun-2016-91.dat-full.root","READ"); //lead 6X0, 0cm away, 32GeV, after turning off pixel telescope
  // TFile *inputfile = TFile::Open("t1065-jun-2016-90.dat-full.root","READ"); //lead 6X0, 3.4cm away, 32GeV, before turning off pixel telescope
  // TFile *inputfile = TFile::Open("t1065-jun-2016-55To59.root","READ"); //tungsten 6X0, 3.4cm away, 32 GeV
  // TFile *inputfile = TFile::Open("t1065-jun-2016-34To35.root","READ"); //tungsten 6X0, 0cm away, 8 GeV
  //TFile *inputfile = TFile::Open("t1065-jun-2016-32To33.root","READ"); //tungsten 2X0, 0cm away, 32 GeV
  //TFile *inputfile = TFile::Open("t1065-jun-2016-81.dat-full.root","READ"); //lead 6X0, 0cm away, 32GeV, before turning off pixel telescope
  //TFile *inputfile = TFile::Open("t1065-jun-2016-94.dat-full.root","READ"); //lead 6X0, 1cm away, 32GeV, before turning off pixel telescope
  // TFile *inputfile = TFile::Open("t1065-jun-2016-59.dat-full.root","READ"); //tungsten 6X0, 3.4cm away, 32 GeV
  TFile *inputfile = TFile::Open(filename.c_str(),"READ"); //tungsten 6X0, 3.4cm away, 32 GeV

  TTree *tree = (TTree*)inputfile->Get("pulse");

  // get the variables from the ntuple
  float amp[36];
  float integral[36];
  float gauspeak[36];
  float linearTime45[36];

  int event;

  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("gauspeak",1);
  tree->SetBranchStatus("linearTime45",1);
  tree->SetBranchStatus("amp",1);
  tree->SetBranchStatus("int",1);

  tree->SetBranchStatus("event",1);
  
  tree->SetBranchAddress("gauspeak",gauspeak);
  tree->SetBranchAddress("linearTime45",linearTime45);
  tree->SetBranchAddress("amp",amp);
  tree->SetBranchAddress("int",integral);

  tree->SetBranchAddress("event",&event);

  TH1F *histTotalCharge;
  TH1F *histTOFCenter;
  TH1F *histTOFFlatAvg;
  TH1F *histTOFChargeWeightedAvg;
  TH1F *histTOFRingOneFlatAvg;
  TH1F *histTOFRingOneChargeWeightedAvg;
  TH1F *histNChannelRingOne;
  TH1F *histChargeRingOne;
  TH1F *histChargeCenterOverTotalCharge;
  TH1F *histChargeRingOneOverTotalCharge;
  TH1F *histDeltaTCombined;

  TH1F *histDeltaT[7];
  for(int j=0; j<7; j++) histDeltaT[j]= new TH1F(Form("histDeltaT_%d",j),"; Time [ns];Number of Events", 250, 4.4, 5.6);

  TH1F *histDeltaT_smear[7];
  for(int j=0; j<7; j++) histDeltaT_smear[j]= new TH1F(Form("histDeltaT_smear_%d",j),"; Time [ns];Number of Events", 250, 4.4, 5.6);
  
  histTotalCharge = new TH1F("histTotalCharge","; Time [ns];Number of Events", 200, 0,100);
  histTOFCenter = new TH1F("histTOFCenter","; Time [ns];Number of Events", 200, -6,-4);
  histTOFFlatAvg = new TH1F("histTOFFlatAvg","; Time [ns];Number of Events", 200, -6,-4);
  histTOFChargeWeightedAvg = new TH1F("histTOFChargeWeightedAvg","; Time [ns];Number of Events", 200, -6,-4);
  histTOFRingOneFlatAvg = new TH1F("histTOFRingOneFlatAvg","; Time [ns];Number of Events", 200, -10,10);
  histTOFRingOneChargeWeightedAvg = new TH1F("histTOFRingOneChargeWeightedAvg","; Time [ns];Number of Events", 200, -10,10);
  histNChannelRingOne = new TH1F("histNChannelRingOne","; Time [ns];Number of Events", 7, -0.5,6.5);
  histChargeRingOne = new TH1F("histChargeRingOne","; Time [ns];Number of Events", 200, 0,100);
  histChargeCenterOverTotalCharge = new TH1F("histChargeCenterOverTotalCharge","; Time [ns];Number of Events", 50, 0,1);
  histChargeRingOneOverTotalCharge = new TH1F("histChargeRingOneOverTotalCharge","; Time [ns];Number of Events", 50, 0,1);
  histDeltaTCombined= new TH1F("histDeltaTCombined","; Time [ns];Number of Events", 100, 4,5.5);

  //read all entries and fill the histograms
  Long64_t nentries = tree->GetEntries();

  // to generate the seed for the TRandom3
  srand(time(NULL));
  int seed = rand();

  // to set the smear time (this is in ns, so 0.05 means a smearing of 50 ps)
  float smear = 0.5;

  // to set the range for the histograms
  // 50 ps smearing
  //float width = 0.5;
  // 500 ps smearing
  float width = 4;

  // The same seed will be used in the first and second event loops to generate the same random sequence with TRandom.
  TRandom3 *r = new TRandom3(seed); 


  std::cout<<"Number of events in Sample: "<<nentries<<std::endl;  
  for (Long64_t iEntry=0;iEntry<nentries;iEntry++) {
    if (iEntry %1000 == 0) cout << "Processing Event " << iEntry << "\n";
    tree->GetEntry(iEntry);    
  
    float photekTimeGauss = gauspeak[0];
    float photekAmp = amp[0];
    float photekCharge = integral[0];

    float centerAmp = amp[1];
    float centerCharge = integral[1];
    float centerTime = linearTime45[1];

    float RingOneCharge = 0;
    int NumberOfChannelsInRingOne = 0;
    float RingOneTimeFlatAvg = 0.0;
    float RingOneTimeChargeWeightedAvg = 0.0;
    float TotalCharge = 0;
    int NumberOfChannels = 0;
    float TimeFlatAvg = 0;
    float TimeChargeWeightedAvg = 0;

    float DeltaT[7] = {-99.};
    float DeltaT_smear[7] = {-99.};
    
    
    //require photek (for electron selection)
    //if( !(photekAmp > 0.05 && photekCharge > 2) ) continue;
    
    //require signal in the central pixel
    //if( !(centerAmp > 0.03 && centerCharge > 2) ) continue;


    // 8GeV cuts, 1 mm
    //if( !(photekAmp > 0.015 && photekCharge > 0.4 ) ) continue;
    //require signal in the central pixel
    //if( !(centerCharge > 2.5 && centerAmp > 0.02 ) ) continue;

    // 8GeV cuts, 32 mm
    //if( !(photekAmp > 0.04 && photekCharge > 1 ) ) continue;
    //require signal in the central pixel
    //if( !(centerCharge > 1 && centerAmp > 0.01 ) ) continue;

    // 16GeV cuts, 1 mm
    //if( !(photekAmp > 0.03 && photekCharge > 0.8 ) ) continue;
    //require signal in the central pixel
    //if( !(centerCharge > 6 && centerAmp > 0.07 ) ) continue;

    // 32GeV cuts, 1 mm
    if( !(photekAmp > 0.1 && photekCharge > 2 ) ) continue;
    //require signal in the central pixel
    if( !(centerCharge > 11 && centerAmp > 0.15 ) ) continue;

    // 32GeV cuts, 10 mm
    //if( !(photekAmp > 0.1 && photekCharge > 2 ) ) continue;
    //require signal in the central pixel
    //if( !(centerCharge > 8 && centerAmp > 0.1 ) ) continue;

    // 32GeV cuts, 32 mm
    //if( !(photekAmp > 0.09 && photekCharge > 2 ) ) continue;
    //require signal in the central pixel
    //if( !(centerCharge > 3 && centerAmp > 0.05 ) ) continue;

    // 32GeV cuts, 75 mm
    //if( !(photekAmp > 0.09 && photekCharge > 2 ) ) continue;
    //require signal in the central pixel
    //if( !(centerCharge > 2 && centerAmp > 0.03 ) ) continue;
    

    for ( int j = 1; j <= 7; j++)
      {
        float smear_number = r->Gaus(0, smear);
        if ( iEntry > 1000 && iEntry < 1015 )
        {
          //cout << "smear number " << smear_number << ":\n" << endl;
        }
        if ( amp[j] > 0.01 && integral[j] > 1 ) 
        {
          DeltaT[j-1] = gauspeak[0] - linearTime45[j];
          DeltaT_smear[j-1] = gauspeak[0] - smear_number - linearTime45[j];
          histDeltaT[j-1]->Fill(DeltaT[j-1]);
          histDeltaT_smear[j-1]->Fill(DeltaT_smear[j-1]);
          NumberOfChannels++;
          TotalCharge += integral[j];
          TimeFlatAvg += linearTime45[j];
          TimeChargeWeightedAvg += linearTime45[j]*integral[j];

          if (j > 1) 
          {
            //cout << j << " : " << linearTime45[j] << " " << integral[j] << " : " << RingOneTimeFlatAvg << " " << RingOneTimeChargeWeightedAvg << " " << TimeFlatAvg << " " << TimeChargeWeightedAvg << "\n";
            NumberOfChannelsInRingOne++;
            RingOneCharge += integral[j];
            RingOneTimeFlatAvg += linearTime45[j];
            RingOneTimeChargeWeightedAvg += linearTime45[j]*integral[j];
          }
        }
      }

    
    /*for(int jj=0; jj<7; jj++) 
    {
      histDeltaT[jj]->Fill(DeltaT[jj]);
    }*/
  
    // do the time averaging
    float DeltaTCombined = -99.;

    
    //if(amp[0]>0.03 && integral[0]>3. && integral[1]>3. && integral[5]>3. && integral[4]>3. && integral[3]>3)
    //DeltaTCombined = gauspeak[0]-(linearTime45[1] + linearTime45[4] + linearTime45[5]+ linearTime45[3])/4.;
    
    histDeltaTCombined->Fill(DeltaTCombined);    
    
    RingOneTimeFlatAvg = RingOneTimeFlatAvg / NumberOfChannelsInRingOne;
    RingOneTimeChargeWeightedAvg = RingOneTimeChargeWeightedAvg / RingOneCharge;
    TimeFlatAvg = TimeFlatAvg / NumberOfChannels;
    TimeChargeWeightedAvg = TimeChargeWeightedAvg / TotalCharge;

    if (NumberOfChannelsInRingOne==3) 
    {
      histTOFRingOneFlatAvg->Fill( RingOneTimeFlatAvg - photekTimeGauss);
      histTOFRingOneChargeWeightedAvg->Fill( RingOneTimeChargeWeightedAvg - photekTimeGauss);      
    }

    // Fill histograms
    histTotalCharge->Fill(TotalCharge);
    histNChannelRingOne->Fill(NumberOfChannelsInRingOne);
    histChargeRingOne->Fill(RingOneCharge);
    histChargeCenterOverTotalCharge->Fill(centerCharge / (centerCharge + RingOneCharge));
    histChargeRingOneOverTotalCharge->Fill(RingOneCharge / (centerCharge + RingOneCharge));

    histTOFCenter->Fill(centerTime - photekTimeGauss);
    histTOFFlatAvg->Fill( TimeFlatAvg - photekTimeGauss);
    histTOFChargeWeightedAvg->Fill( TimeChargeWeightedAvg - photekTimeGauss);
    
  }

  float meanT[7];
  for ( int i = 0; i < 7; i++ )
  {
    meanT[i] = histDeltaT[i]->GetMean();
    std::cout << "MEAN " << i << "-->" << meanT[i] << std::endl;
  }

  float meanT_smear[7];
  for ( int i = 0; i < 7; i++ )
  {
    meanT_smear[i] = histDeltaT_smear[i]->GetMean();
    std::cout << "MEAN " << i << "-->" << meanT_smear[i] << std::endl;
  }


  TH1F *histTOFCenter_C;
  TH1F *histTOFFlatAvg_C;
  TH1F *histTOFChargeWeightedAvg_C;
  TH1F *histTOFRingOneFlatAvg_C;
  TH1F *histTOFRingOneChargeWeightedAvg_C;
  TH1F *histDeltaTCombined_C;

  // this is the histogram that will have the gaussian fit to it
  TH1F *histDeltaT_C[7];
  for(int j=0; j < 7; j++) histDeltaT_C[j]= new TH1F(Form("histDeltaT_C_%d",j),"; #Deltat (ns) ; Entries / (0.01 ns)", 80, -0.4, 0.4);

  TH1F *histDeltaTshifted_C[7];
  for(int j=0; j < 7; j++) histDeltaTshifted_C[j]= new TH1F(Form("histDeltaTshifted_C_%d",j),"; #Deltat (ns) ; Entries / (0.01 ns)", 80, -0.4, 0.4);

  // TOF histograms of each pixel with the smear added
  TH1F *histDeltaTshifted_smear_C[7];
  for(int j=0; j < 7; j++) histDeltaTshifted_smear_C[j]= new TH1F(Form("histDeltaTshifted_smear_C_%d",j),"; #Deltat (ns) ; Entries / (0.01 ns)", 80, -1 * width, width);
  
  // TOF histogram of each pixel with the smear added, for the code that only uses the ring
  TH1F *histDeltaTshifted_smear2_C[6];
  for(int j=0; j < 6; j ++) histDeltaTshifted_smear2_C[j]= new TH1F(Form("histDeltaTshifted_smear2_C_%d",j),"; #Deltat (ns) ; Entries / (0.01 ns)", 80, -1 * width, width);
  
  histTOFCenter_C = new TH1F("histTOFCenter_C","; Time [ns];Number of Events", 200, -6,-4);

  // TOF of largest pixels, by charge weighting. No fit is done on this histogram.
  TH1F *histTOF_largest[7];
  for(int j=0; j < 7; j++) histTOF_largest[j]= new TH1F(Form("histTOF_largest_%d",j),"; #Deltat (ns) ; Entries / (0.01 ns)", 80, -0.3, 0.3);

  // TOF of largest pixels, by charge weighting. The Gaussian fit will be done to this one.
  TH1F *histTOF_largest_C[7];
  for(int j=0; j < 7; j++) histTOF_largest_C[j]= new TH1F(Form("histTOF_largest_C_%d",j),"; #Deltat (ns) ; Entries / (0.01 ns)", 80, -0.3, 0.3);

  // TOF of largest pixels, by equal weighting. The Gaussian fit will be done to this one.
  TH1F *histTOF_largest_equal_C[7];
  for(int j=0; j < 7; j++) histTOF_largest_equal_C[j]= new TH1F(Form("histTOF_largest_equal_C_%d",j),"; #Deltat (ns) ; Entries / (0.01 ns)", 80, -0.3, 0.3);

  // TOF of largest pixels, with the smear added. For this, the times will be added with charge weighting.
  TH1F *histTOF_largest_smear_C[7];
  for(int j=0; j < 7; j++) histTOF_largest_smear_C[j]= new TH1F(Form("histTOF_largest_smear_C_%d",j),"; #Deltat (ns) ; Entries / (0.01 ns)", 80, -1 * width, width);

  // TOF of largest pixels, with the smear added. For these, the times of each pixel will be added with equal weighting
  TH1F *histTOF_largest_smear_equal_C[7];
  for(int j=0; j < 7; j++) histTOF_largest_smear_equal_C[j]= new TH1F(Form("histTOF_largest_smear_equal_C_%d",j),"; #Deltat (ns) ; Entries / (0.01 ns)", 80, -1 * width, width);

  TH1F* histMaxIndex_C = new TH1F("histMaxIndex","; Index;Number of Events", 7, 0.5,7.5);
  TH1F* histPixelsCombined_C = new TH1F("histPixelsCombined","; Number of Pixels;Number of Events", 7, 0.5,7.5);

  TH1F* histEnergyRatio_C = new TH1F("histEnergyRatio","; Energy center / energy 7;Number of Events", 50, -0.1 , 1);
  TH1F* histEnergyRatioPhotekCut_C = new TH1F("histEnergyRatioPhotekCut","; Energy center / energy 7;Number of Events", 50, -0.1 , 1);
  TH1F* histEnergyRatioEventCut_C = new TH1F("histEnergyRatioEventCut","; Energy center / energy 7;Number of Events", 50, -0.1 , 1);
  TH1F* histEnergyRatioAllCuts_C = new TH1F("histEnergyRatioAllCuts","; Energy center / energy 7;Number of Events", 50, -0.1 , 1);

  TH2F* histEnergyPhotekAmp_C = new TH2F("histEnergyVsPhotekAmp","; Amplitude (V);Energy Ratio ;Number of Events", 100, -0.1, 0.6, 100, 0 , 1 );
  TH2F* histEnergyPhotekAmpPhotekCut_C = new TH2F("histEnergyVsPhotekAmpPhotekCut","; Amplitude (V);Energy Ratio ;Number of Events", 100, -0.1, 0.6, 100, 0 , 1 );
  TH2F* histEnergyPhotekAmpEventCut_C = new TH2F("histEnergyVsPhotekAmpEventCut","; Amplitude (V);Energy Ratio ;Number of Events", 100, -0.1, 0.6, 100, 0 , 1 );
  TH2F* histEnergyPhotekAmpAllCuts_C = new TH2F("histEnergyVsPhotekAmpAllCuts","; Amplitude (V);Energy Ratio ;Number of Events", 100, -0.1, 0.6, 100, 0 , 1 );

  TH1F* histTotalChargeContained_C = new TH1F("histTotalChargeContained","; Charge (pC);Number of Events", 100, -10 , 150);
  TH1F* histChargeCenterContained_C = new TH1F("histChargeCenterContained","; Charge (pC);Number of Events", 100, -10, 150);

  TH1F *histChargeContained[7];
  for(int j=0; j < 7; j++) histChargeContained[j]= new TH1F(Form("histChargeContained_%d",j),"; Charge (pC);Number of Events", 100, -10, 150);

  TH2F* histEnergyDeltaT2D_C = new TH2F("histEnergyDeltaT2D","; Energy (pC); #Deltat (ns); Number of Events", 100, 0, 100, 100, -0.15, 0.15);
  TH2F* histEnergyDeltaT2D_adj_C = new TH2F("histEnergyDeltaT2D_adj","; Energy (pC); #Deltat (ns); Number of Events", 100, 0, 100, 100, -0.15, 0.15);

  TH1F* histTOFEnergyCut_C = new TH1F("histTOFEnergyCut","; Time (ns);Number of Events", 80, -0.3, 0.3);
  TH1F* histTOFEnergyCutCenter_C = new TH1F("histTOFEnergyCutCenter","; Time (ns);Number of Events", 80, -0.3, 0.3);
  
  TH2F* histDeltaTCharge_C = new TH2F("histDeltaTCharge","; Charge (pC);Time (ns);Number of Events", 300, -5, 100, 300, -0.3 , 0.3 );
  TH2F* histDeltaTChargeCorr_C = new TH2F("histDeltaTChargeCorr","; Charge (pC);Time (ns);Number of Events", 300, -5, 100, 300, -0.3 , 0.3 );
  TH1F* histTOF_corr_C = new TH1F("histTOF_corr","; Time (ns); Number of Events", 300, -0.3, 0.3);
  TH1F* histTOF_corr_fit_C = new TH1F("histTOF_corr_fit","; Time (ns); Number of Events", 300, -0.3, 0.3);

  // histograms to plot the distibution of the variable being cut on. Cuts change based on the beam energy.
  TH1F* histPhotekAmpCut_C = new TH1F("histPhotekAmpCut","; Photek Amplitude (V);Number of Events", 300, -0.5 , 0.5);
  TH1F* histPhotekChargeCut_C = new TH1F("histPhotekChargeCut","; Photek Charge (pC);Number of Events", 300, -0.5 , 20);
  TH1F* histCenterChargeCut_C = new TH1F("histCenterChargeCut","; Center Charge (pC);Number of Events", 300, -0.5 , 100);
  TH1F* histCenterAmpCut_C = new TH1F("histCenterAmplitudeCut","; Center Amplitude (V);Number of Events", 300, -0.5 , 1);

  histTOFFlatAvg_C = new TH1F("histTOFFlatAvg_C","; Time [ns];Number of Events", 200, -6,-4);
  histTOFChargeWeightedAvg_C = new TH1F("histTOFChargeWeightedAvg_C","; Time [ns];Number of Events", 200, -1,-1);
  histTOFRingOneFlatAvg_C = new TH1F("histTOFRingOneFlatAvg_C","; Time [ns];Number of Events", 200, -10,10);
  histTOFRingOneChargeWeightedAvg_C = new TH1F("histTOFRingOneChargeWeightedAvg_C","; Time [ns];Number of Events", 200, -10,10);
  histDeltaTCombined_C = new TH1F("histDeltaTCombined_C","; Time [ns];Number of Events", 100, 4,5.5);

  TH1F* histRandomNumberTest_C = new TH1F("histRandomNumberTest_C","; number; number of events", 200, -0.5, 0.5);
  TH1F* histLinearTime_C = new TH1F("histLinearTime_C","; time; number of events", 200, -0, 150);


  // Delete memory allocated to r from first for loop. Then use same seed for the TRandom so the same random sequence is used in both event loops.
  delete r;
  r = new TRandom3(seed);

  std::cout<<"Number of events in Sample: "<<nentries<<std::endl; 
  for (Long64_t iEntry=0;iEntry<nentries;iEntry++)
  {
    if (iEntry %1000 == 0) cout << "Processing Event " << iEntry << "\n";
    tree->GetEntry(iEntry);  

    float photekTimeGauss = gauspeak[0];
    float photekAmp = amp[0];
    float photekCharge = integral[0];

    float totalEnergy = 0;
    float Energy = 0;
    float percent = 0;

    float largestIndex = 0;
    float largestIndexIntegral = 0;

    float ratio = 0;
    float ratio1 = 0;
    float ratio2 = 0;
    float ratio3 = 0;

    float centerAmp = amp[1];
    float centerCharge = integral[1];
    float centerTime = linearTime45[1];

    // fill histograms of variable being cut on, plots full range
    histPhotekAmpCut_C->Fill( photekAmp );
    histPhotekChargeCut_C->Fill( photekCharge );
    histCenterAmpCut_C->Fill( centerAmp );
    histCenterChargeCut_C->Fill( centerCharge );


    std::vector< Pixel > vect1;
    Pixel pixel1;
    for ( int j = 1; j <= 7; j++)
    {
      pixel1.index = j;
      if ( j == 1 )
      {
        pixel1.charge = 2*integral[j]; // 2* to account for the 6db attenuation
      }
      else pixel1.charge = integral[j];
      vect1.push_back( pixel1 ) ;
    }

    // energy ratio histogram with no cuts applied, goes through every event before event cuts are applied
    float energy_center1 = vect1[0].charge;
    float energy_total1 = vect1[0].charge + vect1[1].charge + vect1[2].charge + vect1[3].charge + vect1[4].charge + vect1[5].charge + vect1[6].charge;
    ratio1 = energy_center1 / energy_total1;

    histEnergyRatio_C->Fill( ratio1 );
    histEnergyPhotekAmp_C->Fill( photekAmp, ratio1 );

    if ( !(photekAmp > 0.1 && photekCharge > 2 )) continue;

    float energy_center3 = vect1[0].charge;
    float energy_total3 = vect1[0].charge + vect1[1].charge + vect1[2].charge + vect1[3].charge + vect1[4].charge + vect1[5].charge + vect1[6].charge;
    ratio3 = energy_center3 / energy_total3;

    if ( energy_center3 > 0 )
    {
      histEnergyRatioPhotekCut_C->Fill( ratio3 );
      histEnergyPhotekAmpPhotekCut_C->Fill( photekAmp, ratio3 );
    }

    // require photek amp and charge to be above cuts (to select electron events and not background), and require central pixel to be above amp and charge cuts

    // 8GeV cuts, 1 mm 
    //if( !(photekAmp > 0.015 && photekCharge > 0.4 ) ) continue;
    //require signal in the central pixel
    //if( !(centerCharge > 2.5 && centerAmp > 0.02 ) ) continue;

    // 8GeV cuts, 32 mm
    //if( !(photekAmp > 0.04 && photekCharge > 1 ) ) continue;
    //require signal in the central pixel
    //if( !(centerCharge > 1 && centerAmp > 0.01 ) ) continue;

    // 16GeV cuts, 1 mm
    //if( !(photekAmp > 0.03 && photekCharge > 0.8 ) ) continue;
    //require signal in the central pixel
    //if( !(centerCharge > 6 && centerAmp > 0.07 ) ) continue;

    // 32GeV cuts, 1 mm
    if( !(photekAmp > 0.1 && photekCharge > 2 ) ) continue;
    //require signal in the central pixel
    if( !(centerCharge > 11 && centerAmp > 0.15 ) ) continue;

    // 32GeV cuts, 10 mm
    //if( !(photekAmp > 0.1 && photekCharge > 2 ) ) continue;
    //require signal in the central pixel
    //if( !(centerCharge > 8 && centerAmp > 0.1 ) ) continue;

    // 32GeV cuts, 32 mm
    //if( !(photekAmp > 0.09 && photekCharge > 2 ) ) continue;
    //require signal in the central pixel
    //if( !(centerCharge > 3 && centerAmp > 0.05 ) ) continue;

    // 32GeV cuts, 75 mm
    //if( !(photekAmp > 0.09 && photekCharge > 2 ) ) continue;
    //require signal in the central pixel
    //if( !(centerCharge > 2 && centerAmp > 0.03 ) ) continue;

      



    // Generates random numbers from the Gaussian distribution with mean 0 and sigma 0.05 (50 picoseconds) to perform the time smearing for each pixel
    // float smear_number[7];
    // for ( int i = 0; i <= 6; i++)
    // {
    //   Double_t random = r->Gaus(0, 0.05);
    //   histRandomNumberTest_C->Fill( random );
    // }





    std::vector< Pixel > vect;
    Pixel pixel;
    for ( int j = 1; j <= 7; j++)
    {
      // smear number is a random number choosen from a Gaussian distribution with mean 0 and sigma 0.05 (50 ps) to perform the time smearing for each pixel
      float smear_number = r->Gaus(0, smear);

      if ( amp[j] > 0.01 && integral[j] > 1 )
      {
        histDeltaT_C[j-1]->Fill(gauspeak[0] - linearTime45[j] - meanT[j-1]);
        histDeltaTshifted_C[j-1]->Fill(gauspeak[0] - linearTime45[j] - meanT[j-1]);
        histDeltaTshifted_smear_C[j-1]->Fill(gauspeak[0] + smear_number - linearTime45[j] - meanT_smear[j-1]);
      } 
      pixel.index = j;
      if ( j == 1 )
      {
        pixel.charge = 2*integral[j]; // 2* to account for the 6db attenuation
      }
      else pixel.charge = integral[j];
      pixel.time = gauspeak[0] - (linearTime45[j] + meanT[j-1]);
      pixel.time_smear = gauspeak[0] + smear_number - linearTime45[j] - meanT[j-1];
      pixel.amp = amp[j];
      vect.push_back( pixel ) ;
    }

    histLinearTime_C->Fill( linearTime45[1] );

    auto sortPixel = []( Pixel a, Pixel b ) {return a.charge > b.charge ?  true : false;};
    //std::sort( vect.begin(), vect.end(), sortPixel );

    float average;
    float average_equal;
    float average_smear;
    float average_smear_equal;
    float average_corrected;
    float offset_correction;
    float offset_slope;
    float offset_quadratic;
    float offset_cubic;

    // this will find which pixel has the maximum charge and plot the results in the MaxIndex histogram. largestIndex and largestIndexIntegral start as 0
    // channels are 1 indexed for this
    for ( int j = 1; j <= 7; j++)
    {
      if ( integral[j] > largestIndexIntegral )
      {
        largestIndexIntegral = integral[j];
        largestIndex = j;
      }
    }
    histMaxIndex_C->Fill( largestIndex );


    float energy_center = vect[0].charge;
    float energy_total = vect[0].charge + vect[1].charge + vect[2].charge + vect[3].charge + vect[4].charge + vect[5].charge + vect[6].charge;
    float time_center = vect[0].time;
    ratio = energy_center / energy_total;

    histEnergyRatioEventCut_C->Fill( ratio );
    histEnergyPhotekAmpEventCut_C->Fill( photekAmp, ratio );


    // Sort the vect and plot the combination of the largest energy pixels
    std::sort( vect.begin(), vect.end(), sortPixel );

    // These values will be used to determine if the pixel passes the required cuts
    float energy1 = vect[0].charge;           // 6db attenuator on the center pixel already accounted for
    float energy2 = vect[1].charge;
    float energy3 = vect[2].charge;
    float energy4 = vect[3].charge;
    float energy5 = vect[4].charge;
    float energy6 = vect[5].charge;
    float energy7 = vect[6].charge;

    float amp1 = vect[0].amp;
    float amp2 = vect[1].amp;
    float amp3 = vect[2].amp;
    float amp4 = vect[3].amp;
    float amp5 = vect[4].amp;
    float amp6 = vect[5].amp;
    float amp7 = vect[6].amp;

    float weight1 = 0;    // If the pixel passes the cuts, these weights will be changed to the charge of the pixel (and time will be charge weighted) 
    float weight2 = 0;
    float weight3 = 0;
    float weight4 = 0;
    float weight5 = 0;
    float weight6 = 0;
    float weight7 = 0;     

    float count1 = 0;     // To count how many pixels are combined total. Also, if the pixel passes the cuts, this will be changed to 1 and will be used for the equal pixel weighting.
    float count2 = 0;
    float count3 = 0;
    float count4 = 0;
    float count5 = 0;
    float count6 = 0;
    float count7 = 0;
    float Pixels = 0;

    // This makes the TOF histogram with pixels with the largest charges, weighting the pixel time by the pixel charge (if pixel passes cuts). This is done event by event, so the TOF_largest_x will not necissarily have x pixels combined in every event.
    // Also makes histograms of the total energy contained in the pixels used for the TOF histogram.
    // 6db attenuator on the center pixel is already accounted for.
    // For an additional pixel to be added, it must have a charge > 1 pC and an amplitude > 0.01 V. If the pixel does not pass these cuts, then it is added with a weight of 0.
        
    if ( energy1 > 1 && amp1 > 0.01 )
    {
      weight1 = energy1;
      count1 = 1;
    }

    Energy = weight1;
    average = vect[0].time;
    average_smear = vect[0].time_smear;
    histTOF_largest[0]->Fill( average );
    histTOF_largest_C[0]->Fill( average );
    histTOF_largest_equal_C[0]->Fill( average );
    histTOF_largest_smear_C[0]->Fill( average_smear );
    histTOF_largest_smear_equal_C[0]->Fill( average_smear );
    histChargeContained[0]->Fill( Energy );

    if ( energy2 > 1 && amp2 > 0.01 )
    {
      weight2 = energy2;
      count2 = 1;
    }

    Energy = weight1 + weight2;
    average = ( weight1*vect[0].time + weight2*vect[1].time )/( Energy );
    average_equal = ( count1*vect[0].time + count2*vect[1].time )/( count1 + count2 );
    average_smear = ( weight1*vect[0].time_smear + weight2*vect[1].time_smear )/( Energy );
    average_smear_equal = ( count1*vect[0].time_smear + count2*vect[1].time_smear )/ (count1 + count2 );
    histTOF_largest[1]->Fill( average );
    histTOF_largest_C[1]->Fill( average ); 
    histTOF_largest_equal_C[1]->Fill ( average_equal );
    histTOF_largest_smear_C[1]->Fill( average_smear );
    histTOF_largest_smear_equal_C[1]->Fill( average_smear_equal );
    histChargeContained[1]->Fill( Energy );

    if ( energy3 > 1 && amp3 > 0.01 )
    {
      weight3 = energy3;
      count3 = 1;
    }
        
    Energy = weight1 + weight2 + weight3;
    average = ( weight1*vect[0].time + weight2*vect[1].time + weight3*vect[2].time )/( Energy );
    average_equal = ( count1*vect[0].time + count2*vect[1].time + count3*vect[2].time )/( count1 + count2 + count3 );
    average_smear = ( weight1*vect[0].time_smear + weight2*vect[1].time_smear + weight3*vect[2].time_smear )/( Energy );
    average_smear_equal = ( count1*vect[0].time_smear + count2*vect[1].time_smear + count3*vect[2].time_smear )/( count1 + count2 + count3 );
    histTOF_largest[2]->Fill( average );
    histTOF_largest_C[2]->Fill( average );
    histTOF_largest_equal_C[2]->Fill( average_equal );
    histTOF_largest_smear_C[2]->Fill( average_smear );
    histTOF_largest_smear_equal_C[2]->Fill( average_smear_equal );
    histChargeContained[2]->Fill( Energy );

    if ( energy4 > 1 && amp4 > 0.01 )
    {
      weight4 = energy4;
      count4 = 1;
    }
        
    Energy = weight1 + weight2 + weight3 + weight4;
    average = ( weight1*vect[0].time + weight2*vect[1].time + weight3*vect[2].time + weight4*vect[3].time )/( Energy );
    average_equal = ( count1*vect[0].time + count2*vect[1].time + count3*vect[2].time + count4*vect[3].time )/( count1 + count2 + count3 + count4 );
    average_smear = ( weight1*vect[0].time_smear + weight2*vect[1].time_smear + weight3*vect[2].time_smear + weight4*vect[3].time_smear )/( Energy );
    average_smear_equal = ( count1*vect[0].time_smear + count2*vect[1].time_smear + count3*vect[2].time_smear + count4*vect[3].time_smear )/( count1 + count2 + count3 + count4 );
    histTOF_largest[3]->Fill( average );
    histTOF_largest_C[3]->Fill( average );
    histTOF_largest_equal_C[3]->Fill( average_equal );
    histTOF_largest_smear_C[3]->Fill( average_smear ); 
    histTOF_largest_smear_equal_C[3]->Fill( average_smear_equal );    
    histChargeContained[3]->Fill( Energy );

    if ( energy5 > 1 && amp5 > 0.01 )
    {
      weight5 = energy5;
      count5 = 1;
    }
        
    Energy = weight1 + weight2 + weight3 + weight4 + weight5;
    average = ( weight1*vect[0].time + weight2*vect[1].time + weight3*vect[2].time + weight4*vect[3].time + weight5*vect[4].time )/( Energy );
    average_equal = ( count1*vect[0].time + count2*vect[1].time + count3*vect[2].time + count4*vect[3].time + count5*vect[4].time )/( count1 + count2 + count3 + count4 + count5 );
    average_smear = ( weight1*vect[0].time_smear + weight2*vect[1].time_smear + weight3*vect[2].time_smear + weight4*vect[3].time_smear + weight5*vect[4].time_smear )/( Energy );
    average_smear_equal = ( count1*vect[0].time_smear + count2*vect[1].time_smear + count3*vect[2].time_smear + count4*vect[3].time_smear + count5*vect[4].time_smear )/( count1 + count2 + count3 + count4 + count5 );
    histTOF_largest[4]->Fill( average );
    histTOF_largest_C[4]->Fill( average );
    histTOF_largest_equal_C[4]->Fill( average_equal );
    histTOF_largest_smear_C[4]->Fill( average_smear );
    histTOF_largest_smear_equal_C[4]->Fill( average_smear_equal );
    histChargeContained[4]->Fill( Energy );

    if ( energy6 > 1 && amp6 > 0.01 )
    {
      weight6 = energy6;
      count6 = 1;
    }
        
    Energy = weight1 + weight2 + weight3 + weight4 + weight5 + weight6;
    average = ( weight1*vect[0].time + weight2*vect[1].time + weight3*vect[2].time + weight4*vect[3].time + weight5*vect[4].time + weight6*vect[5].time )/( Energy );
    average_equal = ( count1*vect[0].time + count2*vect[1].time + count3*vect[2].time + count4*vect[3].time + count5*vect[4].time + count6*vect[5].time )/( count1 + count2 + count3 + count4 + count5 + count6 );
    average_smear = ( weight1*vect[0].time_smear + weight2*vect[1].time_smear + weight3*vect[2].time_smear + weight4*vect[3].time_smear + weight5*vect[4].time_smear + weight6*vect[5].time_smear )/( Energy );
    average_smear_equal = ( count1*vect[0].time_smear + count2*vect[1].time_smear + count3*vect[2].time_smear + count4*vect[3].time_smear + count5*vect[4].time_smear + count6*vect[5].time_smear )/( count1 + count2 + count3 + count4 + count5 + count6 );
    histTOF_largest[5]->Fill( average );
    histTOF_largest_C[5]->Fill( average );
    histTOF_largest_equal_C[5]->Fill( average_equal );
    histTOF_largest_smear_C[5]->Fill( average_smear );
    histTOF_largest_smear_equal_C[5]->Fill( average_smear_equal );
    histChargeContained[5]->Fill( Energy );

    if ( energy7 > 1 && amp7 > 0.01 )
    {
      weight7 = energy7;
      count7 = 1;
    }
      
    // values of offset are from the cubic fit of the TProfile, pol3 
    offset_correction = 0.06116;
    offset_slope = -0.003860;
    offset_quadratic = 0.00004994;
    offset_cubic = -0.0000000007116;

    Energy = weight1 + weight2 + weight3 + weight4 + weight5 + weight6 + weight7;
    average = ( weight1*vect[0].time + weight2*vect[1].time + weight3*vect[2].time + weight4*vect[3].time + weight5*vect[4].time + weight6*vect[5].time + weight7*vect[6].time )/( Energy );
    average_equal = ( count1*vect[0].time + count2*vect[1].time + count3*vect[2].time + count4*vect[3].time + count5*vect[4].time + count6*vect[5].time + count7*vect[6].time )/( count1 + count2 + count3 + count4 + count5 + count6 + count7 );
    average_smear = ( weight1*vect[0].time_smear + weight2*vect[1].time_smear + weight3*vect[2].time_smear + weight4*vect[3].time_smear + weight5*vect[4].time_smear + weight6*vect[5].time_smear + weight7*vect[6].time_smear )/( Energy );
    average_smear_equal = ( count1*vect[0].time_smear + count2*vect[1].time_smear + count3*vect[2].time_smear + count4*vect[3].time_smear + count5*vect[4].time_smear + count6*vect[5].time_smear + count7*vect[6].time_smear )/( count1 + count2 + count3 + count4 + count5 + count6 + count7 );
    histTOF_largest[6]->Fill( average );
    histTOF_largest_C[6]->Fill( average );
    histTOF_largest_equal_C[6]->Fill( average_equal );
    histTOF_largest_smear_C[6]->Fill( average_smear );
    histTOF_largest_smear_equal_C[6]->Fill( average_smear_equal ); 
    histChargeContained[6]->Fill( Energy );

    totalEnergy = energy1 + energy2 + energy3 + energy4 + energy5 + energy6 + energy7;

    // energy histogram with all pixel cuts applied
    float energy_total2 = weight1 + weight2 + weight3 + weight4 + weight5 + weight6 + weight7;
    ratio2 = energy_center / energy_total2;

    histEnergyRatioAllCuts_C->Fill( ratio2 );
    histEnergyPhotekAmpAllCuts_C->Fill( photekAmp, ratio2 );

    // apply offset correction so there is no dependence of time on charge
    average_corrected = average - offset_correction - Energy * offset_slope - Energy * Energy * offset_quadratic - Energy * Energy * Energy * offset_cubic;

    // Makes 2D histogram with delta t and charge for event, using overall cuts on Photek and center pixel (charge and amplitude). Plots just energy of pixels that pass pixel cuts.
    histDeltaTCharge_C->Fill( Energy , average );
    // Same but with corrected average time (based on fit)
    histDeltaTChargeCorr_C->Fill( Energy, average_corrected );

    // Makes a histogram of the TOF with time correction, has all 7 pixels combined
    histTOF_corr_C->Fill( average_corrected );
    histTOF_corr_fit_C->Fill( average_corrected ); 

    // Makes histogram of charge (or energy) vs. number of events with this charge
    histTotalChargeContained_C->Fill( totalEnergy );

    // Makes histogram of number of pixels combined (number of pixels above noise)
    Pixels = count1 + count2 + count3 + count4 + count5 + count6 + count7;
    histPixelsCombined_C->Fill( Pixels );

    // Makes a histogram with the charge contained in the center pixel
    histChargeCenterContained_C->Fill( energy_center );

    // Makes histogram of time vs number of events (to determine time resolution) using events based on center pixel and photek passing charge and amplitude cuts
    // Apply cuts such that only a specific energy (or charge) range is used for the events plotted
    // Time is weighted by charge in the pixel, same as the histTOF_largest above
    if ( Energy >= 20 && Energy <= 40 )
    {
      histTOFEnergyCut_C->Fill( average_corrected );
    }

    histEnergyDeltaT2D_C->Fill( Energy, average );
    histEnergyDeltaT2D_adj_C->Fill( Energy, average_corrected );

    if ( average >= 0.1 )
    {
      cout << "delta T " << average << ":\n" << endl;
      cout << "Event " << event << ":\n" << endl; 
    }

    // same but only use center pixel
    if ( energy_center >= 30 && energy_center <= 40 )
    {
      average = time_center ;
      histTOFEnergyCutCenter_C->Fill( average );
    }

    //std::cout << "event: " << iEntry << "-->" << vect[0].charge << " " << vect[1].charge << " " << vect[2].charge << std::endl;
  }
  



  // first loop fills histograms and finds their mean to shift to 0
  delete r;
  r = new TRandom3(seed);


  TH1F *histDeltaT_ring[6];
  for(int j=0; j<6; j++) histDeltaT_ring[j]= new TH1F(Form("histDeltaT_ring_%d",j),"; Time [ns];Number of Events", 250, 4.4, 5.6);

  TH1F *histDeltaT_smear_ring[6];
  for(int j=0; j<6; j++) histDeltaT_smear_ring[j]= new TH1F(Form("histDeltaT_smear_ring_%d",j),"; Time [ns];Number of Events", 250, 4.4, 5.6);


  std::cout<<"Number of events in Sample: "<<nentries<<std::endl;  
  for (Long64_t iEntry=0;iEntry<nentries;iEntry++)
  {
    if (iEntry %1000 == 0) cout << "Processing Event " << iEntry << "\n";
    tree->GetEntry(iEntry);    
  
    float photekTimeGauss = gauspeak[0];
    float photekAmp = amp[0];
    float photekCharge = integral[0];

    float PixelAmp = 0;
    float PixelCharge = 0;
    float PixelTime = 0;

    float largestIndexIntegral = 0;
    int largestIndex = 0;

    float DeltaT_ring[6] = {-99.};
    float DeltaT_smear_ring[6] = {-99.};

    for ( int j = 2; j <= 7; j++)
    {
      if ( integral[j] > largestIndexIntegral )
      {
        largestIndexIntegral = integral[j];
        largestIndex = j;
      }
    }

    PixelAmp = amp[largestIndex];
    PixelCharge = integral[largestIndex];
    PixelTime = linearTime45[largestIndex];

    // require photek amp and charge to be above cuts (to select electron events and not background), and require one of the first ring pixels to be above amp and charge cuts
    if( !(photekAmp > 0.1 && photekCharge > 2 ) ) continue;
    //require signal in largest of the ring pixels
    if( !(PixelCharge > 1 && PixelAmp > 0.01 ) ) continue;
    

    for ( int j = 2; j <= 7; j++)
    {
      float smear_number = r->Gaus(0, smear);

      if ( amp[j] > 0.01 && integral[j] > 1 ) 
      {
        DeltaT_ring[j-2] = gauspeak[0] - linearTime45[j];
        DeltaT_smear_ring[j-2] = gauspeak[0] - smear_number - linearTime45[j];
        histDeltaT_ring[j-2]->Fill(DeltaT_ring[j-2]);
        histDeltaT_smear_ring[j-2]->Fill(DeltaT_smear_ring[j-2]);
      }
    }
  }

  // get the means of the delta T ring histograms, both the smear and the non-smeared one. These means will be used to shift the histograms to 0 in the next event loop
  float meanT_ring[6];
  for ( int i = 0; i < 6; i++ )
  {
    meanT_ring[i] = histDeltaT_ring[i]->GetMean();
  }

  float meanT_smear_ring[6];
  for ( int i = 0; i < 6; i++ )
  {
    meanT_smear_ring[i] = histDeltaT_smear_ring[i]->GetMean();
  }


  // this event loop does the analysis
  // TOF of largest pixels, by charge weighting. The Gaussian fit will be done to this one.
  TH1F *histTOF_largest_ring_C[6];
  for(int j=0; j < 6; j++) histTOF_largest_ring_C[j]= new TH1F(Form("histTOF_largest_ring_C_%d",j),"; #Deltat (ns) ; Entries / (0.01 ns)", 80, -0.4, 0.4);

  // TOF of largest pixels, by equal weighting. The Gaussian fit will be done to this one.
  TH1F *histTOF_largest_equal_ring_C[6];
  for(int j=0; j < 6; j++) histTOF_largest_equal_ring_C[j]= new TH1F(Form("histTOF_largest_equal_ring_C_%d",j),"; #Deltat (ns) ; Entries / (0.01 ns)", 80, -0.4, 0.4);

  // TOF of largest pixels, with the smear added. For this, the times will be added with charge weighting.
  TH1F *histTOF_largest_smear_ring_C[6];
  for(int j=0; j < 6; j++) histTOF_largest_smear_ring_C[j]= new TH1F(Form("histTOF_largest_smear_ring_C_%d",j),"; #Deltat (ns) ; Entries / (0.01 ns)", 80, -1 * width, width);

  // TOF of largest pixels, with the smear added. For these, the times of each pixel will be added with equal weighting
  TH1F *histTOF_largest_smear_equal_ring_C[6];
  for(int j=0; j < 6; j++) histTOF_largest_smear_equal_ring_C[j]= new TH1F(Form("histTOF_largest_smear_equal_ring_C_%d",j),"; #Deltat (ns) ; Entries / (0.01 ns)", 80, -1 * width, width);

  TH1F* histPixelsCombinedRing_C = new TH1F("histPixelsCombinedRing","; Number of Pixels;Number of Events", 6, 0.5,6.5);


  // Delete memory allocated to r from first for loop. Then use same seed for the TRandom so the same random sequence is used in both event loops.
  delete r;
  r = new TRandom3(seed);


  // first ring analysis
  std::cout<<"Number of events in Sample: "<<nentries<<std::endl; 
  for (Long64_t iEntry=0;iEntry<nentries;iEntry++)
  {
    if (iEntry %1000 == 0) cout << "Processing Event " << iEntry << "\n";
    tree->GetEntry(iEntry);  

    float photekTimeGauss = gauspeak[0];
    float photekAmp = amp[0];
    float photekCharge = integral[0];

    float centerAmp = amp[1];
    float centerCharge = integral[1];
    float centerTime = linearTime45[1];

    float PixelAmp = 0;
    float PixelCharge = 0;
    float PixelTime = 0;

    float largestIndexIntegral = 0;
    int largestIndex = 0;


    for ( int j = 2; j <= 7; j++)
    {
      if ( integral[j] > largestIndexIntegral )
      {
        largestIndexIntegral = integral[j];
        largestIndex = j;
      }
    }

    PixelAmp = amp[largestIndex];
    PixelCharge = integral[largestIndex];
    PixelTime = linearTime45[largestIndex];



    // require photek amp and charge to be above cuts (to select electron events and not background), and require one of the first ring pixels to be above amp and charge cuts
    //if( !(photekAmp > 0.1 && photekCharge > 2 ) ) continue;



    // 8GeV cuts, 1 mm 
    //if( !(photekAmp > 0.015 && photekCharge > 0.4 ) ) continue;

    // 8GeV cuts, 32 mm
    //if( !(photekAmp > 0.04 && photekCharge > 1 ) ) continue;

    // 16GeV cuts, 1 mm
    //if( !(photekAmp > 0.03 && photekCharge > 0.8 ) ) continue;

    // 32GeV cuts, 1 mm
    if( !(photekAmp > 0.1 && photekCharge > 2 ) ) continue;

    // 32GeV cuts, 10 mm
    //if( !(photekAmp > 0.1 && photekCharge > 2 ) ) continue;

    // 32GeV cuts, 32 mm
    //if( !(photekAmp > 0.09 && photekCharge > 2 ) ) continue;

    // 32GeV cuts, 75 mm
    //if( !(photekAmp > 0.09 && photekCharge > 2 ) ) continue;

    //require signal in largest of the ring pixels
    if( !(PixelCharge > 1 && PixelAmp > 0.01 ) ) continue;


    std::vector< Pixel > vect2;
    Pixel pixel2;
    for ( int j = 2; j <= 7; j++)
    {
      // smear number is a random number choosen from a Gaussian distribution with mean 0 and sigam 0.05 (50 ps) to perform the time smearing for each pixel
      float smear_number = r->Gaus(0, smear);
      if ( iEntry > 1000 && iEntry < 1015 )
      {
        //cout << "smear number " << smear_number << ":\n" << endl;
      }
      if ( amp[j] > 0.01 && integral[j] > 1 )
      {
        histDeltaTshifted_smear2_C[j-2]->Fill(gauspeak[0] + smear_number - linearTime45[j] - meanT_smear_ring[j-2]);
      } 
      pixel2.index = j;
      pixel2.charge = integral[j];
      pixel2.time = gauspeak[0] - (linearTime45[j] + meanT_ring[j-2]);
      pixel2.time_smear = gauspeak[0] + smear_number - linearTime45[j] - meanT_smear_ring[j-2];
      pixel2.amp = amp[j];
      vect2.push_back( pixel2 ) ;
    }

    // Sort the vect and plot the combination of the largest energy pixels
    auto sortPixel = []( Pixel a, Pixel b ) {return a.charge > b.charge ?  true : false;};
    std::sort( vect2.begin(), vect2.end(), sortPixel );

    // These values will be used to determine if the pixel passes the required cuts
    float energy2 = vect2[1].charge;
    float energy3 = vect2[2].charge;
    float energy4 = vect2[3].charge;
    float energy5 = vect2[4].charge;
    float energy6 = vect2[5].charge;
    float energy7 = vect2[6].charge;

    float amp2 = vect2[1].amp;
    float amp3 = vect2[2].amp;
    float amp4 = vect2[3].amp;
    float amp5 = vect2[4].amp;
    float amp6 = vect2[5].amp;
    float amp7 = vect2[6].amp;

    float weight2 = 0;      // If the pixel passes the cuts, these weights will be changed to the charge of the pixel (and time will be charge weighted) 
    float weight3 = 0;
    float weight4 = 0;
    float weight5 = 0;
    float weight6 = 0;
    float weight7 = 0;     

    float count2 = 0;     // To count how many pixels are combined total. Also, if the pixel passes the cuts, this will be changed to 1 and will be used for the equal pixel weighting.
    float count3 = 0;
    float count4 = 0;
    float count5 = 0;
    float count6 = 0;
    float count7 = 0;
    float Pixels = 0;

    float Energy = 0;
    float average = 0;
    float average_equal = 0;
    float average_smear = 0;
    float average_smear_equal = 0;

    // This makes the TOF histogram with pixels with the largest charges, weighting the pixel time by the pixel charge (if pixel passes cuts). This is done event by event, so the TOF_largest_x will not necissarily have x pixels combined in every event.
    // For an additional pixel to be added, it must have a charge > 1 pC and an amplitude > 0.01 V. If the pixel does not pass these cuts, then it is added with a weight of 0.
        
    if ( energy2 > 1 && amp2 > 0.01 )
    {
      weight2 = energy2;
      count2 = 1;
    }

    Energy = weight2;
    average = vect2[0].time;
    average_smear = vect2[0].time_smear;
    histTOF_largest_ring_C[0]->Fill( average );
    histTOF_largest_equal_ring_C[0]->Fill( average );
    histTOF_largest_smear_ring_C[0]->Fill( average_smear );
    histTOF_largest_smear_equal_ring_C[0]->Fill( average_smear );

    if ( energy3 > 1 && amp3 > 0.01 )
    {
      weight3 = energy3;
      count3 = 1;
    }

    Energy = weight2 + weight3;
    average = ( weight2*vect2[0].time + weight3*vect2[1].time )/( Energy );
    average_equal = ( count2*vect2[0].time + count3*vect2[1].time )/( count2 + count3 );
    average_smear = ( weight2*vect2[0].time_smear + weight3*vect2[1].time_smear )/( Energy );
    average_smear_equal = ( count2*vect2[0].time_smear + count3*vect2[1].time_smear )/ (count2 + count3 );
    histTOF_largest_ring_C[1]->Fill( average ); 
    histTOF_largest_equal_ring_C[1]->Fill ( average_equal );
    histTOF_largest_smear_ring_C[1]->Fill( average_smear );
    histTOF_largest_smear_equal_ring_C[1]->Fill( average_smear_equal );
    
    if ( energy4 > 1 && amp4 > 0.01 )
    {
      weight4 = energy4;
      count4 = 1;
    }

    Energy = weight2 + weight3 + weight4 ;
    average = ( weight2*vect2[0].time + weight3*vect2[1].time + weight4*vect2[2].time )/( Energy );
    average_equal = ( count2*vect2[0].time + count3*vect2[1].time + count4*vect2[2].time )/( count2 + count3 + count4 );
    average_smear = ( weight2*vect2[0].time_smear + weight3*vect2[1].time_smear + weight4*vect2[2].time_smear )/( Energy );
    average_smear_equal = ( count2*vect2[0].time_smear + count3*vect2[1].time_smear + count4*vect2[2].time_smear )/( count2 + count3 + count4 );
    histTOF_largest_ring_C[2]->Fill( average );
    histTOF_largest_equal_ring_C[2]->Fill( average_equal );
    histTOF_largest_smear_ring_C[2]->Fill( average_smear );
    histTOF_largest_smear_equal_ring_C[2]->Fill( average_smear_equal );

    if ( energy5 > 1 && amp5 > 0.01 )
    {
      weight5 = energy5;
      count5 = 1;
    }

    Energy = weight2 + weight3 + weight4 + weight5 ;
    average = ( weight2*vect2[0].time + weight3*vect2[1].time + weight4*vect2[2].time + weight5*vect2[3].time )/( Energy );
    average_equal = ( count2*vect2[0].time + count3*vect2[1].time + count4*vect2[2].time + count5*vect2[3].time )/( count2 + count3 + count4 + count5  );
    average_smear = ( weight2*vect2[0].time_smear + weight3*vect2[1].time_smear + weight4*vect2[2].time_smear + weight5*vect2[3].time_smear )/( Energy );
    average_smear_equal = ( count2*vect2[0].time_smear + count3*vect2[1].time_smear + count4*vect2[2].time_smear + count5*vect2[3].time_smear )/( count2 + count3 + count4 + count5 );
    histTOF_largest_ring_C[3]->Fill( average );
    histTOF_largest_equal_ring_C[3]->Fill( average_equal );
    histTOF_largest_smear_ring_C[3]->Fill( average_smear );
    histTOF_largest_smear_equal_ring_C[3]->Fill( average_smear_equal );

    if ( energy6 > 1 && amp6 > 0.01 )
    {
      weight6 = energy6;
      count6 = 1;
    }

    Energy = weight2 + weight3 + weight4 + weight5 + weight6;
    average = ( weight2*vect2[0].time + weight3*vect2[1].time + weight4*vect2[2].time + weight5*vect2[3].time + weight6*vect2[4].time )/( Energy );
    average_equal = ( count2*vect2[0].time + count3*vect2[1].time + count4*vect2[2].time + count5*vect2[3].time + count6*vect2[4].time )/( count2 + count3 + count4 + count5 + count6 );
    average_smear = ( weight2*vect2[0].time_smear + weight3*vect2[1].time_smear + weight4*vect2[2].time_smear + weight5*vect2[3].time_smear + weight6*vect2[4].time_smear )/( Energy );
    average_smear_equal = ( count2*vect2[0].time_smear + count3*vect2[1].time_smear + count4*vect2[2].time_smear + count5*vect2[3].time_smear + count6*vect2[4].time_smear )/( count2 + count3 + count4 + count5 + count6 );
    histTOF_largest_ring_C[4]->Fill( average );
    histTOF_largest_equal_ring_C[4]->Fill( average_equal );
    histTOF_largest_smear_ring_C[4]->Fill( average_smear );
    histTOF_largest_smear_equal_ring_C[4]->Fill( average_smear_equal ); 
    
    if ( energy7 > 1 && amp7 > 0.01 )
    {
      weight7 = energy7;
      count7 = 1;
    }

    Energy = weight2 + weight3 + weight4 + weight5 + weight6 + weight7;
    average = ( weight2*vect2[0].time + weight3*vect2[1].time + weight4*vect2[2].time + weight5*vect2[3].time + weight6*vect2[4].time + weight7*vect2[5].time )/( Energy );
    average_equal = ( count2*vect2[0].time + count3*vect2[1].time + count4*vect2[2].time + count5*vect2[3].time + count6*vect2[4].time + count7*vect2[5].time )/( count2 + count3 + count4 + count5 + count6 + count7 );
    average_smear = ( weight2*vect2[0].time_smear + weight3*vect2[1].time_smear + weight4*vect2[2].time_smear + weight5*vect2[3].time_smear + weight6*vect2[4].time_smear + weight7*vect2[5].time_smear )/( Energy );
    average_smear_equal = ( count2*vect2[0].time_smear + count3*vect2[1].time_smear + count4*vect2[2].time_smear + count5*vect2[3].time_smear + count6*vect2[4].time_smear + count7*vect2[5].time_smear )/( count2 + count3 + count4 + count5 + count6 + count7 );
    histTOF_largest_ring_C[5]->Fill( average );
    histTOF_largest_equal_ring_C[5]->Fill( average_equal );
    histTOF_largest_smear_ring_C[5]->Fill( average_smear );
    histTOF_largest_smear_equal_ring_C[5]->Fill( average_smear_equal ); 

    // Makes histogram of number of pixels combined (number of pixels above noise)
    Pixels = count2 + count3 + count4 + count5 + count6 + count7;
    histPixelsCombinedRing_C->Fill( Pixels );

    //std::cout << "event: " << iEntry << "-->" << vect[0].charge << " " << vect[1].charge << " " << vect[2].charge << std::endl;
  }


  //Normalize hists
  histTotalCharge = NormalizeHist(histTotalCharge);
  histChargeRingOne = NormalizeHist(histChargeRingOne);
  histNChannelRingOne = NormalizeHist(histNChannelRingOne);
  histChargeCenterOverTotalCharge = NormalizeHist(histChargeCenterOverTotalCharge);
  histChargeRingOneOverTotalCharge = NormalizeHist(histChargeRingOneOverTotalCharge);

  // All fits use log likelihood method (QMLES)
  // Do Gaussian fit of delta T distributions for each pixel alone
  TF1* f1_g1[7];
  for(int j=0; j<7; j++)
  {
    double mean = histDeltaT_C[j]->GetMean();
    double rms = histDeltaT_C[j]->GetRMS();
    double xmin = mean-2.0*rms;
    double xmax = mean+2.0*rms;
    f1_g1[j] = new TF1( Form("g_fit_%d",j), "gaus(0)", xmin, xmax);
    //cout << "\nFitting Channel #" << j << ":\n" << endl;
    histDeltaT_C[j]->Fit(Form("g_fit_%d",j),"QMLES","",xmin,xmax);
  }

  // Do Gaussian fit of delta T distributions for each pixel alone with the 50 ps smear applied
  TF1* f1_g2[7];
  for(int j=0; j<7; j++)
  {
    double mean = histDeltaTshifted_smear_C[j]->GetMean();
    double rms = histDeltaTshifted_smear_C[j]->GetRMS();
    double xmin = mean-2.0*rms;
    double xmax = mean+2.0*rms;
    f1_g2[j] = new TF1( Form("g_fit_%d",j), "gaus(0)", xmin, xmax);
    //cout << "\nFitting Channel #" << j << ":\n" << endl;
    histDeltaTshifted_smear_C[j]->Fit(Form("g_fit_%d",j),"QMLES","",xmin,xmax);
  }

  // Do Gaussian fit of time resolution largest pixel distributions, with charge weighting
  TF1* f1_g3[7];
  for(int j=0; j<7; j++) 
  {
    double mean = histTOF_largest_C[j]->GetMean();
    double rms = histTOF_largest_C[j]->GetRMS();
    double xmin = mean-2.0*rms;
    double xmax = mean+2.0*rms;
    f1_g3[j] = new TF1( Form("g_fit_%d",j), "gaus(0)", xmin, xmax); 
    //cout << "\nFitting Channel #" << j << ":\n" << endl;
    histTOF_largest_C[j]->Fit(Form("g_fit_%d",j),"QMLES","", xmin,xmax);
  }

  // Do Gaussian fit of time resolution largest pixel distributions, with equal pixel weighting
  TF1* f1_g4[7];
  for(int j=0; j<7; j++) 
  {
    double mean = histTOF_largest_equal_C[j]->GetMean();
    double rms = histTOF_largest_equal_C[j]->GetRMS();
    double xmin = mean-2.0*rms;
    double xmax = mean+2.0*rms;
    f1_g4[j] = new TF1( Form("g_fit_%d",j), "gaus(0)", xmin, xmax); 
    //cout << "\nFitting Channel #" << j << ":\n" << endl;
    histTOF_largest_equal_C[j]->Fit(Form("g_fit_%d",j),"QMLES","", xmin,xmax);
  }

  // Do Gaussian fit of the time resolution largest pixel distributions with the 50 ps smear added (charge weighting)
  TF1* f1_g5[7];
  for(int j=0; j<7; j++) 
  {
    double mean = histTOF_largest_smear_C[j]->GetMean();
    double rms = histTOF_largest_smear_C[j]->GetRMS();
    double xmin = mean-2.0*rms;
    double xmax = mean+2.0*rms;
    f1_g5[j] = new TF1( Form("g_fit_%d",j), "gaus(0)", xmin, xmax); 
    //cout << "\nFitting Channel #" << j << ":\n" << endl;
    histTOF_largest_smear_C[j]->Fit(Form("g_fit_%d",j),"QMLES","", xmin,xmax);
  }

  // Do Gaussian fit of the time resolution largest pixel distributions with the 50 ps smear added (equal weighting)
  TF1* f1_g6[7];
  for(int j=0; j<7; j++) 
  {
    double mean = histTOF_largest_smear_equal_C[j]->GetMean();
    double rms = histTOF_largest_smear_equal_C[j]->GetRMS();
    double xmin = mean-2.0*rms;
    double xmax = mean+2.0*rms;
    f1_g6[j] = new TF1( Form("g_fit_%d",j), "gaus(0)", xmin, xmax); 
    //cout << "\nFitting Channel #" << j << ":\n" << endl;
    histTOF_largest_smear_equal_C[j]->Fit(Form("g_fit_%d",j),"QMLES","", xmin,xmax);
  }

  // Do Gaussian fit of time resolution for 7 pixels combined time resolution, with time-charge correction performed.
  TF1* f1_g7;
  double mean = histTOF_corr_fit_C->GetMean();
  double rms = histTOF_corr_fit_C->GetRMS();
  double xmin = mean - 2.0 * rms;
  double xmax = mean + 2.0 * rms;
  f1_g7 = new TF1(Form("g_fit"), "gaus(0)", xmin, xmax);
  //cout << "\nFitting Corrected Time Resolution" << endl;
  histTOF_corr_fit_C->Fit(Form("g_fit"),"QMLES","", xmin, xmax);

  // Do Gaussian fit of time resolution for 7 pixels combined time resolution, with charge cuts performed.
  TF1* f1_g8;
  double mean_2 = histTOFEnergyCut_C->GetMean();
  double rms_2 = histTOFEnergyCut_C->GetRMS();
  double xmin_2 = mean_2 - 2.0 * rms_2;
  double xmax_2 = mean_2 + 2.0 * rms_2;
  f1_g8 = new TF1(Form("g_fit"), "gaus(0)", xmin_2, xmax_2);
  //cout << "\nFitting Corrected Time Resolution" << endl;
  histTOFEnergyCut_C->Fit(Form("g_fit"),"QMLES","", xmin_2, xmax_2);

  // Do Gaussian fit of time resolution for ring largest pixels with charge weighing
  TF1* f1_g9[6];
  for(int j=0; j<6; j++) 
  {
    double mean = histTOF_largest_ring_C[j]->GetMean();
    double rms = histTOF_largest_ring_C[j]->GetRMS();
    double xmin = mean-2.0*rms;
    double xmax = mean+2.0*rms;
    f1_g9[j] = new TF1( Form("g_fit_%d",j), "gaus(0)", xmin, xmax); 
    //cout << "\nFitting Channel #" << j << ":\n" << endl;
    histTOF_largest_ring_C[j]->Fit(Form("g_fit_%d",j),"QMLES","", xmin,xmax);
  }

  // Do Gaussian fit of time resolution for ring largest pixels with equal weighting
  TF1* f1_g10[6];
  for(int j=0; j<6; j++) 
  {
    double mean = histTOF_largest_equal_ring_C[j]->GetMean();
    double rms = histTOF_largest_equal_ring_C[j]->GetRMS();
    double xmin = mean-2.0*rms;
    double xmax = mean+2.0*rms;
    f1_g10[j] = new TF1( Form("g_fit_%d",j), "gaus(0)", xmin, xmax); 
    //cout << "\nFitting Channel #" << j << ":\n" << endl;
    histTOF_largest_equal_ring_C[j]->Fit(Form("g_fit_%d",j),"QMLES","", xmin,xmax);
  }

  // Do Gaussian fit of time resolution for ring largest pixels with charge weighting, smearing applied
  TF1* f1_g11[6];
  for(int j=0; j<6; j++) 
  {
    double mean = histTOF_largest_smear_ring_C[j]->GetMean();
    double rms = histTOF_largest_smear_ring_C[j]->GetRMS();
    double xmin = mean-2.0*rms;
    double xmax = mean+2.0*rms;
    f1_g11[j] = new TF1( Form("g_fit_%d",j), "gaus(0)", xmin, xmax); 
    //cout << "\nFitting Channel #" << j << ":\n" << endl;
    histTOF_largest_smear_ring_C[j]->Fit(Form("g_fit_%d",j),"QMLES","", xmin,xmax);
  }

  // Do Gaussian fit of time resolution for ring largest pixels with equal weighting, smearing applied
  TF1* f1_g12[6];
  for(int j=0; j<6; j++) 
  {
    double mean = histTOF_largest_smear_equal_ring_C[j]->GetMean();
    double rms = histTOF_largest_smear_equal_ring_C[j]->GetRMS();
    double xmin = mean-2.0*rms;
    double xmax = mean+2.0*rms;
    f1_g12[j] = new TF1( Form("g_fit_%d",j), "gaus(0)", xmin, xmax); 
    //cout << "\nFitting Channel #" << j << ":\n" << endl;
    histTOF_largest_smear_equal_ring_C[j]->Fit(Form("g_fit_%d",j),"QMLES","", xmin,xmax);
  }

  // creates root file
  TFile *file = TFile::Open(outputFilename.c_str(), "RECREATE");
  file->cd();
  file->WriteTObject(histTotalCharge,"histTotalCharge", "WriteDelete");
  file->WriteTObject(histTOFCenter,"histTOFCenter", "WriteDelete");
  file->WriteTObject(histTOFFlatAvg,"histTOFFlatAvg", "WriteDelete");
  file->WriteTObject(histTOFChargeWeightedAvg,"histTOFChargeWeightedAvg", "WriteDelete");
  file->WriteTObject(histTOFRingOneFlatAvg,"histTOFRingOneFlatAvg", "WriteDelete");
  file->WriteTObject(histTOFRingOneChargeWeightedAvg,"histTOFRingOneChargeWeightedAvg", "WriteDelete");
  file->WriteTObject(histNChannelRingOne,"histNChannelRingOne", "WriteDelete");
  file->WriteTObject(histChargeRingOne,"histChargeRingOne", "WriteDelete");
  file->WriteTObject(histChargeCenterOverTotalCharge,"histChargeCenterOverTotalCharge", "WriteDelete");
  file->WriteTObject(histChargeRingOneOverTotalCharge,"histChargeRingOneOverTotalCharge", "WriteDelete");
  file->WriteTObject(histDeltaTCombined,"histDeltaTCombined", "WriteDelete");


  for( int i = 0; i < 7; i++ )
  {
    histTOF_largest[i]->Write( Form("TOF_largest_%d", i+1) );
    histTOF_largest_C[i]->Write( Form("TOF_largest_%d_fit", i+1) );
    histTOF_largest_equal_C[i]->Write( Form("TOF_largest_equal_%d", i+1) );
    histTOF_largest_smear_C[i]->Write( Form("TOF_largest_smear_%d", i+1) );
    histTOF_largest_smear_equal_C[i]->Write( Form("TOF_largest_smear_equal_%d", i+1) );
  }

  for ( int i = 0; i < 7; i++ )
  {
    histChargeContained[i]->Write( Form("Charge_Contained_%d", i+1) );
  }

  file->WriteTObject(histMaxIndex_C,"Max Index", "WriteDelete");
  file->WriteTObject(histPixelsCombined_C,"Pixels Combined", "WriteDelete");
  file->WriteTObject(histPixelsCombinedRing_C,"Pixels Combined Ring", "WriteDelete");

  file->WriteTObject(histEnergyPhotekAmp_C,"Energy Ratio vs Photek Amp");
  file->WriteTObject(histEnergyPhotekAmpPhotekCut_C,"Energy Ratio vs Photek Amp, Photek cut");
  file->WriteTObject(histEnergyPhotekAmpEventCut_C,"Energy Ratio vs Photek Amp, event cut");
  file->WriteTObject(histEnergyPhotekAmpAllCuts_C,"Energy Ratio vs Photek Amp, all cuts");

  file->WriteTObject(histEnergyRatio_C,"E Center / E 7 Ratio", "WriteDelete");
  file->WriteTObject(histEnergyRatioPhotekCut_C,"E Center / E 7 Ratio, Photek cut", "WriteDelete");
  file->WriteTObject(histEnergyRatioEventCut_C,"E Center / E 7 Ratio, event cuts", "WriteDelete");
  file->WriteTObject(histEnergyRatioAllCuts_C,"E Center / E 7 Ratio, all cuts", "WriteDelete");

  file->WriteTObject(histTotalChargeContained_C,"Total Charge Contained", "WriteDelete");

  file->WriteTObject(histChargeCenterContained_C,"ChargeCenterContained", "WriteDelete");
  file->WriteTObject(histEnergyDeltaT2D_C,"EnergyDeltaT","WriteDelete");
  file->WriteTObject(histEnergyDeltaT2D_adj_C,"EnergyDeltaT_adj","WriteDelete");

  file->WriteTObject(histTOFEnergyCut_C,"TOF with energy cut", "WriteDelete");
  file->WriteTObject(histTOFEnergyCutCenter_C,"TOF with energy cut for center pixel", "WriteDelete");

  file->WriteTObject(histDeltaTCharge_C,"Delta T vs Charge", "WriteDelete");
  file->WriteTObject(histDeltaTChargeCorr_C,"Delta T vs Charge Correction", "WriteDelete");
  file->WriteTObject(histTOF_corr_C, "TOF with time-charge correction", "WriteDelete");
  file->WriteTObject(histTOF_corr_fit_C, "TOF with time-charge correction fitted", "WriteDelete");

  file->WriteTObject(histPhotekAmpCut_C,"Photek Amplitude Cut", "WriteDelete");
  file->WriteTObject(histPhotekChargeCut_C,"Photek Charge Cut", "WriteDelete");
  file->WriteTObject(histCenterChargeCut_C,"Center Charge Cut", "WriteDelete");
  file->WriteTObject(histCenterAmpCut_C,"Center Amplitude Cut", "WriteDelete");

  file->WriteTObject(histRandomNumberTest_C,"Random number", "WriteDelete");
  file->WriteTObject(histLinearTime_C,"linearTime45","WriteDelete");

  for( int i = 0; i < 7; i++ )
  {
    histDeltaT[i]->Write( Form("deltaT_%d", i+1) );
    histDeltaTshifted_C[i]->Write( Form("deltaT_shifted_%d", i+1) );
    histDeltaT_C[i]->Write( Form("deltaT_shifted_%d_fit", i+1) );
    histDeltaTshifted_smear_C[i]->Write( Form("deltaTshifted_smear_%d_fit", i+1) );
  }

  for (int i = 0; i < 6; i++ )
  {
    histDeltaTshifted_smear2_C[i]->Write(Form("deltaTshifted_smear2_%d", i+1) );
    histTOF_largest_ring_C[i]->Write( Form("TOF_largest_ring_%d", i+1) );
    histTOF_largest_equal_ring_C[i]->Write( Form("TOF_largest_equal_ring_%d", i+1) );
    histTOF_largest_smear_ring_C[i]->Write( Form("TOF_largest_smear_ring_%d", i+1) );
    histTOF_largest_smear_equal_ring_C[i]->Write( Form("TOF_largest_smear_equal_ring_%d", i+1) );
  }




  // plot fits and sigmas on pdf files
  TCanvas *c = new TCanvas ("c","c",800, 600); 
  TLatex *tex = new TLatex();
  tex->SetNDC(); // Sets coords such that (0,0) is bottom left & (1,1) is top right.
  tex->SetTextSize(0.060);
  tex->SetTextFont(42);
  tex->SetTextColor(kBlack);

  c->cd();

  histMaxIndex_C->Draw();
  c->SaveAs( Form("MaxIndex.pdf"));

  histPixelsCombined_C->Draw();
  c->SaveAs( Form("PixelsCombined.pdf"));

  histPixelsCombinedRing_C->Draw();
  c->SaveAs( Form("PixelsCombinedRing.pdf"));



  // plots of energy contained and time resolution with energy cuts on event selection
  histTOFEnergyCut_C->Draw();
  gStyle->SetOptFit(1);
  gStyle->SetOptStat(1);
  c->SaveAs( Form("time_res_energy.pdf"));

  // plot of total energy contained
  histTotalChargeContained_C->Draw();
  gStyle->SetOptFit(1);
  gStyle->SetOptStat(1);
  c->SaveAs( Form("total_charge_contained.pdf"));

/*
  // all pixel plots
  for (int i = 0; i < 7; i++ )
  {
    histTOF_largest_C[i]->Draw();
    gStyle->SetOptFit(1);
    gStyle->SetOptStat(1);
    c->SaveAs( Form("all_charge_%d.pdf", i+1) );
  }

  for (int i = 0; i < 7; i++ )
  {
    histTOF_largest_equal_C[i]->Draw();
    gStyle->SetOptFit(1);
    gStyle->SetOptStat(1);
    c->SaveAs( Form("all_equal_%d.pdf", i+1) );
  }

  for (int i = 0; i < 7; i++ )
  {
    histTOF_largest_smear_C[i]->Draw();
    gStyle->SetOptFit(1);
    gStyle->SetOptStat(1);
    c->SaveAs( Form("all_smear_charge_%d.pdf", i+1) );
  }

  for (int i = 0; i < 7; i++ )
  {
    histTOF_largest_smear_equal_C[i]->Draw();
    gStyle->SetOptFit(1);
    gStyle->SetOptStat(1);
    c->SaveAs( Form("all_smear_equal_%d.pdf", i+1) );
  }

  // ring pixels plots
  for (int i = 0; i < 6; i++ )
  {
    histTOF_largest_ring_C[i]->Draw();
    gStyle->SetOptFit(1);
    gStyle->SetOptStat(1);
    c->SaveAs( Form("ring_charge_%d.pdf", i+1) );
  }

  for (int i = 0; i < 6; i++ )
  {
    histTOF_largest_equal_ring_C[i]->Draw();
    gStyle->SetOptFit(1);
    gStyle->SetOptStat(1);
    c->SaveAs( Form("ring_equal_%d.pdf", i+1) );
  }

  for (int i = 0; i < 6; i++ )
  {
    histTOF_largest_smear_ring_C[i]->Draw();
    gStyle->SetOptFit(1);
    gStyle->SetOptStat(1);
    c->SaveAs( Form("ring_smear_charge_%d.pdf", i+1) );
  }

  for (int i = 0; i < 6; i++ )
  {
    histTOF_largest_smear_equal_ring_C[i]->Draw();
    gStyle->SetOptFit(1);
    gStyle->SetOptStat(1);
    c->SaveAs( Form("ring_smear_equal_%d.pdf", i+1) );
  }
*/

  file->Close();
  delete file;
}


void MultiChannelStudy_TimingMethod1()
{
  // DoMultiChannelStudy("t1065-jun-2016-90.dat-full.root","output.90.root");
  // DoMultiChannelStudy("t1065-jun-2016-94.dat-full.root","output.94.root");
  // DoMultiChannelStudy("t1065-jun-2016-81.dat-full.root","output.81.root");
  //DoMultiChannelStudy("../../raw/combine_32gev_1cm.root","output_32gev_1cm.root");
  //DoMultiChannelStudy("../../raw/combine_16gev_1mm.root","output_16gev_1mm.root");
  //DoMultiChannelStudy("../../raw/combine_8gev_1mm.root","output_8gev_1mm.root");
  //DoMultiChannelStudy("../../raw/combine_8gev_32mm.root","output_8gev_32mm.root");
  //DoMultiChannelStudy("../../raw/combine_32gev_1mm.root","output_32gev_1mm.root");
  //DoMultiChannelStudy("../../raw/combine_32gev_10mm.root","output_32gev_10mm.root");
  //DoMultiChannelStudy("../../raw/combine_32gev_32mm.root","output_32gev_32mm.root");
  //DoMultiChannelStudy("../../raw/combine_32gev_75mm.root","output_32gev_75mm.root");
  DoMultiChannelStudy("../../raw/combine_32gev_1mm_fixed.root","output_32gev_1mm_fix.root");
  //DoMultiChannelStudy("../../raw/combine_32gev_32mm_fixed.root","output_32gev_32mm_fix.root");
  //DoMultiChannelStudy("../../raw/combine_16gev_1mm_fixed.root","output_16gev_1mm_fix.root");
  //DoMultiChannelStudy("../../raw/combine_8gev_1mm_fixed.root","output_8gev_1mm_fix.root");

  //DoMultiChannelStudy("../../raw/144-155/t1065-jun-2016-148.dat-full.root","output_148.root");
  //DoMultiChannelStudy("../../raw/combine_120GeV_1mm.root","output_120gev_1mm.root");
}