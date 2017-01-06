#include <TFile.h>
#include <TTree.h>
#include <TROOT.h>
#include <TRandom.h>
#include <TApplication.h>
#include <TCanvas.h>
#include <TPostScript.h>
#include <TAxis.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TF1.h>
#include <TLine.h>
#include <TStyle.h>
#include <TProfile.h>
#include <TMapFile.h>
#include <TPaveStats.h>

#include <fstream>
#include <iomanip>
#include <string>
#include <iostream>

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>  
#include <string.h>
#include <unistd.h>
#include <sys/mman.h>
#include <errno.h>
#include <ctype.h>
#include <math.h>
#include <time.h>

int graphic_init();

TStyle* style;

int 
main(int argc, char **argv){

  graphic_init();
  TGraphErrors* gr[4][9];

  TH2F* h2 = new TH2F("h2", "h2", 200, 0, 1024, 400, 0, 4096);  
  h2->SetBit(TH1::kNoStats);

  TH1F* h_tc = new TH1F("h_tc", "Trigger cells", 1024, 0, 1023);
  
  TCanvas* c1 = new TCanvas("c1", "c1", 10, 10, 700, 500);
  c1->SetFillStyle(4000);
  c1->Divide(4, 2);
  
  char title[200];  
  sprintf( title, "%s-%s.ps", argv[0], argv[1]);
  TPostScript* psf1 = new TPostScript( title, 112);

  uint event_header;
  uint temp[3];
  uint samples[9][1024];

  // loop over root files
  sprintf( title, "data/%s.dat", argv[1]);

  FILE* fpin = fopen( title, "r");

  for( int event = 0; event < 2; event++){
    printf("---- event  %5d\n", event);

    fread( &event_header, sizeof(uint), 1, fpin);  
    printf("event header 0 = 0x%08x\n", event_header);
    fread( &event_header, sizeof(uint), 1, fpin);  
    printf("event_header 1 = 0x%08x\n", event_header);
    fread( &event_header, sizeof(uint), 1, fpin);  
    printf("event header 2 = 0x%08x\n", event_header);
    fread( &event_header, sizeof(uint), 1, fpin);  
    printf("event_header 3 = 0x%08x\n", event_header);

    double index[1024], err[1024], amplitude[4][9][1024];
    for( int i = 0; i < 1024; i++){
      index[i] = i;
      err[i] = 0;
    }

    for( int group = 0; group < 4; group++){
      fread( &event_header, sizeof(uint), 1, fpin);  

      printf("group %d description = 0x%08x\n", group, event_header);

      uint tc = (event_header >> 20) & 0xfff;

      h_tc->Fill( tc);

      // if( 800 < tc)
      //	printf("event = %5d,  gr%d  header = 0x%08x,  tc = %4d\n", event, group, event_header, tc);
    
      int nsample = (event_header & 0xfff)/3;
      printf("nsample = %d\n", nsample);

      for(int i = 0; i < nsample; i++){
	fread( &temp, sizeof(uint), 3, fpin);  
	samples[0][i] =  temp[0] & 0xfff;
	samples[1][i] = (temp[0] >> 12) & 0xfff;
	samples[2][i] = (temp[0] >> 24) | ((temp[1] & 0xf) << 8);
	samples[3][i] = (temp[1] >>  4) & 0xfff;
	samples[4][i] = (temp[1] >> 16) & 0xfff;
	samples[5][i] = (temp[1] >> 28) | ((temp[2] & 0xff) << 4);
	samples[6][i] = (temp[2] >>  8) & 0xfff;
	samples[7][i] =  temp[2] >> 20;
      }

      for(int j = 0; j < nsample/8; j++){
	fread( &temp, sizeof(uint), 3, fpin);  
	samples[8][j*8+0] =  temp[0] & 0xfff;
	samples[8][j*8+1] = (temp[0] >> 12) & 0xfff;
	samples[8][j*8+2] = (temp[0] >> 24) | ((temp[1] & 0xf) << 8);
	samples[8][j*8+3] = (temp[1] >>  4) & 0xfff;
	samples[8][j*8+4] = (temp[1] >> 16) & 0xfff;
	samples[8][j*8+5] = (temp[1] >> 28) | ((temp[2] & 0xff) << 4);
	samples[8][j*8+6] = (temp[2] >>  8) & 0xfff;
	samples[8][j*8+7] =  temp[2] >> 20;
      }
      
      for( int i = 0; i < 9; i++){
	for( int j = 0; j < 1024; j++)
	  amplitude[group][i][j] = samples[i][j];	
      }	

      for( int i = 0; i < 9; i++)
	gr[group][i] = new TGraphErrors( 1024, index, amplitude[group][i], err, err);

      fread( &event_header, sizeof(uint), 1, fpin);  
    }
  
    char stitle[200];

    psf1->NewPage();
    for( int i = 0; i < 4; i++){
      c1->cd(i+1);
      h2->SetAxisRange(1500,2500,"y");
      h2->Draw();
      sprintf(stitle, "Ch #%d", i*8);
      h2->SetTitle(stitle);
      gr[i][0]->Draw("p");
      c1->Update();
    }

    for( int i = 0; i < 4; i++){
      c1->cd(i+5);
      h2->Draw();
      sprintf(stitle, "Digitized trigger: Group #%d", i);
      h2->SetTitle(stitle);
      gr[i][8]->Draw("p");
      c1->Update();
    }
  }
  
  psf1->NewPage();
  c1->cd(1);
  h_tc->Draw();
  c1->Update();

  psf1->Close();

  fclose(fpin);

  return 0;
}



int
graphic_init( void){

  style = new TStyle("style", "style");
  
  style->SetLabelFont(132,"X");
  style->SetLabelFont(132,"Y");
  style->SetTitleFont(132,"X");
  style->SetTitleFont(132,"Y");
  style->SetTitleFont(132,"");
  style->SetTitleFontSize( 0.07);
  style->SetStatFont(132);
  style->GetAttDate()->SetTextFont(132);

  style->SetStatW(0.20);
  style->SetStatH(0.23);

  style->SetFuncColor(2);
  style->SetFuncWidth(2);
  style->SetLineWidth(2);
  
  style->SetOptFile(0);
  style->SetOptTitle(1);

  style->SetFrameBorderMode(0);
  style->SetCanvasBorderMode(0);
  style->SetPadBorderMode(0);
  style->SetTitleStyle(4000);
  style->SetPadColor(0);
  style->SetCanvasColor(0);

  style->SetTitleFillColor(0);

  style->SetTitleBorderSize(0);
  //  style->SetTitleX(0.3);
  //  style->SetTitleY(0.06);
  style->SetStatColor(0);
  style->SetStatBorderSize(1);

  style->SetOptStat("emri");
  // style->SetOptStat(1);
  style->SetOptFit(1);
  style->SetTitleOffset( 1.0,"Y");

  style->SetMarkerStyle(20);
  style->SetMarkerSize( 0.3);
  style->SetMarkerColor(4);

  // style->SetOptDate(21);

  //  style->SetPadGridX(1);
  //  style->SetPadGridY(1);

  style->cd();

  return 0;
}
