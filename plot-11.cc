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
TFile*  file;
TTree*  tree;


int 
main(int argc, char **argv){

  double off_mean[4][9][1024];

  FILE* fp1;
  char stitle[200];
  int dummy;

  for( int i = 0; i < 4; i++){
    sprintf( stitle, "v1740_bd0_group_%d_offset.txt", i);

    fp1 = fopen( stitle, "r");
    printf("offset data : %s\n", stitle);

    for( int k = 0; k < 1024; k++)      
      for( int j = 0; j < 9; j++){      
	dummy = fscanf( fp1, "%lf ", &off_mean[i][j][k]);       
	if( k < 2 && 0)
	  printf("%5d  %8.4f\n", j, off_mean[i][j][k]);
      }
  
    fclose(fp1);
  }

  graphic_init();
  TGraphErrors* gr[4][9];

  char hid[200], htitle[200];
  sprintf( hid, "%s", "h2");
  sprintf( htitle, "%s", "h2");

  TH2F* h2 = new TH2F( hid, htitle, 1024, 0, 1024, 400, -2048, 2048);  
  h2->SetBit(TH1::kNoStats);

  TH1F* h_tc = new TH1F("h_tc", "Trigger cells", 1024, 0, 1023);

  TCanvas* c1 = new TCanvas("c1", "c1", 10, 10, 700, 500);
  c1->SetFillStyle(4000);
  c1->Divide(4, 2);
  
  char title[200];  
  sprintf( title, "%s-%s.ps", argv[0], argv[1]);
  TPostScript* psf1 = new TPostScript( title, 112);

  int event;
  ushort  b_c[4][9][1024], tc[4]; 

  // loop over root files
  sprintf( title, "%s.root", argv[1]);

  file = new TFile( title);
  printf("%s data processed.\n", title);

  tree = (TTree*)file->Get("pulse");

  TBranch* t_event = tree->GetBranch("event");
  TBranch* t_tc    = tree->GetBranch("tc");
  TBranch* t_b_c   = tree->GetBranch("b_c");

  t_event->SetAddress( &event);
  t_tc->SetAddress( tc);
  t_b_c->SetAddress( b_c);

  double index[1024], err[1024];
  for( int i = 0; i < 1024; i++){
    index[i] = i;
    err[i] = 0;
  }

  double amplitude[4][9][1024];

  for( int eventn = 0; eventn < tree->GetEntries(); eventn++){
    tree->GetEntry(eventn);

    if(eventn%100 == 0) printf("---- event  %5d,  %5d\n", eventn, event);

    for( int group = 0; group < 4; group++){

      h_tc->Fill( tc[group]);

      for( int i = 0; i < 9; i++)
	for( int j = 0; j < 1024; j++){
	  amplitude[group][i][j] = (double)b_c[group][i][j] - off_mean[group][i][(j+tc[group])%1024];  
	  amplitude[group][i][j] += 235;
	}

      
      for( int i = 0; i < 9; i++)
	gr[group][i] = new TGraphErrors( 1024, index, amplitude[group][i], err, err);
    }

    
    if( eventn%100 == 0 ){
      char stitle[200];	
      psf1->NewPage();
      for( int i = 0; i < 4; i++){
	c1->cd(i+1);
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
      
      
  }

  psf1->NewPage();
  c1->cd(1);
  h_tc->Draw();
  c1->Update();

  psf1->Close();


  file->Close();


  return 0;
  printf("dummy = %d\n", dummy);
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
