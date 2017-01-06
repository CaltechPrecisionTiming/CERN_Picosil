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

  TH2F* h2 = new TH2F("h2", "h2", 200, 0, 1024, 400, 0, 4096);  
  h2->SetBit(TH1::kNoStats);

  TH1F* h_tc = new TH1F("h_tc", "Trigger cells", 1024, 0, 1023);
  
  TH1F* h_offset[4][9][1024];
  TH1F* h_offmean[4][9];
  TH1F* h_offrms[4][9];

  char hid[200], htitle[200];

  for( int i = 0; i < 4; i++)
    for( int j = 0; j < 9; j++)
      for( int k = 0; k < 1024; k++){
     
	sprintf( hid, "h_offset_%d_%d_%d", i, j, k);
	sprintf( htitle, "Offset_%d_%d_%d", i, j, k);

	h_offset[i][j][k] = new TH1F( hid, htitle, 250, 1900, 2400);
      }


  for( int i = 0; i < 4; i++)
    for( int j = 0; j < 9; j++){
      sprintf( hid, "h_offmean_%d_%d", i, j);
      sprintf( htitle, "Offmean_%d_%d", i, j);

      h_offmean[i][j] = new TH1F( hid, htitle, 250, 1900, 2400);

      sprintf( hid, "h_offrms_%d_%d", i, j);
      sprintf( htitle, "Offrms_%d_%d", i, j);

      h_offrms[i][j] = new TH1F( hid, htitle, 100, 0., 5.);
    }

  uint event_header;
  uint temp[3];
  uint samples[9][1024];
  uint dummy;

  char title[200];  


  // loop over root files
  sprintf( title, "/kdrive/data1/caen/2015-11/11-25/%s.dat", argv[1]);

  FILE* fpin = fopen( title, "r");

  for( int event = 0; event < 1000; event++){
    printf("---- event  %5d\n", event);

    dummy = fread( &event_header, sizeof(uint), 1, fpin);  
    printf("event header 0 = 0x%08x\n", event_header);
    dummy = fread( &event_header, sizeof(uint), 1, fpin);  
    printf("event header 1 = 0x%08x\n", event_header);
    dummy = fread( &event_header, sizeof(uint), 1, fpin);  
    printf("event header 2 = 0x%08x\n", event_header);
    dummy = fread( &event_header, sizeof(uint), 1, fpin);  
    printf("event header 3 = 0x%08x\n", event_header);

    for( int group = 0; group < 4; group++){
      dummy = fread( &event_header, sizeof(uint), 1, fpin);  

      uint tc = (event_header >> 20) & 0xfff;

      h_tc->Fill( tc);

      int nsample = (event_header & 0xfff)/3;
      printf("nsample = %d\n", nsample);

      for(int i = 0; i < nsample; i++){
	dummy = fread( &temp, sizeof(uint), 3, fpin);  
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
      

      for( int i = 0; i < 9; i++)
	for( int j = 0; j < 1024; j++){
	  h_offset[group][i][(tc+j)%1024]->Fill( samples[i][j]);
	}
        
      dummy = fread( &event_header, sizeof(uint), 1, fpin);  
    }
  }

  TCanvas* c1 = new TCanvas("c1", "c1", 10, 10, 700, 500);
  c1->SetFillStyle(4000);
  c1->Divide(5, 2);
  
  sprintf( title, "%s-%s.ps", argv[0], argv[1]);
  TPostScript* psf1 = new TPostScript( title, 112);

  double off_mean[4][9][1024], off_rms[4][9][1024];

  for( int i =0; i < 4; i++)
    for( int j =0; j < 9; j++)
      for( int k =0; k < 1024; k++){
	h_offmean[i][j]->Fill( h_offset[i][j][k]->GetMean());	
	h_offrms[i][j]->Fill( h_offset[i][j][k]->GetRMS());	
	off_mean[i][j][k] =  h_offset[i][j][k]->GetMean();
	off_rms[i][j][k] =  h_offset[i][j][k]->GetRMS();

	if( h_offset[i][j][k]->GetRMS() > 5){
	  printf("i = %4d, j = %4d, k = %4d  =>  %8.2f\n", i, j, k, off_rms[i][j][k]);
	}
      }

  for( int i =0; i < 4; i++){
    psf1->NewPage();

    for( int j =0; j < 9; j++){
      c1->cd(1+j);
      h_offmean[i][j]->Draw();
      c1->Update();
    }
  }

  for( int i =0; i < 4; i++){
    psf1->NewPage();

    for( int j =0; j < 9; j++){
      c1->cd(1+j);
      h_offrms[i][j]->Draw();
      c1->Update();
    }
  }

  for( int i =0; i < 4; i++)
    for( int k =0; k < 1; k++){
      psf1->NewPage();

      for( int j =0; j < 9; j++){
	c1->cd(1+j);
	h_offset[i][j][k]->Draw();
	c1->Update();
      }
    }

  for( int i =0; i < 4; i++)
    for( int k =1023; k < 1024; k++){
      psf1->NewPage();

      for( int j =0; j < 9; j++){
	c1->cd(1+j);
	h_offset[i][j][k]->Draw();
	c1->Update();
      }
    }

  psf1->NewPage();
  c1->cd(1);
  h_offset[1][3][1009]->Draw();
  c1->Update();

  c1->cd(2);
  h_offset[1][3][1010]->Draw();
  c1->Update();

  c1->cd(3);
  h_offset[1][3][1011]->Draw();
  c1->Update();

  c1->cd(4);
  h_offset[1][3][1012]->Draw();
  c1->Update();


  psf1->Close();

  fclose(fpin);

  for( int i = 0; i < 4; i++){
    sprintf( title, "v1740_bd1_group_%d_offset.txt", i);

    fpin = fopen( title, "w");

    for( int k = 0; k < 1024; k++){      
      for( int j = 0; j < 9; j++){      
	fprintf( fpin, "%8.2lf ", off_mean[i][j][k]);       
	//	printf("%5d  %8.4f\n", j, off_mean[0][0][j]);
      }
      fprintf( fpin, "\n");       
    }

    fclose(fpin);
  }


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
