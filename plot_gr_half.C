int 
plot_gr_half(char *datafile){

  double off_mean[4][9][1024];
  double tcal_dV[4][1024];

  FILE* fp1;
  char stitle[200];
  int dummy;
  double fdummy;

  for( int i = 0; i < 4; i++){
    sprintf( stitle, "v1740_bd1_group_%d_offset.txt", i);

    fp1 = fopen( stitle, "r");
    printf("Amplitude offset data : %s\n", stitle);

    for( int k = 0; k < 1024; k++)      
      for( int j = 0; j < 9; j++){      
	dummy = fscanf( fp1, "%lf ", &off_mean[i][j][k]);       
	if( k < 2 && 0)
	  printf("%5d  %8.4f\n", j, off_mean[i][j][k]);
      }
  
    fclose(fp1);
  }

  for( int i = 0; i < 4; i++){
    sprintf( stitle, "v1740_bd1_group_%d_dV.txt", i);

    fp1 = fopen( stitle, "r");
    printf("Time calibration data : %s\n", stitle);

    for( int k = 0; k < 1024; k++)      
      dummy = fscanf( fp1, "%lf %lf %lf %lf %lf ", 
		      &fdummy, &fdummy, &fdummy, &fdummy, &tcal_dV[i][k]);       
    fclose(fp1);
  }

  double dV_sum[4] = {0, 0, 0, 0};
  for( int i = 0; i < 4; i++)
    for( int j = 0; j < 1024; j++)
      dV_sum[i] += tcal_dV[i][j];


  double tcal[4][1024];
  for( int i = 0; i < 4; i++)
    for( int j = 0; j < 1024; j++)
      tcal[i][j] = tcal_dV[i][j]/dV_sum[i]*200.0;

  TGraphErrors* gr[4][9];

  TH2F* h2 = new TH2F("h2", "h2", 200, 0, 200, 400, 0, 4096);  
  h2->SetBit(TH1::kNoStats);

  TH1F* h_tc = new TH1F("h_tc", "Trigger cells", 1024, 0, 1023);
  
  TCanvas* c1 = new TCanvas("c1", "c1", 10, 10, 700, 500);
  c1->SetFillStyle(4000);
  c1->Divide(2, 2);
  
  char title[200];  
  sprintf( title, "plot-gr-half-%s.ps", datafile);
  TPostScript* psf1 = new TPostScript( title, 112);

  unsigned int event_header;
  unsigned int temp[3];
  unsigned int samples[9][1024];

  sprintf( title, "%s.dat", datafile);

  FILE* fpin = fopen( title, "r");

  // loop over event
  for( int event = 0; event < 20; event++){
    printf("---- event  %5d\n", event);

    dummy = fread( &event_header, sizeof(unsigned int), 1, fpin);  
    dummy = fread( &event_header, sizeof(unsigned int), 1, fpin);  
    dummy = fread( &event_header, sizeof(unsigned int), 1, fpin);  
    dummy = fread( &event_header, sizeof(unsigned int), 1, fpin);  

    double index[4][1024], err[1024], amplitude[4][9][1024];
    for( int i = 0; i < 1024; i++)
      err[i] = 0;
    
    for( int group = 0; group < 2; group++){

      dummy = fread( &event_header, sizeof(unsigned int), 1, fpin);        
      unsigned int tc = (event_header >> 20) & 0xfff;

      for( int i = 0; i < 1024; i++){
	if(i==0) 
	  index[group][i] = 0;
	else
	  index[group][i] = (tcal[group][(i-1+tc)%1024] + index[group][i-1]);
      }

      h_tc->Fill( tc);

      int nsample = (event_header & 0xfff)/3;

      for(int i = 0; i < nsample; i++){
	dummy = fread( &temp, sizeof(unsigned int), 3, fpin);  
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
	dummy = fread( &temp, sizeof(unsigned int), 3, fpin);  
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
	for( int j = 0; j < 1024; j++){
	  samples[i][j] -= off_mean[group][i][(j+tc)%1024];  
	  samples[i][j] += 2100;
	  amplitude[group][i][j] = samples[i][j];	
	}	
      }

      for( int i = 0; i < 9; i++){
	gr[group][i] = new TGraphErrors( 1024, index[group], amplitude[group][i], err, err);
	gr[group][i]->SetMarkerColor(4);
      }
      dummy = fread( &event_header, sizeof(unsigned int), 1, fpin);  
    }
  
    char stitle[200];
      
    h2->SetAxisRange(0, 2500, "y");
    h2->GetXaxis()->SetTitle("Time (ns)");
    h2->GetYaxis()->SetTitle("ADC count");

    psf1->NewPage();
    for( int i = 0; i < 1; i++){
      c1->cd(i+1);
      h2->Draw();
      sprintf(stitle, "Evenet %d :  Gr#0 Ch #%d", event, i );
      h2->SetTitle(stitle);
      gr[0][0]->Draw("p");
      c1->Update();
    }

    c1->cd(2);
    h2->Draw();
    sprintf(stitle, "Evenet %d :  Gr#0 Ch #%d ", event, 1);
    h2->SetTitle(stitle);
    gr[0][1]->Draw("p");
    c1->Update();

    c1->cd(3);
    h2->Draw();
    sprintf(stitle, "Evenet %d :  Gr#0 Ch #%d ", event, 2);
    h2->SetTitle(stitle);
    gr[0][2]->Draw("p");
    c1->Update();

    c1->cd(4);
    h2->Draw();
    sprintf(stitle, "Evenet %d :  Group Trigger 0 ", event);
    h2->SetTitle(stitle);
    gr[0][8]->Draw("p");
    c1->Update();
  }

  psf1->Close();

  fclose(fpin);

  return 0;
}

