#include "Aux.hh"
#include <math.h>
#include <string>

//*********************************************************
// Get amplification factor used for the silicon sensor
//*********************************************************
double GetAmplificationFactor ( double measuredAmplitude ) {
  
  int index_firstBinAboveInput = -1;
  for (int i=0; i < nPoints; ++i) {
    index_firstBinAboveInput = i;
    if (measuredAmplitude < outputAmplitude[i]) break;
  }
  
  double answer = 0; 

  if (measuredAmplitude > outputAmplitude[21]) answer =amplificationFactor[21];
  else if (index_firstBinAboveInput == 0) answer = amplificationFactor[0];
  else {
    
    //cout << "index_firstBinAboveInput = " << index_firstBinAboveInput << " : "
    //	 << amplificationFactor[index_firstBinAboveInput-1] << " " << outputAmplitude[index_firstBinAboveInput]
    //	 << "\n";
    double x = measuredAmplitude - outputAmplitude[index_firstBinAboveInput-1];
    double y = amplificationFactor[index_firstBinAboveInput-1] + x * (amplificationFactor[index_firstBinAboveInput] - amplificationFactor[index_firstBinAboveInput-1]) / (outputAmplitude[index_firstBinAboveInput] - outputAmplitude[index_firstBinAboveInput-1]);
    //cout << "x = " << x << " , y = " << y << "\n";
    answer = y;
  }

  //cout << measuredAmplitude << " " << answer << "\n";

  return answer;
  
}


TGraphErrors* GetTGraph(  float* channel, float* time )
{		
  //Setting Errors
  float errorX[1024], errorY[1024], channelFloat[1024];
  float _errorY = 0.00; //5%error on Y
  for ( int i = 0; i < 1024; i++ )
    {
      errorX[i]       = .0;
      errorY[i]       = _errorY*channel[i];
      channelFloat[i] = -channel[i];
    }
  TGraphErrors* tg = new TGraphErrors( 1024, time, channelFloat, errorX, errorY );
  return tg;
};


TGraphErrors GetTGraph(  short* channel, float* time )
{		
  //Setting Errors
  float errorX[1024], errorY[1024], channelFloat[1024];
  float _errorY = 0.00; //5%error on Y
  for ( int i = 0; i < 1024; i++ )
    {
      errorX[i]       = .0;
      errorY[i]       = _errorY*channel[i];
      channelFloat[i] = -channel[i];
    }
  //TGraphErrors* tg = new TGraphErrors( 1024, time, channelFloat, errorX, errorY );
  TGraphErrors tg( 1024, time, channelFloat, errorX, errorY );
  return tg;
};

////////////////////////////////////////////
// find minimum of the pulse
////////////////////////////////////////////
int FindMin( int n, short *a) {
  return FindRealMin(n,a);
}

int FindMinAbsolute( int n, short *a) {
  
  if (n <= 0 || !a) return -1;
  float xmin = a[5];
  int loc = 0;
  for  (int i = 5; i < n-10; i++) {
    if (xmin > a[i] && a[i+1] < 0.5*a[i] && a[i] < -40. )  
      {
	//std::cout << i << " " << a[i] << std::endl;
	xmin = a[i];
	loc = i;
	//if ( a[i+5]>a[i] && a[i+10]>a[i+5] ) {
	//break;
      }
  }
  //std::cout << "loc0: " << loc << std::endl;
  return loc;
}

int FindRealMin( int n, short *a) {  
  if (n <= 0 || !a) return -1;
  float xmin = a[5];
  int loc = 0;
  
  float noise = 0;
  
  for ( int i = 5; i < 100; i++)
    {
      if( abs(a[i]) > noise ) 
	{
	  noise = abs(a[i]);
	}
    }

  for  (int i = 5; i < n-10; i++) {
    if (xmin > a[i] && a[i+1] < 0.5*a[i] && a[i] < -3*noise && a[i] < -50.)  
      {
	//std::cout << a[i] << std::endl;
	xmin = a[i];
	loc = i;
	//if ( a[i+5]>a[i] && a[i+10]>a[i+5] ) {
	//break;
      }
  }
 
  float xmin_init = xmin;
  float xmin_new = a[5];
  int loc_new = loc;

  bool stop = false;
  while( !stop )
    {
      for ( int i = 5; i < loc_new -25; i++ )
        {
          if ( a[i] < xmin_new && 0.5*a[i] > a[i+1] && a[i] < 0.15* xmin_init )
            {
              xmin_new = a[i];
              loc_new = i;
            }
        }

      xmin_init = xmin_new;

      if( loc_new == loc ) break;
      //std::cout << "new peak @ " << loc_new << ", ymin: " << xmin_new << std::endl;
      if ( xmin_new > -2*noise || xmin_new > -40 ) loc_new = 0;
      xmin_new = a[5];
      loc = loc_new;
    }

  //std::cout << "LOC2: " << loc << std::endl;                                                                                                                               
  /*                                                                
  while ( xmin_init != xmin_new ) {
    for (int i = 5; i < loc - 50; i++) {
      if (xmin_new > a[i] && a[i+1] < 0.5*a[i] && a[i] < xmin_init*2/3 )  {
        xmin_new = a[i];
        loc = i;
      }
    }
    xmin_init = xmin_new
    xmin_new = a[5]
  }
  */
  return loc_new;
}

int FindMinFirstPeakAboveNoise( int n, short *a) {
  
  const float noise = 50;
  if (n <= 0 || !a) return -1;
  int loc = 0;
  
  for  (int i = 10; i < n-10; i++) {   
    if ( abs(a[i]) > noise 
	&& 
	(a[i] < a[i-1] && a[i] < a[i-2] && a[i] < a[i-3] && a[i] < a[i-4] && a[i] < a[i-5] )
	&&
	(a[i] < a[i+1] && a[i] < a[i+2] && a[i] < a[i+3] && a[i] < a[i+4] && a[i] < a[i+5] )
	)  {
      loc = i;
      break;
    }
  }  
  return loc;
}


// find the mean time from gaus fit
float GausFit_MeanTime(TGraphErrors* pulse, const float index_first, const float index_last)
{
  TF1* fpeak = new TF1("fpeak","gaus", index_first, index_last);
  pulse->Fit("fpeak","Q","", index_first, index_last);
  
  float timepeak = fpeak->GetParameter(1);
  delete fpeak;
  
  return timepeak;
};

float GausFit_MeanTime(TGraphErrors* pulse, const float index_first, const float index_last, TString fname)
{
  TF1* fpeak = new TF1("fpeak","gaus", index_first, index_last);
  //float max = pulse->GetMaximum();
  double max = -9999;
  double* yy = pulse->GetY();
  for ( int i = 0; i < 1024; i++ )
    {
      if ( yy[i] > max ) max = yy[i];
    }
  //std::cout << "max: " << max << std::endl;
  if( max < 42 || index_first < 10 || index_last > 1010 ) return -99999;
  pulse->Fit("fpeak","Q","", index_first, index_last);
  
  TCanvas* c = new TCanvas("canvas","canvas",800,400) ;
  float timepeak = fpeak->GetParameter(1);
  pulse->GetXaxis()->SetLimits( timepeak-10, timepeak+10);
  pulse->SetMarkerSize(0.5);
  pulse->SetMarkerStyle(20);
  pulse->Draw("AP");
  c->SaveAs(fname+".pdf");
  delete fpeak;
  
  return timepeak;
};

float RisingEdgeFitTime(TGraphErrors * pulse, const float index_min, const float constantFraction, TString fname, bool makePlot )
{
  double x_low, x_high, y, dummy;
  //pulse->GetPoint(index_min-7, x_low, y);
  //pulse->GetPoint(index_min-2, x_high, y);
  pulse->GetPoint(index_min-12, x_low, y);
  pulse->GetPoint(index_min-7, x_high, y);
  pulse->GetPoint(index_min, dummy, y);  
  TF1* flinear = new TF1("flinear","[0]*x+[1]", x_low, x_high );
  
  pulse->Fit("flinear","Q","", x_low, x_high );
  double slope = flinear->GetParameter(0);
  double b     = flinear->GetParameter(1);
  
  if ( makePlot )
    {
      std::cout << "make plot" << std::endl;
      TCanvas* c = new TCanvas("canvas","canvas",800,400) ;
      pulse->GetXaxis()->SetLimits(x_low-3, x_high+3);
      pulse->SetMarkerSize(0.3);
      pulse->SetMarkerStyle(20);
      pulse->Draw("AP");
      c->SaveAs(fname+"LinearFit.pdf");
      //delete c;
    }
  delete flinear;
  
  return (constantFraction*y-b)/slope;
};


void RisingEdgeFitTime(TGraphErrors * pulse, const float index_min, float* tstamp, int event, TString fname, bool makePlot )
{
  double x_low, x_high, y, dummy;
  double ymax;
  pulse->GetPoint(index_min, x_low, ymax);
  for ( int i = 1; i < 100; i++ )
    {
      pulse->GetPoint(index_min-i, x_low, y);
      if ( y < 0.2*ymax ) break;
    }
  for ( int i = 1; i < 100; i++ )
    {
      pulse->GetPoint(index_min-i, x_high, y);
      if ( y < 0.9*ymax ) break;
    }
  //pulse->GetPoint(index_min-8, x_low, y);
  //pulse->GetPoint(index_min-3, x_high, y);


  //pulse->GetPoint(index_min-12, x_low, y);
  //pulse->GetPoint(index_min-7, x_high, y);
  pulse->GetPoint(index_min, dummy, y);
  
  TF1* flinear = new TF1("flinear","[0]*x+[1]", x_low, x_high );
  float max = -9999;
  double* yy = pulse->GetY();
  
  for ( int i = 0; i < 1024; i++ )
    {
      if ( yy[i] > max ) max = yy[i];
    }
  //std::cout << "max: " << max << std::endl;

  /*if( max < 10 || index_min < 0 || index_min > 1023 )
    {
      std::cout << "DEB: skipping event--> " << event << std::endl;
      return;
    }
  */
  pulse->Fit("flinear","Q","", x_low, x_high );
  double slope = flinear->GetParameter(0);
  double b     = flinear->GetParameter(1);
  
  if ( makePlot )
    {
      std::cout << "make plot" << std::endl;
      TCanvas* c = new TCanvas("canvas","canvas",800,400) ;
      pulse->GetXaxis()->SetLimits(x_low-10, x_high+10);
      pulse->SetMarkerSize(0.3);
      pulse->SetMarkerStyle(20);
      pulse->Draw("AP");
      c->SaveAs(fname+"LinearFit.pdf");
      //delete c;
    }
  tstamp[0] = (0.0*y-b)/slope;
  //std::cout << "----" << tstamp[0]  << std::endl;
  tstamp[1] = (0.15*y-b)/slope;
  tstamp[2] = (0.30*y-b)/slope;
  tstamp[3] = (0.45*y-b)/slope;
  tstamp[4] = (0.60*y-b)/slope;
  
  delete flinear;
};


double GetGaussTime( TGraphErrors* pulse )
{
  return 0;
};


float GetBaseline(TGraphErrors * pulse, int i_low, int i_high, TString fname )
{
  double x_low, x_high, y, dummy;
  pulse->GetPoint(i_low, x_low, y);
  pulse->GetPoint(i_high, x_high, y);
  
  TF1* flinear = new TF1("flinear","[0]", x_low, x_high );
  
  pulse->Fit("flinear","RQ","", x_low, x_high );
  
  /* std::cout << "make plot" << std::endl;
  std::cout << x_low << x_high << fname << std::endl;
  TCanvas* c = new TCanvas("canvas","canvas",800,400) ;
  pulse->GetXaxis()->SetLimits(x_low-3, x_high+3);
  pulse->SetMarkerSize(1);
  pulse->SetMarkerStyle(20);
  pulse->Draw("AP");
  c->SaveAs(fname+"LinearFit.pdf"); */
  
  float a = flinear->GetParameter(0);
  delete flinear;
  
  return a;
}

float GetBaseline( int peak, short *a ) {

  float tmpsum = 0;
  float tmpcount = 0;
  //std::cout << "GGG\n";
  if (peak < 300) {
    for  (int i = peak + 200; i < 1000; i++) {
      // std::cout << i << " : " << a[i] << "\n";
      tmpsum += a[i];
      tmpcount += 1.0;
    }
  } else {
    for  (int i = 5; i < peak-200; i++) {
      // std::cout << i << " : " << a[i] << "\n";
      tmpsum += a[i];
      tmpcount += 1.0;
    }
  }
  // std::cout << tmpsum / tmpcount << "\n";

  return tmpsum / tmpcount;
}


float GetPulseIntegral(int peak, short *a, std::string option) 
{
  float integral = 0.;

  if (option == "full") {
    for (int i=5; i < 1020; i++) {
      integral += a[i] * 0.2 * 1e-9 * (1.0/4096.0) * (1.0/50.0) * 1e12; //in units of pC, for 50Ohm termination
    }
  }
  else {
    for (int i=peak-20; i < peak+25; i++) {
      integral += a[i] * 0.2 * 1e-9 * (1.0/4096.0) * (1.0/50.0) * 1e12; //in units of pC, for 50Ohm termination
    }
  }
  return -1.0 * integral;

}

TGraphErrors* GetTGraphFilter( short* channel, float* time, TString pulseName, bool makePlot )
{
  float Gauss[1024];
  //Setting Errors
  float errorX[1024], errorY[1024], channelFloat[1024];
  float _errorY = 0.00; //5%error on Y
  
  for ( int i = 0; i < 1024; i++ )
    {
      
      errorX[i]       = .0;
      errorY[i]       = _errorY*channel[i];
      channelFloat[i] = -channel[i];
    }
  
  TF1 *fb = new TF1("fb","gaus(0)", 0.0, 204.6);
  fb->SetParameter(1, 100);
  float sigma =1.0;
  fb->SetParameter(2, sigma);
  fb->SetParameter(0, 1/(sqrt(3.1415*2.0)*sigma) );
  //eval Gaussian
  float step = 0.2;//200ps
  for ( int i = 0; i < 1024; i++ )
    {
      Gauss[i] = fb->Eval( float(i)*step );
    }
  
  float channelFloatFiltered[2048];
  for ( int i = 0; i < 2048; i++ )
    {
      float convolvedPoint = 0;
      for ( int j = 0; j <= i; j++ )
	{
	  if ( i < 1024 )
	    {
	      convolvedPoint += channelFloat[i-j]*Gauss[1023-j];
	    }
	  else
	    {
	      if ( 1023-(i-1023)-j >= 0 ) convolvedPoint += channelFloat[1023-j]*Gauss[1023-(i-1023)-j];
	    }
	}
      //if ( i < 1024 ) channelFloatFiltered[i] = convolvedPoint;
      channelFloatFiltered[i] = convolvedPoint;
    }
  
  float channelFloatFilteredFix[1024];
  for ( int i = 0; i < 1024; i++ )
    {
      channelFloatFilteredFix[i] = 0.2*channelFloatFiltered[i+523];
    }
  
  TGraphErrors* tg = new TGraphErrors( 1024, time, channelFloat, errorX, errorY );
  TGraphErrors* tg2 = new TGraphErrors( 1024, time, channelFloatFilteredFix, errorX, errorY );

  if (makePlot) {
    TCanvas* c = new TCanvas("canvas","canvas",800,400) ;         
    tg2->GetXaxis()->SetLimits(50, 70);
    tg->GetXaxis()->SetLimits(50, 70);
    //tg2->Fit("fb","","", 0.0, 204.6 );
    tg2->SetMarkerSize(0.5);
    tg->SetMarkerSize(0.5);
    tg2->SetMarkerStyle(20);
    tg->SetMarkerStyle(20);
    tg2->Draw("AP");
    tg2->SetMarkerColor(kBlue);
    tg->Draw("sameP");
    c->SaveAs(pulseName + "GausPulse.pdf");
  }
  return tg;
};

