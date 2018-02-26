#include "functions.h"
#include "TSpectrum.h"



float DummyFunc()
{
     return 1;
}


std::vector<std::pair<unsigned int,unsigned int> > CheckDarkCounts(float* vals)
{
  std::vector<std::pair<unsigned int,unsigned int> > ret;
  
  
  //-------------
  // count pulses
  int nPulses = 0;
  for(unsigned int sIt = 0; sIt < nS; ++sIt)
  {
    // begin of pulse
    if( vals[sIt] < 0.85 )
    {
      ++nPulses;
      std::pair<unsigned int,unsigned int> dummy;
      dummy.first = sIt;
      
      for(unsigned int sIt2 = sIt; sIt2 < nS; ++sIt2)
      {
        // end of pulse
        if( vals[sIt2] > 0.85 )
        {
          sIt = sIt2;
          dummy.second = sIt2;
          ret.push_back(dummy);
          break;
        } // end of pulse
      }
    } // begin of pulse
  }
  
  return ret;
}
/*
void CountDarkCounts (float* vals)
{


   Int_t npeaks = 10;
   int ampNS = 1024;
   TH1F * hPulse = new TH1F ("hPulse", "hPulse", ampNS, 0, ampNS);

   for (int iS = 0; iS < ampNS; iS++) hPulse->SetBinContent(iS+1, (4100-vals[iS])/1000);

   TSpectrum *s = new TSpectrum(2*npeaks);
   Int_t nfound = s->Search(hPulse, 2,"",0.10);
   printf("Found %d candidate peaks to fit\n",nfound);
   
   //Estimate background using TSpectrum::Background
   TH1 *hb = s->Background(hPulse,20,"same");
//   if (hb) c1->Update();
//   if (np <0) return;
   return;

}*/

void AnalyzeDarkCounts(float* vals, unsigned int& pedSMin, unsigned int& pedSMax, unsigned int& ampSMin, unsigned int& ampNS)
{
  std::vector<std::pair<unsigned int,unsigned int> > darkCounts = CheckDarkCounts(vals);
  
  if( darkCounts.size() > 1 )
  {
    unsigned int maxDurationIt = 0;
    float maxDuration = -999999999;
    for(unsigned int dcIt = 0; dcIt < darkCounts.size(); ++dcIt)
    {
      std::pair<unsigned int,unsigned int> dummy = darkCounts.at(dcIt);
      if( (dummy.second-dummy.first) > maxDuration  )
      {
        maxDuration = dummy.second-dummy.first;
        maxDurationIt = dcIt;
      }
    }
    
    if( maxDurationIt > 0 )
    {
      if( (darkCounts.at(maxDurationIt).first-darkCounts.at(maxDurationIt-1).second) > 30 )
      {
        pedSMin = darkCounts.at(maxDurationIt-1).second + 15;
        pedSMax = darkCounts.at(maxDurationIt).first - 15;
      }
      else if( (darkCounts.at(maxDurationIt).first-darkCounts.at(maxDurationIt-1).second) > 20 )
      {
        pedSMin = darkCounts.at(maxDurationIt-1).second + 10;
        pedSMax = darkCounts.at(maxDurationIt).first - 10;
      }
      else if( (darkCounts.at(maxDurationIt).first-darkCounts.at(maxDurationIt-1).second) > 10 )
      {
        pedSMin = darkCounts.at(maxDurationIt-1).second + 5;
        pedSMax = darkCounts.at(maxDurationIt).first - 5;
      }
      else
      {
        pedSMin = darkCounts.at(maxDurationIt-1).second;
        pedSMax = darkCounts.at(maxDurationIt).first;
      }
      
      if( (darkCounts.at(maxDurationIt).second-darkCounts.at(maxDurationIt).first) > 30 )
      {
        ampSMin = darkCounts.at(maxDurationIt).first + 15;
        ampNS = darkCounts.at(maxDurationIt).second - darkCounts.at(maxDurationIt).first - 30;
      }
      else if( (darkCounts.at(maxDurationIt).second-darkCounts.at(maxDurationIt).first) > 20 )
      {
        ampSMin = darkCounts.at(maxDurationIt).first + 10;
        ampNS = darkCounts.at(maxDurationIt).second - darkCounts.at(maxDurationIt).first - 20;
      }
      else if( (darkCounts.at(maxDurationIt).second-darkCounts.at(maxDurationIt).first) > 10 )
      {
        ampSMin = darkCounts.at(maxDurationIt).first + 5;
        ampNS = darkCounts.at(maxDurationIt).second - darkCounts.at(maxDurationIt).first - 10;
      }
      else
      {
        ampSMin = darkCounts.at(maxDurationIt).first;
        ampNS = darkCounts.at(maxDurationIt).second - darkCounts.at(maxDurationIt).second;
      }
    }
    else
    {
      if( darkCounts.at(0).first > 0 )
      {
        pedSMin = 0;
        pedSMax = std::max(5,std::max(0,int(darkCounts.at(0).first-15)));
        
        if( (darkCounts.at(0).second-darkCounts.at(0).first) > 30 )
        {
          ampSMin = darkCounts.at(0).first + 15;
          ampNS = darkCounts.at(0).second - darkCounts.at(0).first - 30;
        }
        else if( (darkCounts.at(0).second-darkCounts.at(0).first) > 20 )
        {
          ampSMin = darkCounts.at(0).first + 10;
          ampNS = darkCounts.at(0).second - darkCounts.at(0).first - 20;
        }
        else if( (darkCounts.at(0).second-darkCounts.at(0).first) > 10 )
        {
          ampSMin = darkCounts.at(0).first + 5;
          ampNS = darkCounts.at(0).second - darkCounts.at(0).first - 10;
        }
        else
        {
          ampSMin = darkCounts.at(0).first;
          ampNS = darkCounts.at(0).second - darkCounts.at(0).second;
        }
      }
      else
      {
        pedSMin = darkCounts.at(0).second + 15;
        pedSMax = darkCounts.at(0).second + 15 + 50;
        
        ampSMin = 0;
        ampNS = std::max(int(darkCounts.at(0).second-15),1);
      }
    }        
  }
  else if( darkCounts.size() == 1 )
  {
    if( darkCounts.at(0).first > 0 )
    {
      pedSMin = 0;
      pedSMax = std::max(5,std::max(0,int(darkCounts.at(0).first-15)));
      
      if( (darkCounts.at(0).second-darkCounts.at(0).first) > 30 )
      {
        ampSMin = darkCounts.at(0).first + 15;
        ampNS = darkCounts.at(0).second - darkCounts.at(0).first - 30;
      }
      else if( (darkCounts.at(0).second-darkCounts.at(0).first) > 20 )
      {
        ampSMin = darkCounts.at(0).first + 10;
        ampNS = darkCounts.at(0).second - darkCounts.at(0).first - 20;
      }
      else if( (darkCounts.at(0).second-darkCounts.at(0).first) > 10 )
      {
        ampSMin = darkCounts.at(0).first + 5;
        ampNS = darkCounts.at(0).second - darkCounts.at(0).first - 10;
      }
      else
      {
        ampSMin = darkCounts.at(0).first;
        ampNS = darkCounts.at(0).second - darkCounts.at(0).second;
      }
    }
    else
    {
      pedSMin = darkCounts.at(0).second + 15;
      pedSMax = darkCounts.at(0).second + 15 + 50;
      
      ampSMin = 0;
      ampNS = std::max(int(darkCounts.at(0).second-15),1);
    }
  }
}



float GetPedestal(float* vals, const unsigned int& s1, const unsigned int& s2)
{
  float ped = 0.;
  
  for(unsigned int sIt = s1; sIt < s2; ++sIt)
  {
    ped += vals[sIt];
  }
  
  return 1. * ped / (s2-s1);
}



std::pair<std::pair<float,std::pair<float,float> >,float> GetAmplitudePulse(float* vals, const unsigned int& s1, const unsigned int& s2, const bool& isNegative)
{
  float ped = GetPedestal(vals,s1,s2);
  
  float comp = -999999999.;
  if( isNegative ) comp = 999999999;
  
  float timeAmp = -1;
  
  float integral = 0.;
  for(unsigned int sIt = s2; sIt < nS; ++sIt)
  {
    if( !isNegative )
    {
      integral += vals[sIt] - ped;
      
      if( vals[sIt] > comp )
      {
        comp = vals[sIt];
        timeAmp = sIt * 1. / GS_s;
      }
    }
    else
    {
      integral += -1. * (vals[sIt] - ped);
      
      if( vals[sIt] < comp )
      {
        comp = vals[sIt];
        timeAmp = sIt * 1. / GS_s;
      }
    }
  }

//  std::cout << " comp = " << comp << " :: ped = " << ped << " :: integral = " << integral << std::endl;
  
  std::pair<float,float> ret(fabs(comp-ped),integral);
  std::pair<float,std::pair<float,float> > ret2(ped,ret);
  std::pair<std::pair<float,std::pair<float,float> >,float> ret3(ret2,timeAmp);
  return ret3;
}

std::pair<std::pair<float,std::pair<float,float> >,float> GetAmplitudeIntegral(float* vals, const unsigned int& s1, const unsigned int& s2, const bool& isNegative)
{
  float ped = GetPedestal(vals,s1,s2);
  
  float timeAmp = -1;
  
  float integral_half  = 0.;
  float integral_short = 0.;

  for(unsigned int sIt = s2; sIt < (s2+200); ++sIt)
  {
    if( !isNegative )
    {
      if (sIt < (s2+200))  integral_half  += vals[sIt] - ped;      
      if (sIt < (s2+100))  integral_short += vals[sIt] - ped;      
    }
    else
    {
      if (sIt < (s2+200))  integral_half  += -1.*(vals[sIt] - ped);      
      if (sIt < (s2+100))  integral_short += -1.*(vals[sIt] - ped);          
    }
  }

//  std::cout << " comp = " << comp << " :: ped = " << ped << " :: integral = " << integral << std::endl;
  
  std::pair<float,float> ret(integral_half,integral_short);
  std::pair<float,std::pair<float,float> > ret2(ped,ret);
  std::pair<std::pair<float,std::pair<float,float> >,float> ret3(ret2,timeAmp);
  return ret3;
}



std::pair<float,float> GetAmplitudeSquare(float* vals, const unsigned int& s1, const unsigned int& s2, const float& SqTH, const unsigned int& nS, const bool& isNegative)
{
  float ped = GetPedestal(vals,s1,s2);
  float threshold = SqTH*1000;
  int sMin = 0;
  
//   threshold = 100;
  
  for(unsigned int sIt = 0; sIt < 1024; ++sIt)
  {
     if( isNegative)
     {
//        std::cout << "vals = " << vals[sIt] << " :: ped = "  << ped << " :: threshold = " << threshold << std::endl;
       if( fabs(vals[sIt]-ped) > threshold )
       {
         sMin = sIt;         
         break;
       }
     }
  }
  
  sMin+=10;
  
  
  
  float amp = GetPedestal(vals,sMin,sMin+nS);
  //std::cout << "sMin = " << sMin << " :: amp = " << amp << std::endl;
  
  std::pair<float,float> ret(ped,fabs(amp-ped));
  return ret;
}



std::pair<std::pair<float,float>,std::pair<TF1*,TF1*> > GetTimeTh(float* vals, const unsigned int& s1, const float& ped, const float& amp,
                                                                  const float& th, const unsigned int& nSBef, const unsigned int& nSAft,
                                                                  const bool& isFraction, const bool& isNegative, const bool& deleteFit)
{
  float threshold = th;
  if( isFraction )
    threshold = amp * th;
  
  unsigned int sBef1 = 0;
  unsigned int sBef2 = 0;
  for(unsigned int sIt1 = s1; sIt1 < nS; ++sIt1)
  {
     if( isNegative )
     {
       if( vals[sIt1] < (ped-threshold) )
       {
         sBef1 = std::max(0,int(sIt1-1));
         
         for(unsigned int sIt2 = sIt1; sIt2 < nS; ++sIt2)
         {
           if( vals[sIt2] > (ped-threshold) )
           {
             sBef2 = sIt2 - 1;
             break;
           }
         }
         
         break;
       }
     }
     
     else
     {
       if( vals[sIt1] > (ped+threshold) )
       {
         sBef1 = std::max(0,int(sIt1-1));
         
         for(unsigned int sIt2 = sIt1; sIt2 < nS; ++sIt2)
         {
           if( vals[sIt2] < (ped+threshold) )
           {
             sBef2 = sIt2 - 1;
             break;
           }
         }
         
         break;
       }
     }
   }
//    std::cout << "sBef1 = " << sBef1 << " :: sBef2 = " << sBef2 << std::endl;
  
  std::pair<float,float> ret1(0,0);
  std::pair<TF1*,TF1*> ret2(NULL,NULL);
  std::pair<std::pair<float,float>,std::pair<TF1*,TF1*> > ret(ret1,ret2);
  
  
  if( sBef1 > 0 )
  {
    TGraph* g = PointsAroundTh(vals,sBef1,nSBef,nSAft);
    g -> SetMarkerColor(kBlue+1);
    g -> SetMarkerStyle(21);
    g -> GetXaxis()->SetLimits(0,1024);
    g -> GetYaxis()->SetLimits(0,1);
    g -> Fit("pol1","Q");
    //g -> Draw("ALPE");
    
    TF1* fitFunc = (TF1*)( g->GetFunction("pol1") );
    fitFunc -> SetLineColor(kMagenta);
    fitFunc -> SetLineWidth(1);
    
    float time = -1.;
    if( isNegative ) time = ( (ped-threshold) - fitFunc->GetParameter(0) ) / fitFunc->GetParameter(1);
    else             time = ( (ped+threshold) - fitFunc->GetParameter(0) ) / fitFunc->GetParameter(1);
    
    ret.first.first = time;
    ret.second.first = fitFunc;
    
    if( deleteFit )
    {
      delete fitFunc;
      delete g;
    }
  }
  
  if( sBef2 > 0 )
  {
    TGraph* g = PointsAroundTh(vals,sBef2,nSAft,nSBef);
    g -> Fit("pol1","Q");
    
    TF1* fitFunc = (TF1*)( g->GetFunction("pol1") );
    fitFunc -> SetLineColor(kMagenta);
    fitFunc -> SetLineWidth(1);
    
    float time = -1.;
    if( isNegative ) time = ( (ped-threshold) - fitFunc->GetParameter(0) ) / fitFunc->GetParameter(1);
    else             time = ( (ped+threshold) - fitFunc->GetParameter(0) ) / fitFunc->GetParameter(1);
    
    ret.first.second = time;
    ret.second.second = fitFunc;
    
    if( deleteFit )
    {
      delete fitFunc;
      delete g;
    }
  }
  
  return ret;
}

std::pair<std::pair<float,float>,std::pair<TF1*,TF1*> > GetTimeThWeighed(float* vals, float* times, const unsigned int& s1, const float& ped, const float& amp,
                                                                  const float& th, const unsigned int& nSBef, const unsigned int& nSAft,
                                                                  const bool& isFraction, const bool& isNegative, const bool& deleteFit)
{
  float threshold = th;
  if( isFraction )
    threshold = amp * th;
  
  unsigned int sBef1 = 0;
  unsigned int sBef2 = 0;
  for(unsigned int sIt1 = s1; sIt1 < nS; ++sIt1)
  {
     if( isNegative )
     {
       if( vals[sIt1] < (ped-threshold) )
       {
         sBef1 = std::max(0,int(sIt1-1));
         
         for(unsigned int sIt2 = sIt1; sIt2 < nS; ++sIt2)
         {
           if( vals[sIt2] > (ped-threshold) )
           {
             sBef2 = sIt2 - 1;
             break;
           }
         }
         
         break;
       }
     }
     
     else
     {
       if( vals[sIt1] > (ped+threshold) )
       {
         sBef1 = std::max(0,int(sIt1-1));
         
         for(unsigned int sIt2 = sIt1; sIt2 < nS; ++sIt2)
         {
           if( vals[sIt2] < (ped+threshold) )
           {
             sBef2 = sIt2 - 1;
             break;
           }
         }
         
         break;
       }
     }
   }
//    std::cout << "sBef1 = " << sBef1 << " :: sBef2 = " << sBef2 << std::endl;
  
  std::pair<float,float> ret1(0,0);
  std::pair<TF1*,TF1*> ret2(NULL,NULL);
  std::pair<std::pair<float,float>,std::pair<TF1*,TF1*> > ret(ret1,ret2);
  
  
  if( sBef1 > 0 )
  {
    TGraph* g = PointsAroundThWeighed(vals,times,sBef1,nSBef,nSAft);
    g -> SetMarkerColor(kBlue+1);
    g -> SetMarkerStyle(21);
    g -> GetXaxis()->SetLimits(0,1024);
    g -> GetYaxis()->SetLimits(0,1);
    g -> Fit("pol1","Q");
    //g -> Draw("ALPE");
    
    TF1* fitFunc = (TF1*)( g->GetFunction("pol1") );
    fitFunc -> SetLineColor(kMagenta);
    fitFunc -> SetLineWidth(1);
    
    float time = -1.;
    if( isNegative ) time = ( (ped-threshold) - fitFunc->GetParameter(0) ) / fitFunc->GetParameter(1);
    else             time = ( (ped+threshold) - fitFunc->GetParameter(0) ) / fitFunc->GetParameter(1);
    
    ret.first.first = time;
    ret.second.first = fitFunc;
    
    if( deleteFit )
    {
      delete fitFunc;
      delete g;
    }
  }
  
  if( sBef2 > 0 )
  {
    TGraph* g = PointsAroundThWeighed(vals,times,sBef2,nSAft,nSBef);
    g -> Fit("pol1","Q");
    
    TF1* fitFunc = (TF1*)( g->GetFunction("pol1") );
    fitFunc -> SetLineColor(kMagenta);
    fitFunc -> SetLineWidth(1);
    
    float time = -1.;
    if( isNegative ) time = ( (ped-threshold) - fitFunc->GetParameter(0) ) / fitFunc->GetParameter(1);
    else             time = ( (ped+threshold) - fitFunc->GetParameter(0) ) / fitFunc->GetParameter(1);
    
    ret.first.second = time;
    ret.second.second = fitFunc;
    
    if( deleteFit )
    {
      delete fitFunc;
      delete g;
    }
  }
  
  return ret;
}


std::vector <int>  GetSpikes(float* vals, const float& spikeAmp)
{
    
    float ped = GetPedestal(vals,10,100);
    int spikeId = 0;
    std::vector<int> spikeIt;    
    std::vector<float> spikeMax;
//     std::cout << "start spikes hunting --> spikeAmp = " << spikeAmp << " :: ped = " << ped << std::endl;
  
    for(unsigned int sIt = 0; sIt < nS; ++sIt)
    {
//         std::cout << "sIt = " << sIt << " :: vals = " << vals[sIt] << std::endl;
        
        if( vals[sIt] > (ped + spikeAmp*0.5*4000)  || vals[sIt] < (ped - spikeAmp*0.5*4000)  )
        {
            
            spikeIt.push_back(sIt);
            spikeMax.push_back(vals[sIt]-ped);
            std::cout << " spikeId = " << spikeId << " :: spikeIt = " << spikeIt.at(spikeId) << " :: spikeMax = " << spikeMax.at(spikeId) << std::endl;
            
            spikeId++;

        }
    }

    
  
    // return the number of spikes found and their position
//     std::vector <int> ret(spikeId,spikeIt);
    return spikeIt;
}



/*
std::pair<std::pair<float,float>,TF1*> GetTimeTemplateFit(float* vals, const TProfile* p_template, const float& ped, const bool& draw)
{
  std::pair<float,float> ret1(0,0);
  std::pair<std::pair<float,float>,TF1*> ret(ret1,NULL);
  
  
  TGraphErrors* graph = new TGraphErrors();
  for(unsigned int sIt = 0; sIt < nS; ++sIt)
   {
     graph -> SetPoint(sIt,sIt/GS_s,vals[sIt]);
     graph -> SetPointError(sIt,0.,0.0005);
   }
  
  
  histoFunc* templateHistoFunc = new histoFunc((TH1F*)(p_template));
  TF1* fitFunc = new TF1("fitFunc",templateHistoFunc,0.,300.,4,"histoFunc");
  fitFunc -> SetNpx(10000);
  fitFunc -> SetParLimits(0,0.,0.5);
  fitFunc -> SetParLimits(1,0.,50.);
  fitFunc -> SetParLimits(3,0.,100.);
  fitFunc -> SetParameter(0,ped);
  fitFunc -> SetParameter(1,10.);
  fitFunc -> FixParameter(2,1.);
  fitFunc -> SetParameter(3,30.);
  fitFunc -> SetLineWidth(2);
  fitFunc -> SetLineColor(kRed);
  fitFunc -> SetParName(0,"c");
  fitFunc -> SetParName(1,"N");
  fitFunc -> SetParName(2,"k");
  fitFunc -> SetParName(3,"t_{0}");
  graph -> Fit("fitFunc","QNS+");
  
  ret.first.first = fitFunc -> GetParameter(1);
  ret.first.second = fitFunc -> GetParameter(3);
  ret.second = fitFunc;
  
  if( draw )
  {
    graph -> SetMarkerSize(0.5);
    graph -> Draw("P,same");
    fitFunc -> Draw("same");
  }
  else
  {
    delete graph;
  }
  return ret;
}
*/

TGraph* PointsAroundTh(float* vals, const unsigned int& sBef, const unsigned int& nSBef, const unsigned int& nSAft)
{
  // n-point interpolation: nSBef points before threshold, BSAft points after threshold
  TGraph* g = new TGraph();

  int point = 0;
  for(unsigned int iS = 0; iS < nSBef; ++iS)
  {
    unsigned int sample = sBef - (nSBef-iS-1);
    float time = sample * 1. / GS_s;
    
    g -> SetPoint(point,time,vals[sample]);
    ++point;
  }
  for(unsigned int iS = 0; iS < nSAft; ++iS)
  {
    unsigned int sample = sBef + 1 +iS;
    float time = sample * 1. / GS_s;
    
    g -> SetPoint(point,time,vals[sample]);
    ++point;
  }
  
  return g;
}

TGraph* PointsAroundThWeighed(float* vals, float* times, const unsigned int& sBef, const unsigned int& nSBef, const unsigned int& nSAft)
{
  // n-point interpolation: nSBef points before threshold, BSAft points after threshold
  TGraph* g = new TGraph();

  int point = 0;
  for(unsigned int iS = 0; iS < nSBef; ++iS)
  {
    unsigned int sample = sBef - (nSBef-iS-1);
    g -> SetPoint(point,times[sample],vals[sample]);
    ++point;
  }
  for(unsigned int iS = 0; iS < nSAft; ++iS)
  {
    unsigned int sample = sBef + 1 +iS;
    g -> SetPoint(point,times[sample],vals[sample]);
    ++point;
  }
  
  return g;
}
