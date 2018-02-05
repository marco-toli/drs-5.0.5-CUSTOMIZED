#ifndef __functions_h__
#define __functions_h__

// #include "histoFunc.h"

#include <iostream>

#include "TF1.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TProfile.h"
#include "TSpectrum.h"

const unsigned int nS = 1024;
const float GS_s = 5.12;



//void CountDarkCounts(float* vals);

float DummyFunc();

std::vector<std::pair<unsigned int,unsigned int> > CheckDarkCounts(float* vals);
void AnalyzeDarkCounts(float* vals, unsigned int& pedSMin, unsigned int& pedSMax, unsigned int& ampSMin, unsigned int& ampNS);

float GetPedestal(float* vals, const unsigned int& s1, const unsigned int& s2);
std::pair<std::pair<float,std::pair<float,float> >,float> GetAmplitudePulse(float* vals, const unsigned int& s1, const unsigned int& s2, const bool& isNegative);
std::pair<std::pair<float,std::pair<float,float> >,float> GetAmplitudeIntegral(float* vals, const unsigned int& s1, const unsigned int& s2, const bool& isNegative);

std::pair<float,float> GetAmplitudeSquare(float* vals, const unsigned int& s1, const unsigned int& s2, const unsigned int& sMin, const unsigned int& nS, const bool& isNegative);

std::pair<std::pair<float,float>,std::pair<TF1*,TF1*> > GetTimeTh(float* vals, const unsigned int& s1, const float& ped, const float& amp,
                                                                  const float& th, const unsigned int& nSBef, const unsigned int& nSAft,
                                                                  const bool& isFraction, const bool& isNegative, const bool& deleteFit = true);
// std::pair<std::pair<float,float>,TF1*> GetTimeTemplateFit(float* vals, const TProfile* p_template, const float& ped, const bool& draw = false);

TGraph* PointsAroundTh(float* vals, const unsigned int& sBef, const unsigned int& nSBef, const unsigned int& nSAft);

#endif
