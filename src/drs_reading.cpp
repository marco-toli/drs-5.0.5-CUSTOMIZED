/*******************************************************************************************\
  Name:         drs_reading.cpp
  Created by:   Stefan Ritt
  Modifed by:	Zo� Favier (zoe.favier@etu.unistra.fr)
  For:			Stefan Gundacker (stefan.gundacker@cern.ch)
  Contents:     Application to read out a DRS4 evaluation board and to save data in a text file
  Used for:		SiPM characterization - CERN EP-CMX-DA
  Date:			08/2017
  For Windows use VS C++ 2010 project
  For Unix, compiled with: g++ drs_reading.cpp DRS.cpp averager.cpp mxml.c -o drs_reading `root-config --cflags --libs`
\********************************************************************************************/


// #ifdef _MSC_VER
// #include <windows.h>
// #endif 


#include <math.h>
#include <time.h>
#include <algorithm>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "TFile.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TKey.h"
#include "TString.h"
#include "TTree.h"
#include "TBranch.h"
#include "TH1F.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TLine.h"
#include "TProfile.h"
#include "TSpline.h"
#include "TCanvas.h"
// #include "TDirectory.h." 
#include "TObject.h"
#include "TRandom3.h"
#include "TApplication.h"

#include "strlcpy.h"
#include "DRS.h"
#include "mxml.h"
#include "functions.h"

#include <iostream> // for std::cout
#include <string> // for std::string
#include <fstream>
#include <istream> 
#include <sstream> // for std::istringstream
#include <unistd.h>


using namespace std;

auto const NB_REG = 4; // this value is fixed and shouldn't be changed
                       // it coressponds to the ideal number of points used for the linear regression 

const unsigned int pedSMin = 0;
const unsigned int pedSMax = 40;
const unsigned int nSBef = 2;
const unsigned int nSAft = 2;
unsigned int ampSMin;
const float spikeAmp = 0.003;
const int NCH = 8;
bool isPositive[NCH];
const float SqTH = 0.01;
const unsigned int ampNS = 40;



int main(int argc, char *argv[])
{

	/* reading input parameters */ 
        
        
	if(argc < 2)
        {
		std::cout << "ERROR: too few arguments\nusage: " << argv[0] << " output_file_name (NEVENTS) (DEBUG) (SAVEWF) (sampling speed) (TrigLevel) (TrigDelay) " << std::endl;
		return -1;
	}

	
	TApplication *theApp = new TApplication("App", &argc, argv);
	
	
        
        std::string name_out_file = "../drs-data/";
        
        
        bool SAVEWF   = false;
        bool DEBUG    = false;
        bool MULTI    = true; // if using multiple boards in a daisy chain
        
	int NEVENTS        = 10;
	double GSPS        = 5.12;
        double TrigLevel   = 0.03;
        int TrigDelay      = 0; //number of samples or ns? from right to left
//         double BiasVoltage = -1;
        
        if (argc > 1)  name_out_file += argv[1];     
        if (argc > 2)  NEVENTS = atoi(argv[2]);
//         if (argc > 3)  DEBUG  = (bool) argv[3];
//         if (argc > 4)  SAVEWF = (bool) argv[4];
//         DEBUG = false;
//         SAVEWF = true;
//         MULTI = true;
        
        
        std::string name_out_root = name_out_file + ".root";
        std::string name_out_text = name_out_file + ".txt";
        

        //setting positive of channels for analysis
        isPositive[0] = true; 
        isPositive[1] = true; 
        isPositive[2] = true; 
        isPositive[3] = true;
        
        isPositive[4] = false;
        isPositive[5] = false;
        isPositive[6] = false;
        isPositive[7] = false;
        
//         for (int i= 0; i< NCH; i++) isPositive[i] = false;            

        DRS *drs;
	DRSBoard *b, *mb;
        int nBoards;
                
	
	/* do initial scan - searching for connected boards */
	drs = new DRS();
        if (MULTI) drs->SortBoards();
   
   
	/* show any found board(s) */
        if (!DEBUG)
        {
            for (int i=0 ; i<drs->GetNumberOfBoards() ; i++)
            {
		b = drs->GetBoard(i);
		printf("Found DRS4 evaluation board, serial #%d, firmware revision %d\n", 
		b->GetBoardSerialNumber(), b->GetFirmwareVersion());
            }
        

            /* exit if no board found */
            nBoards = drs->GetNumberOfBoards();
            if (nBoards == 0)
            {
		perror("No DRS4 evaluation board found\n");
		return -1;
            }
	
            std::cout << "initializing board parameters..." << std::endl;
            
            /* use first board with highest serial number as the master board */
            if (MULTI) mb = drs->GetBoard(0);
        
            /* continue working with first board only */
            for (int i=0 ; i<drs->GetNumberOfBoards() ; i++) 
            {
                b = drs->GetBoard(i);

                /* initialize board */
                b->Init();
                
                /* select external reference clock for slave modules */
                /* NOTE: this only works if the clock chain is connected */
                if (i > 0) 
                {
                    if (b->GetFirmwareVersion() >= 21260) // this only works with recent firmware versions
                    { 
                        if (b->GetScaler(5) > 300000)        // check if external clock is connected
                            b->SetRefclk(true);               // switch to external reference clock
                    }
                }

                /* set sampling frequency */
                b->SetFrequency(GSPS, true);

                /* enable transparent mode needed for analog trigger */
//                 b->SetTranspMode(0);

                /* set input range to -0.5V ... +0.5V */
                b->SetInputRange(0);

                /* use following lines to enable hardware trigger*/
                b->EnableTrigger(1, 1);           
//                     b->EnableTrigger(1, 1);    
                
                if (b->GetBoardType() >= 8) // Evaluation Board V4&5
                {

                    if (i == 0) 
                    {
                        /* master board: enable hardware trigger */                            
                        b->SetTranspMode(1);
//                         b->SetTriggerSource(1<<0);        // set CH1 as source                                                                  
                        b->SetTriggerSource(0x3); // 
                                                  //COMBINED TRIGGERS
                                                  // hexadecimal 0x9 means binary 0000 0000 1001 --> bit 8 and 11 ON, means the trigger is on CH1 and CH4  --> coincidence trigger!
                                                  // hexadecimal 0x5 means binary 0000 0000 0101 --> to trigger on ch1 && ch3
                                                  // hexadecimal 0xA means binary 0000 0000 1010 --> to trigger on ch2 && ch4
                                                  // hexadecimal 0x3 means binary 0000 0000 0011 --> to trigger on ch1 && ch2
                                                  // hexadecimal 0x7 means binary 0000 0000 0111 --> to trigger on ch1 && ch2 && ch3
                                                  // hexadecimal 0xC means binary 0000 0000 0011 --> to trigger on ch3 && ch4
                                                  
                        
                                                  //SINGLE CHANNEL TRIGGER
                                                  // hexadecimal 0x1 means binary 0000 0000 0001 -->trigger on single channel 1
                                                  // hexadecimal 0x2 means binary 0000 0000 0010 -->trigger on single channel 2
                                                  // hexadecimal 0x4 means binary 0000 0000 0100 -->trigger on single channel 3
                                                  // hexadecimal 0x8 means binary 0000 0000 1000 -->trigger on single channel 4
                        
                        b->SetTriggerLevel    (TrigLevel);
                        b->SetTriggerDelayNs  (TrigDelay);
                        b->SetTriggerPolarity (true);      // if the argument is false trigger on positive edge
                        
                    } 
                    else 
                    {
                        /* slave boards: enable hardware trigger on Trigger IN */
                        b->SetTriggerSource(1<<4);        // set Trigger IN as source
                        b->SetTriggerPolarity(false);     // isPositive edge
                    } 
      
                    

                }
            }
            std::cout << "board have been initialized" << std::endl;
        }
        
                
                

        //count running time
	float time;
        clock_t t1, t2;
        t1 = clock();
        
        
        //open output file    
        TFile *f_root_out = TFile::Open(name_out_root.c_str(),"RECREATE");
        f_root_out->cd();
	
	
        //defining variables
        float time_array[8][1024];
        float wave_array[8][1024];
        float ave_array [8][1024];
        
	/* set tree parameters */
	std::string tree_name = "ntu";
        
	TTree *tree = new TTree(tree_name.c_str(), tree_name.c_str());
        
        float *sample_time  = new float[8192];
        float *sample_value = new float[8192];
        float *t_time   = new float[8];
        float *t_time_c = new float[8];
        float *t_ped    = new float[8];
        float *t_amp    = new float[8];
        float *t_int    = new float[8];
        
        
        if (SAVEWF)
        {
             tree->Branch("sample_time",  sample_time,  "sample_time[8192]/F");
             tree->Branch("sample_value", sample_value, "sample_value[8192]/F");
        }
        
        tree->Branch("t_time", t_time, "t_time[8]/F");
        tree->Branch("t_time_c", t_time_c, "t_time_c[8]/F");//calibrated samples time
        tree->Branch("t_ped",  t_ped,  "t_ped[8]/F");       
        tree->Branch("t_amp",  t_amp,  "t_amp[8]/F");       
        tree->Branch("t_int",  t_int,  "t_int[8]/F");

        
        //defining histograms and DQM plots
        int NBINS = 4000;
        
        TCanvas * cRawCTR = new TCanvas ("cRawCTR", "cRawCTR", 500, 500);
        TH1F * hRawCTR = new TH1F("hRawCTR", "hRawCTR", 2000, -10, 10);
        TH1F * hRawCTR_calib = new TH1F("hRawCTR_calib", "hRawCTR_calib", 2000, -10, 10);
        hRawCTR_calib->SetLineColor(kRed+1);
        
        TCanvas * cAmplitudes = new TCanvas ("cAmplitudes", "cAmplitudes", 1000,1000);
        cAmplitudes->Divide(2,4);
        
        TCanvas * cPulses = new TCanvas ("cPulses", "cPulses", 1000,1000);
        cPulses->Divide(2,4);
        
        TGraphErrors* gPulse[NCH];        
        TGraphErrors* gPulseAve[NCH];
        
        
        TH1F * hPed[NCH];
        TH1F * hAmp[NCH];
        TH1F * hInt[NCH];
        
        for (int iCh = 0; iCh < NCH; iCh++)
        {
            gPulse[iCh]    = new TGraphErrors ();
            gPulseAve[iCh] = new TGraphErrors ();            
            
            hPed[iCh] = new TH1F (Form("hPed_%d", iCh), Form("hPed_%d", iCh), NBINS, 0, 1);
            hAmp[iCh] = new TH1F (Form("hAmp_%d", iCh), Form("hAmp_%d", iCh), NBINS, 0, 1);
            hInt[iCh] = new TH1F (Form("hInt_%d", iCh), Form("hInt_%d", iCh), NBINS, 0, 1000);
            
        }
        
        //DQM graphical settings
        for (int iCh = 0; iCh < NCH; iCh++)
        {
//             gPulse[iCh]     ->SetMarkerStyle(20);s
//             gPulse[iCh]     ->SetMarkerSize(0.5);
        }
        
        TLine * triggerLine = new TLine(0, TrigLevel, 1024, TrigLevel);
        triggerLine->SetLineColor(kRed+1);
        triggerLine->SetLineStyle(7);
        triggerLine->SetLineWidth(2);
        
        
        
        
        std::cout << "histograms have been initialized" << std::endl;
        
//         ofstream file_writing(name_out_text.c_str(), ios::out | ios::trunc); // open the file in writing, clear the open file if existing
        unsigned int microseconds = 0.02e6;
        float max_amp = -999;
        float min_amp = 999;
        float max_ct  = -999;
        float min_ct  = 999;
        
        
        
        for (int iEv=0 ; iEv<NEVENTS ; iEv++)
	{
            
            if (DEBUG) usleep(microseconds);
            
            //read waveforms from all boards
            if (!DEBUG)     
            {
                /* start boards (activate domino wave), master is last */
                if (MULTI)
                {
                    for (int j=drs->GetNumberOfBoards()-1 ; j>=0 ; j--)
                    {
                        drs->GetBoard(j)->StartDomino();
                    }
                    /* wait for trigger on master board */
                    std::cout << "Waiting for trigger..." << std::endl;
                    fflush(stdout);
                    while (mb->IsBusy());

//                     fprintf(f, "Event #%d =====================================================\n", j);
                    for (int iBoard=0 ; iBoard<drs->GetNumberOfBoards() ; iBoard++) 
                    {
//                         std::cout << "reading board: " << iBoard << std::endl;
                        b = drs->GetBoard(iBoard);
                        if (b->IsBusy()) 
                        {
                            iEv--; /* skip that event, must be some fake trigger */
                            break;
                        }
                        
                        /* read all waveforms from all boards */
                        b->TransferWaves(0, 8);
			
                        // Note: On the evaluation board input #1 is connected to channel 0 and 1 of
                        //the DRS chip, input #2 is connected to channel 2 and 3 and so on. So to
                        //get the input #2 we have to read DRS channel #2, not #1.
              
                        for (int iCh = 0; iCh < 4; iCh++)
                        {
//                             std::cout << "reading channel: " << iCh << " - board: " << iBoard << std::endl;
                
                
                            b->GetTime(0, iCh*2, b->GetTriggerCell(0), time_array[iCh + iBoard*4]);  // read time (X) array of first channel in ns [could use GetTriggerCell(6) too]
                            b->GetWave(0, iCh*2, wave_array[iCh + iBoard*4]); // decode waveform (Y) array of first channel in mV 
                      
                            //saving all waveforms in root tree
                            if (SAVEWF)
                            {
                                for (int iSample = 0; iSample < 1024; iSample++)
                                {
                                    sample_time[(iCh + iBoard*4)*iSample]  = time_array[iCh + iBoard*4][iSample];
                                    sample_value[(iCh + iBoard*4)*iSample] = wave_array[iCh + iBoard*4][iSample];
                                }
                            }
                        }
                    }//end of reading board cycle
                }
                
//                 else
//                 {
//                     
//                                         
//                 /*
//                 }
//                 else
//                 {
// 
// 
// 		// start board (domino wave) 
//                     b = drs->GetBoard(0);
//                     b->StartDomino();
// 		// wait for trigger
//                     std::cout << "Waiting for trigger..." << std::endl;
//                     fflush(stdout);
//                     while (b->IsBusy());
//                     b->TransferWaves(0, 8);
//                 }*/
// 
//                 }
            }
            
            //generate random pulses
            else if (DEBUG) 
            {
                for (int iCh = 0; iCh < NCH; iCh++)
                {
                    
                    int my_pol;
                    float offset;
                    if (isPositive[iCh]) 
                    {
                        my_pol = 1;
                        offset = 0;
                    }
                    else
                    {
                        my_pol = -1;
                        offset = 0;
                    }
                    std::cout << " my_pol[" << iCh << "] = " << my_pol << std::endl;
                
                    //random pulse generator
//                    int t0 = (int) gRandom->Gaus(0, 0.) + 100;
                    int t0 = 100;
                
                    
                    for (int iSample = 0; iSample < 1024; iSample++)
                    {
//                         TRandom *rand = new TRandom3();
// //                         float my_ped = (rand->Rndm()-0.5)/10;
                        float my_ped = gRandom->Uniform(0,10);
//                        std::cout << " my_ped[" << iCh << "][" << iSample << "] = " << my_ped << std::endl;
                                            
                        wave_array[iCh][iSample] = my_ped + offset;
                    
                        if (iCh!= 3 && iCh !=4)
                        {
                            if      (iSample == t0+1)  wave_array[iCh][iSample] += 100*my_pol;
                            else if (iSample == t0+2)  wave_array[iCh][iSample] += 200*my_pol;
                            else if (iSample == t0+3)  wave_array[iCh][iSample] += 300*my_pol;
                            else if (iSample == t0+4)  wave_array[iCh][iSample] += 400*my_pol;
                            else if (iSample  > t0+4 
                                    && iSample  < t0+300) wave_array[iCh][iSample] += 500*my_pol;
                        }
// *
                    }
                }
            }


            
            //search for spikes using ch3 and ch4 as references
            std::vector <int> spikesIt_mb;
            std::vector <int> spikesIt_b2;            
            if(!DEBUG)
            {
                spikesIt_mb = GetSpikes(wave_array[3],spikeAmp);
                spikesIt_b2 = GetSpikes(wave_array[4],spikeAmp);
            }
            
            
            //pulse shape analysis of all channels
            for (int iCh = 0; iCh < NCH; iCh++)
            {
                if (!DEBUG)
                {
                    //clean up spikes
                    //master board
                    if (iCh<4)
                    {
                        for (int iSpike = 0; iSpike < spikesIt_mb.size(); iSpike++) 
                        {                        
                            if (spikesIt_mb.at(iSpike)>0) wave_array[iCh][spikesIt_mb.at(iSpike)] = wave_array[iCh][spikesIt_mb.at(iSpike)-1];
                        }
                    }             
                    //slave board
                    else
                    {
                        for (int iSpike = 0; iSpike < spikesIt_b2.size(); iSpike++) 
                        {                        
                            if (spikesIt_b2.at(iSpike)>0) wave_array[iCh][spikesIt_b2.at(iSpike)] = wave_array[iCh][spikesIt_b2.at(iSpike)-1];
                        }
                    }
                }
                
                //filling graphs    
                for (int iSample = 0; iSample <1024; iSample++)
                {
                    gPulse[iCh]->SetPoint(iSample, iSample, wave_array[iCh][iSample]/1000);
                    ave_array[iCh][iSample] += wave_array[iCh][iSample]/1000;
                }
                
                cPulses->cd(iCh+1);
                if (isPositive[iCh])  // for isPositive
                {
                    std::pair<std::pair<float,std::pair<float,float> >,float >  valsAmp = GetAmplitudePulse(wave_array[iCh],pedSMin,pedSMax,!isPositive[iCh]);
                //      std::pair<std::pair<float,float>,float> valsAmp = GetAmplitudePulse(chVal[iCh],pedSMin,pedSMax,false);
                    t_ped[iCh] = valsAmp.first.first;
                    t_amp[iCh] = valsAmp.first.second.first;
                    t_int[iCh] = valsAmp.first.second.second;
                    
                    std::pair<std::pair<float,float>,std::pair<TF1*,TF1*> > valsTime = GetTimeTh(wave_array[iCh],pedSMax,t_ped[iCh],t_amp[iCh],0.5,nSBef,nSAft,true,false,false);
                    t_time[iCh] = valsTime.first.first;
//                     t_tot[iCh]  = (valsTime.first.second - t_time[iCh]) * t_amp[iCh]; // (time2-time1) * amplitude of square wave = duration*amplitude
                    
                    std::pair<std::pair<float,float>,std::pair<TF1*,TF1*> > valsTimeW = GetTimeThWeighed(wave_array[iCh], time_array[iCh], pedSMax,t_ped[iCh],t_amp[iCh],0.5,nSBef,nSAft,true,false,false);
                    t_time_c[iCh] = valsTimeW.first.first;

                    
                }
                else    //for negative
                {
                    std::pair<float,float> valsAmp = GetAmplitudeSquare(wave_array[iCh],pedSMin,pedSMax,SqTH,ampNS,!isPositive[iCh]);
                    t_ped[iCh] = valsAmp.first;
                    t_amp[iCh] = valsAmp.second; // amplitude as linear fit across the square wave plateau
                    
                    std::pair<std::pair<float,float>,std::pair<TF1*,TF1*> > valsTime = GetTimeTh(wave_array[iCh],pedSMax,t_ped[iCh],t_amp[iCh],0.5,nSBef,nSAft,true,true,false);
                    t_time[iCh] = valsTime.first.first;
                    t_int[iCh]  = (valsTime.first.second - t_time[iCh]) * t_amp[iCh]; // (time2-time1) * amplitude of square wave = duration*amplitude
                    
                    std::pair<std::pair<float,float>,std::pair<TF1*,TF1*> > valsTimeW = GetTimeThWeighed(wave_array[iCh], time_array[iCh], pedSMax,t_ped[iCh],t_amp[iCh],0.5,nSBef,nSAft,true,true,false);
                    t_time_c[iCh] = valsTimeW.first.first;
                    
                    std::cout << "time[ " << iCh << "] = " << t_time[iCh] << " :: time_calib = " << t_time_c[iCh] << std::endl;
                }                                       
                
    
                //normalizing to Volts
                t_ped[iCh]/= 1000;
                t_amp[iCh]/= 1000;
                t_int[iCh]/= 1000;
                
                std::cout << "t_amp[" << iCh << "] = " << t_amp[iCh] << " :: t_int = " << t_int[iCh] << " :: t_ped = " << t_ped[iCh] << " :: time = " << t_time[iCh] << std::endl;
                
                gPulse[iCh]->Draw("ALPE");
                if(iCh<3)   triggerLine->Draw("same");
                
                
                int my_pol;
                if (isPositive[iCh]) my_pol = 1;
                else                 my_pol = -1;
                
                TLine * maxAmpLine = new TLine(0, t_amp[iCh]*my_pol + t_ped[iCh], 1024, t_amp[iCh]*my_pol + t_ped[iCh]);
                maxAmpLine->SetLineColor(kGreen+1);
                maxAmpLine->SetLineStyle(7);
                maxAmpLine->SetLineWidth(2);
                maxAmpLine->Draw("same");
                gPad -> Update();
                    
                    
            }
                    /* write the maximum amplitude of CH1 ; the maximum amplitude of CH4 ; 
                     *			the integrale of CH1 ; the integrale of CH4 ;
                     *			the delay between CH2 and CH3 at the threshold in a text file */
                    
                    
                    
// 		if (TEXTFILE) file_writing << amplitude_maxCH1 << " " << amplitude_maxCH4 << " " << integCH1 << " " << integCH4 << " " 	<< th2 - th1 << endl;

            tree->Fill();
            std::cout << "filling tree..." << std::endl;
        
            //filling rough histos
            float CT = t_time[5]-t_time[6];
            float CT_calib = t_time_c[5]-t_time_c[6];
//                std::cout << "t0 = " <<  t_time[0] <<  " :: t1 = " <<  t_time[1] << " :: t2 = " <<  t_time[2] << " :: t3 = " <<  t_time[3] << " :: CT = " << CT << std::endl;
            hRawCTR -> Fill(CT);
            hRawCTR_calib -> Fill(CT_calib);
            
            if (CT > max_ct) max_ct = CT;
            if (CT < min_ct) min_ct = CT;
            
            for (int iCh = 0; iCh < NCH; iCh++) 
            {
                hPed[iCh] -> Fill (t_ped[iCh]);
                hAmp[iCh] -> Fill (t_amp[iCh]);
                hInt[iCh] -> Fill (t_int[iCh]);
                    
                if (t_amp[iCh] >max_amp) max_amp = t_amp[iCh];
                if (t_amp[iCh] <min_amp) min_amp = t_amp[iCh];
            }
                
            //updating histo drawing
            cRawCTR->cd();
            hRawCTR->Draw();                
            hRawCTR_calib->Draw("same");    
            //hRawCTR->GetXaxis()->SetRangeUser(min_ct*1.5, max_ct*1.5);
            gPad -> Update();
            
            for (int iCh = 0; iCh < NCH; iCh++)
            {
                cAmplitudes->cd(iCh+1);
                hAmp[iCh]->Draw();                
                hAmp[iCh]->GetXaxis()->SetRangeUser(min_amp*0.95, max_amp*1.05);
            }                
            gPad -> Update();
            

                
                
            std::cout << "\rEvent #%d read successfully\n" << iEv << std::endl;
	}
	
	
	for (int iCh = 0; iCh < NCH; iCh++)
        {
            for (int iSample = 0; iSample <1024; iSample++)
            {
                gPulseAve[iCh]->SetPoint(iSample, iSample/GSPS, ave_array[iCh][iSample]/NEVENTS);
            }
            gPulseAve[iCh]->Write();
        }
	
	
        int bytes = f_root_out-> Write();
        std::cout << " nr of Kb written. " << int (bytes/1024) << std::endl;
        f_root_out->Close();


	/* delete DRS object -> close USB connection */
	delete drs;

	t2 = clock();
        time = (float)(t2-t1)/CLOCKS_PER_SEC;
        printf("\ntime = %0.1f seconds \n\n", time);


        theApp->Run();
        return 0;
	
}



	



