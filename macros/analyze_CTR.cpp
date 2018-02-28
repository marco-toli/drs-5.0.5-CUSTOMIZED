// g++ -Wall -o analyze_CTR.exe analyze_CTR.cpp `root-config --cflags --glibs`


// #include "VectorSmallestInterval.h"
// #include "VectorSmallestInterval.cc"
// #include "fitUtils.h"
// #include "setTDRStyle.h"
// #include "ConfigParser.h"

#include <iostream>
#include <map>

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TProfile2D.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TApplication.h"
#include "TFormula.h"



int main(int argc, char** argv)
{
    
    
    //reading input file
    
//     std::string filename = argv[1];
    std::string filename = "../../drs-data/test_70_56.root";
    TFile * RunFile = new TFile(filename.c_str(),"READ");
    std::cout << "reading file: " <<  filename.c_str() << std::endl;
    TApplication* theApp = new TApplication("App", &argc, argv);
    

    
    
    // defining histos
    int NBINS = 2000;
    TH1F * hCTR           = new TH1F ("hCTR", "hCTR",           NBINS, -10, 10 );
    TH1F * hCTRCalib      = new TH1F ("hCTRCalib", "hCTRCalib", NBINS, -10, 10 );
    TH1F * hCTR_corr      = new TH1F ("hCTR_corr", "hCTR_corr", NBINS, -10, 10 );
    TH1F * hCTRCalib_corr = new TH1F ("hCTRCalib_corr", "hCTRCalib_corr", NBINS, -10, 10 );
    
    TProfile* pAmpWalk      = new TProfile ("pAmpWalk", "pAmpWalk", NBINS, 0, 5);
    TProfile* pAmpWalkCalib = new TProfile ("pAmpWalkCalib", "pAmpWalkCalib", NBINS, 0, 5);
    
    int NCH = 8;
    TH1F * hAmp[NCH];
    for (int iCh = 0; iCh<NCH; iCh++)
    {
        hAmp[iCh] = new TH1F (Form("hAmp_%d", iCh), Form("hAmp_%d", iCh), NBINS, 0, 1);
    }
    
    
    //setting threshold?
    float tl1 = 0.065, th1 = 0.09, tl2 = 0.1, th2 = 0.14;
    
    //reading input tree
    float t_time[8];
    float t_time_c[8];
    float t_ped[8];
    float t_amp[8];
    float t_int[8];
        
    /* 
    float *sample_time  = new float[8192];
    float *sample_value = new float[8192];
    
    if (SAVEWF)
    {
            tree->Branch("sample_time",  sample_time,  "sample_time[8192]/F");
            tree->Branch("sample_value", sample_value, "sample_value[8192]/F");
    }*/
    
    
    TTree* tree = (TTree*) RunFile->Get("ntu");
    
    tree->SetBranchAddress("t_time",   &t_time);
    tree->SetBranchAddress("t_time_c", &t_time_c);//calibrated samples time
    tree->SetBranchAddress("t_ped",    &t_ped);       
    tree->SetBranchAddress("t_amp",    &t_amp);       
    tree->SetBranchAddress("t_int",    &t_int);
    
    
    
    
    int NEVENTS = tree->GetEntries();
    std::cout << "nEvents = " << NEVENTS << std::endl;
    
    
    
    
    
    
    
    std::cout << "(1) looping over events to get amp walk corrections..." << std::endl;
    
    for (Int_t iEvt= 0; iEvt < NEVENTS; iEvt++) 
    {
        tree->GetEntry(iEvt);
        std::cout << "iEvent: " << iEvt << std::endl;
        
        for (int iCh = 0; iCh < NCH; iCh++)
        {
            std::cout << "t_time[" << iCh << "] = " << t_time[iCh] << std::endl;
            hAmp[iCh]->Fill(t_amp[iCh]);
        }
        
        pAmpWalk    ->Fill(t_amp[0]/t_amp[1], t_time[5] - t_time[6] );
        pAmpWalk    ->Fill(t_amp[0]/t_amp[1], t_time_c[5] - t_time_c[6] );
        
        if (t_amp[0]>tl1 && t_amp[0] <th1 &&
            t_amp[1]>tl2 && t_amp[1] <th2
        )
        {            
            hCTR        ->Fill(t_time[5] - t_time[6]);
            hCTRCalib   ->Fill(t_time_c[5] - t_time_c[6]);
            
        }
                
    }//end of loop over events
        
    
    
    TCanvas * cAmps = new TCanvas ("cAmps", "cAmps", 1000,500);
    cAmps->Divide(4,2);
    for (int iCh = 0; iCh<NCH; iCh++)
    {
        cAmps->cd(iCh+1);
        hAmp[iCh]->Draw();
    }
    
    TCanvas * cCTR = new TCanvas ("cCTR", "cCTR", 500, 500);
    hCTR->Draw();
    hCTRCalib->Draw("same");
    hCTRCalib->SetLineColor(kGreen+2);
    
    TCanvas *cAmpWalk = new TCanvas ("cAmpWalk", "cAmpWalk", 500, 500);
    pAmpWalk->Draw();
    pAmpWalk->SetMarkerStyle(20);
    pAmpWalk->SetMarkerColor(kBlack);
    pAmpWalk->SetLineColor(kBlack);
    pAmpWalk->GetYaxis()->SetRangeUser(-1, 1);
    
    float minRatio = 1.0, maxRatio = 2.5;
    TF1* fitCorrRatio = new TF1("fitCorrRatio","[0]*log([1]*x)+[2]",minRatio,maxRatio);
    fitCorrRatio->SetParameters(-0.34,1e-65, -55);
    pAmpWalk->Rebin(8);
    for (int i = 0; i< 3; i++) pAmpWalk->Fit(fitCorrRatio, "QR");
    /*
    TF1* fitCorrCalibRatio = new TF1("fitCorrCalibRatio","[0]*log([1]*x)+[2]",minRatio,maxRatio);
    fitCorrCalibRatio->SetParameters(-0.22,1e-65, -14);
    pAmpWalkCalib->Fit(fitCorrCalibRatio, "QR");
        */
    
    
    
    
    
    
    std::cout << "(2) looping over events to apply amp walk corrections..." << std::endl;    
    for (Int_t iEvt= 0; iEvt < NEVENTS; iEvt++) 
    {
        tree->GetEntry(iEvt);
//         std::cout << "iEvent: " << iEvt << std::endl;
        
        
//         pAmpWalk    ->Fill(t_amp[0]/t_amp[1], t_time[5] - t_time[6] );
        
        if (t_amp[0]>tl1 && t_amp[0] <th1 &&
            t_amp[1]>tl2 && t_amp[1] <th2
        )
        {            
            hCTR_corr        ->Fill(t_time[5] - t_time[6] + fitCorrRatio->Eval(t_amp[0]/t_amp[1]));
//             hCTRCalib_corr   ->Fill(t_time_c[5] - t_time_c[6] + fitCorrCalibRatio->Eval(t_amp[0]/t_amp[1]));
            
        }
                
    }//end of loop over events
    
    
    cCTR->cd();
    hCTR_corr     ->SetLineColor(kBlue+1);
    hCTR_corr     ->Draw("same");
//     hCTRCalib_corr->SetLineColor(kRed+1);
//     hCTRCalib_corr->Draw("same");
    
    
    theApp->Run();
    
    
}