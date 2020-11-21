//-----------------------------------------------------------------------------------//
//--  RENO_test_EnergyS_ntuple.C - By M.A. Acero O., A.A. Aguilar-A. - 2020-06-19  --//
//-----------------------------------------------------------------------------------//
//-- This macro computes the change of the number of events given a change on the  --//
//-- energy resolution uncertainty. The result is stored as a function giving the  --//
//-- corresponding derivative to be used in the minumization process.              --//
//-----------------------------------------------------------------------------------//
#include "constants.h"
#include <math.h>
#include <iostream>
#include <string>
//---*****************************************************************************---//
//------------------------ CONSTANTS ------------------------------------------------//
//---*****************************************************************************---//
//---*****************************************************---//
//-- Global variables and quantities ------------------------//
//---*****************************************************---//
//Efficiencies and DAQ time (PRL128 (2018))
double emuem[nDet] ={0.7647,0.7647};
double daqTime[nDet] = {1807.88,2193.04};
//---*****************************************************---//
//IBD rate per day w/o oscillations
double noOsc_IBDrate_perday[nDet];
//---*****************************************************---//
double s2th_13; //oscillation parameter to be fitted
double dm2_ee;   //oscillation parameter to be fitted
double e; //Escale Pull term
double spc[nDet][NB];
double NoscTot[nDet];
int iAD;
int iNR;
double xbins[NB+1];
//Number of energy scale factors to test
const int de = 6;
//---*****************************************************---//
TH1F *nosc_spect_hist[nDet];
TH1F *nosc_spect_hist_e[nDet][de];
TH1F *nosc_spect_hist_diff[nDet][de];
TH1F *nosc_spect_hist_deri[nDet][de];
//---*****************************************************---//

int RENO_EScaleDeriv_ntuple()
{
    //------------- Style --------------
    //gROOT->SetStyle("Plain");
    gStyle->SetOptStat(0);
    //------------- Style --------------

    cout << "Let's begin..." << endl;
    //-- File to get noOsc normalizations
    TString filePath = dirName;
    ifstream IBDrates_file(filePath + "/files/RENO_noOsc_IBDrates_perday.txt");
    cout << "Reading noOsc normalizations file ..." << endl;
    for (int i=0; i< nDet; i ++){
        IBDrates_file >> noOsc_IBDrate_perday[i];
        cout << "noOscIBD_rates_perday " << i << ": " << noOsc_IBDrate_perday[i] << endl;
    }//for

    //-- File with the predicted spectra (for different oscillation paramenters)
    cout << "Reading file - Loop in progress..." << endl;
    TFile *fntuple = new TFile(filePath + "/files_root/RENO-ntuple_BFosc.root","READ");
    TTree *T = (TTree*)fntuple->Get("T");
    TCut cutES;

    //-- Binning definition
    double delta_bins2 = (5.6 - 1.2)/22; // 0.2 MeV/bin
    
    for (int i = 0 ; i < (NB-3) ; i++){
        xbins[i] = 1.2 + delta_bins2*i;
    }
    xbins[23] = xbins[22] + 0.4;
    xbins[24] = xbins[23] + 0.4;
    xbins[25] = 8.4 - 1.4;
    xbins[26] = 8.4;

    //-- Declaration of Spectra
    for(int iAD = 0 ; iAD < nDet ; iAD++) {
        //- Original spectra
        nosc_spect_hist[iAD] = new TH1F(Form("nosc_spect_hist_%d",iAD),"",NB,xbins);
        nosc_spect_hist[iAD]->GetXaxis()->SetTitle("Prompt Reconstructed Energy (MeV)");
        nosc_spect_hist[iAD]->GetYaxis()->SetTitle("Events");
        nosc_spect_hist[iAD]->SetLineColor(kBlack);
        nosc_spect_hist[iAD]->SetLineWidth(2);
        for (int j = 0 ; j < de ; j++){
            //- Corrected spectra (energy scale)
            nosc_spect_hist_e[iAD][j] = new TH1F(Form("nosc_spect_hist_e_%d_%d",iAD,j),"",NB,xbins);
            nosc_spect_hist_e[iAD][j]->SetLineColor(kBlue+j);
        }
    }
    
    double e;
    double sesc       = 0.0015;          //from PRL128 (2018)
    double lowe       = -5*sesc;
    double highe      = +5*sesc;
    double deltae     = (highe - lowe)/(de-1);
    double f_ePos[NB] = {0.0};
    double f_eNeg[NB] = {0.0};
    double spcNew[nDet][NB];
    double Eold_i,Enew_i,Enew_j;
    
    for (int iAD = 0 ; iAD < nDet ; iAD++)
    {
        //-- Filling the original/unmodified histogram
        //cutES = Form("(id==%d)",iAD);
        cutES = Form("(1.0/(1.0 - %e*((sin( 1.267 * %e * Ln/En ))**2) - ((cos(0.5 * asin(sqrt(%e))))**4) * %e * (sin( 1.267 * %e * Ln/En ))**2))*(id==%d)",ssq2th13RENO,dmsqeeRENO,ssq2th13RENO,ssq2th12RENO,dmsq21RENO,iAD);
        T->Draw(Form("Ep >> nosc_spect_hist_%d",iAD),cutES,"");
        
        //-- Filling the corrected/modified histograms
        for (int ee = 0 ; ee < de ; ee++) {
            e = lowe + ee*deltae;
            //cutES = Form("(1+%f)*(id==%d)",e,iAD);
            //T->Draw(Form("Ep >> nosc_spect_hist_e_%d_%d",iAD,ee),cutES,"");
            T->Draw(Form("(1+%f)*Ep >> nosc_spect_hist_e_%d_%d",e,iAD,ee),cutES);
        }//End for(ee)
    }// End for(iAD)

    //-- Computing the difference produced by the energy scale
    //- (Original - Corrected) spectra
    for (int ee = 0 ; ee < de ; ee++) {
        e = lowe + ee*deltae;
        //- Near Detector
        nosc_spect_hist_diff[0][ee] = new TH1F(*nosc_spect_hist[0]);
        nosc_spect_hist_deri[0][ee] = new TH1F(*nosc_spect_hist[0]);
        nosc_spect_hist_diff[0][ee]->GetYaxis()->SetTitle("Events(Original - Corrected)");
        //- Far Detector
        nosc_spect_hist_diff[1][ee] = new TH1F(*nosc_spect_hist[1]);
        nosc_spect_hist_deri[1][ee] = new TH1F(*nosc_spect_hist[1]);
        nosc_spect_hist_diff[1][ee]->GetYaxis()->SetTitle("Events(Original - Corrected)");
        for (int i = 1 ; i <= NB ; i++) {
            double binCont = nosc_spect_hist[0]->GetBinContent(i) - nosc_spect_hist_e[0][ee]->GetBinContent(i);
            nosc_spect_hist_diff[0][ee]->SetBinContent(i,binCont);
            nosc_spect_hist_deri[0][ee]->SetBinContent(i,binCont/(e));
            
            binCont = nosc_spect_hist[1]->GetBinContent(i) - nosc_spect_hist_e[1][ee]->GetBinContent(i);
            nosc_spect_hist_diff[1][ee]->SetBinContent(i,binCont);
            nosc_spect_hist_deri[1][ee]->SetBinContent(i,binCont/(e));
        }
        nosc_spect_hist_diff[0][ee]->SetLineColor(kAzure+ee);
        nosc_spect_hist_diff[1][ee]->SetLineColor(kAzure+ee);
        nosc_spect_hist_deri[0][ee]->SetLineColor(kAzure+ee);
        nosc_spect_hist_deri[1][ee]->SetLineColor(kAzure+ee);
    }

    TF1 *fFit4[nDet];
    TF1 *fFit7[nDet];
    for (int i = 0 ; i < nDet ; i++) {
        fFit4[i] = new TF1(Form("fFit4_%d",i),"[0]+[1]*x+[2]*(x**2)+[3]*(x**3)+[4]*(x**4)",xbins[0],xbins[NB]);
        fFit4[i]->SetParameters(0.0,0.0,0.0,0.0,0.0);
        nosc_spect_hist_deri[i][1]->Fit(Form("fFit4_%d",i),"r");
        fFit4[i]->SetLineStyle(2);

        fFit7[i] = new TF1(Form("fFit7_%d",i),"[0]+[1]*x+[2]*(x**2)+[3]*(x**3)+[4]*(x**4)+[5]*(x**5)+[6]*(x**6)+[7]*(x**7)",xbins[0],xbins[NB]);
        fFit7[i]->SetParameters(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0);
        nosc_spect_hist_deri[i][1]->Fit(Form("fFit7_%d",i),"r");
    }
    
    TFile *fEscale = new TFile(filePath + "/files_root/RENO_EScaleDerivative.root","recreate");
    for (int i = 0 ; i < nDet ; i++) {
        fFit4[i]->Write();
        fFit7[i]->Write();
    }
    
    //-- Plotting section -----------------------//
    //-- Defining histogram frames --//
    //-- ND Histograms frames
    double max_ND = 0.0, min_ND = 0.0;
    for (int ee = 0 ; ee < de ; ee++) {
        double max_temp = TMath::Max(nosc_spect_hist[0]->GetMaximum(),nosc_spect_hist_e[0][ee]->GetMaximum());
        if (max_temp >= max_ND)
            max_ND = max_temp;
        double min_temp = TMath::Min(nosc_spect_hist[0]->GetMinimum(),nosc_spect_hist_e[0][ee]->GetMinimum());
        if (min_temp <= min_ND)
            min_ND = min_temp;
    }
    TH2F *frameND = new TH2F("frameND","",NB,xbins,10,min_ND,1.01*max_ND);
    
    double max_NDdif = 0.0, min_NDdif = 0.0;
    for (int ee = 0 ; ee < de ; ee++) {
        double max_temp = TMath::Max(max_NDdif,nosc_spect_hist_diff[0][ee]->GetMaximum());
        if (max_temp >= max_NDdif)
            max_NDdif = max_temp;
        double min_temp = TMath::Min(min_NDdif,nosc_spect_hist_diff[0][ee]->GetMinimum());
        if (min_temp <= min_NDdif)
            min_NDdif = min_temp;
    }
    TH2F *frameNDdif = new TH2F("frameNDdif","",NB,xbins,10,1.1*min_NDdif,1.1*max_NDdif);
    
    double max_NDder = 0.0, min_NDder = 0.0;
    for (int ee = 0 ; ee < de ; ee++) {
        double max_temp = TMath::Max(max_NDder,nosc_spect_hist_deri[0][ee]->GetMaximum());
        if (max_temp >= max_NDder)
            max_NDder = max_temp;
        double min_temp = TMath::Min(min_NDder,nosc_spect_hist_deri[0][ee]->GetMinimum());
        if (min_temp <= min_NDder)
            min_NDder = min_temp;
    }
    TH2F *frameNDder = new TH2F("frameNDder","",NB,xbins,10,1.1*min_NDder,1.1*max_NDder);

    //-- FD Histograms frames
    double max_FD = 0.0, min_FD = 0.0;
    for (int ee = 0 ; ee < de ; ee++) {
        double max_temp = TMath::Max(nosc_spect_hist[1]->GetMaximum(),nosc_spect_hist_e[1][ee]->GetMaximum());
        if (max_temp >= max_FD)
            max_FD = max_temp;
        double min_temp = TMath::Min(nosc_spect_hist[1]->GetMinimum(),nosc_spect_hist_e[1][ee]->GetMinimum());
            if (min_temp <= min_FD)
                min_FD = min_temp;
    }
    TH2F *frameFD = new TH2F("frameFD","",NB,xbins,10,min_FD,1.01*max_FD);
    
    double max_FDdif = 0.0, min_FDdif = 0.0;
    for (int ee = 0 ; ee < de ; ee++) {
        double max_temp = TMath::Max(max_FDdif,nosc_spect_hist_diff[1][ee]->GetMaximum());
        if (max_temp >= max_FDdif)
            max_FDdif = max_temp;
        double min_temp = TMath::Min(min_FDdif,nosc_spect_hist_diff[1][ee]->GetMinimum());
        if (min_temp <= min_FDdif)
            min_FDdif = min_temp;
    }
    TH2F *frameFDdif = new TH2F("frameFDdif","",NB,xbins,10,1.1*min_FDdif,1.1*max_FDdif);
    
    double max_FDder = 0.0, min_FDder = 0.0;
    for (int ee = 0 ; ee < de ; ee++) {
        double max_temp = TMath::Max(max_FDder,nosc_spect_hist_deri[1][ee]->GetMaximum());
        if (max_temp >= max_FDder)
            max_FDder = max_temp;
        double min_temp = TMath::Min(min_FDder,nosc_spect_hist_deri[1][ee]->GetMinimum());
        if (min_temp <= min_FDder)
            min_FDder = min_temp;
    }
    TH2F *frameFDder = new TH2F("frameFDder","",NB,xbins,10,1.1*min_FDder,1.1*max_FDder);
    //-- END - Defining histogram frames

    TCanvas *canv1 = new TCanvas("canv1", "Near and Far Detectors", 2*600,2*400);
    canv1->Divide(2,2);
    //Near Detector
    canv1->cd(1);
    frameND->Draw();
    nosc_spect_hist[0]->Draw("same hist");
    for (int ee = 0 ; ee < de ; ee++)
        nosc_spect_hist_e[0][ee]->Draw("same hist");
    canv1->cd(2);
    frameNDdif->Draw();
    for (int ee = 0 ; ee < de ; ee++)
        nosc_spect_hist_diff[0][ee]->Draw("same hist");

    //Far Detector
    canv1->cd(3);
    frameFD->Draw();
    nosc_spect_hist[1]->Draw("same hist");
    for (int ee = 0 ; ee < de ; ee++)
        nosc_spect_hist_e[1][ee]->Draw("same hist");
    canv1->cd(4);
    frameFDdif->Draw();
    for (int ee = 0 ;ee < de ; ee++)
        nosc_spect_hist_diff[1][ee]->Draw("same hist");
    
    //-- Derivatives --//
    TCanvas *canv2 = new TCanvas("canv2", "Near and Far Detectors Derivatives", 2*600,2*400);
    canv2->Divide(2,2);
    //Near Detector
    canv2->cd(1);
    frameND->Draw();
    nosc_spect_hist[0]->Draw("same hist");
    for (int ee = 0 ; ee < de ; ee++)
        nosc_spect_hist_e[0][ee]->Draw("same hist");
    canv2->cd(2);
    frameNDder->Draw();
    for (int ee = 0 ; ee < de ; ee++)
        nosc_spect_hist_deri[0][ee]->Draw("same hist");
    fFit4[0]->Draw("same");
    fFit7[0]->Draw("same");

    //Far Detector
    canv2->cd(3);
    frameFD->Draw();
    nosc_spect_hist[1]->Draw("same hist");
    for (int ee = 0 ; ee < de ; ee++)
        nosc_spect_hist_e[1][ee]->Draw("same hist");
    canv2->cd(4);
    frameFDder->Draw();
    for (int ee = 0 ;ee < de ; ee++)
        nosc_spect_hist_deri[1][ee]->Draw("same hist");
    fFit4[1]->Draw("same");
    fFit7[1]->Draw("same");

    std::cout << "Succesful run!!" << endl;

    //canv1->Print(Form("e_%.4f.pdf",e));
    
//    TCanvas *canv2 = new TCanvas("canv2", "", 2*600,2*400);
//    //Near Detector
//    canv2->cd();
//    frameND->Draw();
//    nosc_spect_hist[0]->Draw("same hist");
//    for (int ee = 0 ; ee < de ; ee++)
//        //nosc_spect_hist_e[0][ee]->Draw("same hist");
//        //nosc_spect_hist_diff[0][ee]->Draw("same hist");
//        nosc_spect_hist_deri[0][ee]->Draw("same hist");


    return 0;

}//-END
