//-----------------------------------------------------------------------------------//
//--  MINOS_EnergyS_ntuple.C - By M.A. Acero O., A.A. Aguilar-A. - 2021-11-05  --//
//-----------------------------------------------------------------------------------//
//-- This macro computes the change of the number of events given a change on the  --//
//-- energy resolution uncertainty. The result is stored as a function giving the  --//
//-- corresponding derivative to be used in the minumization process.              --//
//-----------------------------------------------------------------------------------//
#include "EnergyS_constants.h" //CAUTION - To check
//#include "constants.h"
#include <math.h>
#include <iostream>
#include <string>
//---*****************************************************************************---//
//------------------------ CONSTANTS ------------------------------------------------//
//---*****************************************************************************---//
//---*****************************************************---//
//-- Global variables and quantities ------------------------//
//---*****************************************************---//
//---*****************************************************---//
double e;       //Escale Pull term
// histogram binning for the data (PRL 110, 251801 (2013))
//const int NB_numu110    = 23;
//double    lo110         =  0.5;
//double    hi110         = 14.0;
double xbins_numu110[NB_numu110+1];
//Number of energy scale factors to test
const int de = 6;
//---*****************************************************---//
TH1F *nosc_spect_hist;
TH1F *nosc_spect_hist_e1[de];
TH1F *nosc_spect_hist_diff1[de];
TH1F *nosc_spect_hist_deri1[de];
TH1F *nosc_spect_hist_e2[de];
TH1F *nosc_spect_hist_diff2[de];
TH1F *nosc_spect_hist_deri2[de];

TH1F *nosc_spect_hist_2nd_diff[de];
TH1F *nosc_spect_hist_2nd_deri[de];
//---*****************************************************---//

int MINOS_EnergyS_ntuple()
{
    //------------- Style --------------
    //gROOT->SetStyle("Plain");
    gStyle->SetOptStat(0);
    //------------- Style --------------

    cout << "Let's begin..." << endl;

    //-- File with the predicted spectra (for different oscillation paramenters)
    TString filePath = dirName;
    cout << "Reading file - Loop in progress..." << endl;
    TFile *fntuple = new TFile(filePath + "/data/minos_ntuple.root","READ");
    TTree *Tnumu = (TTree*)fntuple->Get("Tnumu");
    TCut cutES;

    //-- Binning definition
    double delta_bins110 = (9.0 - 0.5)/17.0; // 0.5 GeV/bin
    for (int i = 0 ; i < (NB_numu110 - 5) ; i++)
        xbins_numu110[i] = lo110 + i*delta_bins110;
    xbins_numu110[18] = xbins_numu110[17] + 0.75;
    xbins_numu110[19] = xbins_numu110[18] + 0.75;
    xbins_numu110[20] = xbins_numu110[19] + 0.75;
    xbins_numu110[21] = xbins_numu110[20] + 0.75;
    xbins_numu110[22] = xbins_numu110[21] + 1.0;
    xbins_numu110[23] = hi110;

    //-- Declaration of Spectra
    //- Original spectra
    nosc_spect_hist = new TH1F("nosc_spect_hist","",NB_numu110,xbins_numu110);
    nosc_spect_hist->GetXaxis()->SetTitle("Reconstructed Energy (GeV)");
    nosc_spect_hist->GetYaxis()->SetTitle("Events");
    nosc_spect_hist->SetLineColor(kBlack);
    nosc_spect_hist->SetLineWidth(2);
    //- Corrected spectra (energy scale)
    for (int j = 0 ; j < de ; j++){
        nosc_spect_hist_e1[j] = new TH1F(Form("nosc_spect_hist_e1_%d",j),"",NB_numu110,xbins_numu110);
        nosc_spect_hist_e1[j]->SetLineColor(kBlue-j);
        nosc_spect_hist_e2[j] = new TH1F(Form("nosc_spect_hist_e2_%d",j),"",NB_numu110,xbins_numu110);
        nosc_spect_hist_e2[j]->SetLineColor(kBlue-j);
    }
    
    double e;
    double sesc       = 0.06;          //from J. Mitchell PhD. Thesis
    double lowe       = -1*sesc;
    double highe      = +1*sesc;
    double deltae     = (highe - lowe)/(de-1);
    double f_ePos[NB_numu110] = {0.0};
    double f_eNeg[NB_numu110] = {0.0};
    double spcNew[NB_numu110];
    double Eold_i,Enew_i,Enew_j;
    
    //-- Filling the original/unmodified histogram (no oscillated spectrum)
    Tnumu->Draw("Ereco >> nosc_spect_hist");

    double fdel = 0.7;        
    //-- Filling the corrected/modified histograms
    for (int ee = 0 ; ee < de ; ee++) {
        e = lowe + ee*deltae;
        Tnumu->Draw(Form("(1+%f)*Ereco >> nosc_spect_hist_e1_%d",e,ee));
        Tnumu->Draw(Form("(1+%f)*Ereco >> nosc_spect_hist_e2_%d",e+fdel*deltae,ee)); //scaling twice 
    }//End for(ee)

    //-- Computing the difference produced by the energy scale
    //- (Scaled - Original) spectra
    for (int ee = 0 ; ee < de ; ee++) {
        e = lowe + ee*deltae;
        nosc_spect_hist_diff1[ee] = new TH1F(*nosc_spect_hist);
        nosc_spect_hist_diff1[ee] ->SetName(Form("nosc_spec_hist_diff1_%d",ee));
        nosc_spect_hist_deri1[ee] = new TH1F(*nosc_spect_hist);
        nosc_spect_hist_deri1[ee] ->SetName(Form("nosc_spec_hist_deri1_%d",ee));
        nosc_spect_hist_diff1[ee]->GetYaxis()->SetTitle("Events(Scaled - Original) 1 esc");

        nosc_spect_hist_diff2[ee] = new TH1F(*nosc_spect_hist);
        nosc_spect_hist_diff2[ee] ->SetName(Form("nosc_spec_hist_diff2_%d",ee));
        nosc_spect_hist_deri2[ee] = new TH1F(*nosc_spect_hist);
        nosc_spect_hist_deri2[ee] ->SetName(Form("nosc_spec_hist_deri2_%d",ee));
        nosc_spect_hist_diff2[ee]->GetYaxis()->SetTitle("Events(Scaled - Original) 2 esc");

        nosc_spect_hist_2nd_diff[ee] = new TH1F(*nosc_spect_hist);
        nosc_spect_hist_2nd_diff[ee] ->SetName(Form("nosc_spec_hist_2nd_diff_%d",ee));
        nosc_spect_hist_2nd_deri[ee] = new TH1F(*nosc_spect_hist);
        nosc_spect_hist_2nd_deri[ee] ->SetName(Form("nosc_spec_hist_2nd_deri_%d",ee));

        for (int i = 1 ; i <= NB_numu110 ; i++) {
            double binCont1 = nosc_spect_hist_e1[ee]->GetBinContent(i) - nosc_spect_hist->GetBinContent(i);
            nosc_spect_hist_diff1[ee]->SetBinContent(i,binCont1);
            nosc_spect_hist_deri1[ee]->SetBinContent(i,binCont1/(e));
            double binCont2 = nosc_spect_hist_e2[ee]->GetBinContent(i) - nosc_spect_hist->GetBinContent(i);
            nosc_spect_hist_diff2[ee]->SetBinContent(i,binCont2);
            nosc_spect_hist_deri2[ee]->SetBinContent(i,binCont2/(e+fdel*deltae));
	}

        //for (int i = 1 ; i <= NB_numu110 ; i++) {
	//    //double binCont3 = nosc_spect_hist_deri2[ee]->GetBinContent(i) - nosc_spect_hist_deri1[ee]->GetBinContent(i);
	//    double binCont3 = nosc_spect_hist_deri2[ee]->GetBinContent(i) - nosc_spect_hist_deri1[3]->GetBinContent(i);
        //    nosc_spect_hist_2nd_diff[ee]->SetBinContent(i,binCont3);
        //    nosc_spect_hist_2nd_deri[ee]->SetBinContent(i,binCont3/(e));	    
	// }
        nosc_spect_hist_diff1[ee]->SetLineColor(kAzure+ee);
        nosc_spect_hist_deri1[ee]->SetLineColor(kAzure+ee);
        nosc_spect_hist_diff2[ee]->SetLineColor(kGreen+ee);
        nosc_spect_hist_deri2[ee]->SetLineColor(kGreen+ee);

        //nosc_spect_hist_2nd_diff[ee]->SetLineColor(kRed-ee);
        //nosc_spect_hist_2nd_deri[ee]->SetLineColor(kRed-ee);
    }


    for (int ee = 0 ; ee < de ; ee++) {
        e = lowe + ee*deltae;

        for (int i = 1 ; i <= NB_numu110 ; i++) {
	    double binCont3 = nosc_spect_hist_deri2[ee]->GetBinContent(i) - nosc_spect_hist_deri1[ee]->GetBinContent(i);
	    //double binCont3 = nosc_spect_hist_deri2[ee]->GetBinContent(i) - nosc_spect_hist_deri1[3]->GetBinContent(i);
            nosc_spect_hist_2nd_diff[ee]->SetBinContent(i,binCont3);
            nosc_spect_hist_2nd_deri[ee]->SetBinContent(i,binCont3/(fdel*deltae));
        }

        nosc_spect_hist_2nd_diff[ee]->SetLineColor(kRed-ee);
        nosc_spect_hist_2nd_deri[ee]->SetLineColor(kRed-ee);

    }


    TF1 *fFit4_1e;
    fFit4_1e = new TF1("fFit4_1e","[0]+[1]*x+[2]*(x**2)+[3]*(x**3)+[4]*(x**4)",xbins_numu110[0],xbins_numu110[NB_numu110]);
    fFit4_1e->SetParameters(0.0,0.0,0.0,0.0,0.0);
    nosc_spect_hist_deri1[1]->Fit("fFit4_1e","r");
    fFit4_1e->SetLineStyle(2);

    TF1 *fFit4_2e;
    fFit4_2e = new TF1("fFit4_2e","[0]+[1]*x+[2]*(x**2)+[3]*(x**3)+[4]*(x**4)",xbins_numu110[0],xbins_numu110[NB_numu110]);
    fFit4_2e->SetParameters(0.0,0.0,0.0,0.0,0.0);
    nosc_spect_hist_deri2[1]->Fit("fFit4_2e","r");
    fFit4_2e->SetLineStyle(2);

    TF1 *fFit7_1e;
    fFit7_1e = new TF1("fFit7_1e","[0]+[1]*x+[2]*(x**2)+[3]*(x**3)+[4]*(x**4)+[5]*(x**5)+[6]*(x**6)+[7]*(x**7)",xbins_numu110[0],xbins_numu110[NB_numu110]);
    fFit7_1e->SetParameters(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0);
    nosc_spect_hist_deri1[1]->Fit("fFit7_1e","r");

    TF1 *fFit7_2e;
    fFit7_2e = new TF1("fFit7_2e","[0]+[1]*x+[2]*(x**2)+[3]*(x**3)+[4]*(x**4)+[5]*(x**5)+[6]*(x**6)+[7]*(x**7)",xbins_numu110[0],xbins_numu110[NB_numu110]);
    fFit7_2e->SetParameters(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0);
    nosc_spect_hist_deri2[1]->Fit("fFit7_2e","r");

    TF1 *fFitGD_1e;
    fFitGD_1e = new TF1("fFitGD_1e","[0]+[1]*exp(-0.5*pow((x-[2])/[3],2))+[4]*exp(-0.5*pow((x-[5])/[6],2))",xbins_numu110[0],xbins_numu110[NB_numu110]);
    fFitGD_1e->SetParameters(0.0,-2e6,3.0, 2.0,1e6,5.0,2.0);
    fFitGD_1e->SetLineColor(kBlack);
    nosc_spect_hist_deri1[4]->Fit("fFitGD_1e","r");

    TF1 *fFitGD_2e;
    fFitGD_2e = new TF1("fFitGD_2e","[0]+[1]*exp(-0.5*pow((x-[2])/[3],2))+[4]*exp(-0.5*pow((x-[5])/[6],2))",xbins_numu110[0],xbins_numu110[NB_numu110]);
    fFitGD_2e->SetParameters(0.0,-2e6,3.0, 2.0,1e6,5.0,2.0);
    fFitGD_2e->SetLineColor(kBlack);
    nosc_spect_hist_deri2[4]->Fit("fFitGD_2e","r");

    TF1 *fFitGD_2nd;
    fFitGD_2nd = new TF1("fFitGD_2nd","[0]+[1]*exp(-0.5*pow((x-[2])/[3],2))+[4]*exp(-0.5*pow((x-[5])/[6],2))",xbins_numu110[0],xbins_numu110[NB_numu110]);
    fFitGD_2nd->SetParameters(0.0,-2e6,3.0, 2.0,1e6,5.0,2.0);
    fFitGD_2nd->SetLineColor(kBlack);
    nosc_spect_hist_2nd_deri[4]->Fit("fFitGD_2nd","r");
    
    TFile *fEscale = new TFile(filePath + "/data/minos_EScaleDerivative.root","recreate");
    fFit4_1e->Write();
    fFit7_1e->Write();
    fFitGD_1e->Write();
    fFit4_2e->Write();
    fFit7_2e->Write();
    fFitGD_2e->Write();

    fFitGD_2nd->Write();

    nosc_spect_hist->Write();
    for (int i = 0; i < de; i++){
      nosc_spect_hist_e1[i]->Write();
      nosc_spect_hist_diff1[i]->Write();
      nosc_spect_hist_deri1[i]->Write();
      nosc_spect_hist_e2[i]->Write();
      nosc_spect_hist_diff2[i]->Write();
      nosc_spect_hist_deri2[i]->Write();

      nosc_spect_hist_2nd_diff[i]->Write();
      nosc_spect_hist_2nd_deri[i]->Write();
    }

    //-- Plotting section -----------------------//
    //-- Defining histogram frames --//
    double max_numu = 0.0, min_numu = 0.0;
    for (int ee = 0 ; ee < de ; ee++) {
        double max_temp = TMath::Max(nosc_spect_hist->GetMaximum(),nosc_spect_hist_e1[ee]->GetMaximum());
        if (max_temp >= max_numu)
            max_numu = max_temp;
        double min_temp = TMath::Min(nosc_spect_hist->GetMinimum(),nosc_spect_hist_e1[ee]->GetMinimum());
        if (min_temp <= min_numu)
            min_numu = min_temp;
    }
    TH2F *framenumu = new TH2F("framenumu","",NB_numu110,xbins_numu110,10,min_numu,1.01*max_numu);
    
    double max_numudif = 0.0, min_numudif = 0.0;
    for (int ee = 0 ; ee < de ; ee++) {
        double max_temp = TMath::Max(max_numudif,nosc_spect_hist_diff1[ee]->GetMaximum());
        if (max_temp >= max_numudif)
            max_numudif = max_temp;
        double min_temp = TMath::Min(min_numudif,nosc_spect_hist_diff1[ee]->GetMinimum());
        if (min_temp <= min_numudif)
            min_numudif = min_temp;
    }
    TH2F *framenumudif = new TH2F("framenumudif","",NB_numu110,xbins_numu110,10,1.1*min_numudif,1.1*max_numudif);
    
    double max_numuder = 0.0, min_numuder = 0.0;
    for (int ee = 0 ; ee < de ; ee++) {
        double max_temp = TMath::Max(max_numuder,nosc_spect_hist_deri1[ee]->GetMaximum());
        if (max_temp >= max_numuder)
            max_numuder = max_temp;
        double min_temp = TMath::Min(min_numuder,nosc_spect_hist_deri1[ee]->GetMinimum());
        if (min_temp <= min_numuder)
            min_numuder = min_temp;
    }
    TH2F *framenumuder = new TH2F("framenumuder","",NB_numu110,xbins_numu110,10,1.1*min_numuder,1.1*max_numuder);
    //-- END - Defining histogram frames

    TCanvas *canv1 = new TCanvas("canv1", "MINOS", 1*600,2*400);
    canv1->Divide(1,2);
    //Muon neutrinos
    canv1->cd(1);
    framenumu->Draw();
    nosc_spect_hist->Draw("same hist");
    for (int ee = 0 ; ee < de ; ee++){
        nosc_spect_hist_e1[ee]->Draw("same hist");
        nosc_spect_hist_e2[ee]->Draw("same hist");}
 
   canv1->cd(2);
    framenumudif->Draw();
    for (int ee = 0 ; ee < de ; ee++){
      nosc_spect_hist_diff1[ee]->Draw("same hist");
      nosc_spect_hist_diff1[ee]->Draw("same hist");}
    
    //-- Derivatives --//
    TCanvas *canv2 = new TCanvas("canv2", "MINOS Derivatives", 1*600,2*400);
    canv2->Divide(1,2);
    //Muon Neutrinos
    canv2->cd(1);
    framenumu->Draw();
    nosc_spect_hist->Draw("same hist");
    for (int ee = 0 ; ee < de ; ee++){
        nosc_spect_hist_e1[ee]->Draw("same hist");
        nosc_spect_hist_e2[ee]->Draw("same hist");}
    canv2->cd(2);
    framenumuder->Draw();
    for (int ee = 0 ; ee < de ; ee++){
        nosc_spect_hist_deri1[ee]->Draw("same hist");
        nosc_spect_hist_deri2[ee]->Draw("same hist");}
    fFit4_1e->Draw("same");
    fFit7_1e->Draw("same");
    fFitGD_1e->Draw("same");

    fFit4_2e->Draw("same");
    fFit7_2e->Draw("same");
    fFitGD_2e->Draw("same");


    std::cout << "Succesful run!!" << endl;

    return 0;

}//-END
