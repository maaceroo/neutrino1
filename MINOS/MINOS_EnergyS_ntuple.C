//-----------------------------------------------------------------------------------//
//--  MINOS_EnergyS_ntuple.C - By M.A. Acero O., A.A. Aguilar-A. - 2021-11-05  --//
//-----------------------------------------------------------------------------------//
//-- This macro computes the change of the number of events given a change on the  --//
//-- energy resolution uncertainty. The result is stored as a function giving the  --//
//-- corresponding derivative to be used in the minumization process.              --//
//-----------------------------------------------------------------------------------//
#include "EnergyS_constants.h" 
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
double xbins_numu110   [NB_numu110     +1];
double xbins_numub110  [NB_numubar110  +1];
double xbins_numubWS110[NB_numubarWS110+1];
//Number of energy scale factors to test
const int de = 6;
//---*****************************************************---//
//neutrinos
TH1F *numu_nosc_spect_hist;
TH1F *numu_nosc_spect_hist_e1[de];
TH1F *numu_nosc_spect_hist_diff1[de];
TH1F *numu_nosc_spect_hist_deri1[de];
TH1F *numu_nosc_spect_hist_e2[de];
TH1F *numu_nosc_spect_hist_diff2[de];
TH1F *numu_nosc_spect_hist_deri2[de];

TH1F *numu_nosc_spect_hist_2nd_diff[de];
TH1F *numu_nosc_spect_hist_2nd_deri[de];
//anti neutrinos
TH1F *numub_nosc_spect_hist;
TH1F *numub_nosc_spect_hist_e1[de];
TH1F *numub_nosc_spect_hist_diff1[de];
TH1F *numub_nosc_spect_hist_deri1[de];
TH1F *numub_nosc_spect_hist_e2[de];
TH1F *numub_nosc_spect_hist_diff2[de];
TH1F *numub_nosc_spect_hist_deri2[de];

TH1F *numub_nosc_spect_hist_2nd_diff[de];
TH1F *numub_nosc_spect_hist_2nd_deri[de];
//anti neutrinos WS
TH1F *numubWS_nosc_spect_hist;
TH1F *numubWS_nosc_spect_hist_e1[de];
TH1F *numubWS_nosc_spect_hist_diff1[de];
TH1F *numubWS_nosc_spect_hist_deri1[de];
TH1F *numubWS_nosc_spect_hist_e2[de];
TH1F *numubWS_nosc_spect_hist_diff2[de];
TH1F *numubWS_nosc_spect_hist_deri2[de];

TH1F *numubWS_nosc_spect_hist_2nd_diff[de];
TH1F *numubWS_nosc_spect_hist_2nd_deri[de];
//---*****************************************************---//

int MINOS_EnergyS_ntuple()
{
  //------------- Style --------------
  //gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  //------------- Style --------------
  
  cout << "Let's begin..." << endl;
  
  //-- File with the predicted spectra (for different oscillation parameters)
  TString filePath = dirName;
  cout << "Reading file - Loop in progress..." << endl;
  TFile *fntuple  = new TFile(filePath + "/data/minos_ntuple.root","READ");
  TTree *Tnumu    = (TTree*)fntuple->Get("Tnumu");    //-- neutrinos
  TTree *Tnumub   = (TTree*)fntuple->Get("Tnumub");   //-- anti neutrinos
  TTree *TnumubWS = (TTree*)fntuple->Get("TnumubWS"); //-- anti neutrinos WS
  TCut cutES;
  
  //-- Binning definition
  //numu
  double delta_bins110 = (9.0 - 0.5)/17.0; // 0.5 GeV/bin
  for (int i = 0 ; i < (NB_numu110 - 5) ; i++)
    xbins_numu110[i] = lo110 + i*delta_bins110;
  xbins_numu110[18] = xbins_numu110[17] + 0.75;
  xbins_numu110[19] = xbins_numu110[18] + 0.75;
  xbins_numu110[20] = xbins_numu110[19] + 0.75;
  xbins_numu110[21] = xbins_numu110[20] + 0.75;
  xbins_numu110[22] = xbins_numu110[21] + 1.0;
  xbins_numu110[23] = hi110;
  //numub
  double delta_binsb110 = (12.0 - lob110)/(NB_numubar110-1); // 1.0 GeV/bin
  for (int i = 0 ; i < (NB_numubar110) ; i++)
    xbins_numub110[i] = lob110 + delta_binsb110*i;
  xbins_numub110[12] = xbins_numub110[11] + 2.0;
  //numub WS
  double delta_binsbWS110 = (20.0 - lobWS110)/(NB_numubarWS110-1); // 2.0 GeV/bin
  for (int i = 0 ; i < (NB_numubarWS110) ; i++)
    xbins_numubWS110[i] = lobWS110 + delta_binsbWS110*i;
  xbins_numubWS110[11] = xbins_numubWS110[10] + 5.0;
  
  //-- Declaration of Spectra
  //- Original spectra
  //numu
  numu_nosc_spect_hist = new TH1F("numu_nosc_spect_hist","",NB_numu110,xbins_numu110);
  numu_nosc_spect_hist->GetXaxis()->SetTitle("Reconstructed Energy (GeV)");
  numu_nosc_spect_hist->GetYaxis()->SetTitle("Events");
  numu_nosc_spect_hist->SetLineColor(kBlack);
  numu_nosc_spect_hist->SetLineWidth(2);
  //numub
  numub_nosc_spect_hist = new TH1F("numub_nosc_spect_hist","",NB_numubar110,xbins_numub110);
  numub_nosc_spect_hist->GetXaxis()->SetTitle("Reconstructed Energy (GeV)");
  numub_nosc_spect_hist->GetYaxis()->SetTitle("Events");
  numub_nosc_spect_hist->SetLineColor(kBlack);
  numub_nosc_spect_hist->SetLineWidth(2);
  //numub WS
  numubWS_nosc_spect_hist = new TH1F("numubWS_nosc_spect_hist","",NB_numubarWS110,xbins_numubWS110);
  numubWS_nosc_spect_hist->GetXaxis()->SetTitle("Reconstructed Energy (GeV)");
  numubWS_nosc_spect_hist->GetYaxis()->SetTitle("Events");
  numubWS_nosc_spect_hist->SetLineColor(kBlack);
  numubWS_nosc_spect_hist->SetLineWidth(2);
  //- Corrected spectra (energy scale)
  for (int j = 0 ; j < de ; j++){
    //numu
    numu_nosc_spect_hist_e1[j] = new TH1F(Form("numu_nosc_spect_hist_e1_%d",j),"",NB_numu110,xbins_numu110);
    numu_nosc_spect_hist_e1[j]->SetLineColor(kBlue-j);
    numu_nosc_spect_hist_e2[j] = new TH1F(Form("numu_nosc_spect_hist_e2_%d",j),"",NB_numu110,xbins_numu110);
    numu_nosc_spect_hist_e2[j]->SetLineColor(kBlue-j);
    //numub
    numub_nosc_spect_hist_e1[j] = new TH1F(Form("numub_nosc_spect_hist_e1_%d",j),"",NB_numubar110,xbins_numub110);
    numub_nosc_spect_hist_e1[j]->SetLineColor(kBlue-j);
    numub_nosc_spect_hist_e2[j] = new TH1F(Form("numub_nosc_spect_hist_e2_%d",j),"",NB_numubar110,xbins_numub110);
    numub_nosc_spect_hist_e2[j]->SetLineColor(kBlue-j);
    //numub WS
    numubWS_nosc_spect_hist_e1[j] = new TH1F(Form("numubWS_nosc_spect_hist_e1_%d",j),"",NB_numubarWS110,xbins_numubWS110);
    numubWS_nosc_spect_hist_e1[j]->SetLineColor(kBlue-j);
    numubWS_nosc_spect_hist_e2[j] = new TH1F(Form("numubWS_nosc_spect_hist_e2_%d",j),"",NB_numubarWS110,xbins_numubWS110);
    numubWS_nosc_spect_hist_e2[j]->SetLineColor(kBlue-j);
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
  //numu
  Tnumu   ->Draw("Ereco    >> numu_nosc_spect_hist");
  //numub
  Tnumub  ->Draw("Erecob   >> numub_nosc_spect_hist");
  //numub WS
  TnumubWS->Draw("ErecobWS >> numubWS_nosc_spect_hist");
  
  double fdel = 0.7;
  //-- Filling the corrected/modified histograms
  for (int ee = 0 ; ee < de ; ee++) {
    e = lowe + ee*deltae;
    //numu
    Tnumu->Draw(Form("(1+%f)*Ereco >> numu_nosc_spect_hist_e1_%d",e,ee));
    Tnumu->Draw(Form("(1+%f)*Ereco >> numu_nosc_spect_hist_e2_%d",e+fdel*deltae,ee)); //scaling twice 
    //numub
    Tnumub->Draw(Form("(1+%f)*Erecob >> numub_nosc_spect_hist_e1_%d",e,ee));
    Tnumub->Draw(Form("(1+%f)*Erecob >> numub_nosc_spect_hist_e2_%d",e+fdel*deltae,ee)); //scaling twice 
    //numub WS
    TnumubWS->Draw(Form("(1+%f)*ErecobWS >> numubWS_nosc_spect_hist_e1_%d",e,ee));
    TnumubWS->Draw(Form("(1+%f)*ErecobWS >> numubWS_nosc_spect_hist_e2_%d",e+fdel*deltae,ee)); //scaling twice 
  }//End for(ee)
  
  //-- Computing the difference produced by the energy scale
  //- (Scaled - Original) spectra
  for (int ee = 0 ; ee < de ; ee++) {
    e = lowe + ee*deltae;
    //numu ----------------------------------------------------------
    numu_nosc_spect_hist_diff1[ee] = new TH1F(*numu_nosc_spect_hist);
    numu_nosc_spect_hist_diff1[ee] ->SetName(Form("numu_nosc_spec_hist_diff1_%d",ee));
    numu_nosc_spect_hist_deri1[ee] = new TH1F(*numu_nosc_spect_hist);
    numu_nosc_spect_hist_deri1[ee] ->SetName(Form("numu_nosc_spec_hist_deri1_%d",ee));
    numu_nosc_spect_hist_diff1[ee]->GetYaxis()->SetTitle("Events(Scaled - Original) 1 esc");
    
    numu_nosc_spect_hist_diff2[ee] = new TH1F(*numu_nosc_spect_hist);
    numu_nosc_spect_hist_diff2[ee] ->SetName(Form("numu_nosc_spec_hist_diff2_%d",ee));
    numu_nosc_spect_hist_deri2[ee] = new TH1F(*numu_nosc_spect_hist);
    numu_nosc_spect_hist_deri2[ee] ->SetName(Form("numu_nosc_spec_hist_deri2_%d",ee));
    numu_nosc_spect_hist_diff2[ee]->GetYaxis()->SetTitle("Events(Scaled - Original) 2 esc");
    
    numu_nosc_spect_hist_2nd_diff[ee] = new TH1F(*numu_nosc_spect_hist);
    numu_nosc_spect_hist_2nd_diff[ee] ->SetName(Form("numu_nosc_spec_hist_2nd_diff_%d",ee));
    numu_nosc_spect_hist_2nd_deri[ee] = new TH1F(*numu_nosc_spect_hist);
    numu_nosc_spect_hist_2nd_deri[ee] ->SetName(Form("numu_nosc_spec_hist_2nd_deri_%d",ee));
    
    //numub ----------------------------------------------------------
    numub_nosc_spect_hist_diff1[ee] = new TH1F(*numub_nosc_spect_hist);
    numub_nosc_spect_hist_diff1[ee] ->SetName(Form("numub_nosc_spec_hist_diff1_%d",ee));
    numub_nosc_spect_hist_deri1[ee] = new TH1F(*numub_nosc_spect_hist);
    numub_nosc_spect_hist_deri1[ee] ->SetName(Form("numub_nosc_spec_hist_deri1_%d",ee));
    numub_nosc_spect_hist_diff1[ee]->GetYaxis()->SetTitle("Events(Scaled - Original) 1 esc");
    
    numub_nosc_spect_hist_diff2[ee] = new TH1F(*numub_nosc_spect_hist);
    numub_nosc_spect_hist_diff2[ee] ->SetName(Form("numub_nosc_spec_hist_diff2_%d",ee));
    numub_nosc_spect_hist_deri2[ee] = new TH1F(*numub_nosc_spect_hist);
    numub_nosc_spect_hist_deri2[ee] ->SetName(Form("numub_nosc_spec_hist_deri2_%d",ee));
    numub_nosc_spect_hist_diff2[ee]->GetYaxis()->SetTitle("Events(Scaled - Original) 2 esc");
    
    numub_nosc_spect_hist_2nd_diff[ee] = new TH1F(*numub_nosc_spect_hist);
    numub_nosc_spect_hist_2nd_diff[ee] ->SetName(Form("numub_nosc_spec_hist_2nd_diff_%d",ee));
    numub_nosc_spect_hist_2nd_deri[ee] = new TH1F(*numub_nosc_spect_hist);
    numub_nosc_spect_hist_2nd_deri[ee] ->SetName(Form("numub_nosc_spec_hist_2nd_deri_%d",ee));
    
    //numub WS ----------------------------------------------------------
    numubWS_nosc_spect_hist_diff1[ee] = new TH1F(*numubWS_nosc_spect_hist);
    numubWS_nosc_spect_hist_diff1[ee] ->SetName(Form("numubWS_nosc_spec_hist_diff1_%d",ee));
    numubWS_nosc_spect_hist_deri1[ee] = new TH1F(*numubWS_nosc_spect_hist);
    numubWS_nosc_spect_hist_deri1[ee] ->SetName(Form("numubWS_nosc_spec_hist_deri1_%d",ee));
    numubWS_nosc_spect_hist_diff1[ee]->GetYaxis()->SetTitle("Events(Scaled - Original) 1 esc");
    
    numubWS_nosc_spect_hist_diff2[ee] = new TH1F(*numubWS_nosc_spect_hist);
    numubWS_nosc_spect_hist_diff2[ee] ->SetName(Form("numubWS_nosc_spec_hist_diff2_%d",ee));
    numubWS_nosc_spect_hist_deri2[ee] = new TH1F(*numubWS_nosc_spect_hist);
    numubWS_nosc_spect_hist_deri2[ee] ->SetName(Form("numubWS_nosc_spec_hist_deri2_%d",ee));
    numubWS_nosc_spect_hist_diff2[ee]->GetYaxis()->SetTitle("Events(Scaled - Original) 2 esc");
    
    numubWS_nosc_spect_hist_2nd_diff[ee] = new TH1F(*numubWS_nosc_spect_hist);
    numubWS_nosc_spect_hist_2nd_diff[ee] ->SetName(Form("numubWS_nosc_spec_hist_2nd_diff_%d",ee));
    numubWS_nosc_spect_hist_2nd_deri[ee] = new TH1F(*numubWS_nosc_spect_hist);
    numubWS_nosc_spect_hist_2nd_deri[ee] ->SetName(Form("numubWS_nosc_spec_hist_2nd_deri_%d",ee));
    
    for (int i = 1 ; i <= NB_numu110 ; i++) {
      //numu
      double binCont1 = numu_nosc_spect_hist_e1[ee]->GetBinContent(i) - numu_nosc_spect_hist->GetBinContent(i);
      numu_nosc_spect_hist_diff1[ee]->SetBinContent(i,binCont1);
      numu_nosc_spect_hist_deri1[ee]->SetBinContent(i,binCont1/(e));
      double binCont2 = numu_nosc_spect_hist_e2[ee]->GetBinContent(i) - numu_nosc_spect_hist->GetBinContent(i);
      numu_nosc_spect_hist_diff2[ee]->SetBinContent(i,binCont2);
      numu_nosc_spect_hist_deri2[ee]->SetBinContent(i,binCont2/(e+fdel*deltae));
      if(i <= NB_numubar110) {
	//numub
	double binCont1b = numub_nosc_spect_hist_e1[ee]->GetBinContent(i) - numub_nosc_spect_hist->GetBinContent(i);
	numub_nosc_spect_hist_diff1[ee]->SetBinContent(i,binCont1b);
	numub_nosc_spect_hist_deri1[ee]->SetBinContent(i,binCont1b/(e));
	double binCont2b = numub_nosc_spect_hist_e2[ee]->GetBinContent(i) - numub_nosc_spect_hist->GetBinContent(i);
	numub_nosc_spect_hist_diff2[ee]->SetBinContent(i,binCont2b);
	numub_nosc_spect_hist_deri2[ee]->SetBinContent(i,binCont2b/(e+fdel*deltae));
      }
      if(i <= NB_numubarWS110) {
	//numub
	double binCont1bWS = numubWS_nosc_spect_hist_e1[ee]->GetBinContent(i) - numubWS_nosc_spect_hist->GetBinContent(i);
	numubWS_nosc_spect_hist_diff1[ee]->SetBinContent(i,binCont1bWS);
	numubWS_nosc_spect_hist_deri1[ee]->SetBinContent(i,binCont1bWS/(e));
	double binCont2bWS = numubWS_nosc_spect_hist_e2[ee]->GetBinContent(i) - numubWS_nosc_spect_hist->GetBinContent(i);
	numubWS_nosc_spect_hist_diff2[ee]->SetBinContent(i,binCont2bWS);
	numubWS_nosc_spect_hist_deri2[ee]->SetBinContent(i,binCont2bWS/(e+fdel*deltae));
      }
    }
    
    //for (int i = 1 ; i <= NB_numu110 ; i++) {
    //    //double binCont3 = nosc_spect_hist_deri2[ee]->GetBinContent(i) - nosc_spect_hist_deri1[ee]->GetBinContent(i);
    //    double binCont3 = nosc_spect_hist_deri2[ee]->GetBinContent(i) - nosc_spect_hist_deri1[3]->GetBinContent(i);
    //    nosc_spect_hist_2nd_diff[ee]->SetBinContent(i,binCont3);
    //    nosc_spect_hist_2nd_deri[ee]->SetBinContent(i,binCont3/(e));	    
    // }
    //numu -------------------------------------------------
    numu_nosc_spect_hist_diff1[ee]->SetLineColor(kAzure+ee);
    numu_nosc_spect_hist_deri1[ee]->SetLineColor(kAzure+ee);
    numu_nosc_spect_hist_diff2[ee]->SetLineColor(kGreen+ee);
    numu_nosc_spect_hist_deri2[ee]->SetLineColor(kGreen+ee);
    //numub -------------------------------------------------
    numub_nosc_spect_hist_diff1[ee]->SetLineColor(kAzure+ee);
    numub_nosc_spect_hist_deri1[ee]->SetLineColor(kAzure+ee);
    numub_nosc_spect_hist_diff2[ee]->SetLineColor(kGreen+ee);
    numub_nosc_spect_hist_deri2[ee]->SetLineColor(kGreen+ee);
    //numub WS ------------------------------------------------
    numubWS_nosc_spect_hist_diff1[ee]->SetLineColor(kAzure+ee);
    numubWS_nosc_spect_hist_deri1[ee]->SetLineColor(kAzure+ee);
    numubWS_nosc_spect_hist_diff2[ee]->SetLineColor(kGreen+ee);
    numubWS_nosc_spect_hist_deri2[ee]->SetLineColor(kGreen+ee);
    
  }
  
  
  for (int ee = 0 ; ee < de ; ee++) {
    e = lowe + ee*deltae;
    
    for (int i = 1 ; i <= NB_numu110 ; i++) {
      //numu
      double binCont3 = numu_nosc_spect_hist_deri2[ee]->GetBinContent(i) - numu_nosc_spect_hist_deri1[ee]->GetBinContent(i);
      numu_nosc_spect_hist_2nd_diff[ee]->SetBinContent(i,binCont3);
      numu_nosc_spect_hist_2nd_deri[ee]->SetBinContent(i,binCont3/(fdel*deltae));
      if(i <= NB_numubar110) {
	//numub
	double binCont3b = numub_nosc_spect_hist_deri2[ee]->GetBinContent(i) - numub_nosc_spect_hist_deri1[ee]->GetBinContent(i);
	numub_nosc_spect_hist_2nd_diff[ee]->SetBinContent(i,binCont3b);
	numub_nosc_spect_hist_2nd_deri[ee]->SetBinContent(i,binCont3b/(fdel*deltae));
      }
      if(i <= NB_numubarWS110) {
	//numub WS
	double binCont3bWS = numubWS_nosc_spect_hist_deri2[ee]->GetBinContent(i) - numubWS_nosc_spect_hist_deri1[ee]->GetBinContent(i);
	numubWS_nosc_spect_hist_2nd_diff[ee]->SetBinContent(i,binCont3bWS);
	numubWS_nosc_spect_hist_2nd_deri[ee]->SetBinContent(i,binCont3bWS/(fdel*deltae));
      }
    }
    
    //numu
    numu_nosc_spect_hist_2nd_diff[ee]->SetLineColor(kRed-ee);
    numu_nosc_spect_hist_2nd_deri[ee]->SetLineColor(kRed-ee);
    //numub
    numub_nosc_spect_hist_2nd_diff[ee]->SetLineColor(kRed-ee);
    numub_nosc_spect_hist_2nd_deri[ee]->SetLineColor(kRed-ee);
    //numub WS
    numubWS_nosc_spect_hist_2nd_diff[ee]->SetLineColor(kRed-ee);
    numubWS_nosc_spect_hist_2nd_deri[ee]->SetLineColor(kRed-ee);
    
  }
  
  //numu
  TF1 *fFit4_1e;
  fFit4_1e = new TF1("fFit4_1e","[0]+[1]*x+[2]*(x**2)+[3]*(x**3)+[4]*(x**4)",xbins_numu110[0],xbins_numu110[NB_numu110]);
  fFit4_1e->SetParameters(0.0,0.0,0.0,0.0,0.0);
  numu_nosc_spect_hist_deri1[1]->Fit("fFit4_1e","r");
  fFit4_1e->SetLineStyle(2);
  
  TF1 *fFit4_2e;
  fFit4_2e = new TF1("fFit4_2e","[0]+[1]*x+[2]*(x**2)+[3]*(x**3)+[4]*(x**4)",xbins_numu110[0],xbins_numu110[NB_numu110]);
  fFit4_2e->SetParameters(0.0,0.0,0.0,0.0,0.0);
  numu_nosc_spect_hist_deri2[1]->Fit("fFit4_2e","r");
  fFit4_2e->SetLineStyle(2);
  
  TF1 *fFit7_1e;
  fFit7_1e = new TF1("fFit7_1e","[0]+[1]*x+[2]*(x**2)+[3]*(x**3)+[4]*(x**4)+[5]*(x**5)+[6]*(x**6)+[7]*(x**7)",xbins_numu110[0],xbins_numu110[NB_numu110]);
  fFit7_1e->SetParameters(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0);
  numu_nosc_spect_hist_deri1[1]->Fit("fFit7_1e","r");
  
  TF1 *fFit7_2e;
  fFit7_2e = new TF1("fFit7_2e","[0]+[1]*x+[2]*(x**2)+[3]*(x**3)+[4]*(x**4)+[5]*(x**5)+[6]*(x**6)+[7]*(x**7)",xbins_numu110[0],xbins_numu110[NB_numu110]);
  fFit7_2e->SetParameters(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0);
  numu_nosc_spect_hist_deri2[1]->Fit("fFit7_2e","r");
  
  TF1 *fFitGD_1e;
  fFitGD_1e = new TF1("fFitGD_1e","[0]+[1]*exp(-0.5*pow((x-[2])/[3],2))+[4]*exp(-0.5*pow((x-[5])/[6],2))",xbins_numu110[0],xbins_numu110[NB_numu110]);
  fFitGD_1e->SetParameters(0.0,-2e6,3.0, 2.0,1e6,5.0,2.0);
  fFitGD_1e->SetLineColor(kBlack);
  numu_nosc_spect_hist_deri1[5]->Fit("fFitGD_1e","r");
  
  TF1 *fFitGD_2e;
  fFitGD_2e = new TF1("fFitGD_2e","[0]+[1]*exp(-0.5*pow((x-[2])/[3],2))+[4]*exp(-0.5*pow((x-[5])/[6],2))",xbins_numu110[0],xbins_numu110[NB_numu110]);
  fFitGD_2e->SetParameters(0.0,-2e6,3.0, 2.0,1e6,5.0,2.0);
  fFitGD_2e->SetLineColor(kBlack);
  numu_nosc_spect_hist_deri2[5]->Fit("fFitGD_2e","r");
  
  TF1 *fFitGD_2nd;
  fFitGD_2nd = new TF1("fFitGD_2nd","[0]+[1]*exp(-0.5*pow((x-[2])/[3],2))+[4]*exp(-0.5*pow((x-[5])/[6],2))",xbins_numu110[0],xbins_numu110[NB_numu110]);
  fFitGD_2nd->SetParameters(0.0,-2e6,3.0, 2.0,1e6,5.0,2.0);
  fFitGD_2nd->SetLineColor(kBlack);
  numu_nosc_spect_hist_2nd_deri[5]->Fit("fFitGD_2nd","r");
  
  //numub
  TF1 *fbFitGD_1e;
  fbFitGD_1e = new TF1("fbFitGD_1e","[0]+[1]*exp(-0.5*pow((x-[2])/[3],2))+[4]*exp(-0.5*pow((x-[5])/[6],2))",xbins_numub110[0],xbins_numub110[NB_numubar110]);
  fbFitGD_1e->SetParameters(0.0,-2e6,3.0, 2.0,1e6,5.0,2.0);
  fbFitGD_1e->SetLineColor(kBlack);
  numub_nosc_spect_hist_deri1[4]->Fit("fbFitGD_1e","r");
  
  TF1 *fbFitGD_2e;
  fbFitGD_2e = new TF1("fbFitGD_2e","[0]+[1]*exp(-0.5*pow((x-[2])/[3],2))+[4]*exp(-0.5*pow((x-[5])/[6],2))",xbins_numub110[0],xbins_numub110[NB_numubar110]);
  fbFitGD_2e->SetParameters(0.0,-2e6,3.0, 2.0,1e6,5.0,2.0);
  fbFitGD_2e->SetLineColor(kBlack);
  numub_nosc_spect_hist_deri2[4]->Fit("fbFitGD_2e","r");
  
  TF1 *fbFitGD_2nd;
  fbFitGD_2nd = new TF1("fbFitGD_2nd","[0]+[1]*exp(-0.5*pow((x-[2])/[3],2))+[4]*exp(-0.5*pow((x-[5])/[6],2))",xbins_numub110[0],xbins_numub110[NB_numubar110]);
  fbFitGD_2nd->SetParameters(0.0,2e6,3.0, 2.0,-1e6,5.0,2.0);
  fbFitGD_2nd->SetLineColor(kBlack);
  numub_nosc_spect_hist_2nd_deri[4]->Fit("fbFitGD_2nd","r");
  
  //numub WS
  TF1 *fbWSFitGD_1e;
  fbWSFitGD_1e = new TF1("fbWSFitGD_1e","[0]+[1]*exp(-0.5*pow((x-[2])/[3],2))+[4]*exp(-0.5*pow((x-[5])/[6],2))",xbins_numubWS110[0],xbins_numubWS110[NB_numubarWS110]);
  fbWSFitGD_1e->SetParameters(0.0,-2e6,3.0, 2.0,1e6,5.0,2.0);
  fbWSFitGD_1e->SetLineColor(kBlack);
  numubWS_nosc_spect_hist_deri1[4]->Fit("fbWSFitGD_1e","r");
  
  TF1 *fbWSFitGD_2e;
  fbWSFitGD_2e = new TF1("fbWSFitGD_2e","[0]+[1]*exp(-0.5*pow((x-[2])/[3],2))+[4]*exp(-0.5*pow((x-[5])/[6],2))",xbins_numubWS110[0],xbins_numubWS110[NB_numubarWS110]);
  fbWSFitGD_2e->SetParameters(0.0,-2e6,3.0, 2.0,1e6,5.0,2.0);
  fbWSFitGD_2e->SetLineColor(kBlack);
  numubWS_nosc_spect_hist_deri2[4]->Fit("fbWSFitGD_2e","r");
  
  TF1 *fbWSFitGD_2nd;
  fbWSFitGD_2nd = new TF1("fbWSFitGD_2nd","[0]+[1]*exp(-0.5*pow((x-[2])/[3],2))+[4]*exp(-0.5*pow((x-[5])/[6],2))",xbins_numubWS110[0],xbins_numubWS110[NB_numubarWS110]);
  fbWSFitGD_2nd->SetParameters(0.0, 2e6,2.0,5.0,-1e5,7.0,6.0);
  fbWSFitGD_2nd->SetLineColor(kBlack);
  numubWS_nosc_spect_hist_2nd_deri[5]->Fit("fbWSFitGD_2nd","r");
  
  TFile *fEscale = new TFile(filePath + "/data/minos_EScaleDerivative.root","recreate");
  //numu
  //fFit4_1e    ->Write();
  //fFit7_1e    ->Write();
  fFitGD_1e    ->Write();
  //fFit4_2e    ->Write();
  //fFit7_2e    ->Write();
  fFitGD_2e    ->Write();
  fFitGD_2nd   ->Write();
  //numub
  fbFitGD_1e   ->Write();
  fbFitGD_2e   ->Write();
  fbFitGD_2nd  ->Write();
  //numub WS
  fbWSFitGD_1e ->Write();
  fbWSFitGD_2e ->Write();
  fbWSFitGD_2nd->Write();
  
  numu_nosc_spect_hist   ->Write();
  numub_nosc_spect_hist  ->Write();
  numubWS_nosc_spect_hist->Write();
  for (int i = 0; i < de; i++){
    //numu
    numu_nosc_spect_hist_e1[i]   ->Write();
    numu_nosc_spect_hist_diff1[i]->Write();
    numu_nosc_spect_hist_deri1[i]->Write();
    numu_nosc_spect_hist_e2[i]   ->Write();
    numu_nosc_spect_hist_diff2[i]->Write();
    numu_nosc_spect_hist_deri2[i]->Write();
    
    numu_nosc_spect_hist_2nd_diff[i]->Write();
    numu_nosc_spect_hist_2nd_deri[i]->Write();
    
    //numub
    numub_nosc_spect_hist_e1[i]   ->Write();
    numub_nosc_spect_hist_diff1[i]->Write();
    numub_nosc_spect_hist_deri1[i]->Write();
    numub_nosc_spect_hist_e2[i]   ->Write();
    numub_nosc_spect_hist_diff2[i]->Write();
    numub_nosc_spect_hist_deri2[i]->Write();
    
    numub_nosc_spect_hist_2nd_diff[i]->Write();
    numub_nosc_spect_hist_2nd_deri[i]->Write();
    
    //numub WS
    numubWS_nosc_spect_hist_e1[i]   ->Write();
    numubWS_nosc_spect_hist_diff1[i]->Write();
    numubWS_nosc_spect_hist_deri1[i]->Write();
    numubWS_nosc_spect_hist_e2[i]   ->Write();
    numubWS_nosc_spect_hist_diff2[i]->Write();
    numubWS_nosc_spect_hist_deri2[i]->Write();
    
    numubWS_nosc_spect_hist_2nd_diff[i]->Write();
    numubWS_nosc_spect_hist_2nd_deri[i]->Write();
  }
  
  //-- Plotting section -----------------------//
  //-- Defining histogram frames --//
  double max_numu = 0.0, min_numu = 0.0;
  double max_numub = 0.0, min_numub = 0.0;
  double max_numubWS = 0.0, min_numubWS = 0.0;
  for (int ee = 0 ; ee < de ; ee++) {
    //numu
    double max_temp = TMath::Max(numu_nosc_spect_hist->GetMaximum(),numu_nosc_spect_hist_e1[ee]->GetMaximum());
    if (max_temp >= max_numu)
      max_numu = max_temp;
    double min_temp = TMath::Min(numu_nosc_spect_hist->GetMinimum(),numu_nosc_spect_hist_e1[ee]->GetMinimum());
    if (min_temp <= min_numu)
      min_numu = min_temp;
    //numub
    max_temp = TMath::Max(numub_nosc_spect_hist->GetMaximum(),numub_nosc_spect_hist_e1[ee]->GetMaximum());
    if (max_temp >= max_numub)
      max_numub = max_temp;
    min_temp = TMath::Min(numub_nosc_spect_hist->GetMinimum(),numub_nosc_spect_hist_e1[ee]->GetMinimum());
    if (min_temp <= min_numub)
      min_numub = min_temp;
    //numub WS
    max_temp = TMath::Max(numubWS_nosc_spect_hist->GetMaximum(),numubWS_nosc_spect_hist_e1[ee]->GetMaximum());
    if (max_temp >= max_numubWS)
      max_numubWS = max_temp;
    min_temp = TMath::Min(numubWS_nosc_spect_hist->GetMinimum(),numubWS_nosc_spect_hist_e1[ee]->GetMinimum());
    if (min_temp <= min_numubWS)
      min_numubWS = min_temp;
  }
  TH2F *framenumu    = new TH2F("framenumu",   "",NB_numu110,     xbins_numu110,   10,min_numu,   1.01*max_numu);
  TH2F *framenumub   = new TH2F("framenumub",  "",NB_numubar110,  xbins_numub110,  10,min_numub,  1.01*max_numub);
  TH2F *framenumubWS = new TH2F("framenumubWS","",NB_numubarWS110,xbins_numubWS110,10,min_numubWS,1.01*max_numubWS);
  
  double max_numudif    = 0.0, min_numudif    = 0.0;
  double max_numubdif   = 0.0, min_numubdif   = 0.0;
  double max_numubWSdif = 0.0, min_numubWSdif = 0.0;
  for (int ee = 0 ; ee < de ; ee++) {
    //numu
    double max_temp = TMath::Max(max_numudif,numu_nosc_spect_hist_diff1[ee]->GetMaximum());
    if (max_temp >= max_numudif)
      max_numudif = max_temp;
    double min_temp = TMath::Min(min_numudif,numu_nosc_spect_hist_diff1[ee]->GetMinimum());
    if (min_temp <= min_numudif)
      min_numudif = min_temp;
    //numub
    max_temp = TMath::Max(max_numubdif,numub_nosc_spect_hist_diff1[ee]->GetMaximum());
    if (max_temp >= max_numubdif)
      max_numubdif = max_temp;
    min_temp = TMath::Min(min_numubdif,numub_nosc_spect_hist_diff1[ee]->GetMinimum());
    if (min_temp <= min_numubdif)
      min_numubdif = min_temp;
    //numub WS
    max_temp = TMath::Max(max_numubWSdif,numubWS_nosc_spect_hist_diff1[ee]->GetMaximum());
    if (max_temp >= max_numubWSdif)
      max_numubWSdif = max_temp;
    min_temp = TMath::Min(min_numubWSdif,numubWS_nosc_spect_hist_diff1[ee]->GetMinimum());
    if (min_temp <= min_numubWSdif)
      min_numubWSdif = min_temp;
  }
  TH2F *framenumudif    = new TH2F("framenumudif",   "",NB_numu110,     xbins_numu110,   10,1.1*min_numudif,   1.1*max_numudif);
  TH2F *framenumubdif   = new TH2F("framenumubdif",  "",NB_numubar110,  xbins_numub110,  10,1.1*min_numubdif,  1.1*max_numubdif);
  TH2F *framenumubWSdif = new TH2F("framenumubWSdif","",NB_numubarWS110,xbins_numubWS110,10,1.1*min_numubWSdif,1.1*max_numubWSdif);
  
  double max_numuder    = 0.0, min_numuder    = 0.0;
  double max_numubder   = 0.0, min_numubder   = 0.0;
  double max_numubWSder = 0.0, min_numubWSder = 0.0;
  for (int ee = 0 ; ee < de ; ee++) {
    //numu
    double max_temp = TMath::Max(max_numuder,numu_nosc_spect_hist_deri1[ee]->GetMaximum());
    if (max_temp >= max_numuder)
      max_numuder = max_temp;
    double min_temp = TMath::Min(min_numuder,numu_nosc_spect_hist_deri1[ee]->GetMinimum());
    if (min_temp <= min_numuder)
      min_numuder = min_temp;
    //numub
    max_temp = TMath::Max(max_numubder,numub_nosc_spect_hist_deri1[ee]->GetMaximum());
    if (max_temp >= max_numubder)
      max_numubder = max_temp;
    min_temp = TMath::Min(min_numubder,numub_nosc_spect_hist_deri1[ee]->GetMinimum());
    if (min_temp <= min_numubder)
      min_numubder = min_temp;
    //numub WS
    max_temp = TMath::Max(max_numubWSder,numubWS_nosc_spect_hist_deri1[ee]->GetMaximum());
    if (max_temp >= max_numubWSder)
      max_numubWSder = max_temp;
    min_temp = TMath::Min(min_numubWSder,numubWS_nosc_spect_hist_deri1[ee]->GetMinimum());
    if (min_temp <= min_numubWSder)
      min_numubWSder = min_temp;
  }
  TH2F *framenumuder    = new TH2F("framenumuder",   "",NB_numu110,     xbins_numu110,   10,1.1*min_numuder,   1.1*max_numuder);
  TH2F *framenumubder   = new TH2F("framenumubder",  "",NB_numubar110,  xbins_numub110,  10,1.1*min_numubder,  1.1*max_numubder);
  TH2F *framenumubWSder = new TH2F("framenumubWSder","",NB_numubarWS110,xbins_numubWS110,10,1.1*min_numubWSder,1.1*max_numubWSder);
  //-- END - Defining histogram frames
  
  //numu
  TCanvas *canv1 = new TCanvas("canv1", "MINOS", 1*600,2*400);
  canv1->Divide(1,2);
  //Muon neutrinos
  canv1->cd(1);
  framenumu->Draw();
  numu_nosc_spect_hist->Draw("same hist");
  for (int ee = 0 ; ee < de ; ee++){
    numu_nosc_spect_hist_e1[ee]->Draw("same hist");
    numu_nosc_spect_hist_e2[ee]->Draw("same hist");}
  
  canv1->cd(2);
  framenumudif->Draw();
  for (int ee = 0 ; ee < de ; ee++){
    numu_nosc_spect_hist_diff1[ee]->Draw("same hist");
    numu_nosc_spect_hist_diff1[ee]->Draw("same hist");}
  
  //-- Derivatives --//
  TCanvas *canv2 = new TCanvas("canv2", "MINOS Derivatives", 1*600,2*400);
  canv2->Divide(1,2);
  //Muon Neutrinos
  canv2->cd(1);
  framenumu->Draw();
  numu_nosc_spect_hist->Draw("same hist");
  for (int ee = 0 ; ee < de ; ee++){
    numu_nosc_spect_hist_e1[ee]->Draw("same hist");
    numu_nosc_spect_hist_e2[ee]->Draw("same hist");}
  canv2->cd(2);
  framenumuder->Draw();
  for (int ee = 0 ; ee < de ; ee++){
    numu_nosc_spect_hist_deri1[ee]->Draw("same hist");
    numu_nosc_spect_hist_deri2[ee]->Draw("same hist");}
  fFit4_1e->Draw("same");
  fFit7_1e->Draw("same");
  fFitGD_1e->Draw("same");
  
  fFit4_2e->Draw("same");
  fFit7_2e->Draw("same");
  fFitGD_2e->Draw("same");
  
  //numub
  TCanvas *canv11 = new TCanvas("canv11", "MINOS", 1*600,2*400);
  canv11->Divide(1,2);
  //Muon anti neutrinos
  canv11->cd(1);
  framenumub->Draw();
  numub_nosc_spect_hist->Draw("same hist");
  for (int ee = 0 ; ee < de ; ee++){
    numub_nosc_spect_hist_e1[ee]->Draw("same hist");
    numub_nosc_spect_hist_e2[ee]->Draw("same hist");}
  
  canv11->cd(2);
  framenumubdif->Draw();
  for (int ee = 0 ; ee < de ; ee++){
    numub_nosc_spect_hist_diff1[ee]->Draw("same hist");
    numub_nosc_spect_hist_diff1[ee]->Draw("same hist");}
  
  //-- Derivatives --//
  TCanvas *canv22 = new TCanvas("canv22", "MINOS Derivatives", 1*600,2*400);
  canv22->Divide(1,2);
  //Muon Anti Neutrinos
  canv22->cd(1);
  framenumub->Draw();
  numub_nosc_spect_hist->Draw("same hist");
  for (int ee = 0 ; ee < de ; ee++){
    numub_nosc_spect_hist_e1[ee]->Draw("same hist");
    numub_nosc_spect_hist_e2[ee]->Draw("same hist");}
  canv2->cd(2);
  framenumubder->Draw();
  for (int ee = 0 ; ee < de ; ee++){
    numub_nosc_spect_hist_deri1[ee]->Draw("same hist");
    numub_nosc_spect_hist_deri2[ee]->Draw("same hist");}
  fbFitGD_1e->Draw("same");
  fbFitGD_2e->Draw("same");
  
  //numub WS
  TCanvas *canv111 = new TCanvas("canv111", "MINOS", 1*600,2*400);
  canv111->Divide(1,2);
  //Muon anti neutrinos WS
  canv111->cd(1);
  framenumubWS->Draw();
  numubWS_nosc_spect_hist->Draw("same hist");
  for (int ee = 0 ; ee < de ; ee++){
    numubWS_nosc_spect_hist_e1[ee]->Draw("same hist");
    numubWS_nosc_spect_hist_e2[ee]->Draw("same hist");}
  
  canv111->cd(2);
  framenumubWSdif->Draw();
  for (int ee = 0 ; ee < de ; ee++){
    numubWS_nosc_spect_hist_diff1[ee]->Draw("same hist");
    numubWS_nosc_spect_hist_diff1[ee]->Draw("same hist");}
  
  //-- Derivatives --//
  TCanvas *canv222 = new TCanvas("canv222", "MINOS Derivatives", 1*600,2*400);
  canv222->Divide(1,2);
  //Muon Anti Neutrinos WS
  canv222->cd(1);
  framenumubWS->Draw();
  numubWS_nosc_spect_hist->Draw("same hist");
  for (int ee = 0 ; ee < de ; ee++){
    numubWS_nosc_spect_hist_e1[ee]->Draw("same hist");
    numubWS_nosc_spect_hist_e2[ee]->Draw("same hist");}
  canv2->cd(2);
  framenumubWSder->Draw();
  for (int ee = 0 ; ee < de ; ee++){
    numubWS_nosc_spect_hist_deri1[ee]->Draw("same hist");
    numubWS_nosc_spect_hist_deri2[ee]->Draw("same hist");}
  fbWSFitGD_1e->Draw("same");
  fbWSFitGD_2e->Draw("same");
  
  std::cout << "Succesful run!!" << endl;
  
  return 0;
  
}//-END
