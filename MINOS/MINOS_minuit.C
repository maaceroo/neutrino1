//----------------------------------------------------------//
//-- db_minuit_spec.C                          2021-10-22 --//
//-- Authors: M.A. Acero O. & A.A. Aguilar A              --//
//----------------------------------------------------------//
//----------------------------------------------------------//
//-- Analysis of the MINOS data from PRL110,251801 (2013) --//
//-- This includes muon neutrinos and antineutrinos from  --//
//-- the NuMI beam.                                       --//
//----------------------------------------------------------//
//-- This macro can be executed under ROOT typing         --//
//-- "root[0] .x MINOS_minuit.C"                          --//
//-- The results are printed in a file, with the following--//
//-- information:  sin^2(2th23) delta(m32)^2 chi^2(min)   --//
//----------------------------------------------------------//
//--             Last update: 2023 - 08 - 15              --//
//---****************************************************---//
//---------------------- CONSTANTS -------------------------//
//---****************************************************---//
#include "constants.h"
// histogram binning for the simulated data
//#define  NB 26
//#define  lo 0.7
//#define  hi 12.0
//For the sin^2(2th_23) loop
//#define N_s2t  5                           //number of points in the grid
//#define lo_s2t 0.75                         //sin^2(2th_13) min
//#define lo_s2t 0.0                          //sin^2(2th_13) min
//#define hi_s2t 1.0                          //sin^2(2th_13) max
//For the delta(m32)^2 loop
//#define N_dm2  5                           //number of points in the grid
//#define lo_dm2 1.0e-3                       //delta(m31)^2 min
//#define hi_dm2 4.0e-3                       //delta(m31)^2 max
//---*****************************************************---//

//---*****************************************************---//
//-- Global variables and quantities ------------------------//
//---*****************************************************---//
//---*****************************************************---//
double s2th_23;         //oscillation parameter to be fitted
double dm2_32;          //oscillation parameter to be fitted
double eps;             //normalization factor
double epsNC;           //NC background pull term
double ensc;            //Energy sclae pull Term

double spcNumu[NB_numu110];     //Energy scale numu
double spcNumuNoOsc[NB_numu110];     //Energy scale numu
double spcNumuB[NB_numubar110]; //Energy scale numu bar
double spcNumuBNoOsc[NB_numubar110]; //Energy scale numu bar
double spcNumuBWS[NB_numubarWS110]; //Energy scale numu bar
double spcNumuBWSNoOsc[NB_numubarWS110]; //Energy scale numu bar
double NoscTot;
double noNoscTot;
double NoscTotB;
double noNoscTotB;
double NoscTotBWS;
double noNoscTotBWS;
double totalBgd = 8.92;  // integral of bkgd spectrum
double totalBgdB = 2.055;  // integral of bkgd spectrum
double totalBgdBWS = 2.055;  // integral of bkgd spectrum **** TO BE CORRECTED ****

double xbins_numu110[NB_numu110+1];
double xbins_numubar110[NB_numubar110+1];
double xbins_numubarWS110[NB_numubarWS110+1];

//Table I PRL110.251801(2013)
//{No Osc (Simulated), Oscilated (BF), Observed}
double NuMu_Events[3]    = {3201.0, 2543.0, 2579.0};
double NuMuB_Events[3]   = {313.0,   227.0,  226.0};
double NuMuBWS_Events[3] = {363.0,   324.0,  312.0};
//---*****************************************************---//
TH1F *numu_data_spect_histo;
TH1F *numu_bkgd_spect_histo;
TH1F *numub_data_spect_histo;
TH1F *numub_bkgd_spect_histo;
TH1F *numubWS_data_spect_histo;
TH1F *numubWS_bkgd_spect_histo;
//---*****************************************************---//
TH1F *numu_noosc_spect_hist;
TH1F *numu_nu_nosc_spect_hist;
TH1F *numu_wosc_spect_hist;
TH1F *numu_nu_wosc_spect_hist;
TH1F *numu_bfit_spect_hist;

TH1F *numub_noosc_spect_hist;
TH1F *numub_nu_nosc_spect_hist;
TH1F *numub_wosc_spect_hist;
TH1F *numub_nu_wosc_spect_hist;
TH1F *numub_bfit_spect_hist;

TH1F *numubWS_noosc_spect_hist;
TH1F *numubWS_nu_nosc_spect_hist;
TH1F *numubWS_wosc_spect_hist;
TH1F *numubWS_nu_wosc_spect_hist;
TH1F *numubWS_bfit_spect_hist;
//-------------------
//TF1 *fFit4;
//TF1 *fFit7;
TF1 *fFitGD_1st_der;
TF1 *fFitGD_2nd_der;

TF1 *fbFitGD_1st_der;
TF1 *fbFitGD_2nd_der;

TF1 *fbWSFitGD_1st_der;
TF1 *fbWSFitGD_2nd_der;

int SELECTOR; // 0-numu, 1-numub, 2-numubWS

//---*****************************************************---//
//-- Chi-square function to be minimized --------------------//
//-- It has 3 pull parameters and 2 oscillation parameter  --//
//-- which are to be fitted.                                -//
//-- Likelihood function from: J. Mitchell, University of  --//
//-- Cmabridge, PhD. Thesis (2011)                         --//
//---*****************************************************---//
double chi2(const double *xx)
{
  //Pull parameters considered in the chi² function
  eps   = xx[0];
  epsNC = xx[1];
  ensc  = xx[2];
  //---*****************************************************---//
  double SurvPavg     = 0.0;
  double nll          = 0.0; //nll = -2ln(L)
  double nllnu        = 0.0; //nll = -2ln(L)
  double nllnub       = 0.0; //nll = -2ln(L)  
  double nllnubWS     = 0.0; //nll = -2ln(L)  
  double Nd           = 0.0;
  double Nmc          = 0.0;
  //Section 6.2, J. Mitchell PhD Thesis https://minos-docdb.fnal.gov/cgi-bin/sso/RetrieveFile?docid=8704&filename=jessmitchellthesis.pdf&version=1
  //See also https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.106.181801
  double sigeps       = 0.016;
  double sigepsNC     = 0.20;
  double sigensc      = 0.06;
  
  double spcNumuNew[NB_numu110];  
  double spcNumuBNew[NB_numubar110];
  double spcNumuBWSNew[NB_numubarWS110];  
    
  //--- numu section
  SurvPavg = NoscTot/noNoscTot;
  for (int iBIN = 0 ; iBIN < NB_numu110 ; iBIN++)
    {
      //-- Measured events of the Far Detector (background is substracted)
      double numu_data_sigplusbg = (numu_data_spect_histo->GetBinContent(iBIN+1))*NuMu_Events[2];
      double numu_simu_bg        = (numu_bkgd_spect_histo->GetBinContent(iBIN+1))*totalBgd;
      Nd = ( numu_data_sigplusbg - (1-epsNC)*numu_simu_bg );
      
      //-- Predicted Number of events from neutrino oscillations of the Far Detector
      
      //-- Energy resolution uncertainty implementation
      double binCenter, delta_spc, avgSurvProb_bin;
      binCenter = 0.5*(xbins_numu110[iBIN] + xbins_numu110[iBIN+1]);
      //delta_spc = ensc*(fFitGD->Eval(binCenter));
      delta_spc = ensc*(fFitGD_1st_der->Eval(binCenter)) + 0.5*pow(ensc,2)*(fFitGD_2nd_der->Eval(binCenter));
      avgSurvProb_bin = spcNumu[iBIN]/spcNumuNoOsc[iBIN];
      spcNumuNew[iBIN] = spcNumu[iBIN] + 2.0*delta_spc*avgSurvProb_bin;
      
      Nmc = ( (spcNumuNew[iBIN]/NoscTot)*(NuMu_Events[0]*SurvPavg) )*(1-eps);
      Nmc = fudge2*Nmc;   //2022-04-01 - MAAO
      //cout << "iBin = " << iBIN << "  Nmc = " << Nmc << " Nd = " << Nd << endl;
      
      nllnu += 2*(Nmc - Nd + Nd*log(Nd/Nmc));
    }
  
    nllnu = nllnu + pow(eps,2)/(pow(sigeps,2)) + pow(epsNC,2)/(pow(sigepsNC,2)) + pow(ensc,2)/(pow(sigensc,2));
  //nllnu = nllnu + pow(eps,2)/(2*pow(sigeps,2)) + pow(epsNC,2)/(2*pow(sigepsNC,2)) + pow(ensc,2)/(2*pow(sigensc,2)); 

 //-- numubar section
  SurvPavg = NoscTotB/noNoscTotB;
  for (int iBIN = 0 ; iBIN < NB_numubar110 ; iBIN++)
    {
      //-- Measured events of the Far Detector (background is substracted)
      double numub_data_sigplusbg = (numub_data_spect_histo->GetBinContent(iBIN+1))*NuMuB_Events[2];
      double numub_simu_bg        = (numub_bkgd_spect_histo->GetBinContent(iBIN+1))*totalBgdB;
      Nd = ( numub_data_sigplusbg - (1-epsNC)*numub_simu_bg );
      
      //-- Predicted Number of events from neutrino oscillations of the Far Detector
      
      //-- Energy resolution uncertainty implementation
      double binCenter, delta_spc, avgSurvProb_bin;
      binCenter = 0.5*(xbins_numubar110[iBIN] + xbins_numubar110[iBIN+1]);
      //delta_spc = ensc*(fbFitGD->Eval(binCenter));
      delta_spc = ensc*(fbFitGD_1st_der->Eval(binCenter)) + 0.5*pow(ensc,2)*(fbFitGD_2nd_der->Eval(binCenter));
      avgSurvProb_bin = spcNumuB[iBIN]/spcNumuBNoOsc[iBIN];
      spcNumuBNew[iBIN] = spcNumuB[iBIN] + 2.0*delta_spc*avgSurvProb_bin;
     

      Nmc = ( (spcNumuBNew[iBIN]/NoscTotB)*(NuMuB_Events[0]*SurvPavg) )*(1-eps);
      Nmc = fudgeb2*Nmc;   //2023-03-10 - AAAA
      //cout << "iBin = " << iBIN << "  Nmc = " << Nmc << " Nd = " << Nd << endl;
      
      nllnub += 2*(Nmc - Nd + Nd*log(Nd/Nmc));
    }
  
  nllnub = nllnub + pow(eps,2)/(pow(sigeps,2)) + pow(epsNC,2)/(pow(sigepsNC,2)) + pow(ensc,2)/(pow(sigensc,2));
  //nllnub = nllnub + pow(eps,2)/(2*pow(sigeps,2)) + pow(epsNC,2)/(2*pow(sigepsNC,2)) + pow(ensc,2)/(2*pow(sigensc,2));

 //-- numubarWS section
  SurvPavg = NoscTotBWS/noNoscTotBWS;
  for (int iBIN = 0 ; iBIN < NB_numubarWS110 ; iBIN++)
    {
      //-- Measured events of the Far Detector (background is substracted)
      double numubWS_data_sigplusbg = (numubWS_data_spect_histo->GetBinContent(iBIN+1))*NuMuBWS_Events[2];
      double numubWS_simu_bg        = (numubWS_bkgd_spect_histo->GetBinContent(iBIN+1))*totalBgdBWS;
      Nd = ( numubWS_data_sigplusbg - (1-epsNC)*numubWS_simu_bg );
      
      //-- Predicted Number of events from neutrino oscillations of the Far Detector
      
      //-- Energy resolution uncertainty implementation
      double binCenter, delta_spc, avgSurvProb_bin;
      binCenter = 0.5*(xbins_numubarWS110[iBIN] + xbins_numubarWS110[iBIN+1]);
      //delta_spc = ensc*(fbFitGD->Eval(binCenter));
      delta_spc = ensc*(fbWSFitGD_1st_der->Eval(binCenter)) + 0.5*pow(ensc,2)*(fbWSFitGD_2nd_der->Eval(binCenter));
      avgSurvProb_bin = spcNumuBWS[iBIN]/spcNumuBWSNoOsc[iBIN];
      spcNumuBWSNew[iBIN] = spcNumuBWS[iBIN] + 2.0*delta_spc*avgSurvProb_bin;
     

      Nmc = ( (spcNumuBWSNew[iBIN]/NoscTotBWS)*(NuMuBWS_Events[0]*SurvPavg) )*(1-eps);
      Nmc = fudgebWS2*Nmc;
      //cout << "iBin = " << iBIN << "  Nmc = " << Nmc << " Nd = " << Nd << endl;
      
      nllnubWS += 2*(Nmc - Nd + Nd*log(Nd/Nmc));
    }
  
  nllnubWS = nllnubWS + pow(eps,2)/(pow(sigeps,2)) + pow(epsNC,2)/(pow(sigepsNC,2)) + pow(ensc,2)/(pow(sigensc,2));
  //nllnubWS = nllnub + pow(eps,2)/(2*pow(sigeps,2)) + pow(epsNC,2)/(2*pow(sigepsNC,2)) + pow(ensc,2)/(2*pow(sigensc,2));
  
  if      (SELECTOR == 0) nll = nllnu;
  else if (SELECTOR == 1) nll = nllnub;
  else if (SELECTOR == 2) nll = nllnubWS;

  return nll;
}
//---*****************************************************---//

int MINOS_minuit(const char * minName = "Minuit",
		 const char *algoName = "" ,
		 int randomSeed = +10)
//int randomSeed = -1)
{
  cout << "Let's the minimization begins..." << endl;
  
  TString filePath = dirName;
  //-------------------
  //-- Energy scale derivative function
  //-------------------
  TFile *Esc_File = new TFile(filePath + "/data/minos_EScaleDerivative.root","READ");
  //fFit4  = (TF1*)(Esc_File->Get("fFit4_1e"));
  //fFit7  = (TF1*)(Esc_File->Get("fFit7_1e"));
  fFitGD_1st_der = (TF1*)(Esc_File->Get("fFitGD_2e"));
  fFitGD_2nd_der = (TF1*)(Esc_File->Get("fFitGD_2nd"));

  fbFitGD_1st_der = (TF1*)(Esc_File->Get("fbFitGD_2e"));
  fbFitGD_2nd_der = (TF1*)(Esc_File->Get("fbFitGD_2nd"));

  fbWSFitGD_1st_der = (TF1*)(Esc_File->Get("fbWSFitGD_2e"));
  fbWSFitGD_2nd_der = (TF1*)(Esc_File->Get("fbWSFitGD_2nd"));

  //-------------------
  // Energy Histograms
  //-------------------
  TFile *fenergy = new TFile("./MINOS_spectra_PRL108-PRL110.root","read");
  double bfactor, dfactor;

  //Data-numu
  numu_data_spect_histo = (TH1F*) fenergy->Get("numu110_data_histo");
  dfactor = 1.0/numu_data_spect_histo->Integral();
  numu_data_spect_histo->Scale(dfactor);
  //Data-numubar
  numub_data_spect_histo = (TH1F*) fenergy->Get("numub110_data_histo");
  dfactor = 1.0/numub_data_spect_histo->Integral();
  numub_data_spect_histo->Scale(dfactor);
  //Data-numubarWS
  numubWS_data_spect_histo = (TH1F*) fenergy->Get("numubWS110_data_histo");
  dfactor = 1.0/numubWS_data_spect_histo->Integral();
  numubWS_data_spect_histo->Scale(dfactor);

  //Background-numu
  numu_bkgd_spect_histo = (TH1F*) fenergy->Get("numu110_bkgd_histo");
  bfactor = 1.0/numu_bkgd_spect_histo->Integral();
  numu_bkgd_spect_histo->Scale(bfactor);
  //Background-numubar
  numub_bkgd_spect_histo = (TH1F*) fenergy->Get("numub110_bkgd_histo");
  bfactor = 1.0/numub_bkgd_spect_histo->Integral();
  numub_bkgd_spect_histo->Scale(bfactor);
  //Background-numubar
  numubWS_bkgd_spect_histo = (TH1F*) fenergy->Get("numubWS110_bkgd_histo");
  bfactor = 1.0/numubWS_bkgd_spect_histo->Integral();
  numubWS_bkgd_spect_histo->Scale(bfactor);
    
  // histogram binning
  double delta_bins110 = (9.0 - 0.5)/17.0; // 0.5 GeV/bin
  for (int i = 0 ; i < (NB_numu110 - 5) ; i++)
  xbins_numu110[i]  = lo110 + i*delta_bins110;
  xbins_numu110[18] = xbins_numu110[17] + 0.75;
  xbins_numu110[19] = xbins_numu110[18] + 0.75;
  xbins_numu110[20] = xbins_numu110[19] + 0.75;
  xbins_numu110[21] = xbins_numu110[20] + 0.75;
  xbins_numu110[22] = xbins_numu110[21] + 1.0;
  xbins_numu110[23] = hi110;

  double delta_binsb110 = (12.0 - lob110)/(NB_numubar110-1); // 1.0 GeV/bin
  for (int i = 0 ; i < (NB_numubar110) ; i++)
    xbins_numubar110[i] = lob110 + delta_binsb110*i;
  xbins_numubar110[12] = xbins_numubar110[11] + 2.0;

  double delta_binsbWS110 = (20.0 - lobWS110)/(NB_numubarWS110-1); // 2.0 GeV/bin
  for (int i = 0 ; i < (NB_numubarWS110) ; i++)
    xbins_numubarWS110[i] = lobWS110 + delta_binsbWS110*i;
  xbins_numubarWS110[11] = xbins_numubarWS110[10] + 5.0;
  
  numu_noosc_spect_hist   = new TH1F("numu_noosc_spect_hist",  "",NB_numu110,xbins_numu110);
  numu_nu_nosc_spect_hist = new TH1F("numu_nu_nosc_spect_hist","",NB_numu110,xbins_numu110);
  numu_wosc_spect_hist    = new TH1F("numu_wosc_spect_hist",   "",NB_numu110,xbins_numu110);
  numu_nu_wosc_spect_hist = new TH1F("numu_nu_wosc_spect_hist","",NB_numu110,xbins_numu110);
  numu_bfit_spect_hist    = new TH1F("numu_bfit_spect_hist",   "",NB_numu110,xbins_numu110);
  
  numub_noosc_spect_hist   = new TH1F("numub_noosc_spect_hist",  "",NB_numubar110,xbins_numubar110);
  numub_nu_nosc_spect_hist = new TH1F("numub_nu_nosc_spect_hist","",NB_numubar110,xbins_numubar110);
  numub_wosc_spect_hist    = new TH1F("numub_wosc_spect_hist",   "",NB_numubar110,xbins_numubar110);
  numub_nu_wosc_spect_hist = new TH1F("numub_nu_wosc_spect_hist","",NB_numubar110,xbins_numubar110);
  numub_bfit_spect_hist    = new TH1F("numub_bfit_spect_hist",   "",NB_numubar110,xbins_numubar110);

  numubWS_noosc_spect_hist   = new TH1F("numubWS_noosc_spect_hist",  "",NB_numubarWS110,xbins_numubar110);
  numubWS_nu_nosc_spect_hist = new TH1F("numubWS_nu_nosc_spect_hist","",NB_numubarWS110,xbins_numubar110);
  numubWS_wosc_spect_hist    = new TH1F("numubWS_wosc_spect_hist",   "",NB_numubarWS110,xbins_numubar110);
  numubWS_nu_wosc_spect_hist = new TH1F("numubWS_nu_wosc_spect_hist","",NB_numubarWS110,xbins_numubar110);
  numubWS_bfit_spect_hist    = new TH1F("numubWS_bfit_spect_hist",   "",NB_numubarWS110,xbins_numubar110);

  cout << "binning  -> Done..." << endl;
  
  cout << "Minimization settings..." << endl;
  ROOT::Math::Minimizer* min1 = ROOT::Math::Factory::CreateMinimizer(minName, algoName);
  
  //-- Set tolerance , etc...
  min1->SetMaxFunctionCalls(1000000); // for Minuit/Minuit2
  min1->SetMaxIterations(10000);      // for GSL
  min1->SetTolerance(0.001);
  min1->SetPrintLevel(-1);

  ROOT::Math::Minimizer* min2 = ROOT::Math::Factory::CreateMinimizer(minName, algoName);
  min2->SetMaxFunctionCalls(1000000); // for Minuit/Minuit2
  min2->SetMaxIterations(10000);      // for GSL
  min2->SetTolerance(0.001);
  min2->SetPrintLevel(-1);
  
  //-- File to print oscillation parameters and chi2 values
  string chi2surfname;
  //numu
  ofstream numu_chi2Surface_file;
  chi2surfname = "/data/numu_chi2_s2t-dm2_surface.txt";  //(sin^2(2th13), dm2, chi^2_min)
  numu_chi2Surface_file.open(filePath+(chi2surfname).c_str());
  //numubar
  ofstream numub_chi2Surface_file;
  chi2surfname = "/data/numub_chi2_s2t-dm2_surface.txt";  //(sin^2(2th13), dm2, chi^2_min)
  numub_chi2Surface_file.open(filePath+(chi2surfname).c_str());
  //numubarWD
  ofstream numubWS_chi2Surface_file;
  chi2surfname = "/data/numubWS_chi2_s2t-dm2_surface.txt";  //(sin^2(2th13), dm2, chi^2_min)
  numubWS_chi2Surface_file.open(filePath+(chi2surfname).c_str());
  
  //-- File to print minimized pull parameters
  string PullTfname;
  //-numu
  ofstream numu_minimPullT_file;
  PullTfname = "/data/numu_chi2_pullT_surface.txt";
  numu_minimPullT_file.open(filePath + (PullTfname).c_str());
  //-numubar
  ofstream numub_minimPullT_file;
  PullTfname = "/data/numub_chi2_pullT_surface.txt";
  numub_minimPullT_file.open(filePath + (PullTfname).c_str());
  //-numubar
  ofstream numubWS_minimPullT_file;
  PullTfname = "/data/numubWS_chi2_pullT_surface.txt";
  numubWS_minimPullT_file.open(filePath + (PullTfname).c_str());

  //--------------------------------------------------------------
  //---- LOOP over the oscillated spectra for numu and mnumubar 
  //--------------------------------------------------------------
  
  //-- Files with oscillated spectra in grid
  string numu_grid_spectra = "/data/minosNuMu_gridOscSpectra.txt";
  ifstream numu_grid_file;
  numu_grid_file.open(filePath + (numu_grid_spectra).c_str());

  string numub_grid_spectra = "/data/minosNuMuB_gridOscSpectra.txt";
  ifstream numub_grid_file;
  numub_grid_file.open(filePath + (numub_grid_spectra).c_str());

  string numubWS_grid_spectra = "/data/minosNuMuBWS_gridOscSpectra.txt";
  ifstream numubWS_grid_file;
  numubWS_grid_file.open(filePath + (numubWS_grid_spectra).c_str());

  cout << "Reading numu & numubar files - Loop in progress..." << endl;

  int first = 1;
  std::cout << "Is the numu file open?       " << numu_grid_file.is_open()    << "\n" << std::endl;  
  std::cout << "Is the numubar file open?    " << numub_grid_file.is_open()   << "\n" << std::endl;  
  std::cout << "Is the numubarWS file open?  " << numubWS_grid_file.is_open() << "\n" << std::endl;  

  while (
         numu_grid_file >> s2th_23 >> dm2_32 >> spcNumu[0] >> spcNumu[1] >> spcNumu[2] >> spcNumu[3] >> spcNumu[4] >> spcNumu[5] >> spcNumu[6] >> spcNumu[7] >> spcNumu[8] >> spcNumu[9] >> spcNumu[10] >> spcNumu[11] >> spcNumu[12] >> spcNumu[13] >> spcNumu[14] >> spcNumu[15] >> spcNumu[16] >> spcNumu[17] >> spcNumu[18] >> spcNumu[19] >> spcNumu[20] >> spcNumu[21] >> spcNumu[22] >> NoscTot 
         && 
         numub_grid_file >> s2th_23 >> dm2_32 >> spcNumuB[0] >> spcNumuB[1] >> spcNumuB[2] >> spcNumuB[3] >> spcNumuB[4] >> spcNumuB[5] >> spcNumuB[6] >> spcNumuB[7] >> spcNumuB[8] >> spcNumuB[9] >> spcNumuB[10] >> spcNumuB[11] >> NoscTotB 
         && 
         numubWS_grid_file >> s2th_23 >> dm2_32 >> spcNumuBWS[0] >> spcNumuBWS[1] >> spcNumuBWS[2] >> spcNumuBWS[3] >> spcNumuBWS[4] >> spcNumuBWS[5] >> spcNumuBWS[6] >> spcNumuBWS[7] >> spcNumuBWS[8] >> spcNumuBWS[9] >> spcNumuBWS[10] >> NoscTotBWS  )
    {//file loop
      
      if(first <= 1)
        {
          // numu
	  for(int ibin = 1 ; ibin <= 23 ; ibin++)
            {
	      spcNumuNoOsc[ibin-1] = spcNumu[ibin-1];
	      double bincont = spcNumu[ibin-1]*(NuMu_Events[0]/NoscTot);
	      numu_noosc_spect_hist->SetBinContent(ibin,bincont);
            }
          // numubar
	  for(int ibin = 1 ; ibin <= 12 ; ibin++)
            {
	      spcNumuBNoOsc[ibin-1] = spcNumuB[ibin-1];
	      double bincont = spcNumuB[ibin-1]*(NuMuB_Events[0]/NoscTotB);
	      numub_noosc_spect_hist->SetBinContent(ibin,bincont);
            }
          // numubarWS
	  for(int ibin = 1 ; ibin <= 11 ; ibin++)
            {
	      spcNumuBWSNoOsc[ibin-1] = spcNumuBWS[ibin-1];
	      double bincont = spcNumuBWS[ibin-1]*(NuMuBWS_Events[0]/NoscTotBWS);
	      numubWS_noosc_spect_hist->SetBinContent(ibin,bincont);
            }

	  cout << "first = " << first << endl;

	  cout << "Number of events      (numu): " << numu_noosc_spect_hist->Integral()    << endl;
	  noNoscTot = NoscTot;
	  cout << "Number of events   (numubar): " << numub_noosc_spect_hist->Integral()   << endl;
	  noNoscTotB = NoscTotB;
	  cout << "Number of events (numubarWS): " << numubWS_noosc_spect_hist->Integral() << endl;
	  noNoscTotBWS = NoscTotBWS;

	  first++;
        }//if first loop END
      
      const int N_params = 3; //-- Number of parameter of the chi² function --//

      //-- Steps
      double stp = 1.0e-4;
      double step[N_params] = {stp,stp,stp};
      //-- Initial parameter values
      double st = 0.0;
      double start[N_params] = {st,st,st};
      double lim;
      const double *xs;
      double chi2Min;

      //------------------------
      // numu section
      //------------------------
      
      SELECTOR = 0;
      ROOT::Math::Functor f1(&chi2,N_params); //-- Setting the function to be minimized by using Minuit --//

     //-- Set function and errors
      min1->SetFunction(f1);
      min1->SetErrorDef(2.3);

      //-- Set variable limits
      lim = 1.0e-1;
      min1->SetLimitedVariable(0,  "eps",   start[0],  step[0],  -lim, lim);
      min1->SetLimitedVariable(1,  "epsNC", start[1],  step[1],  -lim, lim);
      min1->SetLimitedVariable(2,  "ensc",  start[2],  step[2],  -lim, lim);
      //min1->SetFixedVariable(0,  "eps",   start[0]);
      //min1->SetFixedVariable(1,  "epsNC", start[1]);
      //min1->SetFixedVariable(2,  "ensc",  start[2]);     

      //-- Calling Minuit minimization
      min1->Minimize();
      xs = min1->X();
      chi2Min = min1->MinValue();
      numu_chi2Surface_file << s2th_23 << "\t" << dm2_32 << "\t" << chi2Min << endl;
      
      if (chi2Min < 0.0) {
	cout << "Critical error: chi2Min is negative!  " << chi2Min << endl;
	break;
      }
      numu_minimPullT_file  << s2th_23 << "\t" << dm2_32 << "\t" << xs[0] << "\t" << xs[1] << "\t" << xs[2] << "\t" << chi2Min << endl;

      if (dm2_32 == hi_dm2)
        {
	  numu_chi2Surface_file << endl;
	  cout << s2th_23 << "\t" << dm2_32 << "\t" << chi2Min << endl;
        }

      //------------------------
      // numubar section
      //------------------------

      SELECTOR = 1;
      ROOT::Math::Functor f2(&chi2,N_params); //-- Setting the function to be minimized by using Minuit --//
     
      //-- set functions and errors
      min2->SetFunction(f2);
      min2->SetErrorDef(2.3);
      
      //-- Set variable limits
      lim = 1.0e-1;
      min2->SetLimitedVariable(0,  "eps",   start[0],  step[0],  -lim, lim);
      min2->SetLimitedVariable(1,  "epsNC", start[1],  step[1],  -lim, lim);
      min2->SetLimitedVariable(2,  "ensc",  start[2],  step[2],  -lim, lim);
      //min2->SetFixedVariable(0,  "eps",   start[0]);
      //min2->SetFixedVariable(1,  "epsNC", start[1]);
      //min2->SetFixedVariable(2,  "ensc",  start[2]);
      
      //-- Calling Minuit minimization
      min2->Minimize();
      xs = min2->X();
      
      chi2Min = min2->MinValue();
      numub_chi2Surface_file << s2th_23 << "\t" << dm2_32 << "\t" << chi2Min << endl;
      
      if (chi2Min < 0.0) {
	cout << "Critical error: chi2Min is negative!  " << chi2Min << endl;
	break;
      }
      
      numub_minimPullT_file  << s2th_23 << "\t" << dm2_32 << "\t" << xs[0] << "\t" << xs[1] << "\t" << xs[2] << "\t" << chi2Min << endl;

      if (dm2_32 == hi_dm2) {
	  numub_chi2Surface_file << endl;
	  cout << s2th_23 << "\t" << dm2_32 << "\t" << chi2Min << endl;
      }

      //------------------------
      // numubarWS section
      //------------------------

      SELECTOR = 2;
      ROOT::Math::Functor f3(&chi2,N_params); //-- Setting the function to be minimized by using Minuit --//
     
      //-- set functions and errors
      min3->SetFunction(f3);
      min3->SetErrorDef(2.3);
      
      //-- Set variable limits
      lim = 1.0e-1;
      min3->SetLimitedVariable(0,  "eps",   start[0],  step[0],  -lim, lim);
      min3->SetLimitedVariable(1,  "epsNC", start[1],  step[1],  -lim, lim);
      min3->SetLimitedVariable(2,  "ensc",  start[2],  step[2],  -lim, lim);
      //min3->SetFixedVariable(0,  "eps",   start[0]);
      //min3->SetFixedVariable(1,  "epsNC", start[1]);
      //min3->SetFixedVariable(2,  "ensc",  start[2]);
      
      //-- Calling Minuit minimization
      min3->Minimize();
      xs = min3->X();
      
      chi2Min = min3->MinValue();
      numubWS_chi2Surface_file << s2th_23 << "\t" << dm2_32 << "\t" << chi2Min << endl;
      
      if (chi2Min < 0.0) {
	cout << "Critical error: chi2Min is negative!  " << chi2Min << endl;
	break;
      }
      
      numubWS_minimPullT_file  << s2th_23 << "\t" << dm2_32 << "\t" << xs[0] << "\t" << xs[1] << "\t" << xs[2] << "\t" << chi2Min << endl;

      if (dm2_32 == hi_dm2) {
	numubWS_chi2Surface_file << endl;
	cout << s2th_23 << "\t" << dm2_32 << "\t" << chi2Min << endl;
      }
      

    }//file loop END


  std::cout << "Successful run!!" << endl;
  numu_minimPullT_file.close();
  numu_chi2Surface_file.close();
  numub_minimPullT_file.close();
  numub_chi2Surface_file.close();
  numubWS_minimPullT_file.close();
  numubWS_chi2Surface_file.close();
  
  //------------------------------------------------------------------------------
  //------------------------------------------------------------------------------
  return 0;
}

