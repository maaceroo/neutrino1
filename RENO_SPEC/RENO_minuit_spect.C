//-----------------------------------------------------------------------------------//
//--  RENO_minuit.C - By M.A. Acero O., A.A. Aguilar-A. D.J. Polo T. - 2018-23-09  --//
//-----------------------------------------------------------------------------------//
//-------- Using 'NumericalMinimization.C' macro example to minimize a chi2 ---------//
//-------- function. It uses a Minimizer class in ROOT.                   -------- --//
//-------- input : minimizer name + algorithm name                         ----------//
//-------- randomSeed: = <0 : fixed value:                                 ----------//
//--------                0 random with seed 0;                            ----------//
//--------               >0 random with given seed.                       -------- --//
//--  Analysis of the RENO data from  S. H. Seo An et al., 1610.04326 (Aug 2017)   --//
//-----------------------------------------------------------------------------------//
//-------- This macro can be executed under ROOT typing                    ----------//
//-------- "root[0] .x RENO_minuit2.C"                                     -----------//
//-------- In the simplest case (only minimization) the script will print  ----------//
//-------- the message "Minimum: f(x[i]): chi^2(min)"                      ----------//
//-------- where x[i] are the values of the parameters which minimize the  ----------//
//-------- chi^2-function and chi^2(min) is the minimum value of the chi^2 ----------//
//-------- function.                                                       ----------//
//-------- In the complete case, the results are printed in a file, with   ----------//
//-------- the following information:                                      ----------//
//--------    sin^2(2th_13)  a   chi^2(min)                                ----------//
//-------- In this case, chi^2(min) is the minimum chi^2 for the pull      ----------//
//-------- terms.                                                          ----------//
//-----------------------------------------------------------------------------------//
#include "constants.h"
#include <math.h>
//---*****************************************************************************---//
//------------------------ CONSTANTS ------------------------------------------------//
//---*****************************************************************************---//
// histogram binning for the simulated data
//#define  NB 27
//#define  lo 1.2
//#define  hi 8.4
//Number of Antineutrino Detectors
#define  nAD 2
//Number of Nuclear Reactors
#define nNR 6
//Fixed neutrino oscillations parameters
#define dm2_21 7.53e-5 //eV^2,                  // 1610.0432v5 (2017)
//#define dm2_ee 2.62e-3 //eV^2,                  // 1610.0432v5 (2017)
#define s22th_12 0.846
//For the sin^2(2th_13) loop
//#define N_s2t  100                            //number of points in the grid
//#define N_dm2  100	                           //number of points in the grid
//#define hi_dm2 3.5e-3                       //dm2 max
double ran = N_dm2 * N_s2t;
//---*****************************************************---//

//---*****************************************************---//
//-- Global variables and quantities ------------------------//
//---*****************************************************---//
//(IBD candidates)/(DAQ live time -days-) from PRL 1610.0432v5 (2017)
//double IBDrate_data[nAD][2] = { {634.20,1.18},{64.38,0.36} };
double IBDrate_data[nAD][2] = { {616.67,1.14},{61.24,0.42} };
//IBD rate (per day), total background and efficiencies (1610.0432v5 (2017))
double totalBgd[nAD][2] = { {17.54,0.83},{3.14,0.23} };
double emuem[nAD] ={0.7644,0.7644};
double daqTime[nAD] = {458.49,489.93};
//---*****************************************************---//
// Information obtained by executing the script "RENO_osc_rate.C"
//IBD rate per day w/o oscillations
double noOsc_IBDrate_perday[nAD] = { 622.68,  65.39};
//<sin^2(1.267 dm2_21 L/E)> for each AD
//double avgSinDelta21[nAD] = { 0.000110395, 0.001236837};
//<sin^2(1.267 dm2_ee L/E)> for each AD
//double avgSinDeltaee[nAD] = { 0.109914046, 0.717994798};
//---*****************************************************---//
double s2th_13; //oscillation parameter to be fitted
double dm2_ee;   //oscillation parameter to be fitted
double a; //absolute normalization factor to be fitted
double e; //Escale Pull term
double epsilon; //Efficiency Pull term
double wrd_array[nAD][nNR];
//double fr[nNR];
//double eps_d[nAD];
double b_d[nAD];
double spc[nAD][NB];
double NoscTot[nAD];
double noNoscTot[nAD];
double ratio[NB];
int iAD;
int iNR;
//---*****************************************************---//
//const int nd = nAD;
TH1F *ratio_histo[nAD];
TH1F *data_spect_histo[nAD];
//TH1F *BkGd_spect_histo[nAD];
//---*****************************************************---//
TH1F *nosc_spect_hist[nAD];
TH1F *nosc_spect_hist_1[nAD];
TH1F *nosc_spect_hist_bf[nAD];
//---*****************************************************---//
//-- Chi-square function to be minimized --------------------//
//-- It has 3 pull parameters, 2 oscillation parameters and--//
//-- a normalization factor.                               --//
//---*****************************************************---//
double chi2(const double *xx)
{
  //Pull parameters considered in the chi² function
  epsilon = xx[0];
  e       = xx[1];
  b_d[0]  = xx[2];
  b_d[1]  = xx[3];
  //fr[0] = xx[4];
  //fr[1] = xx[5];
  //fr[2] = xx[6];
  //fr[3] = xx[7];
  //fr[4] = xx[8];
  //fr[5] = xx[9];
  
  
  //---*****************************************************---//
  
  int iAD;
  int iNR;
  double SurvPavg1        = 0.0;
  double SurvPavg2        = 0.0;
  double Nexp1            = 0.0;
  double Nexp2            = 0.0;
  double sqr_chi          = 0.0;
  double sqrerror         = 0.0;
  double Nobs1            = 0.0;
  double Nobs             = 0.0;
  double Nobs2            = 0.0;
  double Bd               = 0.0;
  double OFN              = 0.0;
  double TFN              = 0.0;
  double sB               = 0.087;
  double seps             = 0.002;
  double sfr_r            = 0.009;
  double sesc             = 0.0015;
  double wrd              = 0.0;
  
  SurvPavg1 = NoscTot[0]/noNoscTot[0];
  SurvPavg2 = NoscTot[1]/noNoscTot[1];
  
  for (int iBIN = 0 ; iBIN < NB ; iBIN++)
    {
      //  int index = iAD*NB + iBIN;
      //-- Measured IDB events of the dth Antineutrino Detector (background is substracted)
      Nobs1 = (data_spect_histo[0]->GetBinContent(iBIN+1))*IBDrate_data[0][0]*0.7644*daqTime[0];
      Nobs2 = (data_spect_histo[1]->GetBinContent(iBIN+1))*IBDrate_data[1][0]*0.7644*daqTime[1];
      
      sqrerror = ( Nobs2/(pow(Nobs1,2)) ) + ( (pow(Nobs2,2))/(pow(Nobs1,3)) );

      // Number expected Events from the simulation
      Nexp1 = spc[0][iBIN]*(SurvPavg1*noOsc_IBDrate_perday[0]/NoscTot[0])*0.7644*daqTime[0];
      Nexp2 = spc[1][iBIN]*(SurvPavg2*noOsc_IBDrate_perday[1]/NoscTot[1])*0.7644*daqTime[1];

      // Compute of the ratios of Data and expect spectra
      
      OFN = Nobs2/Nobs1;
      TFN = ( Nexp2 + b_d[1] )/( Nexp1 + b_d[0]);
      
      //cout << " OFN = " << OFN << " Nobs = " << Nobs <<endl;

      // Chi^2 funtion
      sqr_chi += pow( OFN - (TFN*(1.0 + epsilon + e ) ) ,2 )/sqrerror;
      
    }
  
  for (iAD = 0 ; iAD < nAD ; iAD++)
    {
      //-- Background error of the dth Antineutrino Detector
      sB = totalBgd[iAD][1]*0.7644*daqTime[iAD];
      sqr_chi +=  pow(b_d[iAD]/sB,2) ;
    }
  
  sqr_chi += pow(epsilon/seps,2) + pow(e/sesc,2) ;
   
  return sqr_chi;
  
}
//break;
//---*****************************************************---//

int RENO_minuit_spect(const char * minName = "Minuit",
			   const char *algoName = "" ,
			   int randomSeed = 10)
//int randomSeed = -1)
{
  cout << "Let's begin..." << endl;
  /* 
     TFile *wrd_File = new TFile("files_root/ldist_RENO_2x6.root","READ");
     TH1F *wrd_histo = ((TH1F*)(wrd_File->Get("histo_ldist_RENO_2x6")));;
     double mm;
     for (int blid = 0 ; blid < nAD*nNR ; blid++)
     {
     int id = (blid/nNR);
     int ir = (blid - id*nNR);
     
     wrd_array[id][ir] = wrd_histo->GetBinContent(blid+1);
     
     //cout << blid << "  " << (1/wrd_array[id][ir])   << endl;
     //cout << blid << "  " << wrd_array[id][ir]  << endl;
     }
     
  */
  
  //-------------------
  // Energy Histograms
  //-------------------
  TFile *fenergy = new TFile("files_root/RENOplots.root","read");
  //The histogram of near and far data spectra
  for(int n = 0 ; n < nAD ; n++){
    
    data_spect_histo[n] = (TH1F*) fenergy->Get(Form("data_spect_histo_%d",n));
    double dfactor = 1.0/data_spect_histo[n]->Integral();
    data_spect_histo[n]->Scale(dfactor);
    
    ratio_histo[n] = (TH1F*) fenergy->Get(Form("ratio_histo_%d",n));
    double rfactor = 1.0/ratio_histo[n]->Integral();
    ratio_histo[n]->Scale(rfactor);
    
    
  }

  // define number of bins //
  //const int  NB = 27; 
  //const double lo = 1.2;  
  //const double hi = 8.4;
  //const int nd = 2;
  double xbins[NB+1];
  //xbins[0] = 1.2;
  //cout << xbins[0] <<endl;
  double delta_bins2 = (6.0 - 1.2)/24; // 0.2 MeV/bin
  
  for (int i = 0 ; i < (NB-2) ; i++)
    {
      xbins[i] = 1.2 + delta_bins2*i;
    }
  xbins[25] = xbins[24] + 0.4;
  xbins[26] = 8.4 - 1.4;
  xbins[27] = 8.4;
  
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  for(int iAD = 0 ; iAD < nAD ; iAD++)
    {
      nosc_spect_hist[iAD] = new TH1F(Form("nosc_spect_hist_%d",iAD),"",NB,xbins);
      nosc_spect_hist_1[iAD] = new TH1F(Form("nosc_spect_hist_1_%d",iAD),"",NB,xbins);
      nosc_spect_hist_bf[iAD] = new TH1F(Form("nosc_spect_hist_bf_%d",iAD),"",NB,xbins);
    }
  
  cout << "Minimization settings..." << endl;
  ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer(minName, algoName);
  
  //-- Set tolerance , etc...
  min->SetMaxFunctionCalls(1000000); // for Minuit/Minuit2
  min->SetMaxIterations(10000);      // for GSL
  min->SetTolerance(0.001);
  min->SetPrintLevel(-1);
  
  //-- File to print oscillation parameters and chi2 values
  ofstream chi2Surface_file;
  string s2t_dm2 = "files/chi2_s2t-dm2_surface_spect.txt";  //(sin^2(2th13), a, chi^2_min)
  ofstream minimPullT_file;
  string pterms = "files/chi2_pullTerms_spect.txt";
  //string s2t_eps = "chi2_s2t_curve.txt"; //(sin^2(2th13), chi^2_min)
  chi2Surface_file.open((s2t_dm2).c_str());
  minimPullT_file.open((pterms).c_str());
  
  
  ifstream file("./files/RENO_gridOscSpectra_test.txt"); //50x50 parameter space-grid
  //ifstream file("files/RENO_gridOscSpectra.txt");
  cout << "Reading file - Loop in progress..." << endl;
  int iad = 0;
  //int inr = 0;
  int first2 = 1;
  
  while (file >> iAD >> s2th_13 >> dm2_ee  >> spc[iad][0] >> spc[iad][1] >> spc[iad][2] >> spc[iad][3] >> spc[iad][4] >> spc[iad][5] >> spc[iad][6] >> spc[iad][7] >> spc[iad][8] >> spc[iad][9] >> spc[iad][10] >> spc[iad][11] >> spc[iad][12] >> spc[iad][13] >> spc[iad][14] >> spc[iad][15] >> spc[iad][16] >> spc[iad][17] >> spc[iad][18] >> spc[iad][19] >> spc[iad][20] >> spc[iad][21] >> spc[iad][22] >> spc[iad][23] >> spc[iad][24] >> spc[iad][25] >> spc[iad][26] >> NoscTot[iad])
    //    break;
    {//file loop
     // cout << " a = " << spc[0][26] << endl; 
      if(first2 <= 2)
	{
	  for(int ibin = 1 ; ibin <= NB ; ibin++)
	    {
	      double bincont = spc[iad][ibin-1]*(noOsc_IBDrate_perday[iAD-1]/NoscTot[iad])*emuem[iAD-1]*daqTime[iAD-1];
	      double bincontent = bincont/0.2;
	      nosc_spect_hist[iAD-1]->SetBinContent(ibin,bincont);
	      nosc_spect_hist_1[iAD-1]->SetBinContent(ibin,bincontent);
	    }
	  
	  cout << " iAD = " << iAD  /*<< " iNR = " << iNR */ << "  first2 = " << first2 << endl;
	  //nosc_spect_hist[iAD-1]->Print("all"); // ---- 2017-08-07
	  cout << "Number of events: " << nosc_spect_hist[iAD-1]->Integral() << endl; // ---- 2017-08-07
	  //cout << "NoscTot: " << NoscTot[iad][inr] << endl; // ---- 2017-08-07
	  //cout << "NoscTot         : " << NoscTot[iad] << endl; // ---- 2017-08-07
	  noNoscTot[iad] = NoscTot[iad];
	  
	  first2++;
	}//if first2 loop END
      
      iad++;
      
      if(iad == 2)
	{
	  //At this point, we have read the two AD spectra (two lines) for one point (s2th_13,dm2_31) in the grid
	  //if iad BEGIN
	  //cout << endl;
	  iad = 0;
	  
	  //cout << "Start Minimization..." << endl;
	  const int N_params = 4; //-- Number of parameter of the chi² function --//
	  ROOT::Math::Functor f(&chi2,N_params); //-- Setting the function to be minimized by using Minuit --//
	  //-- Steps
	  double stp = 1.0e-3;
	  double step[N_params] = {stp,stp,
				   stp,stp};

	  //-- Initial parameter values
	  double start[N_params] = {0.0,0.0,
				    0.0,0.0};

	  //-- Calling Minuit function setting
	  min->SetFunction(f);
	  
	  //-- Setting variables
	  double lim = 1.0e-3;
	  
	  min->SetLimitedVariable(0,  "epsilon", start[0],  step[0],  -lim, lim);
	  min->SetLimitedVariable(1,  "e",       start[1],  step[1],  -lim, lim);
	  min->SetLimitedVariable(2,  "b_0",     start[2],  step[2],  -lim, lim);
	  min->SetLimitedVariable(3,  "b_1",     start[3],  step[3],  -lim, lim);
	  min->SetErrorDef(2.3);
	  
	/*
	  min->SetFixedVariable(0,  "epsilon", start[0]);
	  min->SetFixedVariable(1,  "e",       start[1]);
	  min->SetFixedVariable(2,  "b_0",     start[2]);
	  min->SetFixedVariable(3,  "b_1",     start[3]);
	*/
	  //-- Calling Minuit minimization
	  min->Minimize();
	  
	  const double *xs = min->X();
	  
	  double chi2Min = min->MinValue();
	  chi2Surface_file << s2th_13 << "\t" << dm2_ee << "\t" << chi2Min << endl;
	  
	  if (chi2Min < 0.0) {
	    cout << "Critical error: chi2Min is negative!  " << chi2Min << endl;
	    break;
	  }
	  
	  //-- Uncomment if you want to print the pull parameters (also Line 205)
	  minimPullT_file  << s2th_13 << "\t" << dm2_ee  << "\t" << xs[0] << "\t" << xs[1] << "\t" << xs[2] << "\t" << xs[3] << "\t" << chi2Min << endl;
	  //}
	  
	  if (dm2_ee == hi_dm2)
	    {
	      chi2Surface_file << endl;
	      cout << s2th_13 << "\t" << dm2_ee << "\t" << min->MinValue() << endl;
	    }
	}  
    }
  //

  std::cout << "Succesful run!!" << endl;
  
  // Drawing section
  TCanvas *c1 = new TCanvas("c1");
  nosc_spect_hist_1[0]->Draw();
  //nosc_spect_hist_1[1]->Draw("same");
  c1->Print("Plots/nosc_near.pdf");
  
  TCanvas *c2 = new TCanvas("c2");
  //nosc_spect_hist_1[0]->Draw();
  nosc_spect_hist_1[1]->Draw();
  c2->Print("Plots/nosc_far.pdf");

  //////////////////////////////////////////////// Chi2 minimun value ///////////////////////////////////////////////////////////////////
    
  ofstream chi2min;
  string chiminima = "files/chi2_minimun_spect.txt";
  chi2min.open((chiminima).c_str());
  int n;
  const int rows = 3;
  const int columns = ran;
  ifstream matriz("files/chi2_s2t-dm2_surface_spect.txt"); 
  double ** matr;  
  double minimo[rows];
  matr = new double*[rows];
  for(int k = 0 ; k < rows ; k++){
    matr[k] = new double[columns];
  }		
  
  for(int l = 0 ; l < columns ; l++){		
    for(int j = 0 ; j < rows ; j++){
      matriz >> matr[j][l];
    }
  }
  
  for(int i = 0 ; i < rows ; i++){
    minimo [i] = matr[i][0];
    for(int ll = 0 ; ll < columns ; ll++){
      if(matr[i][ll] < minimo[i]){
	minimo[i] = matr[i][ll];
	n =ll;
      }
    }
  }  
  cout << "sin2t_13 = " << matr[0][n]   << " dm2_ee = " << matr[1][n]  << " chi2 = " << matr[2][n]  << " n = " << n  << endl;
  chi2min  << matr[0][n]   <<  "\t" << matr[1][n]  << "\t" << matr[2][n]  << endl;
  chi2min <<endl;
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  int m;
  const int rows3 = nAD*ran + nAD;
  const int columns3 = 31;
  ifstream matriz3("files/RENO_gridOscSpectra_test.txt"); 
  double ** matr3;  
  double bfit[rows3];
  matr3 = new double*[rows3];
  for(int k = 0 ; k < rows3 ; k++){
    matr3[k] = new double[columns3];
  }		
  
  for(int j = 0 ; j < rows3 ; j++){
    for(int l = 0 ; l < columns3 ; l++){		
      matriz3 >> matr3[j][l];
    }
  }
  
  for(int i = 0 ; i < rows3 ; i++){
    if(matr3[i][2] == matr[1][n] && matr3[i][1] == matr[0][n] ){
      m = i;
    }
  }
    
  double best_fit[nAD][NB];
  double NoscTot_bf[nAD];
  
  NoscTot_bf[0] = matr3[m-1][30];
  NoscTot_bf[1] = matr3[m][30];
  
  for(int i=0 ; i < NB ; i++ ){
    
    best_fit[0][i] = matr3[m-1][i+3];
    best_fit[1][i] = matr3[m][i+3];
    
    //  cout << "near = " << best_fit[0][i] << " far = " << best_fit[1][i] << endl; 
    
  }
  
  for(int iAD = 0 ; iAD < nAD ; iAD++)
    {
      for(int ibin = 0 ; ibin < NB ; ibin++)
	{
	  double bincont = best_fit[iAD][ibin]*(noOsc_IBDrate_perday[iAD]/NoscTot_bf[iAD])*emuem[iAD]*daqTime[iAD];
	  double bincontent = bincont/xbins[ibin];
	  nosc_spect_hist_bf[iAD]->SetBinContent(ibin,bincontent);
	  //nosc_spect_hist_bf[iAD]->SetBinContent(ibin,bincontent);
	}
    }
  
  TCanvas *c3 = new TCanvas("c3");
  nosc_spect_hist_bf[0]->Draw();
  //nosc_spect_hist_1[1]->Draw("same");
  c3->Print("Plots/nosc_near_bestfit.pdf");
  
  TCanvas *c4 = new TCanvas("c4");
  //frame_spectrafd->Draw();
  nosc_spect_hist_bf[1]->Draw();
  //nosc_spect_hist_1[1]->Draw("same");
  c4->Print("Plots/nosc_far_bestfit.pdf");
  
  return 0;
}
