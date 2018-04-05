//-----------------------------------------------------------------------------------//
//--  RENO_minuit.C - By M.A. Acero O., A.A. Aguilar-A. D.J. Polo T. - 2017-02-02  --//
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

//---*****************************************************************************---//
//------------------------ CONSTANTS ------------------------------------------------//
//---*****************************************************************************---//
// histogram binning for the simulated data
#define  NB 27
#define  lo 1.2
#define  hi 8.4
//Number of Antineutrino Detectors
#define  nAD 2
//Number of Nuclear Reactors
#define nNR 6
//Fixed neutrino oscillations parameters
#define dm2_21 7.53e-5 //eV^2,                  // 1610.0432v5 (2017)
#define dm2_ee 2.62e-3 //eV^2,                  // 1610.0432v5 (2017)
#define s22th_12 0.846
//For the sin^2(2th_13) loop
#define N_s2t  100                            //number of points in the grid
//#define lo_s2t 0.01                             //sin^2(2th_13) min
#define lo_s2t 0.0                            //sin^2(2th_13) min
#define hi_s2t 0.3                            //sin^2(2th_13) max
//For the a loop
#define N_eps  100 	                           //number of points in the grid
#define lo_eps -1.0e-2                       //a min
#define hi_eps +1.0e-2                       //a max
double ran = N_eps * N_s2t;
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
//double noOsc_IBDrate_perday[nAD][nNR] ={ { 630.53, 623.43, 619.85, 620.52, 625.28, 633.21}, {65.66,  65.43,  65.28,  65.21,  65.34,  65.49}};
double noOsc_IBDrate_perday[nAD][nNR] ={ {  43.18,  93.29, 206.74, 170.21,  73.36,  35.92}, {   9.34,  10.64, 11.72,  11.95,  11.43, 10.31}};
//double noOsc_IBDrate_perday[nAD] = { 622.68,  65.39};
//<sin^2(1.267 dm2_21 L/E)> for each AD
//double avgSinDelta21[nAD] = { 0.000110395, 0.001236837};
double avgSinDelta21[nAD][nNR] = { {0.000267,0.000121,0.000055,0.000067,0.000157,0.000327}, {0.001431,0.001251,0.001163,0.001127,0.001189,0.001314}};
//<sin^2(1.267 dm2_ee L/E)> for each AD
//double avgSinDeltaee[nAD] = { 0.109914046, 0.717994798};
double avgSinDeltaee[nAD][nNR] = { {0.250234,0.123595,0.058540,0.070624,0.156896,0.297120}, {0.760081,0.724878,0.701159,0.689812,0.710250,0.733815}};
//---*****************************************************---//
double s22th_13; //oscillation parameter to be fitted
double a; //absolute normalization factor to be fitted
double wrd_array[nAD][nNR];
double fr[nNR];
double eps_d[nAD];
double b_d[nAD];

//---*****************************************************---//
//-- Chi-square function to be minimized --------------------//
//-- It has 10 pull parameters, 1 oscillation parameter and--//
//-- a normalization factor. The last two are to be fitted. -//
//---*****************************************************---//
double chi2(const double *xx)
{
  //Pull parameters considered in the chi² function
  eps_d[0] = xx[0];
  eps_d[1] = xx[1];
  b_d[0] = xx[2];
  b_d[1] = xx[3];
  fr[0] = xx[4];
  fr[1] = xx[5];
  fr[2] = xx[6];
  fr[3] = xx[7];
  fr[4] = xx[8];
  fr[5] = xx[9];
  
//---*****************************************************---//
  int iAD;
  int iNR;
  double SurvP        = 0.0;
  double sqr_chi      = 0.0;
  double sqrerror     = 0.0;
  double Nobs         = 0.0;
  double Nexp         = 0.0;
  double Bd           = 0.0;
  double sB           = 0.0;
  double seps_d       = 0.002;
  double seps_a       = 0.009;
  double wrd = 0.0;
  
  for (iAD = 0 ; iAD < 2 ; iAD++)
    {
      //-- Measured IDB events of the dth Antineutrino Detector (background is substracted)
      Nobs = (IBDrate_data[iAD][0]/*-totalBgd[iAD][0]*/)*0.7644*daqTime[iAD];
      // Error 
      sqrerror = Nobs;
      for (iNR = 0 ; iNR < nNR ; iNR++)
	{
	  //-- Survival Probability equation. Terms depending on dM²_21 and dM²_ee(~dM²_31) are averaged
	  SurvP = 1.0 - (s22th_13*avgSinDeltaee[iAD][iNR]) - (0.25*pow((1.0 + sqrt(1.0 - s22th_13)),2)*s22th_12*avgSinDelta21[iAD][iNR]);
	  //-- Predicted IBD from neutrino oscillations of the dth Antineutrino Detector
	  // cout << "SurP = " << SurvP <<endl;
	  Nexp = (SurvP*noOsc_IBDrate_perday[iAD][iNR])*0.7644*daqTime[iAD];
	  //Nexp = (SurvP*noOsc_IBDrate_perday[iAD]*(area_bajo_curva[iAD]))*0.7644*daqTime[iAD];
	  //-- determined by baselines and reactor fluxes
	  wrd += Nexp*(1+fr[iNR])/*wrd_array[iAD][iNR]*(area_bajo_curva[iAD])*/;
	}	  
      //wrd = Nexp*(1+fr[iNR]);
      //cout << "wrd = " << wrd << "\t Nobs = " << Nobs <<endl;
      sqr_chi += pow( (Nobs + b_d[iAD] - wrd*(1.0 + a + eps_d[iAD])) ,2 )/sqrerror;
      wrd = 0.0;
    }
  //     cout << "chi2 = " << sqr_chi << endl;
  for (iAD = 0 ; iAD < 2 ; iAD++)
    {
      //-- Background error of the dth Antineutrino Detector
      sB = totalBgd[iAD][1]*0.7644*daqTime[iAD];
      sqr_chi += pow(eps_d[iAD]/seps_d,2) + pow(b_d[iAD]/sB,2);
    }
  
  for (iNR = 0 ; iNR < 6 ; iNR++)
    sqr_chi += pow(fr[iNR]/seps_a,2);
  
  return sqr_chi;
}

//---*****************************************************---//

int RENO_minuit2(const char * minName = "Minuit",
		 const char *algoName = "" ,
		 int randomSeed = +10)
//int randomSeed = -1)
{
  cout << "Let's begin..." << endl;
  
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
  
  
  //break;
  //    cout << "wrd array -> Done..." << endl;
  
  cout << "Minimization settings..." << endl;
  ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer(minName, algoName);
    
  //-- Set tolerance , etc...
  min->SetMaxFunctionCalls(1000000); // for Minuit/Minuit2
  min->SetMaxIterations(10000);      // for GSL
  min->SetTolerance(0.001);
  min->SetPrintLevel(-1);
  
  //-- File to print oscillation parameters and chi2 values
  ofstream chi2Surface_file;
  string s2t_eps = "files/chi2_s2t-eps_surface_RATE.txt";  //(sin^2(2th13), a, chi^2_min)
  ofstream minimPullT_file;
  string pterms = "files/chi2_pullTerms_RATE.txt";
  //string s2t_eps = "chi2_s2t_curve.txt"; //(sin^2(2th13), chi^2_min)
  chi2Surface_file.open((s2t_eps).c_str());
  minimPullT_file.open((pterms).c_str());
  
  int is2t;
  double DeltaLog_s2t = (hi_s2t - lo_s2t)/double(N_s2t - 1);  //linear
  
  //-- For the normalization factor loop
  int ieps;
  //double a; //Definied as a global parameter
  double DeltaLin_eps = (hi_eps - lo_eps)/double(N_eps - 1);
  
    cout << "Loop on s2th13 -> Begin..." << endl;
    for (is2t = 0 ; is2t < N_s2t ; is2t++)
      {
	//s2th_13 = pow(10,(log10(lo_s2t) + double(is2t)*DeltaLog_s2t)); //logarithmic
	s22th_13 = lo_s2t + double(is2t)*DeltaLog_s2t; //linear
	
	for (ieps = 0 ; ieps < N_eps ; ieps++)
	  {
	    a = lo_eps + double(ieps)*DeltaLin_eps;
	    
	    const int N_params = 10; //-- Number of parameter of the chi² function --//
	    ROOT::Math::Functor f(&chi2,N_params); //-- Setting the function to be minimized by using Minuit --//
	    //-- Steps
	    double stp = 1.0e-3;
	    //double stp1 = 1.0e1;
	    double step[N_params] = {stp,stp,
	    			     stp,stp,
	    			     stp,stp,stp,stp,stp,stp};
	    //double step[N_params-6] = {stp,stp,
	    //stp,stp};
	    //double step1[N_params-4] = {stp1,stp1,
	    //stp1,stp1,stp1,stp1};
				  
	    //-- Initial parameter values
	    double start[N_params] = {0.0,0.0,
				      0.0,0.0,
				      0.0,0.0,0.0,0.0,0.0,0.0};
	    //double s100 = 100.0;
	    //double start[N_params] = {0.0,0.0,
	    //0.0,0.0,
	    ////100.09,80.84763,69.72081,69.00108,77.50215,95.1999};
	    //s100,s100,s100,s100,s100,s100};
	    
	    //-- Calling Minuit function setting
	    min->SetFunction(f);
	    
	    //-- Setting variables
	    double lim = 1.0e-3;
	    //double lim2 = 1.2e2;
	    double lim2 = 10.0*lim;
	    min->SetLimitedVariable(0,  "b_1", start[0],  step[0],  -lim, lim);
	    min->SetLimitedVariable(1,  "b_2", start[1],  step[1],  -lim, lim);
	    min->SetLimitedVariable(2,  "e_3", start[2],  step[2],  -lim, lim);
	    min->SetLimitedVariable(3,  "e_4", start[3],  step[3],  -lim, lim);
	    min->SetLimitedVariable(4,  "fr_0", start[4],  step[4],  -lim2, lim2);
	    min->SetLimitedVariable(5,  "fr_1", start[5],  step[5],  -lim2, lim2);
	    min->SetLimitedVariable(6,  "fr_2", start[6],  step[6],  -lim2, lim2);
	    min->SetLimitedVariable(7,  "fr_3", start[7],  step[7],  -lim2, lim2);
	    min->SetLimitedVariable(8,  "fr_4", start[8],  step[8],  -lim2, lim2);
	    min->SetLimitedVariable(9,  "fr_5", start[9],  step[9],  -lim2, lim2);
	    //min->SetFixedVariable(0,  "b_1", start[0]);
	    //min->SetFixedVariable(1,  "b_2", start[1]);
	    //min->SetFixedVariable(2,  "e_3", start[2]);
	    //min->SetFixedVariable(3,  "e_4", start[3]);
	    //min->SetFixedVariable(4,  "fr_0", start[4]*0.0);
	    //min->SetFixedVariable(5,  "fr_1", start[5]*0.0);
	    //min->SetFixedVariable(6,  "fr_2", start[6]*0.0);
	    //min->SetFixedVariable(7,  "fr_3", start[7]*0.0);
	    //min->SetFixedVariable(8,  "fr_4", start[8]*0.0);
	    //min->SetFixedVariable(9,  "fr_5", start[9]*0.0);
	    min->SetErrorDef(2.3);
	    
	    //-- Calling Minuit minimization
	    min->Minimize();
		
	    const double *xs = min->X();
	    
	    chi2Surface_file << s22th_13 << "\t" << a << "\t" << min->MinValue() << endl;
	    
	    //-- Uncomment if you want to print the pull parameters (also Line 205)
	    minimPullT_file  << s22th_13 << "\t" << a  << "\t" << /*xs[0] << "\t" << xs[1] << "\t" << xs[2] << "\t" << xs[3] << "\t" <<*/ xs[4] << "\t" << xs[5] << "\t" << xs[6] << "\t" << xs[7] << "\t" << xs[8] << "\t" << xs[9] << "\t" << min->MinValue() << endl;
	      }
	chi2Surface_file << endl;
	//-- Uncomment if you want to print the pull parameters
	minimPullT_file  << endl;
	    
	if (is2t%10 == 0)
	  std::cout << "Succesful run for sin^2(th13) = " << s22th_13 << "!! \t" << min->MinValue() << endl;
      }
    
    std::cout << "Succesful run!!" << endl;
    
    //////////////////////////////////////////////// Chi2 minimun value ////////////////////////////////////////////////////////////////////
    
    ofstream chi2min;
    string chiminima = "files/chi2_minumun_RATE.txt";
    chi2min.open((chiminima).c_str());
    int n;
    const int rows = 3;
    const int columns = ran;
    ifstream matriz("files/chi2_s2t-eps_surface_RATE.txt"); 
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

    cout << "sin2t_13 = " << matr[0][n]   << " a = " << matr[1][n]  << " chi2 = " << matr[2][n]  << " n = " << n  << endl;
    chi2min  << matr[0][n]   <<  "\t" << matr[1][n]  << "\t" << matr[2][n]  << endl;
    chi2min <<endl;
    /////////////////////////////////////////////
    return 0;
}
