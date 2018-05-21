//---------------------------------------------------------------------//
//--  db_minuit.C - By M.A. Acero O. & A.A. Aguilar-A. - 2018-05-12  --//
//---------------------------------------------------------------------//
//-- Using 'NumericalMinimization.C' macro example to minimize a chi2--//
//-- function. It uses a Minimizer class in ROOT.                    --//
//-- input : minimizer name + algorithm name                         --//
//-- randomSeed: = <0 : fixed value:                                 --//
//--                0 random with seed 0;                            --//
//--               >0 random with given seed.                        --//
//-- Analysis of the Daya Bay data from                              --//
//-- F.P. An et al., PRD 95 072006 (2017)                            --//
//---------------------------------------------------------------------//
//-- This macro can be executed under ROOT typing                    --//
//-- "root[0] .x db_minuit.C"                                        --//
//-- In the simplest case (only minimization) the script will print  --//
//-- the message "Minimum: f(x[i]): chi^2(min)"                      --//
//-- where x[i] are the values of the parameters which minimize the  --//
//-- chi^2-function and chi^2(min) is the minimum value of the chi^2 --//
//-- function.                                                       --//
//-- In the complete case, the results are printed in a file, with   --//
//-- the following information:                                      --//
//--    sin^2(2th_13)  epsilon  chi^2(min)                           --//
//-- In this case, chi^2(min) is the minimum chi^2 for the pull      --//
//-- terms.                                                          --//
//---------------------------------------------------------------------//

//---***************************************************************---//
//---------- CONSTANTS ------------------------------------------------//
//---***************************************************************---//
#include "constants.h"
// histogram binning for the simulated data
//#define  NB 35
//#define  lo 0.7
//#define  hi 12.0
//Number of Antineutrino Detectors
//#define nAD 8
//Number of Nuclear Reactors
//#define nNR 6
//Fixed neutrino oscillations parameters
//#define dm2_21 7.53e-5 //eV^2,                  //2017 Review of Particle Physics (On-line, 2018.05.19).
//#define dm2_31 2.45e-3 //eV^2,                  //2017 Review of Particle Physics (On-line, 2018.05.19).
//#define s22th_12 0.861
//For the sin^2(2th_13) loop
//#define N_s2t  60                            //number of points in the grid
//#define lo_s2t 0.01                             //sin^2(2th_13) min
//#define lo_s2t 0.05                            //sin^2(2th_13) min
//#define hi_s2t 0.25                            //sin^2(2th_13) max
//For the epsilon loop
//#define N_eps  60                           //number of points in the grid
//#define lo_eps -1.0e-2                       //epsilon min
//#define hi_eps +1.0e-2                       //epsilon max
//---*****************************************************---//

//---*****************************************************---//
//-- Global variables and quantities ------------------------//
//---*****************************************************---//
//EH1(AD1, AD2),EH2(AD3, AD8),EH3(AD4, AD5, AD6, AD7)
//(IBD candidates)/(DAQ live time -days-) from PRD 95 072006 (2017)
double IBDrate_data[nAD][2] = { {534.93,23.13},{542.75,23.30},{509.00,22.56},{503.83,22.45},{ 72.71, 8.53},{ 72.94, 8.54},{ 72.33, 8.50},{ 72.88, 8.54} };
//IBD rate (per day), total background and efficiencies (PRD 95 072006 (2017))
double totalBgd[nAD][2] = { {11.96, 0.91},{11.96, 0.91},{ 8.79, 0.66},{ 8.69, 0.78},{ 1.59, 0.06},{ 1.50, 0.06},{ 1.51, 0.06},{ 1.28, 0.07} };
double emuem[nAD] ={0.8044,0.8013,0.8365,0.8363,0.9587,0.9585,0.9581,0.9588};
double daqTime[nAD] = {1117.178,1117.178,1114.337,924.933,1106.915,1106.915,1106.915,917.417};
//---*****************************************************---//
// Information obtained by executing the script "db_osc_rate.C"
//IBD rate per day w/o oscillations
double noOsc_IBDrate_perday[nAD] = { 663.03, 675.43, 610.16, 604.31,  79.73,  80.08,  79.46,  80.25};
//<sin^2(1.267 dm2_21 L/E)> for each AD
double avgSinDelta21[nAD] = { 0.000230234, 0.000226747, 0.000251988, 0.000255204, 0.001917360, 0.001915281, 0.001929247, 0.001938521};
//<sin^2(1.267 dm2_31 L/E)> for each AD
double avgSinDelta31[nAD] = { 0.179606632, 0.176403459, 0.204881620, 0.207666334, 0.762535106, 0.759829143, 0.764013835, 0.762621122};
//---*****************************************************---//
double s2th_13; //oscillation parameter to be fitted
double epsilon; //absolute normalization factor to be fitted
double wrd_array[nAD][nNR];
double alpha[nNR];
double eps_d[nAD];
double eta_d[nAD];

//---*****************************************************---//
//-- Chi-square function to be minimized --------------------//
//-- It has 18 pull parameters, 1 oscillation parameter and--//
//-- a normalization factor. The last two are to be fitted. -//
//---*****************************************************---//
double chi2(const double *xx)
{
    //Pull parameters considered in the chi² function
    eps_d[0] = xx[0];
    eps_d[1] = xx[1];
    eps_d[2] = xx[2];
    eps_d[3] = xx[3];
    eps_d[4] = xx[4];
    eps_d[5] = xx[5];
    eta_d[0] = xx[6];
    eta_d[1] = xx[7];
    eta_d[2] = xx[8];
    eta_d[3] = xx[9];
    eta_d[4] = xx[10];
    eta_d[5] = xx[11];
    alpha[0] = xx[12];
    alpha[1] = xx[13];
    alpha[2] = xx[14];
    alpha[3] = xx[15];
    alpha[4] = xx[16];
    alpha[5] = xx[17];

//---*****************************************************---//
    int iAD;
    int iNR;
    double SurvP        = 0.0;
    double sqr_chi      = 0.0;
    double sqrerror     = 0.0;
    double Md           = 0.0;
    double Td           = 0.0;
    double Bd           = 0.0;
    double sB           = 0.0;
    double seps_d       = 1.0*0.002;
    double seps_a       = 1.0*0.008;
    
    for (iAD = 0 ; iAD < 8 ; iAD++)
    {
	//-- Survival Probability equation. Terms depending on dM²_21 and dM²_32(~dM²_31) are averaged
        SurvP = 1.0 - 0.25*pow((1.0 + sqrt(1.0 - s2th_13)),2)*(s22th_12)*(avgSinDelta21[iAD]) - s2th_13*(avgSinDelta31[iAD]);

	//-- Predicted IBD from neutrino oscillations of the dth Antineutrino Detector
        Td = (SurvP * noOsc_IBDrate_perday[iAD])*emuem[iAD]*daqTime[iAD];
	//-- Measured IDB events of the dth Antineutrino Detector (background is substracted)
        Md = (IBDrate_data[iAD][0] - totalBgd[iAD][0]*emuem[iAD])*daqTime[iAD];
	//-- Background of the dth Antineutrino Detector
        Bd = totalBgd[iAD][0]*emuem[iAD]*daqTime[iAD];
        
        sqrerror = Md + Bd;
	
	//-- Fraction of IBD contribution of the rth reactor to the dth AD 
	//-- determined by baselines and reactor fluxes
	double wrd = 0.0;
        for (iNR = 0 ; iNR < nNR ; iNR++)
            wrd += wrd_array[iAD][iNR]*alpha[iNR];
        
        //-- Testing a fake normalization factor
        sqr_chi += pow( (Md - Td*(1.0 + epsilon + eps_d[iAD] + wrd) + eta_d[iAD]) ,2 )/sqrerror;
        
    }
    
    for (iAD = 0 ; iAD < 8 ; iAD++)
    {
        //-- Background error of the dth Antineutrino Detector
        sB = totalBgd[iAD][1]*emuem[iAD]*daqTime[iAD];
        sqr_chi += pow(eps_d[iAD]/seps_d,2) + pow(eta_d[iAD]/sB,2);
    }
    
    for (iNR = 0 ; iNR < 6 ; iNR++)
        sqr_chi += pow(alpha[iNR]/seps_a,2);
    
    return sqr_chi;
}
//---*****************************************************---//

int db_minuit(const char * minName = "Minuit",
              const char *algoName = "" ,
              int randomSeed = +10)
//int randomSeed = -1)
{
    cout << "Let's begin..." << endl;
    
    TFile *wrd_File = new TFile("files_data/daya-bay-ldist.root","READ");
    TH1F *wrd_histo = ((TH1F*)(wrd_File->Get("histo_ldist")));;
    
    for (int blid = 0 ; blid < nAD*nNR ; blid++)
    {
        int id = (blid/nNR);
        int ir = (blid - id*nNR);

        wrd_array[id][ir] = wrd_histo->GetBinContent(blid+1);
        //cout << blid << "  " << wrd_array[id][ir] << endl;
    }
    //break;
    cout << "wrd array -> Done..." << endl;
    
    cout << "Minimization settings..." << endl;
    ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer(minName, algoName);
    
    //-- Set tolerance , etc...
    min->SetMaxFunctionCalls(1000000); // for Minuit/Minuit2
    min->SetMaxIterations(10000);      // for GSL
    min->SetTolerance(0.001);
    min->SetPrintLevel(-1);
    
    //-- File to print oscillation parameters and chi2 values
    ofstream chi2Surface_file;
    string s2t_eps = "files_data/chi2_s2t-eps_surface_RATE.txt";  //(sin^2(2th13), epsilon, chi^2_min)
    //string s2t_eps = "files_data/chi2_s2t_curve.txt"; //(sin^2(2th13), chi^2_min)
    chi2Surface_file.open((s2t_eps).c_str());

    //-- Uncomment if you want to print the pull parameters (also Line 279)
    //-- File to print minimized pull parameters
/*
    ofstream minimPullT_file;
    string PullT = "files_data/chi2_pullT_surface.txt";
    minimPullT_file.open((PullT).c_str());
*/
    
    //-- For the oscillation parameter loop
    int is2t;
    //double s2th_13; //Definied as a global parameter
    //double DeltaLog_s2t = (log10(hi_s2t)-log10(lo_s2t))/double(N_s2t-1); //logarithmic
    double DeltaLog_s2t = (hi_s2t - lo_s2t)/double(N_s2t - 1);  //linear

    //-- For the normalization factor loop
    int ieps;
    //double epsilon; //Definied as a global parameter
    double DeltaLin_eps = (hi_eps - lo_eps)/double(N_eps - 1);

    cout << "Loop on s2th13 -> Begin..." << endl;
    for (is2t = 0 ; is2t < N_s2t ; is2t++)
    {
        //s2th_13 = pow(10,(log10(lo_s2t) + double(is2t)*DeltaLog_s2t)); //logarithmic
        s2th_13 = lo_s2t + double(is2t)*DeltaLog_s2t; //linear

        for (ieps = 0 ; ieps < N_eps ; ieps++)
        {
            epsilon = lo_eps + double(ieps)*DeltaLin_eps;
            
            const int N_params = 18; //-- Number of parameter of the chi² function --//
            ROOT::Math::Functor f(&chi2,N_params); //-- Setting the function to be minimized by using Minuit --//
            //-- Steps
            double stp = 1.0e-3;
            double step[N_params] = {stp,stp,stp,stp,stp,stp,
                                     stp,stp,stp,stp,stp,stp,
                                     stp,stp,stp,stp,stp,stp};
            //-- Initial parameter values
            double start[N_params] = {0.0,0.0,0.0,0.0,0.0,0.0,
                                      0.0,0.0,0.0,0.0,0.0,0.0,
                                      0.0,0.0,0.0,0.0,0.0,0.0};
        
	    //-- Calling Minuit function setting
            min->SetFunction(f);

	    //-- Setting variables
            double lim = 1.0e-3;
            min->SetLimitedVariable(0,  "e_1", start[0],  step[0],  -lim, lim);
            min->SetLimitedVariable(1,  "e_2", start[1],  step[1],  -lim, lim);
            min->SetLimitedVariable(2,  "e_3", start[2],  step[2],  -lim, lim);
            min->SetLimitedVariable(3,  "e_4", start[3],  step[3],  -lim, lim);
            min->SetLimitedVariable(4,  "e_5", start[4],  step[4],  -lim, lim);
            min->SetLimitedVariable(5,  "e_6", start[5],  step[5],  -lim, lim);
            min->SetLimitedVariable(6,  "n_1", start[6],  step[6],  -lim, lim);
            min->SetLimitedVariable(7,  "n_2", start[7],  step[7],  -lim, lim);
            min->SetLimitedVariable(8,  "n_3", start[8],  step[8],  -lim, lim);
            min->SetLimitedVariable(9,  "n_4", start[9],  step[9],  -lim, lim);
            min->SetLimitedVariable(10, "n_5", start[10], step[10], -lim, lim);
            min->SetLimitedVariable(11, "n_6", start[11], step[11], -lim, lim);
            min->SetLimitedVariable(12, "a_1", start[12], step[12], -lim, lim);
            min->SetLimitedVariable(13, "a_2", start[13], step[13], -lim, lim);
            min->SetLimitedVariable(14, "a_3", start[14], step[14], -lim, lim);
            min->SetLimitedVariable(15, "a_4", start[15], step[15], -lim, lim);
            min->SetLimitedVariable(16, "a_5", start[16], step[16], -lim, lim);
            min->SetLimitedVariable(17, "a_6", start[17], step[17], -lim, lim);
            min->SetErrorDef(2.3);
	    
	    //-- Calling Minuit minimization
            min->Minimize();
        
            const double *xs = min->X();
    
            chi2Surface_file << s2th_13 << "\t" << epsilon << "\t" << min->MinValue() << endl;
            
            //-- Uncomment if you want to print the pull parameters (also Line 205)
            //minimPullT_file  << s2th_13 << "\t" << epsilon << "\t" << xs[0] << "\t" << xs[1] << "\t" << xs[2] << "\t" << xs[3] << "\t" << xs[4] << "\t" << xs[5] << "\t" << xs[6] << "\t" << xs[7] << "\t" << xs[8] << "\t" << xs[9] << "\t" << xs[10] << "\t" << xs[11] << "\t" << xs[12] << "\t" << xs[13] << "\t" << xs[14] << "\t" << xs[15] << "\t" << xs[16] << "\t" << xs[17] << "\t" << min->MinValue() << endl;
        }
        chi2Surface_file << endl;
        //-- Uncomment if you want to print the pull parameters
        //minimPullT_file  << endl;

        if (is2t%10 == 0)
            std::cout << "Succesful run for sin^2(th13) = " << s2th_13 << "!! \t" << min->MinValue() << endl;
    }

    std::cout << "Succesful run!!" << endl;
    
    return 0;
}
