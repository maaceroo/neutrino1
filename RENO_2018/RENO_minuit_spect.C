//-----------------------------------------------------------------------------------//
//--        RENO_minuit.C - By M.A. Acero O., A.A. Aguilar-A. - 2020-02-07         --//
//-----------------------------------------------------------------------------------//
//-------- Using 'NumericalMinimization.C' macro example to minimize a chi2  --------//
//-------- function. It uses a Minimizer class in ROOT.                      --------//
//-------- input : minimizer name + algorithm name                           --------//
//-------- randomSeed: = <0 : fixed value:                                   --------//
//--------                0 random with seed 0;                              --------//
//--------               >0 random with given seed.                          --------//
//--      Analysis of the RENO data from G. Bak et al., PRL121, 201801 (2018)      --//
//-----------------------------------------------------------------------------------//
//-------- This macro can be executed under ROOT typing                      --------//
//-------- "root[0] .x RENO_minuit2.C"                                       --------//
//-------- The results are printed in a file, with the following information: -------//
//--------    sin^2(2th_13)   a   chi^2(min)                                 --------//
//-------- In this case, chi^2(min) is the minimum chi^2 for the pull        --------//
//-------- terms.                                                            --------//
//-----------------------------------------------------------------------------------//
#include "constants.h"
#include <math.h>
#include <iostream>
#include <string>
//---*****************************************************************************---//
//------------------------ CONSTANTS ------------------------------------------------//
//---*****************************************************************************---//
//Number of Antineutrino Detectors
//#define nDet 2
//Number of Nuclear Reactors
//#define nRea 6
//-- N_dm2 and N_s2t are defined in the files 'constants.h'
double ran = N_dm2 * N_s2t;
//---*****************************************************---//

//---*****************************************************---//
//-- Global variables and quantities ------------------------//
//---*****************************************************---//
//(IBD candidates)/(DAQ live time -days-) from PRL128 (2018)
//double IBDrate_data[nDet][2] = { {470.53,0.51},{47.06,0.15} };
//double IBDrate_data[nDet][2] = { {461.00,0.58},{44.82,0.18} };
//IBD rate (per day), total background and efficiencies (PRL128 (2018))
double totalBgd[nDet][2] = { {9.53,0.28},{2.24,0.10} };
double emuem[nDet] ={0.7647,0.7647};
double daqTime[nDet] = {1807.88,2193.04};
//---*****************************************************---//
// Information obtained by executing the script "RENO_osc_rate.C"
//IBD rate per day w/o oscillations
double noOsc_IBDrate_perday[nDet];
//---*****************************************************---//
double s2th_13; //oscillation parameter to be fitted
double dm2_ee;   //oscillation parameter to be fitted
double a; //absolute normalization factor to be fitted
double e; //Escale Pull term
double epsilon; //Efficiency Pull term
double wrd_array[nDet][nRea];
double fr[nRea];
double b_d[nDet];
double spc[nDet][NB];
double NoscTot[nDet];
double noNoscTot[nDet];
double ratio[NB];
int iAD;
int iNR;
double xbins[NB+1];
//---*****************************************************---//
//TH1F *ratio_histo[nDet];
TH1F *data_spect_histo[nDet];
//---*****************************************************---//
TH1F *nosc_spect_hist[nDet];
TH1F *nosc_spect_hist_bf[nDet];
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
    fr[0]   = xx[4];
    fr[1]   = xx[5];
    fr[2]   = xx[6];
    fr[3]   = xx[7];
    fr[4]   = xx[8];
    fr[5]   = xx[9];
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
    double sB[2]            = {0.0326,0.0561}; //from PRL128 (2018)
    double seps             = 0.0021;          //from PRL128 (2018)
    double sfr_r            = 0.009;           //from PRL128 (2018)
    double sesc             = 0.0015;          //from PRL128 (2018)
    //double sesc             = 0.0030;          //from PRL128 (2018)
  
    SurvPavg1 = NoscTot[0]/noNoscTot[0];
    SurvPavg2 = NoscTot[1]/noNoscTot[1];

    //- 12.11.2018 -Begins - This part allows to include properly the energy scale pull term
    double f_ePos[NB] = {0.0};
    double f_eNeg[NB] = {0.0};
    double spcNew[nDet][NB];
    double Eold_i,Enew_i,Enew_j;

    if (e < 0) {
        for (int iBIN = 0 ; iBIN < NB ; iBIN++) {
            Eold_i   = xbins[iBIN];
            Enew_i = (1+e)*xbins[iBIN];
            Enew_j = (1+e)*xbins[iBIN+1];
            f_eNeg[iBIN] = (Eold_i - Enew_i)/(Enew_j - Enew_i);
        }
    }//if
    else {
        for (int iBIN = 0 ; iBIN < NB ; iBIN++){
            Eold_i   = xbins[iBIN+1];
            Enew_i = (1+e)*xbins[iBIN+1];
            Enew_j = (1+e)*xbins[iBIN];
            f_ePos[iBIN] = (Enew_i - Eold_i)/(Enew_i - Enew_j);
        }
    } //else
    //- 12.11.2018 -Ends

    //-- Adding artificial noise to the no-osccilated spectra - 19-07-2019 ----
    //TF1 *gauNoise1 = new TF1("gauNoise","exp(-0.5*((x)/[0])^2)",-5.0e1,5.0e1);
    //gauNoise1->SetParameter(0,10.0);
    //TF1 *gauNoise2 = new TF1("gauNoise","exp(-0.5*((x)/[0])^2)",-5.0e-0,5.0e-0);
    //gauNoise2->SetParameter(0,1.0);
    //-------------------------------------------------------------------------
    for (int iBIN = 0 ; iBIN < NB ; iBIN++){
    //-- Measured IDB events of the dth Antineutrino Detector (background is substracted)
	//Corregir: normalizar data_spect_histo para que Nobs tenga unidades de número de eventos
	//data_spect_histo con área 1 y por el ancho del bin (debe tener chipote en 6-7 MeV)
        Nobs1 = (data_spect_histo[0]->GetBinContent(iBIN+1))*IBDrate_data[0][0]*emuem[0]*daqTime[0];
        Nobs2 = (data_spect_histo[1]->GetBinContent(iBIN+1))*IBDrate_data[1][0]*emuem[1]*daqTime[1];
        //Nobs1 = (data_spect_histo[0]->GetBinContent(iBIN+1))*IBDrate_data[0][0]*daqTime[0];
        //Nobs2 = (data_spect_histo[1]->GetBinContent(iBIN+1))*IBDrate_data[1][0]*daqTime[1];

        sqrerror = ( Nobs2/(pow(Nobs1,2)) ) + ( (pow(Nobs2,2))/(pow(Nobs1,3)) );
        //std::cout << "error = " << sqrerror << std::endl;

        //- 12.11.2018 -Begins - This part allows to include properly the energy scale pull term
        if (e < 0) {
            if (iBIN == 0){
                spcNew[0][iBIN] = spc[0][iBIN] + f_eNeg[iBIN+1]*spc[0][iBIN+1];
                spcNew[1][iBIN] = spc[1][iBIN] + f_eNeg[iBIN+1]*spc[1][iBIN+1];
            }
            else if (iBIN<NB-1){
                spcNew[0][iBIN] = spc[0][iBIN] - f_eNeg[iBIN]*spc[0][iBIN] + f_eNeg[iBIN+1]*spc[0][iBIN+1];
                spcNew[1][iBIN] = spc[1][iBIN] - f_eNeg[iBIN]*spc[1][iBIN] + f_eNeg[iBIN+1]*spc[1][iBIN+1];
            }
            else{
                spcNew[0][iBIN] = spc[0][iBIN] - f_eNeg[iBIN]*spc[0][iBIN];
                spcNew[1][iBIN] = spc[1][iBIN] - f_eNeg[iBIN]*spc[1][iBIN];
            }
        }

        if (e >= 0) {
            if (iBIN == 0){
                spcNew[0][iBIN] = spc[0][iBIN] - f_ePos[iBIN+1]*spc[0][iBIN];
                spcNew[1][iBIN] = spc[1][iBIN] - f_ePos[iBIN+1]*spc[1][iBIN];
            }
            else if (iBIN<NB-1){
                spcNew[0][iBIN] = spc[0][iBIN] + f_ePos[iBIN]*spc[0][iBIN-1] - f_ePos[iBIN+1]*spc[0][iBIN];
                spcNew[1][iBIN] = spc[1][iBIN] + f_ePos[iBIN]*spc[1][iBIN-1] - f_ePos[iBIN+1]*spc[1][iBIN];
            }
            else{
                spcNew[0][iBIN] = spc[0][iBIN] + f_ePos[iBIN]*spc[0][iBIN-1];
                spcNew[1][iBIN] = spc[1][iBIN] + f_ePos[iBIN]*spc[1][iBIN-1];
            }
        }
        //- 12.11.2018 -Ends
        
        Nexp1 = 0.0;
        Nexp2 = 0.0;
        for (iNR = 0 ; iNR < nRea ; iNR++){
            Nexp1 += (1+fr[iNR])*wrd_array[0][iNR]*spcNew[0][iBIN]*(SurvPavg1*noOsc_IBDrate_perday[0]/NoscTot[0])*emuem[0]*daqTime[0];

            Nexp2 += (1+fr[iNR])*wrd_array[1][iNR]*spcNew[1][iBIN]*(SurvPavg2*noOsc_IBDrate_perday[1]/NoscTot[1])*emuem[1]*daqTime[1];
        }
        Nexp1 = fudge*Nexp1; //fudge normalization factor defined in 'constants.h'

        // Compute of the ratios of Data and expect spectra
        OFN = Nobs2/Nobs1;
        TFN = ( (1.0 + epsilon)*Nexp2 - b_d[1] )/( (1.0 + epsilon)*Nexp1 - b_d[0]);
        //-- Adding artificial noise to the no-osccilated spectra -----------------------------
        //double noise1 = gauNoise1->GetRandom();
        //double noise2 = gauNoise2->GetRandom();
        //TFN = ( (1.0 + epsilon)*Nexp2 - b_d[1] )/( (1.0 + epsilon)*Nexp1 - b_d[0] ) + noise;
        //TFN = ( (1.0 + epsilon)*Nexp2 - b_d[1]  + noise2)/( (1.0 + epsilon)*Nexp1 - b_d[0]  + noise1);
        //-------------------------------------------------------------------------------------
        //cout << " OFN = " << OFN << " Nobs = " << Nobs <<endl;

        // Chi^2 funtion
        sqr_chi += pow( OFN - TFN ,2 )/sqrerror;
    }//-- for(iBin) - End
  
    for (iAD = 0 ; iAD < nDet ; iAD++){
      //-- Background error of the dth Antineutrino Detector
        sqr_chi +=  pow(b_d[iAD]/sB[iAD],2) ;
    }
  
    for (iNR = 0 ; iNR < 6 ; iNR++)
    {
        sqr_chi += pow(fr[iNR]/sfr_r,2);
    }

    sqr_chi += pow(epsilon/seps,2) + pow(e/sesc,2) ;

  return sqr_chi;
  
}
//---*****************************************************---//

int RENO_minuit_spect(const char * minName = "Minuit",
			   const char *algoName = "" ,
			   int randomSeed = 10)
//int randomSeed = -1)
{
    cout << "Let's begin..." << endl;
   
    TFile *wrd_File = new TFile("files_root/ldist_RENO.root","READ");
    TH1F *wrd_histo0 = ((TH1F*)(wrd_File->Get("histo_ldist_RENO_near")));;
    TH1F *wrd_histo1 = ((TH1F*)(wrd_File->Get("histo_ldist_RENO_far")));;
    double mm;
    for (int iNR = 0 ; iNR < nRea ; iNR++)
    {
        wrd_array[0][iNR] = wrd_histo0->GetBinContent(iNR+1);
        wrd_array[1][iNR] = wrd_histo1->GetBinContent(iNR+1);
         //cout << wrd_array[0][iNR] << "   " << wrd_array[1][iNR]  << endl;
    }
    
    //-------------------
    // Energy Histograms
    //-------------------
    TFile *fenergy = new TFile("files_root/RENOplots.root","read");
    //The histogram of near and far data spectra
    for(int n = 0 ; n < nDet ; n++){
        data_spect_histo[n] = (TH1F*) fenergy->Get(Form("data_spect_histo_%d",n));
        double dfactor = 1.0/data_spect_histo[n]->Integral();
        data_spect_histo[n]->Scale(dfactor);
    
        //ratio_histo[n] = (TH1F*) fenergy->Get(Form("ratio_histo_%d",n));
        //double rfactor = 1.0/ratio_histo[n]->Integral();
        //ratio_histo[n]->Scale(rfactor);
    }

    double delta_bins2 = (5.6 - 1.2)/22; // 0.2 MeV/bin
    
    for (int i = 0 ; i < (NB-3) ; i++){
        xbins[i] = 1.2 + delta_bins2*i;
    }
    xbins[23] = xbins[22] + 0.4;
    xbins[24] = xbins[23] + 0.4;
    xbins[25] = 8.4 - 1.4;
    xbins[26] = 8.4;

    for(int iAD = 0 ; iAD < nDet ; iAD++)
    {
        nosc_spect_hist[iAD] = new TH1F(Form("nosc_spect_hist_%d",iAD),"",NB,xbins);
        nosc_spect_hist_bf[iAD] = new TH1F(Form("nosc_spect_hist_bf_%d",iAD),"",NB,xbins);
    }
  
    cout << "Minimization settings..." << endl;
    ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer(minName, algoName);
  
    //-- Set tolerance , etc...
    min->SetMaxFunctionCalls(1000000); // for Minuit/Minuit2
    min->SetMaxIterations(10000);      // for GSL
    min->SetTolerance(0.001);
    min->SetPrintLevel(-1);
  
    //-- File to get noOsc normalizations
    ifstream IBDrates_file("files/RENO_noOsc_IBDrates_perday.txt");
    cout << "Reading noOsc normalizations file ..." << endl;
    for (int i=0; i< nDet; i ++)
        {
            IBDrates_file >> noOsc_IBDrate_perday[i];
            cout << "noOscIBD_rates_perday " << i << ": " << noOsc_IBDrate_perday[i] << endl;
        }//for

    //-- File to print oscillation parameters and chi2 values
    ofstream chi2Surface_file;
    string s2t_dm2 = "files/chi2_s2t-dm2_surface_spect.txt";  //(sin^2(2th13), a, chi^2_min)
    ofstream minimPullT_file;
    string pterms = "files/chi2_pullTerms_spect.txt";
    //string s2t_eps = "chi2_s2t_curve.txt"; //(sin^2(2th13), chi^2_min)
    chi2Surface_file.open((s2t_dm2).c_str());
    minimPullT_file.open((pterms).c_str());
  
    //ifstream grid_file("./files/RENO_gridOscSpectra_test.txt",ifstream::in);
    cout << "Reading file - Loop in progress..." << endl;
    //string grid_spectra = "./files/RENO_gridOscSpectra_test.txt";
    //std::stringstream line;
    ifstream grid_file;
    grid_file.open("./files/RENO_gridOscSpectra_test.txt");
    std::cout << "Is the file open? \t";
    if (grid_file.is_open() == 1)
        std::cout << "Yes! \n" << std::endl;
    else
        std::cout << "No! Check it out! \n" << std::endl;
    int iad = 0;
    int first2 = 1;
  
    while (grid_file >> iAD >> s2th_13 >> dm2_ee  >> spc[iad][0] >> spc[iad][1] >> spc[iad][2] >> spc[iad][3] >> spc[iad][4] >> spc[iad][5] >> spc[iad][6] >> spc[iad][7] >> spc[iad][8] >> spc[iad][9] >> spc[iad][10] >> spc[iad][11] >> spc[iad][12] >> spc[iad][13] >> spc[iad][14] >> spc[iad][15] >> spc[iad][16] >> spc[iad][17] >> spc[iad][18] >> spc[iad][19] >> spc[iad][20] >> spc[iad][21] >> spc[iad][22] >> spc[iad][23] >> spc[iad][24] >> spc[iad][25] >> NoscTot[iad])
    {//file loop
        //std::cout << "Inside the while loop!" << std::endl;
        if (!grid_file.good()){
            std::cout << "\n \t File is not good! \n" << std::endl;
            break;
        }
        
        if(first2 <= 2)
        {
            for(int ibin = 1 ; ibin <= NB ; ibin++)
            {
                double bincont = spc[iad][ibin-1]*(noOsc_IBDrate_perday[iAD-1]/NoscTot[iad])*emuem[iAD-1]*daqTime[iAD-1];
                //double bincont = spc[iad][ibin-1]*(noOsc_IBDrate_perday[iAD-1]/NoscTot[iad])*1.0*daqTime[iAD-1];
                double bincontent = bincont/0.2;
                nosc_spect_hist[iAD-1]->SetBinContent(ibin,bincont);
            }
	  
            cout << " iAD = " << iAD  /*<< " iNR = " << iNR */ << "  first2 = " << first2 << endl;
            cout << "Number of events: " << nosc_spect_hist[iAD-1]->Integral() << endl;
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
            const int N_params = 10; //-- Number of parameter of the chi² function --//
            ROOT::Math::Functor f(&chi2,N_params); //-- Setting the function to be minimized by using Minuit --//
            //-- Steps
            double stp = 1.0e-4;
            double step[N_params] = {stp,stp,stp,stp,stp,stp,stp,stp,stp,stp};

            //-- Initial parameter values
            double start[N_params] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

            //-- Calling Minuit function setting
            min->SetFunction(f);
            
            //-- Setting variables
            double lim = 5.0e-2;
	  
            min->SetLimitedVariable(0,  "epsilon", start[0],  step[0],  -lim, lim);
            min->SetLimitedVariable(1,  "e",       start[1],  step[1],  -lim, lim);
            min->SetLimitedVariable(2,  "b_0",     start[2],  step[2],  -lim, lim);
            min->SetLimitedVariable(3,  "b_1",     start[3],  step[3],  -lim, lim);
    
            min->SetLimitedVariable(4,  "f_0",     start[4],  step[4],  -lim, lim);
            min->SetLimitedVariable(5,  "f_1",     start[5],  step[5],  -lim, lim);
            min->SetLimitedVariable(6,  "f_2",     start[6],  step[6],  -lim, lim);
            min->SetLimitedVariable(7,  "f_3",     start[7],  step[7],  -lim, lim);
            min->SetLimitedVariable(8,  "f_4",     start[8],  step[8],  -lim, lim);
            min->SetLimitedVariable(9,  "f_5",     start[9],  step[9],  -lim, lim);
	  	
            /*
             min->SetFixedVariable(0,  "epsilon", start[0]);
             min->SetFixedVariable(1,  "e",       start[1]);
             min->SetFixedVariable(2,  "b_0",     start[2]);
             min->SetFixedVariable(3,  "b_1",     start[3]);
             min->SetFixedVariable(1,  "e",       start[1]);
             */
            /*
             min->SetFixedVariable(4,  "f_0",     start[4]);
             min->SetFixedVariable(5,  "f_1",     start[5]);
             min->SetFixedVariable(6,  "f_2",     start[6]);
             min->SetFixedVariable(7,  "f_3",     start[7]);
             min->SetFixedVariable(8,  "f_4",     start[8]);
             min->SetFixedVariable(9,  "f_5",     start[9]);
             */
            min->SetErrorDef(2.3);

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
            minimPullT_file  << s2th_13 << "\t" << dm2_ee  << "\t"
                << xs[0] << "\t" << xs[1] << "\t" << xs[2] << "\t" << xs[3] << "\t" << xs[4] << "\t"
                << xs[5] << "\t" << xs[6] << "\t" << xs[7] << "\t" << xs[8] << "\t" << xs[9] << "\t"
                << chi2Min << endl;
            //}
            
            if (dm2_ee == hi_dm2)
            {
                chi2Surface_file << endl;
                cout << s2th_13 << "\t" << dm2_ee << "\t" << min->MinValue() << endl;
            }
        }
    }
    
    std::cout << "Succesful run!!" << endl;
    grid_file.close();
    
    // Chi2 minimun value //
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
                n = ll;
            }
        }
    }
    cout << "sin2t_13 = " << matr[0][n]   << " dm2_ee = " << matr[1][n]
        << " chi2 = " << matr[2][n]  << " n = " << n  << endl;
    chi2min  << matr[0][n]   <<  "\t" << matr[1][n]  << "\t" << matr[2][n]  << endl;
    chi2min <<endl;
    
  return 0;

}//-END
