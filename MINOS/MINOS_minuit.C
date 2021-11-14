//----------------------------------------------------------//
//-- db_minuit_spec.C  -  By M.A. Acero O.  -  2021-10-22 --//
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
//---****************************************************---//
//---------- CONSTANTS -------------------------------------//
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
double NoscTot;
double noNoscTot;
double totalBgd = 8.92;  // integral of bkgd spectrum
double xbins_numu110[NB_numu110+1];

//Table I PRL110.251801(2013)
//{No Osc (Simulated), Oscilated (BF), Observed}
double NuMu_Events[3]  = {3201.0, 2543.0, 2579.0};
double NuMuB_Events[3] = {313.0,   227.0,  226.0};
//---*****************************************************---//
TH1F *numu_data_spect_histo;
TH1F *numu_bkgd_spect_histo;
TH1F *numub_data_spect_histo;
TH1F *numub_bkgd_spect_histo;
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
//-------------------
TF1 *fFit4;
TF1 *fFit7;
TF1 *fFitGD;

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
    double Nd           = 0.0;
    double Nmc          = 0.0;
    //Section 6.2, J. Mitchell PhD Thesis https://minos-docdb.fnal.gov/cgi-bin/sso/RetrieveFile?docid=8704&filename=jessmitchellthesis.pdf&version=1
    //See also https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.106.181801
    double sigeps       = 0.016;
    double sigepsNC     = 0.20;
    double sigensc      = 0.06;

    double spcNumuNew[NB_numu110];
    
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
        delta_spc = ensc*(fFitGD->Eval(binCenter));
        avgSurvProb_bin = spcNumu[iBIN]/spcNumuNoOsc[iBIN];
        spcNumuNew[iBIN] = spcNumu[iBIN] + delta_spc*avgSurvProb_bin;

        //Nmc = ( (spcNumu[iBIN]/NoscTot)*(NuMu_Events[0]*SurvPavg) - (1-epsNC)*numu_simu_bg )*(1-eps);
        Nmc = ( (spcNumuNew[iBIN]/NoscTot)*(NuMu_Events[0]*SurvPavg) - (1-epsNC)*numu_simu_bg )*(1-eps);
        //cout << "iBin = " << iBIN << "  Nmc = " << Nmc << " Nd = " << Nd << endl;
        
        nll += 2*(Nmc - Nd + Nd*log(Nd/Nmc));
    }
    
    nll = nll + pow(eps,2)/(2*pow(sigeps,2)) + pow(epsNC,2)/(2*pow(sigepsNC,2)) + pow(ensc,2)/(2*pow(sigensc,2));
    
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
    fFit4  = (TF1*)(Esc_File->Get("fFit4"));
    fFit7  = (TF1*)(Esc_File->Get("fFit7"));
    fFitGD = (TF1*)(Esc_File->Get("fFitGD"));
    //-------------------
    // Energy Histograms
    //-------------------
    TFile *fenergy = new TFile("./MINOS_spectra_PRL108-PRL110.root","read");
    double bfactor;
    //Data
    numu_data_spect_histo = (TH1F*) fenergy->Get("numu110_data_histo");
    double dfactor = 1.0/numu_data_spect_histo->Integral();
    numu_data_spect_histo->Scale(dfactor);
    //Background
    numu_bkgd_spect_histo = (TH1F*) fenergy->Get("numu110_bkgd_histo");
    bfactor = 1.0/numu_bkgd_spect_histo->Integral();
    numu_bkgd_spect_histo->Scale(bfactor);
    
    // histogram binning
    double delta_bins110 = (9.0 - 0.5)/17.0; // 0.5 GeV/bin
    for (int i = 0 ; i < (NB_numu110 - 5) ; i++)
        xbins_numu110[i] = lo110 + i*delta_bins110;
    xbins_numu110[18] = xbins_numu110[17] + 0.75;
    xbins_numu110[19] = xbins_numu110[18] + 0.75;
    xbins_numu110[20] = xbins_numu110[19] + 0.75;
    xbins_numu110[21] = xbins_numu110[20] + 0.75;
    xbins_numu110[22] = xbins_numu110[21] + 1.0;
    xbins_numu110[23] = hi110;

    numu_noosc_spect_hist    = new TH1F("numu_noosc_spect_hist",   "",NB_numu110,xbins_numu110);
    numu_nu_nosc_spect_hist = new TH1F("numu_nu_nosc_spect_hist","",NB_numu110,xbins_numu110);
    numu_wosc_spect_hist    = new TH1F("numu_wosc_spect_hist",   "",NB_numu110,xbins_numu110);
    numu_nu_wosc_spect_hist = new TH1F("numu_nu_wosc_spect_hist","",NB_numu110,xbins_numu110);
    numu_bfit_spect_hist    = new TH1F("numu_bfit_spect_hist",   "",NB_numu110,xbins_numu110);
        
    cout << "wrd array -> Done..." << endl;
    
    cout << "Minimization settings..." << endl;
    ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer(minName, algoName);
    
    //-- Set tolerance , etc...
    min->SetMaxFunctionCalls(1000000); // for Minuit/Minuit2
    min->SetMaxIterations(10000);      // for GSL
    min->SetTolerance(0.001);
    min->SetPrintLevel(-1);
    
    //-- File to print oscillation parameters and chi2 values
    ofstream numu_chi2Surface_file;
    string s2t_dm2 = "/data/numu_chi2_s2t-dm2_surface.txt";  //(sin^2(2th13), dm2, chi^2_min)
    numu_chi2Surface_file.open(filePath+(s2t_dm2).c_str());
    
    //-- File to print minimized pull parameters
    ofstream numu_minimPullT_file;
    string PullT = "/data/numu_chi2_pullT_surface.txt";
    numu_minimPullT_file.open(filePath + (PullT).c_str());

    //ifstream file("files_data/db_gridOscSpectra_1M.txt"); //50x50 parameter space-grid
    string numu_grid_spectra = "/data/minosNuMu_gridOscSpectra.txt";
    ifstream numu_grid_file;
    numu_grid_file.open(filePath + (numu_grid_spectra).c_str());
    cout << "Reading file - Loop in progress..." << endl;
    int first = 1;
    std::cout << "Is the file open?  " << numu_grid_file.is_open() << "\n" << std::endl;
    
    while (numu_grid_file >> s2th_23 >> dm2_32 >> spcNumu[0] >> spcNumu[1] >> spcNumu[2] >> spcNumu[3] >> spcNumu[4] >> spcNumu[5] >> spcNumu[6] >> spcNumu[7] >> spcNumu[8] >> spcNumu[9] >> spcNumu[10] >> spcNumu[11] >> spcNumu[12] >> spcNumu[13] >> spcNumu[14] >> spcNumu[15] >> spcNumu[16] >> spcNumu[17] >> spcNumu[18] >> spcNumu[19] >> spcNumu[20] >> spcNumu[21] >> spcNumu[22] >> NoscTot)
    {//file loop

        if(first <= 1)
        {
            for(int ibin = 1 ; ibin <= 23 ; ibin++)
            {
	        spcNumuNoOsc[ibin-1] = spcNumu[ibin-1];
                double bincont = spcNumu[ibin-1]*(NuMu_Events[0]/NoscTot);
                numu_noosc_spect_hist->SetBinContent(ibin,bincont);
            }
            cout << "first = " << first << endl;
            cout << "Number of events: " << numu_noosc_spect_hist->Integral() << endl;
            noNoscTot = NoscTot;
            first++;
        }//if first loop END

        const int N_params = 3; //-- Number of parameter of the chi² function --//
        ROOT::Math::Functor f(&chi2,N_params); //-- Setting the function to be minimized by using Minuit --//
        //-- Steps
        double stp = 1.0e-4;
        double step[N_params] = {stp,stp,stp};
        //-- Initial parameter values
        double st = 0.0;
        double start[N_params] = {st,st,st};

        //-- Calling Minuit function setting
        min->SetFunction(f);
        
        //-- Setting variables
        double lim = 1.0e-1;
        min->SetLimitedVariable(0,  "eps",   start[0],  step[0],  -lim, lim);
        min->SetLimitedVariable(1,  "epsNC", start[1],  step[1],  -lim, lim);
        min->SetLimitedVariable(2,  "ensc",  start[2],  step[2],  -lim, lim);
//        min->SetFixedVariable(0,  "eps",   start[0]);
//        min->SetFixedVariable(1,  "epsNC", start[1]);
//        min->SetFixedVariable(2,  "ensc",  start[2]);
        
        min->SetErrorDef(2.3);
        
        //-- Calling Minuit minimization
        min->Minimize();
        const double *xs = min->X();
        
        double chi2Min = min->MinValue();
        numu_chi2Surface_file << s2th_23 << "\t" << dm2_32 << "\t" << chi2Min << endl;
        
        if (chi2Min < 0.0) {
            cout << "Critical error: chi2Min is negative!  " << chi2Min << endl;
            break;
        }
        
        numu_minimPullT_file  << s2th_23 << "\t" << dm2_32 << "\t" << xs[0] << "\t" << xs[1] << "\t" << xs[2] << "\t" << min->MinValue() << endl;
        //}
        if (dm2_32 == hi_dm2)
        {
            numu_chi2Surface_file << endl;
            cout << s2th_23 << "\t" << dm2_32 << "\t" << min->MinValue() << endl;
        }
        //numu_minimPullT_file  << endl;
        //if (is2t%10 == 0)
        //std::cout << "Succesful run for sin^2(th13) = " << s2th_23 << "!! \t" << min->MinValue() << endl;
        //}
        
    }//file loop END
    std::cout << "Successful run!!" << endl;
    numu_minimPullT_file.close();
    numu_chi2Surface_file.close();

    //------------------------------------------------------------------------------
    //------------------------------------------------------------------------------
    return 0;
}

