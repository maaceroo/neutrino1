//----------------------------------------------------------//
//-- db_minuit_spec.C  -  By M.A. Acero O.  -  2016-04-11 --//
//----------------------------------------------------------//
//-- Using 'NumericalMinimization.C' macro example to     --//
//-- minimize a chi2 function. It uses a Minimizer class  --//
//-- in ROOT.                                             --//
//-- input : minimizer name + algorithm name              --//
//-- randomSeed: = <0 : fixed value: 0 random with seed 0;--//
//-- >0 random with given seed.                           --//
//----------------------------------------------------------//
//-- Analysis of the Daya Bay data from F.P. An et al.,   --//
//-- PRL112 061801 (2014)                                 --//
//----------------------------------------------------------//
//-- This macro can be executed under ROOT typing "root[0] .x db_minuit_spectral.C" --//
//-- In the simplest case (only minimization) the script will print the message     --//
//-- "Minimum: f(x[i]): chi^2(min)" where x[i] are the values of the parameters     --//
//-- which minimize the chi^2-function and chi^2(min) is the minimum value of the   --//
//-- chi^2-function. In the complete case, the results are printed in a  file, with --//
//-- the following information:        sin^2(2th_13)  delta(m31)^2  chi^2(min)      --//
//-- In this case, chi^2(min) is the minimum chi^2 for  the pull terms.             --//
//---*****************************************************---//
//-- FOR THE SPECTRAL ANALYSIS ONLY                        --//
//---*****************************************************---//
//---*****************************************************---//
//---------- CONSTANTS --------------------------------------//
//---*****************************************************---//
#include "constants.h"
// histogram binning for the simulated data
//#define  NB 26
//#define  lo 0.7
//#define  hi 12.0
//Number of Antineutrino Detectors
//#define nAD 6
//Number of Nuclear Reactors
//#define nNR 6
//Fixed neutrino oscillations parameters
//#define dm2_21 7.59e-5 //eV^2,                //PRL 108 171803 (2012)
//#define s22th_12 0.861
//For the sin^2(2th_13) loop
//#define N_s2t  200                           //number of points in the grid
//#define lo_s2t 0.01                         //sin^2(2th_13) min
//#define lo_s2t 0.0                          //sin^2(2th_13) min
//#define hi_s2t 0.3                          //sin^2(2th_13) max
//For the delta(m31)^2 loop
//#define N_dm2  200                           //number of points in the grid
//#define lo_dm2 1.0e-4                       //delta(m31)^2 min
//#define hi_dm2 1.0e-2                       //delta(m31)^2 max
//---*****************************************************---//

//---*****************************************************---//
//-- Global variables and quantities ------------------------//
//---*****************************************************---//
//EH1(AD1, AD2),EH2(AD3),EH3(AD4, AD5, AD6)
//(IBD candidates)/(DAQ live time -days-) from PRL 112 061801 (2014)
double IBDrate_data[nAD][2] = { {530.31,1.67},{536.75,1.68},{489.93,1.61},{ 73.58,0.62},{ 73.21,0.62},{72.35,0.62} };
//IBD rate (per day), total background and efficiencies (PRL 112 061801 (2014))
double totalBgd[nAD][2]     = { {13.20,0.98},{13.01,0.98},{ 9.57,0.71},{ 3.52,0.14},{ 3.48,0.14},{3.43,0.14} };
double emuem[nAD]           = {0.7957,0.7927,0.8282,0.9577,0.9568,0.9566};
double daqTime[nAD]         = {191.001,191.001,189.645,189.779,189.779,189.779};
//---*****************************************************---//
// Information obtained by executing the script "db_osc_rate.C"
//IBD rate per day w/o oscillations
double noOsc_IBDrate_perday[nAD];// = {663.15,673.95,591.86,78.75,78.46,77.58};
//---*****************************************************---//
//const int dim = N_s2t*N_dm2;
double s2th_13;     //oscillation parameter to be fitted
double dm2_31;      //oscillation parameter to be fitted
double wrd_array[nAD][nNR];
double epsilon;     //normalization factor
double alpha[nNR];  //'uncorrelated reactor uncertainty' pull term
double eps_d[nAD];  //'uncorrelated detection uncertainty' pull term
double eta_d[nAD];  //'background' pull term
double spc[nAD][NB];
double NoscTot[nAD];
double noNoscTot[nAD];
int iAD;
//---*****************************************************---//
const int nEH = 3;
TH1F *data_spect_histo[nEH];
TH1F *bkgd_spect_histo[nEH];
//---*****************************************************---//
TH1F *nosc_spect_hist[nAD];
TH1F *nu_nosc_spect_hist[nAD];
TH1F *wosc_spect_hist[nAD];
TH1F *nu_wosc_spect_hist[nAD];
TH1F *bfit_spect_hist[nAD];
//-------------------
// Covariance-matrix Histograms
//-------------------
//** Declarar la matriz de covarianza fraccionaria como un TH2F
TH2F *fracCovaMatrix_hist;    //Se carga dentro de la función principal db_minuot_spec()
//TFile *fracCovaMatrix_File = new TFile("files_data/db_CovaMatrix_6AD_6x26bins.root","READ");
//fracCovaMatrix_histo->(TH1F*)(fracCovaMatrix_File->Get("covaMat_histo_6x6Det"));
int NBx_cov = NB*nAD;
int NBy_cov = NB*nAD;
TMatrixD statErroMatrix_matrix(NBx_cov,NBy_cov);
TMatrixD fracCovaMatrix_matrix(NBx_cov,NBy_cov);  //Se llena en la función principal db_minuot_spec()
TMatrixD fullCovaMatrix_matrix(NBx_cov,NBy_cov);
TMatrixD inv_fullCovaMatrix_matrix(NBx_cov,NBy_cov);
TMatrixD one_matrix(NBx_cov,NBy_cov);
TMatrixD delta_vector(NBx_cov,1);
TMatrixD transp_delta_vector(1,NBy_cov);
TMatrixD predi_vector(NBx_cov,1);

//---*****************************************************---//
//-- Chi-square function to be minimized --------------------//
//-- It has 19 pull parameters and 2 oscillation parameter --//
//-- which are to be fitted.                                -//
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
    epsilon  = xx[18];
    
    //---*****************************************************---//
    int iNR;
    double SurvPavg     = 0.0;
    double sqr_chi      = 0.0;
    double sqrerror     = 0.0;
    double Md           = 0.0;
    double Td           = 0.0;
    double Bd           = 0.0;
    double sB           = 0.0;
    double seps_d       = 1.0*0.002;
    double salph_r      = 1.0*0.008;
    
    for (iAD = 0 ; iAD < nAD ; iAD++){
        SurvPavg = NoscTot[iAD]/noNoscTot[iAD];
        for (int iBIN = 0 ; iBIN < NB ; iBIN++)
        {
            int index = iAD*NB + iBIN;
            
            //-- Predicted IBD from neutrino oscillations of the dth Antineutrino Detector
            Td = spc[iAD][iBIN]*(SurvPavg*noOsc_IBDrate_perday[iAD]/NoscTot[iAD])*emuem[iAD]*daqTime[iAD];
            //-- Measured IDB events of the dth Antineutrino Detector (background is substracted)
            int idx = -1;
            if (iAD < 2) {
                idx = 0;
            }
            else if (iAD == 2) {
                idx = 1;
            } else {
                idx = 2;
            }
            double data_sigplusbg = (data_spect_histo[idx]->GetBinContent(iBIN+1))*IBDrate_data[iAD][0];
            double simu_bg = (bkgd_spect_histo[idx]->GetBinContent(iBIN+1))*totalBgd[iAD][0]*emuem[iAD];
            Md = (data_sigplusbg - simu_bg)*daqTime[iAD];
            //cout << "iAD = " << iAD << "  iBin = " << iBIN << "  Md = " << Md << " Td = " << Td << endl;
            //-- Background of the dth Antineutrino Detector
            Bd = simu_bg;
            
            sqrerror = Md + Bd;
            
            //Statistics Error matrix
            statErroMatrix_matrix(index,index) = sqrerror;
            
            //-- Fraction of IBD contribution of the rth reactor to the dth AD
            //-- determined by baselines and reactor fluxes
            double wrd = 0.0;
            for (iNR = 0 ; iNR < nNR ; iNR++)
                wrd += wrd_array[iAD][iNR]*alpha[iNR];
            
            predi_vector(index,0) = Td*(1.0 + epsilon + eps_d[iAD] + wrd) - eta_d[iAD];
            delta_vector(index,0) = Md - Td*(1.0 + epsilon + eps_d[iAD] + wrd) + eta_d[iAD];
            transp_delta_vector(0,index) = delta_vector(index,0);
            
            //cout << "iAD = " << iAD << "  iBin = " << iBIN << "  Md = " << Md << " Td = " << Td
            //<< " delta = " << delta_vector(index,0) << endl;
            
            //-- Testing a fake normalization factor
            //double test = 0.0;
            //sqr_chi += pow( (Md - Td*(1.0 + epsilon - test + eps_d[iAD] + wrd) + eta_d[iAD]) ,2 )/sqrerror;
            
            //Modificar para que chi2 = (vec.col)^T x Inv.Mat.Cov x (vec.col)
            //1. vec.col = vector columna con la Predicción en cada Bin
            //   Predicción en cada bin = Td*(1.0 + epsilon - test + eps_d[iAD] + wrd) - eta_d[iAD]
            //2. Matriz diagonal con sqrerror en la diagonal (Mat.Stat)
            //3. La Mat.Covarianza debe ser construida (se requiere db_CovaMatrix_6AD_6x26bins.root que contiene Mat.Cov.Fraccionaria_ij):
            //   Mat.Cov_ij = (Mat.Cov.Fraccionaria_ij * NPred_i * NPred_j) + Mat.Stat_ij
            //   Invertir Mat.Cov
            
        }
    }
    //cout << "here the test begins." << endl;
    //-- Tests
    //TMatrixD *test_mat = new TMatrixD(1,1);
    //TMatrixD *part_mat = new TMatrixD(NBx_cov,1);
    //part_mat = fracCovaMatrix_matrix*predi_vector;
    //test_mat = (predi_vector.T())*part_mat;
    //test_mat.Print();
    //cout << "here the test ends." << endl;
    
    for (int i = 0 ; i < NBx_cov ; i++) {
        for (int j = 0 ; j < NBy_cov ; j++) {
            fullCovaMatrix_matrix(i,j) = statErroMatrix_matrix(i,j)
            + predi_vector(i,0)*predi_vector(j,0)*fracCovaMatrix_matrix(i,j);
        }
    }
    //predi_vector.Print();
    //fullCovaMatrix_matrix.Print();
    /*    TVector values;
     TMatrixD vectors = fullCovaMatrix_matrix.EigenVectors(values);
     for (Int_t i = 0; i < values.GetNrows(); ++i) {
     TVectorD vector(TMatrixTColumn_const<double>(vectors, i));
     cout << "eigen-value " << i << " is " << values(i) << " with eigen-vector";
     vector.Print();
     }
     */
    if(fullCovaMatrix_matrix.Determinant() != 0){
        inv_fullCovaMatrix_matrix = fullCovaMatrix_matrix;
        inv_fullCovaMatrix_matrix.Invert();
    }
    else {
        cout << "Determinant = " << fullCovaMatrix_matrix.Determinant() << endl;
        cout << "Singular fullCovaMatrix!" << endl;
        cout << "Execution aborted!" << endl;
        //break;
    }
    
    /*    TMatrixD *a_vector(NBx_cov,1);
     a_vector = inv_fullCovaMatrix_matrix*delta_vector;
     
     TMatrixD *product_matrix(1,1);
     product_matrix = transp_delta_vector*a_vector;*/
    
    for (int i = 0 ; i < NBx_cov ; i++) {
        for (int j = 0 ; j < NBy_cov ; j++) {
            double delta_i = transp_delta_vector(0,i);
            double delta_j = delta_vector(j,0);
            double invMat_ij = inv_fullCovaMatrix_matrix(i,j);
            sqr_chi += delta_i*invMat_ij*delta_j;
        }
    }
    
    //sqr_chi = product_matrix(0,0);
    
    /*    for (iAD = 0 ; iAD < nAD ; iAD++)
     {
     //-- Background error of the dth Antineutrino Detector
     sB = totalBgd[iAD][1]*emuem[iAD]*daqTime[iAD];
     sqr_chi += pow(eps_d[iAD]/seps_d,2) + pow(eta_d[iAD]/sB,2);
     }
     
     for (iNR = 0 ; iNR < nNR ; iNR++)
     sqr_chi += pow(alpha[iNR]/salph_r,2);
     */
    
    return sqr_chi;
}
//---*****************************************************---//

int db_minuit_spec(const char * minName = "Minuit",
                   const char *algoName = "" ,
                   int randomSeed = +10)
//int randomSeed = -1)
{
    cout << "Let's begin..." << endl;
    
    //-------------------
    // Covariance-matrix Histograms
    //-------------------
    TFile *fracCovaMatrix_File = new TFile("files_data/db_CovaMatrix_6AD_6x26bins.root","READ");
    fracCovaMatrix_hist = (TH2F*)(fracCovaMatrix_File->Get("covaMat_histo_6x6Det"));
    //-- Creating the Rebinned Matrix with root tools.
    NBx_cov = NB*nAD;
    NBy_cov = NB*nAD;
    cout << NBx_cov << "  " << NBy_cov << endl;
    //statErroMatrix_matrix = new TMatrixD(NBx_cov,NBy_cov);
    //fracCovaMatrix_matrix = new TMatrixD(NBx_cov,NBy_cov);
    //fullCovaMatrix_matrix = new TMatrixD(NBx_cov,NBy_cov);
    //inv_fullCovaMatrix_matrix = new TMatrixD(NBx_cov,NBy_cov);
    
    //one_matrix = new TMatrixD(NBx_cov,NBy_cov);
    
    //delta_vector        = new TMatrixD(NBx_cov,1);
    //transp_delta_vector = new TMatrixD(1,NBx_cov);
    //predi_vector        = new TMatrixD(NBx_cov,1);
    
    fracCovaMatrix_matrix.Zero();
    statErroMatrix_matrix.Zero();
    fullCovaMatrix_matrix.Zero();
    inv_fullCovaMatrix_matrix.Zero();
    one_matrix.Zero();
    delta_vector.Zero();
    transp_delta_vector.Zero();
    predi_vector.Zero();

    int i_block = 0;
    int j_block = 0;
    cout << "Filling Covariant Matrix..." << endl;
    for (int i = 0 ; i < NBx_cov ; i++) {
        for (int j = 0 ; j < NBy_cov; j++) {
            i_block = i%26;
            j_block = j%26;
            fracCovaMatrix_matrix(i,j) = fracCovaMatrix_hist->GetBinContent(i+1,j+1);
            //if (i_block != j_block) {
            /*
             if ((i_block == j_block + 1) || (j_block == i_block + 1)) {
             fracCovaMatrix_matrix(i,j) *= 0.99;
             }
             else if ((i_block == j_block + 2) || (j_block == i_block + 2)) {
             fracCovaMatrix_matrix(i,j) *= 0.89;
             }
             else if ((i_block == j_block + 3) || (j_block == i_block + 3)) {
             fracCovaMatrix_matrix(i,j) *= 0.84;
             }
             else if ((i_block == j_block + 4) || (j_block == i_block + 4)) {
             fracCovaMatrix_matrix(i,j) *= 0.79;
             }
             else if ((i_block == j_block + 5) || (j_block == i_block + 5)) {
             fracCovaMatrix_matrix(i,j) *= 0.74;
             }
             */
            //fracCovaMatrix_matrix(i,j) *= 0.90;
            //}
        }
    }
    
    cout << "Filling Covariant Matrix... DONE!" << endl;

    //cout << "Done with the matrix!" << endl;
    //fracCovaMatrix_matrix.Print();
    //break;
    
    TFile *wrd_File = new TFile("files_data/daya-bay-ldist6AD.root","READ");
    TH1F *wrd_histo = (TH1F*)(wrd_File->Get("histo_ldist_6Det"));
    
    //-------------------
    // Energy Histograms
    //-------------------
    TFile *fenergy = new TFile("./PRL112_217days_data.root","read");
    //Three sets of histograms one for each Experimental Hall
    double bfactor[nEH];
    for (int i = 0 ; i < nEH ; i++)
    {
        //Data
        data_spect_histo[i] = (TH1F*) fenergy->Get(Form("data_spect_histo_%d",i));
        double dfactor = 1.0/data_spect_histo[i]->Integral();
        data_spect_histo[i]->Scale(dfactor);
        //Background
        bkgd_spect_histo[i] = (TH1F*) fenergy->Get(Form("bkgd_spect_histo_%d",i));
        bfactor[i] = 1.0/bkgd_spect_histo[i]->Integral();
        bkgd_spect_histo[i]->Scale(bfactor[i]);
    }
    
    // histogram binning
    double xbins[27];
    xbins[0] = lo;
    double delta_bins2 = (7.3 - 1.3)/24.;// = 0.25 MeV/bin
    for (int i=0;i<(NB-1);i++)
    {
        xbins[i+1] = 1.3 + delta_bins2*i;
    }
    xbins[26] = hi;
    for(int iAD = 0 ; iAD < nAD ; iAD++)
    {
        nosc_spect_hist[iAD] = new TH1F(Form("nosc_spect_hist_%d",iAD),"",NB,xbins);
        nu_nosc_spect_hist[iAD] = new TH1F(Form("nu_nosc_spect_hist_%d",iAD),"",NB,xbins);
        wosc_spect_hist[iAD] = new TH1F(Form("wosc_spect_hist_%d",iAD),"",NB,xbins);
        nu_wosc_spect_hist[iAD] = new TH1F(Form("nu_wosc_spect_hist_%d",iAD),"",NB,xbins);
        bfit_spect_hist[iAD] = new TH1F(Form("bfit_spect_hist_%d",iAD),"",NB,xbins);
    }
    
    
    for (int blid = 0 ; blid < nAD*nNR ; blid++)
    {
        int id = (blid/nNR);
        int ir = (blid - id*nNR);
        
        wrd_array[id][ir] = wrd_histo->GetBinContent(blid+1);
        //cout << blid << "  " << wrd_array[id][ir] << endl;
    }
    
    cout << "wrd array -> Done..." << endl;
    
    cout << "Minimization settings..." << endl;
    ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer(minName, algoName);
    
    //-- Set tolerance , etc...
    min->SetMaxFunctionCalls(1000000); // for Minuit/Minuit2
    min->SetMaxIterations(10000);      // for GSL
    min->SetTolerance(0.001);
    min->SetPrintLevel(-1);
    
    //-- File to get noOsc normalizations
    ifstream IBDrates_file("files_data/db6AD_noOsc_IBDrates_perday.txt");
    cout << "Reading noOsc normalizations file ..." << endl;
    for (int i=0; i< nAD; i ++)
	{
	    IBDrates_file >> noOsc_IBDrate_perday[i];
      	    cout << "noOscIBD_rates_perday " << i << ": " << noOsc_IBDrate_perday[i] << endl;
    	}//for
 
    //-- File to print oscillation parameters and chi2 values
    ofstream chi2Surface_file;
    string s2t_dm2 = "files_data/chi2_s2t-dm2_surface_SPEC.txt";  //(sin^2(2th13), dm2, chi^2_min)
    chi2Surface_file.open((s2t_dm2).c_str());
    
    //-- File to print minimized pull parameters
    ofstream minimPullT_file;
    string PullT = "files_data/chi2_pullT_surface.txt";
    minimPullT_file.open((PullT).c_str());

    //ifstream file("files_data/db_gridOscSpectra_1M.txt"); //50x50 parameter space-grid
    string grid_spectra = "files_data/db6AD_gridOscSpectra.txt";
    ifstream grid_file;
    grid_file.open((grid_spectra).c_str());
    cout << "Reading file - Loop in progress..." << endl;
    int iad = 0;
    int first6 = 1;
    std::cout << "Is the file open?  " << grid_file.is_open() << "\n" << std::endl;
    
    while (grid_file >> iAD >> s2th_13 >> dm2_31 >> spc[iad][0] >> spc[iad][1] >> spc[iad][2] >> spc[iad][3] >> spc[iad][4] >> spc[iad][5] >> spc[iad][6] >> spc[iad][7] >> spc[iad][8] >> spc[iad][9] >> spc[iad][10] >> spc[iad][11] >> spc[iad][12] >> spc[iad][13] >> spc[iad][14] >> spc[iad][15] >> spc[iad][16] >> spc[iad][17] >> spc[iad][18] >> spc[iad][19] >> spc[iad][20] >> spc[iad][21] >> spc[iad][22] >> spc[iad][23] >> spc[iad][24] >> spc[iad][25] >> NoscTot[iad])
    {//file loop

        if(first6 <= 6) //Must go up to 6 to take the 6 ADs (I think) ---- 2017-08-07
        {//** CHECK HERE the BINING: error in line 286 in spc second index (corrected?)
            for(int ibin = 1 ; ibin <= 26 ; ibin++)
            {
                double bincont = spc[iad][ibin-1]*(noOsc_IBDrate_perday[iAD-1]/NoscTot[iad])*emuem[iAD-1]*daqTime[iAD-1];
                nu_nosc_spect_hist[iAD-1]->SetBinContent(ibin,bincont);
            }
            cout << "iad = " << iad << "   iAD = " << iAD  << "   first6 = " << first6 << endl; // ---- 2017-08-07
            //nosc_spect_hist[iAD-1]->Print("all"); // ---- 2017-08-07
            cout << "Number of events: " << nu_nosc_spect_hist[iAD-1]->Integral() << endl; // ---- 2017-08-07
            //cout << "NoscTot         : " << NoscTot[iad] << endl; // ---- 2017-08-07
            noNoscTot[iad] = NoscTot[iad];
            first6++;
            //continue; // - Commenting out ---- 2017-08-07
        }//if first6 loop END
        
        //cout << iad+1 << ": " << iAD << " " << s2th_13 << " " << dm2_31 << " " ;
        //for (int ib = 0; ib < NB; ib++) {
        //  cout << spc[iad][ib] << " ";
        //}
        //cout << NoscTot[iad] << " " << sp << endl;
        
        iad++;
        /*
         cout << "*************************************" << endl;
         for(int ibin = 1 ; ibin <= 26 ; ibin++)
         {
         cout << nosc_spect_hist[iAD-1]->GetBinContent(ibin) << "  ";
         }
         cout << "\n*************************************" << endl;
         //break;
         */
        
        
        //At this point, we have read the six AD spectra (six lines) for one point (s2th_13,dm2_31) in the grid
        if (iad == 6)
        {//if iad BEGIN
            //cout << endl;
            iad = 0;
            
            //This is executed for every 6 lines of the file
            const int N_params = 19; //-- Number of parameter of the chi² function --//
            ROOT::Math::Functor f(&chi2,N_params); //-- Setting the function to be minimized by using Minuit --//
            //-- Steps
            double stp = 1.0e-3;
            double step[N_params] = {stp,stp,stp,stp,stp,stp,
                stp,stp,stp,stp,stp,stp,
                stp,stp,stp,stp,stp,stp,stp};
            //-- Initial parameter values
            double start[N_params] = {0.0,0.0,0.0,0.0,0.0,0.0,
                0.0,0.0,0.0,0.0,0.0,0.0,
                0.0,0.0,0.0,0.0,0.0,0.0,0.0};
            
            //-- Calling Minuit function setting
            min->SetFunction(f);
            
            //-- Setting variables
            double lim = 1.0e-1;
            /*min->SetLimitedVariable(0,  "e_1", start[0],  step[0],  -lim, lim);
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
             min->SetLimitedVariable(17, "a_6", start[17], step[17], -lim, lim);*/
            min->SetLimitedVariable(18, "eps", start[18], step[18], -lim, lim);
            //min->SetFixedVariable(18, "eps", start[18]);
            min->SetFixedVariable(0,  "e_1", start[0]);
            min->SetFixedVariable(1,  "e_2", start[1]);
            min->SetFixedVariable(2,  "e_3", start[2]);
            min->SetFixedVariable(3,  "e_4", start[3]);
            min->SetFixedVariable(4,  "e_5", start[4]);
            min->SetFixedVariable(5,  "e_6", start[5]);
            min->SetFixedVariable(6,  "n_1", start[6]);
            min->SetFixedVariable(7,  "n_2", start[7]);
            min->SetFixedVariable(8,  "n_3", start[8]);
            min->SetFixedVariable(9,  "n_4", start[9]);
            min->SetFixedVariable(10, "n_5", start[10]);
            min->SetFixedVariable(11, "n_6", start[11]);
            min->SetFixedVariable(12, "a_1", start[12]);
            min->SetFixedVariable(13, "a_2", start[13]);
            min->SetFixedVariable(14, "a_3", start[14]);
            min->SetFixedVariable(15, "a_4", start[15]);
            min->SetFixedVariable(16, "a_5", start[16]);
            min->SetFixedVariable(17, "a_6", start[17]);
            //min->SetFixedVariable(18, "eps", start[18]);
            min->SetErrorDef(2.3);
            
            //-- Calling Minuit minimization
            
            min->Minimize();
            
            const double *xs = min->X();
            
            double chi2Min = min->MinValue();
            chi2Surface_file << s2th_13 << "\t" << dm2_31 << "\t" << chi2Min << endl;
            
            if (chi2Min < 0.0) {
                cout << "Critical error: chi2Min is negative!  " << chi2Min << endl;
                break;
            }
            
            minimPullT_file  << s2th_13 << "\t" << dm2_31 << "\t" << xs[0] << "\t" << xs[1] << "\t" << xs[2] << "\t" << xs[3] << "\t" << xs[4] << "\t" << xs[5] << "\t" << xs[6] << "\t" << xs[7] << "\t" << xs[8] << "\t" << xs[9] << "\t" << xs[10] << "\t" << xs[11] << "\t" << xs[12] << "\t" << xs[13] << "\t" << xs[14] << "\t" << xs[15] << "\t" << xs[16] << "\t" << xs[17] << "\t" << xs[18] << "\t" << min->MinValue() << endl;
            //}
            if (dm2_31 == hi_dm2)
            {
                chi2Surface_file << endl;
                cout << s2th_13 << "\t" << dm2_31 << "\t" << min->MinValue() << endl;
            }
            //minimPullT_file  << endl;
            //if (is2t%10 == 0)
            //std::cout << "Succesful run for sin^2(th13) = " << s2th_13 << "!! \t" << min->MinValue() << endl;
            //}
        }//if iad END
        
    }//file loop END
    std::cout << "Successful run!!" << endl;
    minimPullT_file.close();
    chi2Surface_file.close();

    //break;
    //------------------------------------------------------------------------------
    //------------------------------------------------------------------------------
    //Drawing no-oscillated spectra
    TCanvas *canv0 = new TCanvas("canv0","",3*700,2*350);
    canv0->Divide(3,2);
    for(int iAD = 0 ; iAD < 6 ; iAD++)
    {
        canv0->cd(iAD+1);
        nosc_spect_hist[iAD]->Draw("");
    }
    //canv0->Print("files_plots/canv0.pdf");
    
    TH1F *data_spect_EHhisto[nEH];
    //-- ADs 0, 1 -> EH 0
    data_spect_EHhisto[0] = (TH1F*)data_spect_histo[0]->Clone();
    data_spect_EHhisto[0]->Scale(IBDrate_data[0][0]*daqTime[0]);
    data_spect_EHhisto[0]->Add(data_spect_histo[0],IBDrate_data[1][0]*daqTime[1]);
    //-- AD 2 -> EH 1
    data_spect_EHhisto[1] = (TH1F*)data_spect_histo[1]->Clone();
    data_spect_EHhisto[1]->Scale(IBDrate_data[2][0]*daqTime[2]);
    //-- ADs 3, 4, 5 -> EH 2
    data_spect_EHhisto[2] = (TH1F*)data_spect_histo[2]->Clone();
    data_spect_EHhisto[2]->Scale(IBDrate_data[3][0]*daqTime[3]);
    data_spect_EHhisto[2]->Add(data_spect_histo[2],IBDrate_data[4][0]*daqTime[4]);
    data_spect_EHhisto[2]->Add(data_spect_histo[2],IBDrate_data[5][0]*daqTime[5]);
    
    
    TH1F *nosc_spect_EHhist[nEH];
    //-- ADs 0, 1 -> EH 0
    nosc_spect_EHhist[0] = (TH1F*)nu_nosc_spect_hist[0]->Clone(); // Nus from AD0
    nosc_spect_EHhist[0]->Add(nu_nosc_spect_hist[1],1); // Nus from AD1
    nosc_spect_EHhist[0]->Add(bkgd_spect_histo[0],1.0/bfactor[0]); //Background from EH0
    //-- AD 2 -> EH 1
    nosc_spect_EHhist[1] = (TH1F*)nu_nosc_spect_hist[2]->Clone(); //Nus from AD2
    nosc_spect_EHhist[1]->Add(bkgd_spect_histo[1],1.0/bfactor[1]); //Background from EH1
    //-- ADs 3, 4, 5 -> EH 2
    nosc_spect_EHhist[2] = (TH1F*)nu_nosc_spect_hist[3]->Clone(); //Nus from AD3
    nosc_spect_EHhist[2]->Add(nu_nosc_spect_hist[4],1); //Nus from AD4
    nosc_spect_EHhist[2]->Add(nu_nosc_spect_hist[5],1); //Nus from AD5
    nosc_spect_EHhist[2]->Add(bkgd_spect_histo[2],1.0/bfactor[2]); //Background from EH2
    
    TCanvas *canv1 = new TCanvas("canv1","Events",3*700,1*350);
    canv1->Divide(3,1);
    for(int iEH = 0 ; iEH < nEH ; iEH++)
    {
        canv1->cd(iEH+1);
        nosc_spect_EHhist[iEH]->Draw("h");
        data_spect_EHhisto[iEH]->Draw("p same");
    }
    //canv1->Print("files_plots/canv1.pdf");
    
    //---------------------------------------
    TCanvas *canv2 = new TCanvas("canv2","Events/MeV",3*700,1*350);
    canv2->Divide(3,1);
    TH1F *nosc_spect_EHhist_EvtperMeV[nEH];
    for(int iEH = 0 ; iEH < nEH ; iEH++)
    {
        nosc_spect_EHhist_EvtperMeV[iEH] = new TH1F(Form("nosc_spect_EHhist_EvtperMeV_%d",iEH),"",NB,xbins);
        for(int ibin = 0 ; ibin < NB ; ibin++)
        {
            double binw = nosc_spect_EHhist[iEH]->GetBinWidth(ibin+1);
            double binc = nosc_spect_EHhist[iEH]->GetBinContent(ibin+1);
            //cout << iEH << "  binw = " << binw << "\t binc = " << binc << endl;
            nosc_spect_EHhist_EvtperMeV[iEH]->SetBinContent(ibin+1,binc/binw);
        }
        canv2->cd(iEH+1);
        nosc_spect_EHhist_EvtperMeV[iEH]->Draw("h");
    }
    //canv2->Print("files_plots/canv2.pdf");
    
    //------------------------------------------------------------------------------
    //------------------------------------------------------------------------------
    return 0;
}

