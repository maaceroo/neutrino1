//----------------------------------------------------------//
//-- db_minuit_spec.C  -  By M.A. Acero O.  -  2018-08-18 --//
//----------------------------------------------------------//
//-- Using 'NumericalMinimization.C' macro example to     --//
//-- minimize a chi2 function. It uses a Minimizer class  --//
//-- in ROOT.                                             --//
//-- input : minimizer name + algorithm name              --//
//-- randomSeed: = <0 : fixed value: 0 random with seed 0;--//
//-- >0 random with given seed.                           --//
//----------------------------------------------------------//
//-- Analysis of the Daya Bay data from F.P. An et al.,   --//
//-- PRD 95 072006 (2017)                                 --//
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
//#define N_s2t  5                           //number of points in the grid
//#define lo_s2t 0.01                         //sin^2(2th_13) min
//#define lo_s2t 0.0                          //sin^2(2th_13) min
//#define hi_s2t 0.3                          //sin^2(2th_13) max
//For the delta(m31)^2 loop
//#define N_dm2  5                           //number of points in the grid
//#define lo_dm2 1.0e-4                       //delta(m31)^2 min
//#define hi_dm2 1.0e-2                       //delta(m31)^2 max
//---*****************************************************---//

//---*****************************************************---//
//-- Global variables and quantities ------------------------//
//---*****************************************************---//
//EH1(AD1, AD2),EH2(AD3),EH3(AD4, AD5, AD6)
//(IBD candidates)/(DAQ live time -days-) from PRD 95 072006 (2017)
double daqTime_Total[nAD]         = {1117.178,1117.178,1114.337,924.933,1106.915,1106.915,1106.915,917.417};
double IBDrate_data_Total[nAD][2] = { {534.93,0.69},{542.75,0.70},{509.00,0.68},{503.83,0.74},{ 72.71,0.26},{ 72.94,0.26},{72.33,0.26},{72.88,0.28} };
//IBD rate (per day), total background (8AD) and efficiencies (PRD 95 072006 (2017))
double totalBgd_Total[nAD][2]     = { {11.94,1.07},{11.94,1.07},{ 8.76,0.78},{ 8.69,0.78},{ 1.56,0.07},{ 1.47,0.07},{ 1.48,0.07},{ 1.28,0.07} };
double emuem_Total[nAD]           = {0.8044,0.8013,0.8365,0.8363,0.9587,0.9585,0.9581,0.9588};
// Information obtained by executing the script "db_osc_rate.C"
//IBD rate per day w/o oscillations
double noOsc_IBDrate_perday[nAD]; // = {663.538, 675.961, 610.675, 604.832, 79.5001, 79.8741, 79.2203, 80.021};
//---*****************************************************---//
//-- Information for 6AD
double daqTime_6AD[nAD]           = {191.001, 191.001, 189.645, 0,0,189.779,189.779,189.779,0.0};
double IBDrate_data_6AD[nAD][2]   = { {530.31,1.67},{536.75,1.68},{489.93,1.61},{0,0},        { 73.58,0.62},{ 73.21,0.62},{72.35,0.62},{0,0} };
double totalBgd_6AD[nAD][2]       = { {13.20,0.98}, { 13.01,0.98},{  9.57,0.71},{0,0},        {  3.52,0.14},{ 3.48,0.14},{ 3.43,0.14},{0,0}        };
double emuem_6AD[nAD]             = {0.7957,0.7927,0.8282,0,0,0.9577,0.9568,0.9566,0.0};
//IBD rate per day w/o oscillations
double noOsc_IBDrate_perday_6AD[nAD]//; = {663.15,673.95,591.86,78.75,78.46,77.58}; //-- Modify the macro "db_osc_rate.C"
//---*****************************************************---//
//---*****************************************************---//
//-- Information for 8AD
double daqTime_8AD[nAD];         //--DAQTimeTotal - DAQTime6AD
double IBDrate_data_8AD[nAD][2]; //-- (IBDTotalPerDay*DAQTimeTotal - IBD6ADPerDay*DAQTime6AD) / DAQTime8AD
double totalBgd_8AD[nAD][2];     //-- (TotalPerDay*DAQTimeTotal - 6ADPerDay*DAQTime6AD) / DAQTime8AD
double emuem_8AD[nAD];           //-- = emuem_Total[nAD];
//IBD rate per day w/o oscillations
double noOsc_IBDrate_perday_8AD[nAD]; // Information obtained by executing the script "db_osc_rate.C"
//---*****************************************************---//
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
double e;           //01.03.2019 energy scale parameter
double NoscTot[nAD];
double noNoscTot[nAD];
double xbins[NB+1];
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
int NBx_cov;
int NBy_cov;
TMatrixD *statErroMatrix_matrix;
TMatrixD *fracCovaMatrix_matrix;  //Se llena en la función principal db_minuit_spec()
TMatrixD *fullCovaMatrix_matrix;
TMatrixD *inv_fullCovaMatrix_matrix;
TMatrixD *one_matrix;
TMatrixD *delta_vector;
TMatrixD *transp_delta_vector;
TMatrixD *predi_vector;

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
    eps_d[6] = xx[6];
    eps_d[7] = xx[7];
    eta_d[0] = xx[8];
    eta_d[1] = xx[9];
    eta_d[2] = xx[10];
    eta_d[3] = xx[11];
    eta_d[4] = xx[12];
    eta_d[5] = xx[13];
    eta_d[6] = xx[14];
    eta_d[7] = xx[15];
    alpha[0] = xx[16];
    alpha[1] = xx[17];
    alpha[2] = xx[18];
    alpha[3] = xx[19];
    alpha[4] = xx[20];
    alpha[5] = xx[21];
    epsilon  = xx[22];
    e        = xx[23]; //-- energy scale pull term
  
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
    double sesc         = 0.005; //-- energy scale pull term

    //- 01.03.2019 -Begins - This part allows to include properly the energy scale pull term
    double f_ePos[NB] = {0.0};
    double f_eNeg[NB] = {0.0};
    double spcNew[nAD][NB];
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
    //- 01.03.2019 -Ends

    for (iAD = 0 ; iAD < nAD ; iAD++){
        SurvPavg = NoscTot[iAD]/noNoscTot[iAD];
        for (int iBIN = 0 ; iBIN < NB ; iBIN++)
        {
            //- 01.03.2019 -Begins----------------------------------------------------
            //- This part allows to include properly the energy scale pull term
            if (e < 0) {
                if (iBIN == 0){
                    spcNew[iAD][iBIN] = spc[iAD][iBIN] + f_eNeg[iBIN+1]*spc[iAD][iBIN+1];
                    //spcNew[1][iBIN] = spc[1][iBIN] + f_eNeg[iBIN+1]*spc[1][iBIN+1];
                }
                else if (iBIN<NB-1){
                    spcNew[iAD][iBIN] = spc[iAD][iBIN] - f_eNeg[iBIN]*spc[iAD][iBIN] + f_eNeg[iBIN+1]*spc[iAD][iBIN+1];
                    //spcNew[1][iBIN] = spc[1][iBIN] - f_eNeg[iBIN]*spc[1][iBIN] + f_eNeg[iBIN+1]*spc[1][iBIN+1];
                }
                else{
                    spcNew[iAD][iBIN] = spc[iAD][iBIN] - f_eNeg[iBIN]*spc[iAD][iBIN];
                    //spcNew[1][iBIN] = spc[1][iBIN] - f_eNeg[iBIN]*spc[1][iBIN];
                }
            }
            
            if (e >= 0) {
                if (iBIN == 0){
                    spcNew[iAD][iBIN] = spc[iAD][iBIN] - f_ePos[iBIN+1]*spc[iAD][iBIN];
                    //spcNew[1][iBIN] = spc[1][iBIN] - f_ePos[iBIN+1]*spc[1][iBIN];
                }
                else if (iBIN<NB-1){
                    spcNew[iAD][iBIN] = spc[iAD][iBIN] + f_ePos[iBIN]*spc[iAD][iBIN-1] - f_ePos[iBIN+1]*spc[iAD][iBIN];
                    //spcNew[1][iBIN] = spc[1][iBIN] + f_ePos[iBIN]*spc[1][iBIN-1] - f_ePos[iBIN+1]*spc[1][iBIN];
                }
                else{
                    spcNew[iAD][iBIN] = spc[iAD][iBIN] + f_ePos[iBIN]*spc[iAD][iBIN-1];
                    //spcNew[1][iBIN] = spc[1][iBIN] + f_ePos[iBIN]*spc[1][iBIN-1];
                }
            }
            //- 01.03.2019 -Ends----------------------------------------------------

            int index = iAD*NB + iBIN;
            //-- Survival Probability equation. Terms depending on dM²_21 and dM²_32(~dM²_31) are averaged
            //SurvP = 1.0 - 0.25*pow((1.0 + sqrt(1.0 - s2th_13)),2)*(s22th_12)*(avgSinDelta21[iAD]) - s2th_13*(avgSinDelta31[iAD]);
        
            //-- Predicted IBD from neutrino oscillations of the dth Antineutrino Detector
            //Td = (SurvP * noOsc_IBDrate_perday[iAD])*emuem[iAD]*daqTime[iAD];
            //Td = spc[iAD][iBIN]*(SurvPavg*noOsc_IBDrate_perday[iAD]/NoscTot[iAD])*emuem[iAD]*daqTime[iAD];
            Td = spcNew[iAD][iBIN]*(SurvPavg*noOsc_IBDrate_perday[iAD]/NoscTot[iAD])*emuem[iAD]*daqTime[iAD];
            //-- Measured IDB events of the dth Antineutrino Detector (background is substracted)
            //Md = (IBDrate_data[iAD][0] - totalBgd[iAD][0]*emuem[iAD])*daqTime[iAD];
            int idx = -1;
            if (iAD < 2) {
                idx = 0;
            }
            else if (iAD == 2 || iAD == 3) {
                idx = 1;
            } else {
                idx = 2;
            }
            double data_sigplusbg = (data_spect_histo[idx]->GetBinContent(iBIN+1))*IBDrate_data[iAD][0];
            double simu_bg = (bkgd_spect_histo[idx]->GetBinContent(iBIN+1))*totalBgd[iAD][0]*emuem[iAD];
            Md = (data_sigplusbg - simu_bg)*daqTime[iAD];
            //cout << "iAD = " << iAD << "  iBin = " << iBIN << "  Md = " << Md << " Td = " << Td << endl;
            //-- Background of the dth Antineutrino Detector
            Bd = simu_bg*daqTime[iAD];
	
            sqrerror = Md + Bd;
            
            //Statistics Error matrix
            statErroMatrix_matrix(index,index) = sqrerror;
	
            //-- Fraction of IBD contribution of the rth reactor to the dth AD
            //-- determined by baselines and reactor fluxes
            double wrd = 0.0;
            for (iNR = 0 ; iNR < nNR ; iNR++)
                wrd += wrd_array[iAD][iNR]*alpha[iNR];
        
            predi_vector(index,0) = Td*(1.0 + epsilon + eps_d[iAD] + wrd) - eta_d[iAD];
            //delta_vector(index,0) = Md - Td*(1.0 + epsilon + eps_d[iAD] + wrd) + eta_d[iAD];
            delta_vector(index,0) = Md - Td*(1.0 + epsilon + eps_d[iAD] + wrd) + eta_d[iAD];
            transp_delta_vector(0,index) = delta_vector(index,0);
            
            //cout << "iAD = " << iAD << "  iBin = " << iBIN << "  Md = " << Md << " Td = " << Td
                 //<< " delta = " << delta_vector(index,0) << endl;
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
            fullCovaMatrix_matrix(i,j) = statErroMatrix_matrix(i,j);
            //+ 0.0*predi_vector(i,0)*predi_vector(j,0)*fracCovaMatrix_matrix(i,j); // ** test A.A. 7/sep/18
            //+ predi_vector(i,0)*predi_vector(j,0)*fracCovaMatrix_matrix(i,j);
        }
    }
    //predi_vector->Print();
    //fullCovaMatrix_matrix.Print();

    if(fullCovaMatrix_matrix->Determinant() != 0){
        inv_fullCovaMatrix_matrix = fullCovaMatrix_matrix;
        inv_fullCovaMatrix_matrix->Invert();
    }
    else {
        cout << "Determinant = " << fullCovaMatrix_matrix->Determinant() << endl;
        cout << "Singular fullCovaMatrix!" << endl;
        cout << "Execution aborted!" << endl;
        break;
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

   for (iAD = 0 ; iAD < nAD ; iAD++)
       {
           //-- Background error of the dth Antineutrino Detector
           sB = totalBgd[iAD][1]*emuem[iAD]*daqTime[iAD];
           //sB = totalBgd[iAD][1]/totalBgd[iAD][0];
           sqr_chi += pow(eps_d[iAD]/seps_d,2) + pow(eta_d[iAD]/sB,2);
       }
    
    for (iNR = 0 ; iNR < nNR ; iNR++)
         sqr_chi += pow(alpha[iNR]/salph_r,2);

    sqr_chi += pow(e/sesc,2);
    
    //cout << "AD " << iAD << "  chi2 = " << sqr_chi << endl;
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
    TFile *fracCovaMatrix_File = new TFile("CovaMatrix/db_CovaMatrix_8AD_8x35bins.root","READ");
    fracCovaMatrix_hist = (TH2F*)(fracCovaMatrix_File->Get("covaMat_histo_8x8Det"));
    //-- Creating the Rebinned Matrix with root tools.
    NBx_cov = NB*nAD;
    NBy_cov = NB*nAD;
    cout << NBx_cov << "  " << NBy_cov << endl;
    statErroMatrix_matrix = new TMatrixD(NBx_cov,NBy_cov);
    fracCovaMatrix_matrix = new TMatrixD(NBx_cov,NBy_cov);
    fullCovaMatrix_matrix = new TMatrixD(NBx_cov,NBy_cov);
    inv_fullCovaMatrix_matrix = new TMatrixD(NBx_cov,NBy_cov);
    
    one_matrix = new TMatrixD(NBx_cov,NBy_cov);
    
    delta_vector        = new TMatrixD(NBx_cov,1);
    transp_delta_vector = new TMatrixD(1,NBx_cov);
    predi_vector        = new TMatrixD(NBx_cov,1);
    
    fracCovaMatrix_matrix->Zero();
    statErroMatrix_matrix->Zero();
    fullCovaMatrix_matrix->Zero();
    inv_fullCovaMatrix_matrix->Zero();
    one_matrix->Zero();
    delta_vector->Zero();
    transp_delta_vector->Zero();
    predi_vector->Zero();
    
    int i_block = 0;
    int j_block = 0;
    for (int i = 0 ; i < NBx_cov ; i++) {
        for (int j = 0 ; j < NBy_cov; j++) {
            i_block = i%35;
            j_block = j%35;
            fracCovaMatrix_matrix(i,j) = fracCovaMatrix_hist->GetBinContent(i+1,j+1);

        }
    }


    //cout << "Done with the matrix!" << endl;
    //fracCovaMatrix_matrix->Print();
    //break;

    TFile *wrd_File = new TFile("files_data/daya-bay-ldist.root","READ");
    TH1F *wrd_histo = (TH1F*)(wrd_File->Get("histo_ldist"));

    //-------------------
    // Energy Histograms
    //-------------------
    TFile *fenergy = new TFile("./PRD95_1230days_data.root","read");
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
    //double xbins[36];
    xbins[0] = 0.7;
    double delta_bins2 = (7.9 - 1.3)/33; // 0.20 MeV/bin
    for (int i = 0 ; i < (NB-1) ; i++)
    {
        xbins[i+1] = 1.3 + delta_bins2*i;
    }
    xbins[35] = hi;
    for(int iAD = 0 ; iAD < nAD ; iAD++)
    {
        nosc_spect_hist[iAD]    = new TH1F(Form("nosc_spect_hist_%d",iAD),   "",NB,xbins);
        nu_nosc_spect_hist[iAD] = new TH1F(Form("nu_nosc_spect_hist_%d",iAD),"",NB,xbins);
        wosc_spect_hist[iAD]    = new TH1F(Form("wosc_spect_hist_%d",iAD),   "",NB,xbins);
        nu_wosc_spect_hist[iAD] = new TH1F(Form("nu_wosc_spect_hist_%d",iAD),"",NB,xbins);
        bfit_spect_hist[iAD]    = new TH1F(Form("bfit_spect_hist_%d",iAD),   "",NB,xbins);
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
    ifstream IBDrates_file("files_data/db_noOsc_IBDrates_perday.txt");
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
    ifstream file("files_data/db_gridOscSpectra.txt");
    cout << "Reading file - Loop in progress..." << endl;
    int iad = 0;
    int first8 = 1;
    
    while (file >> iAD >> s2th_13 >> dm2_31 >>
           spc[iad][0]  >> spc[iad][1]  >> spc[iad][2]  >> spc[iad][3]  >> spc[iad][4]  >> spc[iad][5]  >>
           spc[iad][6]  >> spc[iad][7]  >> spc[iad][8]  >> spc[iad][9]  >> spc[iad][10] >> spc[iad][11] >>
           spc[iad][12] >> spc[iad][13] >> spc[iad][14] >> spc[iad][15] >> spc[iad][16] >> spc[iad][17] >>
           spc[iad][18] >> spc[iad][19] >> spc[iad][20] >> spc[iad][21] >> spc[iad][22] >> spc[iad][23] >>
           spc[iad][24] >> spc[iad][25] >> spc[iad][26] >> spc[iad][27] >> spc[iad][28] >> spc[iad][29] >>
           spc[iad][30] >> spc[iad][31] >> spc[iad][32] >> spc[iad][33] >> spc[iad][34] >>
           NoscTot[iad])
        {//file loop
            cout << s2th_13 << "\t" << dm2_31 << endl;
            if(first8 <= 8) //Must go up to 8 to take the 8 ADs (I think) ---- 2017-08-07
            {
                for(int ibin = 1 ; ibin <= NB ; ibin++)
                {
                    double bincont = spc[iad][ibin-1]*(noOsc_IBDrate_perday[iAD-1]/NoscTot[iad])*emuem[iAD-1]*daqTime[iAD-1];
                    nu_nosc_spect_hist[iAD-1]->SetBinContent(ibin,bincont);
                }
                cout << "iad = " << iad << "   iAD = " << iAD  << "   first8 = " << first8 << endl; // ---- 2017-08-07
                //nosc_spect_hist[iAD-1]->Print("all"); // ---- 2017-08-07
                cout << "Number of events: " << nu_nosc_spect_hist[iAD-1]->Integral() << endl; // ---- 2017-08-07
                //cout << "NoscTot         : " << NoscTot[iad] << endl; // ---- 2017-08-07
                noNoscTot[iad] = NoscTot[iad];
                first8++;
                //continue; // - Commenting out ---- 2017-08-07
            }//if first8 loop END
            
            //cout << iad+1 << ": " << iAD << " " << s2th_13 << " " << dm2_31 << " " ;
            //for (int ib = 0; ib < NB; ib++) {
              //  cout << spc[iad][ib] << " ";
            //}
            //cout << NoscTot[iad] << " " << sp << endl;

            iad++;
/*
            cout << "*************************************" << endl;
            for(int ibin = 1 ; ibin <= NB ; ibin++)
            {
                cout << nosc_spect_hist[iAD-1]->GetBinContent(ibin) << "  ";
            }
            cout << "\n*************************************" << endl;
            //break;
*/
            
            
            //At this point, we have read the 8 AD spectra (8 lines) for one point (s2th_13,dm2_31) in the grid
            if (iad == 8)
            {//if iad BEGIN
                //cout << endl;
                iad = 0;
                
                //This is executed for every 6 lines of the file
                //const int N_params = 23; //-- Number of parameter of the chi² function --//
                const int N_params = 24; //-- Number of parameter of the chi² function --//
                ROOT::Math::Functor f(&chi2,N_params); //-- Setting the function to be minimized by using Minuit --//
                //-- Steps
                double stp = 1.0e-3;
                double step[N_params] = {stp,stp,stp,stp,stp,stp,stp,stp,
                                         stp,stp,stp,stp,stp,stp,stp,stp,
                                         stp,stp,stp,stp,stp,stp,stp,stp};
                //-- Initial parameter values
                double start[N_params] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
                                          0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
                                          0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
                
                //-- Calling Minuit function setting
                min->SetFunction(f);
                
                //-- Setting variables
                double lim = 1.0e-1;
                
		min->SetLimitedVariable(0,  "e_1", start[0],  step[0],  -lim, lim);
        min->SetLimitedVariable(1,  "e_2", start[1],  step[1],  -lim, lim);
        min->SetLimitedVariable(2,  "e_3", start[2],  step[2],  -lim, lim);
        min->SetLimitedVariable(3,  "e_4", start[3],  step[3],  -lim, lim);
        min->SetLimitedVariable(4,  "e_5", start[4],  step[4],  -lim, lim);
        min->SetLimitedVariable(5,  "e_6", start[5],  step[5],  -lim, lim);
        min->SetLimitedVariable(6,  "e_7", start[6],  step[6],  -lim, lim);
        min->SetLimitedVariable(7,  "e_8", start[7],  step[7],  -lim, lim);
		
                
		min->SetLimitedVariable(8,  "n_1", start[8],  step[8],  -lim, lim);
        min->SetLimitedVariable(9,  "n_2", start[9],  step[9],  -lim, lim);
        min->SetLimitedVariable(10, "n_3", start[10], step[10], -lim, lim);
        min->SetLimitedVariable(11, "n_4", start[11], step[11], -lim, lim);
        min->SetLimitedVariable(12, "n_5", start[12], step[12], -lim, lim);
        min->SetLimitedVariable(13, "n_6", start[13], step[13], -lim, lim);
        min->SetLimitedVariable(14, "n_7", start[14], step[14], -lim, lim);
        min->SetLimitedVariable(15, "n_8", start[15], step[15], -lim, lim);
		
                
		min->SetLimitedVariable(16, "a_1", start[16], step[16], -lim, lim);
        min->SetLimitedVariable(17, "a_2", start[17], step[17], -lim, lim);
        min->SetLimitedVariable(18, "a_3", start[18], step[18], -lim, lim);
        min->SetLimitedVariable(19, "a_4", start[19], step[19], -lim, lim);
        min->SetLimitedVariable(20, "a_5", start[20], step[20], -lim, lim);
        min->SetLimitedVariable(21, "a_6", start[21], step[21], -lim, lim);
		
        min->SetLimitedVariable(22, "eps", start[22], step[22], -lim, lim);
        
        min->SetLimitedVariable(23, "e",   start[23], step[23], -lim, lim);
		//-----------------------------------------------------------------//
		/*
        min->SetFixedVariable(0,  "e_1", start[0]);
        min->SetFixedVariable(1,  "e_2", start[1]);
        min->SetFixedVariable(2,  "e_3", start[2]);
        min->SetFixedVariable(3,  "e_4", start[3]);
        min->SetFixedVariable(4,  "e_5", start[4]);
        min->SetFixedVariable(5,  "e_6", start[5]);
        min->SetFixedVariable(6,  "e_7", start[6]);
        min->SetFixedVariable(7,  "e_8", start[7]);
		*/
        /*
		min->SetFixedVariable(8,  "n_1", start[8]);
        min->SetFixedVariable(9,  "n_2", start[9]);
        min->SetFixedVariable(10, "n_3", start[10]);
        min->SetFixedVariable(11, "n_4", start[11]);
        min->SetFixedVariable(12, "n_5", start[12]);
        min->SetFixedVariable(13, "n_6", start[13]);
        min->SetFixedVariable(14, "n_7", start[14]);
        min->SetFixedVariable(15, "n_8", start[15]);
		*/
		/*
        min->SetFixedVariable(16, "a_1", start[16]);
        min->SetFixedVariable(17, "a_2", start[17]);
        min->SetFixedVariable(18, "a_3", start[18]);
        min->SetFixedVariable(19, "a_4", start[19]);
        min->SetFixedVariable(20, "a_5", start[20]);
        min->SetFixedVariable(21, "a_6", start[21]);
		*/
        //min->SetFixedVariable(22, "eps", start[22]);
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
                
        minimPullT_file  << s2th_13 << "\t" << dm2_31 << "\t" << xs[0]  << "\t" << xs[1]
                << "\t" << xs[2]  << "\t" << xs[3]  << "\t" << xs[4]  << "\t" << xs[5]
                << "\t" << xs[6]  << "\t" << xs[7]  << "\t" << xs[8]  << "\t" << xs[9]
                << "\t" << xs[10] << "\t" << xs[11] << "\t" << xs[12] << "\t" << xs[13]
                << "\t" << xs[14] << "\t" << xs[15] << "\t" << xs[16] << "\t" << xs[17]
                << "\t" << xs[18] << "\t" << xs[19] << "\t" << xs[20] << "\t" << xs[21]
                << "\t" << xs[22] << "\t" << xs[23] << "\t" << min->MinValue() << endl;
                //}
        if (dm2_31 == hi_dm2)
        {
            chi2Surface_file << endl;
            cout << s2th_13 << "\t" << dm2_31 << "\t" << min->MinValue() << endl;
        }
        cout << s2th_13 << "\t" << dm2_31 << "\t" << min->MinValue() << endl;
        //minimPullT_file  << endl;
        //if (is2t%10 == 0)
        //std::cout << "Succesful run for sin^2(th13) = " << s2th_13 << "!! \t" << min->MinValue() << endl;
        //}
        }//if iad END

    }//file loop END
    std::cout << "Successful run!!" << endl;
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
