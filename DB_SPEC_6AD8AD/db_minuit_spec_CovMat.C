    //--------------------------------------------------------------------------------------//
//--   db_minuit_spec_CovMat.C  -  By M.A. Acero O. & A.A. Aquilar-A. -  2019-07-05   --//
//--------------------------------------------------------------------------------------//
//--     Analysis of the Daya Bay data from F.P. An et al., PRD 95 072006 (2017)      --//
//--------------------------------------------------------------------------------------//
//-- This macro can be executed under ROOT typing "root[0] .x db_minuit_specCovMat.C" --//
//-- In the simplest case (only minimization) the script will print the message       --//
//-- "Minimum: f(x[i]): chi^2(min)" where x[i] are the values of the parameters       --//
//-- which minimize the chi^2-function and chi^2(min) is the minimum value of the     --//
//-- chi^2-function. In the complete case, the results are printed in a  file, with   --//
//-- the following information:        sin^2(2th_13)  delta(m31)^2  chi^2(min)        --//
//-- In this case, chi^2(min) is the minimum chi^2 for  the pull terms.               --//
//---********************************************************************************---//
//--                        FOR THE SPECTRAL ANALYSIS ONLY                            --//
//---********************************************************************************---//
//---********************************************************************************---//
//---- CONSTANTS ----//
#include "constants.h"
//---*****************************************************---//
//-- Global variables and quantities ------------------------//
//---*****************************************************---//
//EH1(AD1, AD2),EH2(AD3),EH3(AD4, AD5, AD6)
//(IBD candidates)/(DAQ live time -days-) from PRD 95 072006 (2017)
double daqTime_Total[nAD]         = {1117.178, 1117.178, 1114.337, 924.933, 1106.915, 1106.915, 1106.915, 917.417};
double IBDrate_data_Total[nAD][2] = { {534.93,0.69},{542.75,0.70},{509.00,0.68},{503.83,0.74},{ 72.71,0.26},{ 72.94,0.26},{72.33,0.26},{72.88,0.28} };
//IBD rate (per day), total background (8AD) and efficiencies (PRD 95 072006 (2017))
double totalBgd_Total[nAD][2]     = { {11.94,1.07},{11.94,1.07},{ 8.76,0.78},{ 8.69,0.78},{ 1.56,0.07},{ 1.47,0.07},{ 1.48,0.07},{ 1.28,0.07} };
double emuem_Total[nAD]           = {0.8044,0.8013,0.8365,0.8363,0.9587,0.9585,0.9581,0.9588};
double noOsc_IBDrate_perday_6AD8AD[nAD]; // Information obtained by executing the script "db_osc_spec.C"
//---*****************************************************---//
//-- Information for 6AD
double daqTime_6AD[nAD]           = {191.001, 191.001, 189.645, 0.0, 189.779, 189.779, 189.779, 0.0};
double IBDrate_data_6AD[nAD][2]   = { {530.31,1.67},{536.75,1.68},{489.93,1.61},{0,0},        { 73.58,0.62},{ 73.21,0.62},{72.35,0.62},{0,0} };
double totalBgd_6AD[nAD][2]       = { {13.20,0.98}, { 13.01,0.98},{  9.57,0.71},{0,0},        {  3.52,0.14},{ 3.48,0.14},{ 3.43,0.14},{0,0}        };
double emuem_6AD[nAD]             = {0.7957,0.7927,0.8282,0.0,0.9577,0.9568,0.9566,0.0};
//IBD rate per day w/o oscillations
//Need to claculate from the ntuple for the 217-days Period - TO BE DONE!
double noOsc_IBDrate_perday_6AD[nAD] = {664.718016,675.560234,593.574780,0.0,79.048929,78.748559,77.848677,0.0};
//---*****************************************************---//
//---*****************************************************---//
//-- Information for 8AD
double daqTime_8AD[nAD];         //--DAQTimeTotal - DAQTime6AD
double IBDrate_data_8AD[nAD][2]; //-- (IBDTotalPerDay*DAQTimeTotal - IBD6ADPerDay*DAQTime6AD) / DAQTime8AD
double totalBgd_8AD[nAD][2];     //-- (TotalPerDay*DAQTimeTotal - 6ADPerDay*DAQTime6AD) / DAQTime8AD
double emuem_8AD[nAD];           //-- = emuem_Total[nAD];
//IBD rate per day w/o oscillations
double noOsc_IBDrate_perday_8AD[nAD] = {0.0};
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
double spc8AD[nAD][NB];
double NoscTot8AD[nAD];
double noNoscTot8AD[nAD];
double spc6AD[nAD][NB];
double NoscTot6AD[nAD];
double noNoscTot6AD[nAD];
double xbins[NB+1];
int iAD;
//---*****************************************************---//
const int nEH = 3;
TH1F *data_spect_histo8AD[nEH];
TH1F *bkgd_spect_histo8AD[nEH];
TH1F *data_spect_histo6AD[nEH];
TH1F *bkgd_spect_histo6AD[nEH];
//---*****************************************************---//
TH1F *nosc_spect_hist8AD[nAD];
TH1F *nu_nosc_spect_hist8AD[nAD];
TH1F *wosc_spect_hist8AD[nAD];
TH1F *nu_wosc_spect_hist8AD[nAD];
TH1F *bfit_spect_hist8AD[nAD];
TH1F *nosc_spect_hist6AD[nAD];
TH1F *nu_nosc_spect_hist6AD[nAD];
TH1F *wosc_spect_hist6AD[nAD];
TH1F *nu_wosc_spect_hist6AD[nAD];
TH1F *bfit_spect_hist6AD[nAD];
//---*****************************************************---//
//-------------------
// Covariance-matrix Histograms
//-------------------
//** Declarar la matriz de covarianza fraccionaria como un TH2F
//TFile *fracCovaMatrix_File = new TFile("files_data/db_CovaMatrix_6AD_6x26bins.root","READ");
TH2F *fracCovaMatrix_hist;    //Se carga dentro de la función principal db_minuot_spec()
int NBx_cov = 2*NB*nAD;
int NBy_cov = 2*NB*nAD;
TMatrixD statErroMatrix_matrix(NBx_cov,NBy_cov);
TMatrixD fracCovaMatrix_matrix(NBx_cov,NBy_cov);  //Se llena en la función principal db_minuit_spec()
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
  
    //---*****************************************************---//
    int iNR;
    double SurvPavg6AD  = 0.0;
    double Md6AD        = 0.0;
    double Td6AD        = 0.0;
    double Bd6AD        = 0.0;
    double sB6AD        = 0.0;
    double sqrerror6AD  = 0.0;
    double SurvPavg8AD  = 0.0;
    double Md8AD        = 0.0;
    double Td8AD        = 0.0;
    double Bd8AD        = 0.0;
    double sB8AD        = 0.0;
    double sqrerror8AD  = 0.0;
    double sqr_chi      = 0.0;
    double seps_d       = 1.0*0.002;
    double salph_r      = 1.0*0.008;
    //cout << "ChiSq function... " << endl;
    for (iAD = 0 ; iAD < nAD ; iAD++){
        SurvPavg6AD = NoscTot6AD[iAD]/noNoscTot6AD[iAD];
        SurvPavg8AD = NoscTot8AD[iAD]/noNoscTot8AD[iAD];
        for (int iBIN = 0 ; iBIN < NB ; iBIN++)
        {
            int index  = iAD*NB + iBIN;
            int index2 = index + nAD*NB;
            //-- Predicted IBD from neutrino oscillations of the dth Antineutrino Detector
            Td6AD = 0;
            if (NoscTot6AD[iAD]!=0)
            Td6AD = spc6AD[iAD][iBIN]*(SurvPavg6AD*noOsc_IBDrate_perday_6AD[iAD]/NoscTot6AD[iAD])*emuem_6AD[iAD]*daqTime_6AD[iAD];
            Td8AD = 0;
            if (NoscTot8AD[iAD]!=0)
            Td8AD = spc8AD[iAD][iBIN]*(SurvPavg8AD*noOsc_IBDrate_perday_8AD[iAD]/NoscTot8AD[iAD])*emuem_8AD[iAD]*daqTime_8AD[iAD];
            //cout << "Td6AD = " << Td6AD << endl;
            //cout << "Td8AD = " << Td8AD << endl;           

            //-- Measured IDB events of the dth Antineutrino Detector (background is substracted)
            int idx = -1;
            if (iAD < 2) {
                idx = 0;
            }
            else if (iAD == 2 || iAD == 3) {
                idx = 1;
            } else {
                idx = 2;
            }
            double data_sigplusbg_6AD = (data_spect_histo6AD[idx]->GetBinContent(iBIN+1))*IBDrate_data_6AD[iAD][0];
            double simu_bg_6AD = (bkgd_spect_histo6AD[idx]->GetBinContent(iBIN+1))*totalBgd_6AD[iAD][0]*emuem_6AD[iAD];
            Md6AD = (data_sigplusbg_6AD - simu_bg_6AD)*daqTime_6AD[iAD];
            double data_sigplusbg_8AD = (data_spect_histo8AD[idx]->GetBinContent(iBIN+1))*IBDrate_data_8AD[iAD][0];
            double simu_bg_8AD = (bkgd_spect_histo8AD[idx]->GetBinContent(iBIN+1))*totalBgd_8AD[iAD][0]*emuem_8AD[iAD];
            Md8AD = (data_sigplusbg_8AD - simu_bg_8AD)*daqTime_8AD[iAD];
            //-- Background of the dth Antineutrino Detector
            Bd6AD = simu_bg_6AD*daqTime_6AD[iAD];
            Bd8AD = simu_bg_8AD*daqTime_8AD[iAD];

            sqrerror6AD = Md6AD + Bd6AD;
            sqrerror8AD = Md8AD + Bd8AD;

            if (sqrerror6AD==0) sqrerror6AD = 1e7;

            //Statistics Error matrix
            statErroMatrix_matrix(index,index)   = sqrerror6AD;
            statErroMatrix_matrix(index2,index2) = sqrerror8AD;

            //-- Fraction of IBD contribution of the rth reactor to the dth AD
            //-- determined by baselines and reactor fluxes
            double wrd = 0.0;
            for (iNR = 0 ; iNR < nNR ; iNR++)
                wrd += wrd_array[iAD][iNR]*alpha[iNR];
        
            predi_vector(index,0)  = Td6AD*(1.0 + epsilon + eps_d[iAD] + wrd) - eta_d[iAD];
            predi_vector(index2,0) = Td8AD*(1.0 + epsilon + eps_d[iAD] + wrd) - eta_d[iAD];
            //cout << "predi_vector(index2,0)" << index2 << " " << predi_vector(index2,0) << endl;
            delta_vector(index,0)  = Md6AD - Td6AD*(1.0 + epsilon + eps_d[iAD] + wrd) + eta_d[iAD];
            delta_vector(index2,0) = Md8AD - Td8AD*(1.0 + epsilon + eps_d[iAD] + wrd) + eta_d[iAD];
            //cout << "delta_vector(index2,0)" << index2 << " " << delta_vector(index2,0) << endl;
            transp_delta_vector(0,index)  = delta_vector(index,0);
            transp_delta_vector(0,index2) = delta_vector(index2,0);
            //cout << "transp_delta_vector(0,index2)" << index2 << " " << transp_delta_vector(0,index2) << endl;
            //cout << endl;
        }
    }

    for (int i = 0 ; i < NBx_cov ; i++) {
        for (int j = 0 ; j < NBy_cov ; j++) {
            fullCovaMatrix_matrix(i,j) = statErroMatrix_matrix(i,j) + 0*predi_vector(i,0)*predi_vector(j,0)*fracCovaMatrix_matrix(i,j);
        }
    }
    predi_vector.Print();
    //fullCovaMatrix_matrix.Print();
    //exit(); //debug

    if(fullCovaMatrix_matrix.Determinant() != 0){
        inv_fullCovaMatrix_matrix = fullCovaMatrix_matrix;
        inv_fullCovaMatrix_matrix.Invert();
    }
    else {
        cout << "Determinant = " << fullCovaMatrix_matrix.Determinant() << endl;
        cout << "Singular fullCovaMatrix!" << endl;
        cout << "Execution aborted!" << endl;
        //break;
        exit(1);
    }
    
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
           sB6AD = totalBgd_6AD[iAD][1]*emuem_6AD[iAD]*daqTime_6AD[iAD];
           sB8AD = totalBgd_8AD[iAD][1]*emuem_8AD[iAD]*daqTime_8AD[iAD];

           sqr_chi += pow(eps_d[iAD]/seps_d,2) + pow(eta_d[iAD]/(sB6AD+sB8AD),2);
       }
    
    for (iNR = 0 ; iNR < nNR ; iNR++)
         sqr_chi += pow(alpha[iNR]/salph_r,2);
    
    //cout << "AD " << iAD << "  chi2 = " << sqr_chi << endl;
    return sqr_chi;
}
//---*****************************************************---//

int db_minuit_spec_CovMat(const char * minName = "Minuit",
                          const char *algoName = "" ,
                          int randomSeed = +10)
    //int randomSeed = -1)
{
    cout << "Let's begin..." << endl;
  
    //--------------------------------------
    // Covariance-matrix Histograms
    //--------------------------------------
    TFile *fracCovaMatrix_File = new TFile("CombCovaMat/db_CovaMatrix_6AD8AD_16x16Det.root","READ");
    fracCovaMatrix_hist = (TH2F*)(fracCovaMatrix_File->Get("covaMat_histo_16x16"));
    //-- Creating the Rebinned Matrix with root tools.
    //NBx_cov = 2*NB*nAD;
    //NBy_cov = 2*NB*nAD;
    cout << NBx_cov << "  " << NBy_cov << endl;
    //statErroMatrix_matrix = TMatrixD(NBx_cov,NBy_cov);
    //fracCovaMatrix_matrix = TMatrixD(NBx_cov,NBy_cov);
    //fullCovaMatrix_matrix = TMatrixD(NBx_cov,NBy_cov);
    //inv_fullCovaMatrix_matrix = TMatrixD(NBx_cov,NBy_cov);
    
    //one_matrix = TMatrixD(NBx_cov,NBy_cov);
    
    //delta_vector        = TMatrixD(NBx_cov,1);
    //transp_delta_vector = TMatrixD(1,NBx_cov);
    //predi_vector        = TMatrixD(NBx_cov,1);
    
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
    for (int i = 0 ; i < NBx_cov ; i++) {
        for (int j = 0 ; j < NBy_cov; j++) {
            i_block = i%35;
            j_block = j%35;
            fracCovaMatrix_matrix(i,j) = fracCovaMatrix_hist->GetBinContent(i+1,j+1);
        }
    }

    cout << "Done with the matrix!" << endl;
    //fracCovaMatrix_matrix.Print();
    //exit(0); //debug

    TFile *wrd_File = new TFile("files_data/daya-bay-ldist.root","READ");
    TH1F *wrd_histo = (TH1F*)(wrd_File->Get("histo_ldist"));

    //-------------------
    // Energy Histograms
    //-------------------
    //TFile *fenergy0217 = new TFile("./PRL112_217days_data.root","read");
    //TFile *fenergy1230 = new TFile("./PRD95_1230days_data.root","read");
    TFile *fenergy = new TFile("./histos_6AD217Days_8AD1013Days_data.root","read");
    //Three sets of histograms one for each Experimental Hall
    double bfactor[nEH];
    for (int i = 0 ; i < nEH ; i++)
        {
            //Data 217 days
            //data_spect_histo6AD[i] = (TH1F*) fenergy0217->Get(Form("data_spect_histo_%d",i));
            data_spect_histo6AD[i] = (TH1F*) fenergy->Get(Form("data_spect_6AD35B_histo_%d",i));
            double dfactor = 1.0/data_spect_histo6AD[i]->Integral();
            data_spect_histo6AD[i]->Scale(dfactor);
            //Background 217 days
            //bkgd_spect_histo6AD[i] = (TH1F*) fenergy0217->Get(Form("bkgd_spect_histo_%d",i));
            bkgd_spect_histo6AD[i] = (TH1F*) fenergy->Get(Form("bkgd_spect_6AD35B_histo_%d",i));
            bfactor[i] = 1.0/bkgd_spect_histo6AD[i]->Integral();
            bkgd_spect_histo6AD[i]->Scale(bfactor[i]);

            //Data 1013 days
            //data_spect_histo8AD[i] = (TH1F*) fenergy1230->Get(Form("data_spect_histo_%d",i));
            data_spect_histo8AD[i] = (TH1F*) fenergy->Get(Form("data_spect_8AD35B_histo_%d",i));
            dfactor = 1.0/data_spect_histo8AD[i]->Integral();
            data_spect_histo8AD[i]->Scale(dfactor);
            //Background 1013 days
            //bkgd_spect_histo8AD[i] = (TH1F*) fenergy1230->Get(Form("bkgd_spect_histo_%d",i));
            bkgd_spect_histo8AD[i] = (TH1F*) fenergy->Get(Form("bkgd_spect_8AD35B_histo_%d",i));
            bfactor[i] = 1.0/bkgd_spect_histo8AD[i]->Integral();
            bkgd_spect_histo8AD[i]->Scale(bfactor[i]);
        }
    
    // histogram binning
    //double xbins[36]; //1230 Days
    xbins[0] = 0.7;
    double delta_bins2 = (7.9 - 1.3)/33; // 0.20 MeV/bin
    for (int i = 0 ; i < (NB-1) ; i++)
    {
        xbins[i+1] = 1.3 + delta_bins2*i;
    }
    xbins[35] = hi;
    for(int iAD = 0 ; iAD < nAD ; iAD++)
    {
        //1013
        nosc_spect_hist8AD[iAD]    = new TH1F(Form("nosc_spect_hist8AD_%d",iAD),   "",NB,xbins);
        nu_nosc_spect_hist8AD[iAD] = new TH1F(Form("nu_nosc_spect_hist8AD_%d",iAD),"",NB,xbins);
        wosc_spect_hist8AD[iAD]    = new TH1F(Form("wosc_spect_hist8AD_%d",iAD),   "",NB,xbins);
        nu_wosc_spect_hist8AD[iAD] = new TH1F(Form("nu_wosc_spect_hist8AD_%d",iAD),"",NB,xbins);
        bfit_spect_hist8AD[iAD]    = new TH1F(Form("bfit_spect_hist8AD_%d",iAD),   "",NB,xbins);
        //217
        nosc_spect_hist6AD[iAD]    = new TH1F(Form("nosc_spect_hist6AD_%d",iAD),   "",NB,xbins);
        nu_nosc_spect_hist6AD[iAD] = new TH1F(Form("nu_nosc_spect_hist6AD_%d",iAD),"",NB,xbins);
        wosc_spect_hist6AD[iAD]    = new TH1F(Form("wosc_spect_hist6AD_%d",iAD),   "",NB,xbins);
        nu_wosc_spect_hist6AD[iAD] = new TH1F(Form("nu_wosc_spect_hist6AD_%d",iAD),"",NB,xbins);
        bfit_spect_hist6AD[iAD]    = new TH1F(Form("bfit_spect_hist6AD_%d",iAD),   "",NB,xbins);
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
    ifstream IBDrates_file1230("files_data/db_noOsc_IBDrates_perday_1230.txt");
    //ifstream IBDrates_file0217("files_data/db_noOsc_IBDrates_perday_217.txt");
    cout << "Reading noOsc normalizations file ..." << endl;
    for (int i=0; i< nAD; i ++)
	{
        IBDrates_file1230 >> noOsc_IBDrate_perday_6AD8AD[i];
        cout << "noOscIBD_rates_perday_6AD8AD " << i << ": " << noOsc_IBDrate_perday_6AD8AD[i] << endl;
        daqTime_8AD[i]         = daqTime_Total[i] - daqTime_6AD[i];
        noOsc_IBDrate_perday_8AD[i] = ((noOsc_IBDrate_perday_6AD8AD[i]*daqTime_Total[i]) - (noOsc_IBDrate_perday_6AD[i]*daqTime_6AD[i]))/daqTime_8AD[i];
        cout << "noOscIBD_rates_perday_6AD " << i << ": " << noOsc_IBDrate_perday_6AD[i] << endl;
        cout << "noOscIBD_rates_perday_8AD " << i << ": " << noOsc_IBDrate_perday_8AD[i] << endl;
    }//for
 
    //-- File to print oscillation parameters and chi2 values
    ofstream chi2Surface_file;
    string s2t_dm2 = "files_data/chi2_s2t-dm2_surface_SPEC.txt";  //(sin^2(2th13), dm2, chi^2_min)
    chi2Surface_file.open((s2t_dm2).c_str());
  
    //-- File to print minimized pull parameters
    ofstream minimPullT_file;
    string PullT = "files_data/chi2_pullT_surface.txt";
    minimPullT_file.open((PullT).c_str());
  
    //ifstream  file("files_data/db_gridOscSpectra_1230.txt");
    string grid_spectra = "files_data/db_gridOscSpectra_1230.txt";
    ifstream file;
    file.open((grid_spectra).c_str());
    cout << "Reading file - Loop in progress..." << endl;
    int iad = 0;
    int first8 = 1;
    std::cout << "Is the file open?  " << file.is_open() << "\n" << std::endl;

    for (int iAD ; iAD < nAD ; iAD++) {
        //daqTime_8AD[iAD]         = daqTime_Total[iAD] - daqTime_6AD[iAD];
        IBDrate_data_8AD[iAD][0] = (IBDrate_data_Total[iAD][0]*daqTime_Total[iAD] - IBDrate_data_6AD[iAD][0]*daqTime_6AD[iAD]) / daqTime_8AD[iAD];
        totalBgd_8AD[iAD][0]     = (totalBgd_Total[iAD][0]*daqTime_Total[iAD]     - totalBgd_6AD[iAD][0]*daqTime_6AD[iAD]) / daqTime_8AD[iAD];
        emuem_8AD[iAD]           = emuem_Total[iAD];
    }
    
    while (file >> iAD >> s2th_13 >> dm2_31 >> spc6AD[iad][0]  >> spc6AD[iad][1]  >> spc6AD[iad][2]  >> spc6AD[iad][3]  >> spc6AD[iad][4]  >> spc6AD[iad][5]  >> spc6AD[iad][6]  >> spc6AD[iad][7]  >> spc6AD[iad][8]  >> spc6AD[iad][9]  >> spc6AD[iad][10] >> spc6AD[iad][11] >> spc6AD[iad][12] >> spc6AD[iad][13] >> spc6AD[iad][14] >> spc6AD[iad][15] >> spc6AD[iad][16] >> spc6AD[iad][17] >> spc6AD[iad][18] >> spc6AD[iad][19] >> spc6AD[iad][20] >> spc6AD[iad][21] >> spc6AD[iad][22] >> spc6AD[iad][23] >> spc6AD[iad][24] >> spc6AD[iad][25] >> spc6AD[iad][26] >> spc6AD[iad][27] >> spc6AD[iad][28] >> spc6AD[iad][29] >> spc6AD[iad][30] >> spc6AD[iad][31] >> spc6AD[iad][32] >> spc6AD[iad][33] >> spc6AD[iad][34] >> NoscTot6AD[iad] >> spc8AD[iad][0]  >> spc8AD[iad][1]  >> spc8AD[iad][2]  >> spc8AD[iad][3]  >> spc8AD[iad][4]  >> spc8AD[iad][5]  >> spc8AD[iad][6]  >> spc8AD[iad][7]  >> spc8AD[iad][8]  >> spc8AD[iad][9]  >> spc8AD[iad][10] >> spc8AD[iad][11] >> spc8AD[iad][12] >> spc8AD[iad][13] >> spc8AD[iad][14] >> spc8AD[iad][15] >> spc8AD[iad][16] >> spc8AD[iad][17] >> spc8AD[iad][18] >> spc8AD[iad][19] >> spc8AD[iad][20] >> spc8AD[iad][21] >> spc8AD[iad][22] >> spc8AD[iad][23] >> spc8AD[iad][24] >> spc8AD[iad][25] >> spc8AD[iad][26] >> spc8AD[iad][27] >> spc8AD[iad][28] >> spc8AD[iad][29] >> spc8AD[iad][30] >> spc8AD[iad][31] >> spc8AD[iad][32] >> spc8AD[iad][33] >> spc8AD[iad][34] >> NoscTot8AD[iad])
    {//file loop
        //cout << s2th_13 << "\t" << dm2_31 << endl;
        if(first8 <= 8)
        {
            double bincont6AD = 0.0;
            double bincont8AD = 0.0;
            for(int ibin = 1 ; ibin <= NB ; ibin++)
            {
                if (iad != 3 && iad != 7)
                    bincont6AD = spc6AD[iad][ibin-1]*(noOsc_IBDrate_perday_6AD[iAD-1]/NoscTot6AD[iad])*emuem_6AD[iAD-1]*daqTime_6AD[iAD-1];
                nu_nosc_spect_hist6AD[iAD-1]->SetBinContent(ibin,bincont6AD);
                
                bincont8AD = spc8AD[iad][ibin-1]*(noOsc_IBDrate_perday_8AD[iAD-1]/NoscTot8AD[iad])*emuem_8AD[iAD-1]*daqTime_8AD[iAD-1];
                nu_nosc_spect_hist8AD[iAD-1]->SetBinContent(ibin,bincont8AD);
            }
            cout << "iad = " << iad << "   iAD = " << iAD  << "   first8 = " << first8 << endl;
            cout << "Number of events 6AD: " << nu_nosc_spect_hist6AD[iAD-1]->Integral() << endl;
            cout << "Number of events 8AD: " << nu_nosc_spect_hist8AD[iAD-1]->Integral() << endl;
            noNoscTot6AD[iad] = NoscTot6AD[iad];
            noNoscTot8AD[iad] = NoscTot8AD[iad];
            first8++;
        }//if first8 loop END
        
        iad++;
        
        //At this point, we have read the 16 AD spectra (8 lines) for one point (s2th_13,dm2_31) in the grid
        if (iad == 8)
        {//if iad BEGIN
            //cout << endl;
            iad = 0;
            
            //cout << "test " << s2th_13 << "\t" << dm2_31 << endl;

            //This is executed for every 6 lines of the file
            const int N_params = 23; //-- Number of parameter of the chi² function --//
            ROOT::Math::Functor f(&chi2,N_params); //-- Setting the function to be minimized by using Minuit --//
            //-- Steps
            double stp = 1.0e-3;
            double step[N_params] = {stp,stp,stp,stp,stp,stp,stp,stp,
                                        stp,stp,stp,stp,stp,stp,stp,stp,
                                        stp,stp,stp,stp,stp,stp,stp};
            //-- Initial parameter values
            double start[N_params] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
                                        0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
                                        0.0,0.0,0.0,0.0,0.0,0.0,0.0};
            
            //-- Calling Minuit function setting
            min->SetFunction(f);
                
            //-- Setting variables
            double lim = 1.0e-1;
            /*
            min->SetLimitedVariable(0,  "e_1", start[0],  step[0],  -lim, lim);
            min->SetLimitedVariable(1,  "e_2", start[1],  step[1],  -lim, lim);
            min->SetLimitedVariable(2,  "e_3", start[2],  step[2],  -lim, lim);
            min->SetLimitedVariable(3,  "e_4", start[3],  step[3],  -lim, lim);
            min->SetLimitedVariable(4,  "e_5", start[4],  step[4],  -lim, lim);
            min->SetLimitedVariable(5,  "e_6", start[5],  step[5],  -lim, lim);
            min->SetLimitedVariable(6,  "e_7", start[6],  step[6],  -lim, lim);
            min->SetLimitedVariable(7,  "e_8", start[7],  step[7],  -lim, lim);
                */
            /*
            min->SetLimitedVariable(8,  "n_1", start[8],  step[8],  -lim, lim);
            min->SetLimitedVariable(9,  "n_2", start[9],  step[9],  -lim, lim);
            min->SetLimitedVariable(10, "n_3", start[10], step[10], -lim, lim);
            min->SetLimitedVariable(11, "n_4", start[11], step[11], -lim, lim);
            min->SetLimitedVariable(12, "n_5", start[12], step[12], -lim, lim);
            min->SetLimitedVariable(13, "n_6", start[13], step[13], -lim, lim);
            min->SetLimitedVariable(14, "n_7", start[14], step[14], -lim, lim);
            min->SetLimitedVariable(15, "n_8", start[15], step[15], -lim, lim);
                */
            /*
            min->SetLimitedVariable(16, "a_1", start[16], step[16], -lim, lim);
            min->SetLimitedVariable(17, "a_2", start[17], step[17], -lim, lim);
            min->SetLimitedVariable(18, "a_3", start[18], step[18], -lim, lim);
            min->SetLimitedVariable(19, "a_4", start[19], step[19], -lim, lim);
            min->SetLimitedVariable(20, "a_5", start[20], step[20], -lim, lim);
            min->SetLimitedVariable(21, "a_6", start[21], step[21], -lim, lim);
            */
            min->SetLimitedVariable(22, "eps", start[22], step[22], -lim, lim);
            //-----------------------------------------------------------------//
		
            min->SetFixedVariable(0,  "e_1", start[0]);
            min->SetFixedVariable(1,  "e_2", start[1]);
            min->SetFixedVariable(2,  "e_3", start[2]);
            min->SetFixedVariable(3,  "e_4", start[3]);
            min->SetFixedVariable(4,  "e_5", start[4]);
            min->SetFixedVariable(5,  "e_6", start[5]);
            min->SetFixedVariable(6,  "e_7", start[6]);
            min->SetFixedVariable(7,  "e_8", start[7]);
                
            min->SetFixedVariable(8,  "n_1", start[8]);
            min->SetFixedVariable(9,  "n_2", start[9]);
            min->SetFixedVariable(10, "n_3", start[10]);
            min->SetFixedVariable(11, "n_4", start[11]);
            min->SetFixedVariable(12, "n_5", start[12]);
            min->SetFixedVariable(13, "n_6", start[13]);
            min->SetFixedVariable(14, "n_7", start[14]);
            min->SetFixedVariable(15, "n_8", start[15]);
                
            min->SetFixedVariable(16, "a_1", start[16]);
            min->SetFixedVariable(17, "a_2", start[17]);
            min->SetFixedVariable(18, "a_3", start[18]);
            min->SetFixedVariable(19, "a_4", start[19]);
            min->SetFixedVariable(20, "a_5", start[20]);
            min->SetFixedVariable(21, "a_6", start[21]);
                
            //min->SetFixedVariable(22, "eps", start[22]);
            min->SetErrorDef(2.3);
            
            //-- Calling Minuit minimization
            //cout << "Calling minimization..." << endl;
            min->Minimize();
          
            const double *xs = min->X();
            
            double chi2Min = min->MinValue();
            chi2Surface_file << s2th_13 << "\t" << dm2_31 << "\t" << chi2Min << endl;
        
            if (chi2Min < 0.0) {
                cout << "Critical error: chi2Min is negative!  " << chi2Min << endl;
                //break;
                exit(1);
            }
                
            minimPullT_file  << s2th_13 << "\t" << dm2_31 << "\t" << xs[0]  << "\t" << xs[1]
                << "\t" << xs[2]  << "\t" << xs[3]  << "\t" << xs[4]  << "\t" << xs[5]
                << "\t" << xs[6]  << "\t" << xs[7]  << "\t" << xs[8]  << "\t" << xs[9]
                << "\t" << xs[10] << "\t" << xs[11] << "\t" << xs[12] << "\t" << xs[13]
                << "\t" << xs[14] << "\t" << xs[15] << "\t" << xs[16] << "\t" << xs[17]
                << "\t" << xs[18] << "\t" << xs[19] << "\t" << xs[20] << "\t" << xs[21]
                << "\t" << xs[22] << "\t" << min->MinValue() << endl;
                //}
            if (dm2_31 == hi_dm2)
            {
                chi2Surface_file << endl;
                cout << s2th_13 << "\t" << dm2_31 << "\t" << min->MinValue() << endl;
            }
        }//if iad END

    }//file loop END
    std::cout << "Successful run!!" << endl;
    //------------------------------------------------------------------------------
    return 0;
}
