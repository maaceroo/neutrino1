//----------------------------------------------------------//
//--  db_ratio.C   -   By M.A. Acero O.   -   2016-03-17  --//
//-- Macro to plot the ratio of measured to expected sig- --//
//-- nal in each detector using the best fit values for   --//
//-- oscillation parameters obtained from using the macro --//
//-- 'db_minuit.C'.                                       --//
//-- Analysis of the Daya Bay data from F.P. An et al.,   --//
//-- PRL112 061801 (2014)                                 --//
//----------------------------------------------------------//

//---*********************************************************---//
//-- Global variables and quantities ----------------------------//
//---*********************************************************---//
//EH1(AD1, AD2),EH2(AD3),EH3(AD4, AD5, AD6)
//(IBD candidates)/(DAQ live time -days-) from PRL 112 061801 (2014)
//double IBDrate_data[nAD][2] = { {530.31,1.67},{536.75,1.68},{489.93,1.61},{ 73.58,0.62},{ 73.21,0.62},{72.35,0.62} };
//IBD rate (per day), total background and efficiencies (PRL 112 061801 (2014))
//double totalBgd[nAD][2] = { {13.20,0.98},{13.01,0.98},{ 9.57,0.71},{ 3.52,0.14},{ 3.48,0.14},{3.43,0.14} };
//double emuem[nAD] ={0.7957,0.7927,0.8282,0.9577,0.9568,0.9566};
//double daqTime[nAD] = {191.001,191.001,189.645,189.779,189.779,189.779};
//---*********************************************************---//
// Information obtained by executing the script "db_osc_rate.C"
//IBD rate per day w/o oscillations
//double noOsc_IBDrate_perday[nAD] = {663.15,673.95,591.86,78.75,78.46,77.58};
//<sin^2(1.267 dm2_21 L/E)> for each AD
//double avgSinDelta21[nAD] = {0.000225897,0.000220442,0.000244761,0.00194541,0.00193977,0.00195247};
//<sin^2(1.267 dm2_31 L/E)> for each AD
//double avgSinDelta31[nAD] = {0.162935,0.159482,0.183415,0.750218,0.75162,0.752462};
//---------------------------------------------------


void db_ratio()
{//begin
    
    //------------- Style --------------
    gROOT->SetStyle("Plain");
    gStyle->SetOptStat(0);
    //------------- Style --------------
    //----------  Text Style  ---------
    //ft = 10 * fontID + precision
    Int_t ft = 10 * 4 + 2;
    Double_t sz = 0.04;
    //----------  Text Style  ---------

    //---------------------------------------------------
    // Open ntuple file to read simulated data
    TFile *fntuple = new TFile("files_data/db-ntuple_5M.root","READ");
    TTree *T = (TTree*)fntuple->Get("T");
    //---------------------------------------------------
    // Baseline histogram
//    TH1F *BaseLine_histo = new TH1F("BaseLine_histo","",10,0,2000);
//    T->Draw("Ln >> BaseLine_histo","Ln>0");
//    double integ = BaseLine_histo->Integral();
//    BaseLine_histo->Scale(1.0/integ);
    //---------------------------------------------------
    //---------------------------------------------------
    // Weighted Baseline and energy histograms
    const int nAD = 6;
    int sel;
    TH1F *wBL_histo[nAD];
    TH1F *wNE_histo[nAD];
    double wBL_avg[nAD];
    double wNE_avg[nAD];
    for (int i = 0 ; i < nAD ; i++) {
        wBL_histo[i] = new TH1F(Form("wBL_histo_%d",i),"",100,0,2000);
        wNE_histo[i] = new TH1F(Form("wNE_histo_%d",i),"",100,0,12);
        if (i < 3) sel = i;
        else if (i >= 3) sel = i + 1;
        T->Draw(Form("Ln >> wBL_histo_%d",i),Form("id==%d",sel));
        T->Draw(Form("En >> wNE_histo_%d",i),Form("id==%d",sel));
    }
    for (int i = 0 ; i < nAD ; i++) {
        wBL_avg[i] = wBL_histo[i]->GetMean();
        wNE_avg[i] = wNE_histo[i]->GetMean();
        cout << "L_avg(" << i << ") = " << wBL_avg[i] << ", E_avg(" << i << ") = " << wNE_avg[i] << endl;
    }
    //---------------------------------------------------
    //IBD rate per day w/o oscillations
    double noOsc_IBDrate_perday[nAD] = {663.15,673.95,591.86,78.75,78.46,77.58};
    //(IBD candidates)/(DAQ live time -days-) from PRL 112 061801 (2014)
    double IBDrate_data[nAD][2] = { {530.31,1.67},{536.75,1.68},{489.93,1.61},{ 73.58,0.62},{ 73.21,0.62},{72.35,0.62} };
    //IBD rate (per day), total background and efficiencies (PRL 112 061801 (2014))
    double totalBgd[nAD][2] = { {13.20,0.98},{13.01,0.98},{ 9.57,0.71},{ 3.52,0.14},{ 3.48,0.14},{3.43,0.14} };
    double emuem[nAD] ={0.7957,0.7927,0.8282,0.9577,0.9568,0.9566};

    //---------------------------------------------------
    //-- File to print survival probability
    ofstream survP_file;
    string sp_file = "survProb.txt";  //(BL(m), survProb)
    survP_file.open((sp_file).c_str());
    //Oscilltion Probability
    double Nexp;
    double Ndet;
    double Posc_BF[nAD];
    double error;
    double s2t = 0.090;
    double dm2 = 2.32e-3;
    double epsBF = 0.000150754;
    for (int i = 0 ; i < nAD ; i++) {
        Posc_BF[i] = 1.0 - s2t*((sin( 1.267 * dm2 * wBL_avg[i]/wNE_avg[i] ))**2) - (0.25*pow(1.0 + sqrt(1.0-s2t),2)) * 0.861 * pow(sin( 1.267 * 7.59e-5 * wBL_avg[i]/wNE_avg[i] ),2);
        //Nexp = (Posc_BF[i] * noOsc_IBDrate_perday[i])*emuem[i]*(1.0 + epsBF);
        Nexp = (noOsc_IBDrate_perday[i])*emuem[i]*(1.0 + epsBF);
        Ndet = (IBDrate_data[i][0] - totalBgd[i][0]*emuem[i]); //- totalBgd[i][0]*emuem[i]
        error = (Ndet/Nexp)*sqrt(pow(IBDrate_data[i][1],2) + pow(totalBgd[i][1],2))/Nexp;
        cout << "SurvProb = " << Posc_BF[i] << ",  RatioCalc = " << Ndet/Nexp
             << ",  Ratio = " << IBDrate_data[i][0]/(noOsc_IBDrate_perday[i]*(1+epsBF))
             << ",  error = " << error << endl;
        survP_file << wBL_avg[i] << "  " << Posc_BF[i] << "  "
                   << Ndet/Nexp  << "  " << error << endl;
    }
    
    //---------------------------------------------------
    // Drawing section
/*    TH2F *frame_h = new TH2F("frame_h","",1000,0,2000,10,0,0.13);
    frame_h->GetXaxis()->SetTitleFont(ft);
    frame_h->GetXaxis()->SetTitleSize(sz);
    frame_h->GetXaxis()->SetLabelFont(ft);
    frame_h->GetXaxis()->SetLabelSize(sz);
    frame_h->GetXaxis()->SetTitle("Baseline (m)");
    frame_h->GetYaxis()->SetTitleFont(ft);
    frame_h->GetYaxis()->SetTitleSize(sz);
    frame_h->GetYaxis()->SetLabelFont(ft);
    frame_h->GetYaxis()->SetLabelSize(sz);
    frame_h->GetYaxis()->SetTitle("a.u.");
*/
//    TCanvas *canv0 = new TCanvas("canv0","canv0",775,500);

//    frame_h->Draw();
//    BaseLine_histo->Draw("same");
    //---------------------------------------------------

}//end
