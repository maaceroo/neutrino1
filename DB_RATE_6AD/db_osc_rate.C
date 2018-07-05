//---------------------------------------------------------------------//
//-- db_osc_rate.C - By M.A. Acero O. & A.A. Aguilar-A. - 2017-01-02 --//
//---------------------------------------------------------------------//
// This macro can be executed under ROOT typing                        //
//"root[0] .x db_osc_rate.C"                                           //
// For the rate-only analysis, using information from F.P. An et al.,  //
// PRL 112 061801 (2014)                                               //
//---------------------------------------------------------------------//
// 2017-01-31                                                          //
// Modifications to perform a rate-only analysis by computing and using//
// average values for the survival probability and taking sin2(2th) as //
// the only free parameter (to bi fitted).                             //
// 2017-02-03                                                          //
// This macro performs a simple chi^2 analysis of the Daya Bay data by //
// comparing the IBD rate at the six AD reported in PRL 112 061801     //
// (2014)                                                              //
//---------------------------------------------------------------------//

#include "constants.h"
#include "rate_SurvProb.C"
#include <iostream>
#include <fstream>
#include <string>

void db_osc_rate()
{ //begin

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
    // Open ntuple file to read simulated data (Ep, En, Ln) for the 6 AD
    TFile *fntuple = new TFile("files_data/db-ntuple_5M.root","READ");
    TTree *T = (TTree*)fntuple->Get("T");
    TCut cutBF;
    //---------------------------------------------------
    // histogram binning for the simulated data
    double    NB = 26;
    double    lo = 0.7;
    double    hi = 12.0;
    double    xbins[27];
    xbins[0] = lo;
    double delta_bins2 = (7.3 - 1.3)/24.;// = 0.25 MeV/bin
    for (int i=0;i<(NB-1);i++)
    {
        xbins[i+1] = 1.3 + delta_bins2*i;
    }
    xbins[26] = hi;
    //---------------------------------------------------
    const int nAD = 6; //Number os Antineutrino detectors
    TH1F *BFit_spect_histo[nAD];
    TH1F *Posc_AD_BF[nAD];
    TH1F *Posc_AD_surv[nAD];
    TH1F *SinDelta21[nAD];
    TH1F *SinDelta31[nAD];
    for (int i = 0 ; i < nAD ; i++)
    {
        //BF-oscillation Ep spectra
        BFit_spect_histo[i] = new TH1F(Form("BFit_spect_histo_%d",i),"",NB,xbins);// info from db-ntuple.root
        BFit_spect_histo[i]->SetLineColor(i+1);

        //Ocillation prpbability at BF - histogram
        Posc_AD_BF[i]       = new TH1F(Form("Posc_AD_BF_%d",i),"",1000,0,1);//to store <POsc(BF)>
        Posc_AD_BF[i]->SetLineColor(i+1);

        //Ocillation prpbability at BF - histogram
        SinDelta21[i]       = new TH1F(Form("SinDelta21_%d",i),"",1000,0,1);//to store <sin^2(1.267 dm2_21 L/E)>
        SinDelta21[i]->SetLineColor(i+1);

        //Ocillation prpbability at BF - histogram
        SinDelta31[i]       = new TH1F(Form("SinDelta31_%d",i),"",1000,0,1);//to store <sin^2(1.267 dm2_31 L/E)>
        SinDelta31[i]->SetLineColor(i+2);

        //Ocillation prpbability at (s2th,dm2) - histogram
        Posc_AD_surv[i]     = new TH1F(Form("Posc_AD_surv_%d",i),"",1000,0,1);//to store <POsc(s2th,dm2)>
    }
    //---------------------------------------------------
    //IBD rate (per day), total background and efficiencies (PRL 112 061801 (2014))
    //EH1(AD1, AD2),EH2(AD3),EH3(AD4, AD5, AD6)
    double IBDrate_perday[nAD][2] =
    {
        {653.30,2.31},{664.15,2.33},
        {581.97,2.33},
        { 73.31,0.66},{ 73.03,0.66},{72.20,0.66}
    };
    //(IBD candidates)/(DAQ live time -days-) from PRL 112 061801 (2014)
    double IBDrate_data[nAD][2] =
    {
        {530.31,23.03},{536.75,23.17},
        {489.93,22.13},
        { 73.58,8.58},{ 73.21,8.56},{72.35,8.51}
    };
    double totalBgd[nAD][2] =
    {
        {13.20,0.98},{13.01,0.98},
        { 9.57,0.71},
        { 3.52,0.14},{ 3.48,0.14},{3.43,0.14}
    };
    double emuem[nAD] =
    {
        0.7957,0.7927,
        0.8282,
        0.9577,0.9568,0.9566
    };

    int sel;
    double avgPosc_AD[nAD]; //<POsc(s2t_BF,dm2_31)>
    double noOsc_IBDrate_perday[nAD]; //IBD rate per day w/o oscillations
    double integ;
    for (int iAD = 0 ; iAD < nAD ; iAD++)
    {
        if (iAD < 3) sel = iAD;
        else if (iAD >= 3) sel = iAD +1;
        //------------------------------------------------
        //Filling Ocillation prpbability at BF - histogram
        T->Draw(Form("(1.0 - 0.090*((sin( 1.267 * 2.32e-3 * Ln/En ))**2) - ((cos(0.5 * asin(sqrt(0.090))))**4) * 0.861 * (sin( 1.267 * 7.59e-5 * Ln/En ))**2) >> Posc_AD_BF_%d",iAD),Form("id==%d",sel));
        integ = Posc_AD_BF[iAD]->Integral();
        Posc_AD_BF[iAD]->Scale(1.0/integ); //Used to plot the Survival Probabilities for each AD with BF parameters (normalized)
        
        //Average oscillation Probability
        avgPosc_AD[iAD] = Posc_AD_BF[iAD]->GetMean();
        //No-oscillation IDB rate (per day)
        noOsc_IBDrate_perday[iAD] = IBDrate_perday[iAD][0]/avgPosc_AD[iAD];
        //Printing results
        cout << "(avgPosc_AD,noOsc_IBDrate_perday)_" << sel << " = (" << avgPosc_AD[iAD]
             << ", " << noOsc_IBDrate_perday[iAD] << ") " << endl;
        //------------------------------------------------

        //------------------------------------------------
        //(BF-oscillation probability) condition to fill BF-oscillation- Ep spectra
        cutBF = Form("(1.0 - 0.090*((sin( 1.267 * 2.32e-3 * Ln/En ))**2) - ((cos(0.5 * asin(sqrt(0.090))))**4) * 0.861 * (sin( 1.267 * 7.59e-5 * Ln/En ))**2)*(id==%d)",iAD);
        //Filling and normalizing BF-oscillation Ep spectra
        T->Draw(Form("Ep >> BFit_spect_histo_%d",sel),cutBF,"");
        integ = BFit_spect_histo[iAD]->Integral();
        BFit_spect_histo[iAD]->Scale(1.0/integ);
        //------------------------------------------------
    }

    //---------------------------------------------------
    // computing <sin^2(1.267 dm2_21 L/E)> and <sin^2(1.267 dm2_31 L/E)>
    
    double avgSinDelta21[nAD]; //<sin^2(1.267 dm2_21 L/E)> for each AD
    double avgSinDelta31[nAD]; //<sin^2(1.267 dm2_31 L/E)> for each AD

    double dm2_21 = 7.59e-5; //eV^2,                  //PRL 108 171803 (2012)
    double dm2_31 = 2.32e-3; //eV^2,                  //PRL 108 171803 (2012)
    int   selad;

    for (int iAD = 0 ; iAD < nAD ; iAD++)
    {
        if (iAD < 3) selad = iAD;
        else if (iAD >= 3) selad = iAD +1;
        
        T->Draw(Form("((sin(1.267 * %e * Ln/En))**2) >> SinDelta21_%d",dm2_21,iAD),Form("id==%d",selad));
        integ = SinDelta21[iAD]->Integral();
        SinDelta21[iAD]->Scale(1.0/integ); //Used to plot
        avgSinDelta21[iAD] = SinDelta21[iAD]->GetMean();
        T->Draw(Form("((sin(1.267 * %e * Ln/En))**2) >> SinDelta31_%d",dm2_31,iAD),Form("id==%d",selad));
        integ = SinDelta31[iAD]->Integral();
        SinDelta31[iAD]->Scale(1.0/integ); //Used to plot
        avgSinDelta31[iAD] = SinDelta31[iAD]->GetMean();
        //cout << "avgSinDelta21_" << iAD << " = " << avgSinDelta21[iAD] << "\t"
             //<< "avgSinDelta31_" << iAD << " = " << avgSinDelta31[iAD] << endl;
    }
    

    //---------------------------------------------------
    //Survival probability for different values of sin2(2th_13)

    // Survival probability function for the rate-only analysis
    gROOT->ProcessLine(".L rate_SurvProb.C");
    TF1 *SurvProb = new TF1("SurvProb",rate_SurvProb,0,1,2);
    SurvProb->SetLineColor(2);
    SurvProb->SetLineStyle(2);
    SurvProb->SetLineWidth(2);
    //---------------------------------------------------
    
    for(int i = 0 ; i < nAD ; i++){
        SurvProb->SetParameters(avgSinDelta21[i],avgSinDelta31[i]);
        double test = SurvProb->Eval(0.090);
        //cout << i << "test = " << test << endl;
    }
    
    
    double s2t_pt;
    
    const int     N_s2t = 1000;
    
//    double       lo_s2t = 0.0;
    double       lo_s2t = 0.01;
    double       hi_s2t = 1.0;
    double DeltaLin_s2t = (hi_s2t - lo_s2t)/double(N_s2t-1);
    double DeltaLog_s2t = (log10(hi_s2t)-log10(lo_s2t))/double(N_s2t-1);
    
    double SurvP;
    double IBDrate_osc[N_s2t];
    
    // Information and quantities to compute chi^2
    double chi2[N_s2t];
    for (int i = 0 ; i < N_s2t; i++) {
        chi2[i] = 0.0;
    }
    double sig_r = 0.008;
    double sig_d = 0.002;
    
    //File to print the chi-square results
    ofstream file;
    string result = "files_data/db_chi2_rate.txt";
    file.open (result.c_str());
    file << setprecision(5);

    for (int is2t = 0 ; is2t < N_s2t ; is2t++)
    {
        for (int iAD = 0 ; iAD < nAD ; iAD++)
        {
            SurvProb->SetParameter(0,avgSinDelta21[iAD]);
            SurvProb->SetParameter(1,avgSinDelta31[iAD]);
            
            s2t_pt = pow(10,(log10(lo_s2t) + double(is2t)*DeltaLog_s2t));
            //s2t_pt = lo_s2t + double(is2t)*DeltaLin_s2t;
            
            SurvP = SurvProb->Eval(s2t_pt);
            //IBDrate_osc[is2t] = (SurvP * noOsc_IBDrate_perday[iAD])*emuem[iAD] + totalBgd[iAD][0];
            IBDrate_osc[is2t] = (SurvP * noOsc_IBDrate_perday[iAD] + totalBgd[iAD][0])*emuem[iAD];
            //cout << IBDrate_osc[is2t] << "   " << IBDrate_data[iAD][0] << "  " << SurvP << "   " << noOsc_IBDrate_perday[iAD] << endl;
            
            double sqrerror = pow(IBDrate_data[iAD][1],2) + pow(totalBgd[iAD][1],2);
            
            chi2[is2t] += (pow((IBDrate_osc[is2t] - IBDrate_data[iAD][0]),2))/sqrerror;
            
        }//for iAD
        file << s2t_pt << "\t" << SurvP << "\t" << chi2[is2t] << endl;
        //cout << endl;
        //chi2 = 0.0;
    }//for is2t
    file << endl;

    file.close();

    double chi2_min = 5e+5;
    int isel;
    for (int jj = 0 ; jj < N_s2t ; jj++)
    {
        if (chi2[jj] < chi2_min)
        {
            isel = jj;
            chi2_min = chi2[isel];
        }
    }
    cout << "chi2_min = " << chi2_min << ", for jj = " << isel << endl;
    
    //---------------------------------------------------
    //---------------------------------------------------

    printf("******************************************************************************************\n");
    printf("* Declaration of arrays needed for db_minuit.C lines 66, 68 and 70. Copy and paste there. *\n");
    printf("******************************************************************************************\n\n");
    //-- Printing out for db_minuit.C
    printf("noOsc_IBDrate_perday[nAD] = {");
    for (int iAD = 0 ; iAD < nAD ; iAD++)
    {
        printf("%7.2f",noOsc_IBDrate_perday[iAD]);
        if(iAD < nAD-1)
            printf(",");
        else
            printf("};\n");
    }
    
    printf("double avgSinDelta21[nAD] = {");
    for (int iAD = 0 ; iAD < nAD ; iAD++)
    {
        printf("%12.9f",avgSinDelta21[iAD]);
        if(iAD < nAD-1)
            printf(",");
        else
            printf("};\n");
    }
    
    printf("double avgSinDelta31[nAD] = {");
    for (int iAD = 0 ; iAD < nAD ; iAD++)
    {
        printf("%12.9f",avgSinDelta31[iAD]);
        if(iAD < nAD-1)
            printf(",");
        else
            printf("};\n\n");
    }
    printf("******************************************************************************************\n");
    
    //---------------------------------------------------

    
    
    //---------------------------------------------------
    // Drawing section

    TH2F *frame_spectra = new TH2F("frame_spectra","",NB,lo,hi,10,0,18);
    frame_spectra->GetXaxis()->SetTitle("Prompt energy (MeV)");
    frame_spectra->GetYaxis()->SetTitle("Events/day (bkgd-subtracted)");

    TCanvas *canv0 = new TCanvas("canv0","canv0",775,500);
    
    TLegend *leg_spect = new TLegend(0.7,0.6,0.9,0.9);
    leg_spect->SetFillColor(0);
    for (int i = 0 ; i < 1 ; i++)
    {
        leg_spect->AddEntry(BFit_spect_histo[i],Form("AD%d",i+1));
        //leg_spect->AddEntry(wosc_spect_histo[i],Form("Spectra %d ",i));
    }
    
    
    //frame_spectra->Draw();
    BFit_spect_histo[0]->Draw("");
//    wosc_spect_histo[0]->Draw("same");
//    for (int i = 0 ; i < 4 ; i++)
//    {
//        wosc_spect_histo[i]->Draw("same");
//    }
    leg_spect->Draw();
    //break;
    
    //cout << endl << "Mean value: " << BFit_spect_histo[0]->GetMean() << endl;

    //canv0->Print("osc_test.pdf");

    //---------------------------------------------------
    // Drawing Survival Probabilities at the six ADs
    TH2F *frame_POscBF = new TH2F("frame_POscBF","",1000,0.89,1.01,10,-0.01,0.13);
    frame_POscBF->GetXaxis()->SetTitleFont(ft);
    frame_POscBF->GetXaxis()->SetTitleSize(sz);
    frame_POscBF->GetXaxis()->SetLabelFont(ft);
    frame_POscBF->GetXaxis()->SetLabelSize(sz);
    frame_POscBF->GetXaxis()->SetTitle("Survival Probability");
    frame_POscBF->GetYaxis()->SetTitleFont(ft);
    frame_POscBF->GetYaxis()->SetTitleSize(sz);
    frame_POscBF->GetYaxis()->SetLabelFont(ft);
    frame_POscBF->GetYaxis()->SetLabelSize(sz);
    //frame_POscBF->GetYaxis()->SetTitle("a.u.");
    
    TLegend *leg_avgs = new TLegend(0.1,0.6,0.35,0.9);
    leg_avgs->SetTextFont(ft);
    leg_avgs->SetTextSize(0.8*sz);
    leg_avgs->SetFillColor(0);
    for (int i = 0 ; i < nAD ; i++)
    {
        leg_avgs->AddEntry(Posc_AD_BF[i],Form("AD%d P_{avg} = %f",i+1,avgPosc_AD[i]));
    }
    
    TCanvas *canv1 = new TCanvas("canv1","canv1",775,500);
    frame_POscBF->Draw();
    for (int j = 0 ; j < nAD ; j++) {
        Posc_AD_BF[j]->Draw("same");
    }
    leg_avgs->Draw();
    
    canv1->Print("files_plots/POsc_avg.pdf");
    //---------------------------------------------------

} //end
