//-- db_osc_spec.C  -  M.A.Acero O. - A.A.Alexis A. - 08.11.2019 --//
//For the 6AD+8AD spectral analysis, using information from
//F.P. An et al., PRD 95 072006 (2017)
#include "constants.h"
#include <iostream>
#include <fstream>
#include <string>

void db_osc_spec()
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
    // Open  file to read simulated data
    //TFile *fntuple = new TFile("files_data/db-ntuple_0.7-1.3_noER_100M.root","READ"); //-No ERF, 1st bin 0.7-1.3 MeV
    TFile *fntuple = new TFile("files_data/db-ntuple_noER_1M.root","READ");         //-No ERF, 1st bin mel-1.3 MeV
    //TFile *fntuple = new TFile("files_data/db-ntuple_unPhys_noER_100M.root","READ");  //-No ERF, 1st bin w/o restrictions
    TTree *T6AD8AD = (TTree*)fntuple->Get("T6AD8AD"); //-- Events for the 6AD+8AD period
    TCut cutBF_6AD8AD; //-- For the 1013-days (6AD+8AD) period
    TCut cutBF_6AD;    //-- For the  217-days (6AD) period
    //---------------------------------------------------
    // histogram binning for the simulated data
    //double    NB = 35;
    //double    lo = 0.7;
    //double    hi = 12.0;
    double xbins[36];
    xbins[0] = lo;
    double delta_bins2 = (7.9 - 1.3)/33; // 0.20 MeV/bin
    for (int i = 0 ; i < (NB-1) ; i++)
    {
        xbins[i+1] = 1.3 + delta_bins2*i;
    }
    xbins[35] = hi;
    //---------------------------------------------------
    //const int nAD = 8; //Number of Antineutrino detectors
    TH1F *nu_nosc_spect_histo_1013[nAD];
    TH1F *BFit_spect_histo_1013[nAD];
    TH1F *Posc_AD_BF_1230[nAD];
    TH1F *Posc_AD_BF_217[nAD];
    TH1F *Posc_AD_surv_1230[nAD];
    TH1F *nu_nosc_spect_histo_217[nAD];
    TH1F *BFit_spect_histo_217[nAD];
    for (int i = 0 ; i < nAD ; i++)
    {
        //no-oscillation Ep spectra - histograms
        nu_nosc_spect_histo_1013[i] = new TH1F(Form("nu_nosc_spect_histo_1013_%d",i),"",NB,xbins);// info from db-ntuple.root
        nu_nosc_spect_histo_1013[i]->SetLineColor(i+3);
        nu_nosc_spect_histo_217[i] = new TH1F(Form("nu_nosc_spect_histo_217_%d",i),"",NB,xbins);// info from db-ntuple.root
        nu_nosc_spect_histo_217[i]->SetLineColor(i+4);
        //BF-oscillation Ep spectra - histograms
        BFit_spect_histo_1013[i] = new TH1F(Form("BFit_spect_histo_1013_%d",i),"",NB,xbins);// info from db-ntuple.root
        BFit_spect_histo_1013[i]->SetLineColor(i+1);
        BFit_spect_histo_217[i] = new TH1F(Form("BFit_spect_histo_217_%d",i),"",NB,xbins);// info from db-ntuple.root
        BFit_spect_histo_217[i]->SetLineColor(i+1);

        //Ocillation prpbability at BF - histograms
        Posc_AD_BF_1230[i]       = new TH1F(Form("Posc_AD_BF_1230_%d",i),"",1000,0,1);//to store <POsc(BF)>
        Posc_AD_BF_1230[i]->SetLineColor(i+1);
        Posc_AD_BF_217[i]       = new TH1F(Form("Posc_AD_BF_217_%d",i),"",1000,0,1);//to store <POsc(BF)> for 217 days
        Posc_AD_BF_217[i]->SetLineColor(i+1);

        //Oscillation prpbability at (s2th,dm2) - histograms
        Posc_AD_surv_1230[i]     = new TH1F(Form("Posc_AD_surv_1230_%d",i),"",1000,0,1);//to store <POsc(s2th,dm2)>
    }
    //---------------------------------------------------
    //IBD rat (per day) (PRD 95 072006 (2017))
    //EH1(AD1, AD2),EH2(AD3,AD8),EH3(AD4, AD5, AD6, AD7)
    double IBDrate_perday_1230[nAD][2] =
    {
        {653.03,1.37},{665.42,1.38},
        {599.71,1.12},{593.82,1.18},
        { 74.25,0.28},{ 74.60,0.28},{73.98,0.28},{74.73,0.30}
    };
    IBDrate_perday_1230[4][0] = 1.0*IBDrate_perday_1230[4][0];
    IBDrate_perday_1230[5][0] = 1.0*IBDrate_perday_1230[5][0];
    IBDrate_perday_1230[6][0] = 1.0*IBDrate_perday_1230[6][0];
    IBDrate_perday_1230[7][0] = 1.0*IBDrate_perday_1230[7][0];
    //IBD rat (per day) (PRL 112 061801 (2014))
    //EH1(AD1, AD2),EH2(AD3),EH3(AD4, AD5, AD6)
    double IBDrate_perday_217[nAD][2] =
      {
        {653.30,2.31},{664.15,2.33},
        {581.97,2.07},{0.0,0.0},
        { 73.31,0.66},{ 73.03,0.66},{72.20,0.66},{0.0,0.0}
      };
    IBDrate_perday_217[4][0] = 1.0*IBDrate_perday_217[4][0];
    IBDrate_perday_217[5][0] = 1.0*IBDrate_perday_217[5][0];
    IBDrate_perday_217[6][0] = 1.0*IBDrate_perday_217[6][0];

    //Computing <POsc(s2t_BF,dm2_31)> for each AD
    // AD1 -> id = 0; AD2 -> id = 1; AD3 -> id = 2; AD4 -> id = 4; AD5 -> id = 5; AD6 -> id = 6

      /* Used data
       s22th13_BF = 0.089;                        //PRL 112 061801 (2014)
       th13 = 0.5*asin(sqrt(s22th13_BF));
       c4th13 = (cos(th13))**4;
       dm2_31 = 2.32e-3; //eV^2,                  //PRL 108 171803 (2012)
       const double dm2_21 = 7.59e-5; //eV^2,     //PRL 108 171803 (2012)
       const double s22th12 = 0.861;              //PRL 108 171803 (2012)
       */

    FILE *file_IBDrates8AD;
    file_IBDrates8AD = fopen("files_data/db_noOsc_IBDrates_perday_1230.txt","w");
    FILE *file_IBDrates6AD;
    file_IBDrates6AD = fopen("files_data/db_noOsc_IBDrates_perday_217.txt","w");

    int sel;
    double TotNosc_1013[nAD];
    double avgPosc_AD_1230[nAD]; //<POsc(s2t_BF,dm2_31)>
    double noOsc_IBDrate_perday_1230[nAD];
    double integ_1230;
    double integ_1013;

    double TotNosc_217[nAD];
    double avgPosc_AD_217[nAD]; //<POsc(s2t_BF,dm2_31)> for the 217 days
    double noOsc_IBDrate_perday_217[nAD];
    double integ_217;
    for (int iAD = 0 ; iAD < nAD ; iAD++)
    {
        sel = iAD;
        //------------------------------------------------
        //Filling Ocillation prpbability at BF - histogram
        T6AD8AD->Draw(Form("(1.0 - 0.0841*((sin( 1.267 * 2.50e-3 * Ln/En ))**2) - ((cos(0.5 * asin(sqrt(0.0841))))**4) * 0.846 * (sin( 1.267 * 7.53e-5 * Ln/En ))**2) >> Posc_AD_BF_1230_%d",iAD),Form("id==%d",sel));
        integ_1230 = Posc_AD_BF_1230[iAD]->Integral();
        Posc_AD_BF_1230[iAD]->Scale(1.0/integ_1230); //Used to plot the Survival Probabilities for each AD with BF parameters (normalized)
        //Average oscillation Probability
        avgPosc_AD_1230[iAD] = Posc_AD_BF_1230[iAD]->GetMean();
        //No-oscillation IDB rate (per day)
        noOsc_IBDrate_perday_1230[iAD] = IBDrate_perday_1230[iAD][0]/avgPosc_AD_1230[iAD];
        //Printing results
        cout << "(avgPosc_AD,noOsc_IBDrate_perday_1230)_" << sel << " = (" << avgPosc_AD_1230[iAD]
        << ", " << noOsc_IBDrate_perday_1230[iAD] << ") " << endl;
        fprintf(file_IBDrates8AD,"%f \n", noOsc_IBDrate_perday_1230[iAD]);
        //------------------------------------------------
        //Filling Ocillation prpbability at BF - histogram (217 days)
        //--NOTE: The BF oscillation parameters are from the analysis of 217 days (ssq2th13= 0.090, dmsqee=2.59e-3 eV^2)
        if (iAD != 3 && iAD != 7){
	  T6AD8AD->Draw(Form("(1.0 - 0.090*((sin( 1.267 * 2.59e-3 * Ln/En ))**2) - ((cos(0.5 * asin(sqrt(0.090))))**4) * 0.857 * (sin( 1.267 * 7.50e-5 * Ln/En ))**2) >> Posc_AD_BF_217_%d",iAD),Form("id==%d && per==1",sel));
	  integ_217 = Posc_AD_BF_217[iAD]->Integral();
	  Posc_AD_BF_217[iAD]->Scale(1.0/integ_217); //Used to plot the Survival Probabilities for each AD with BF parameters (normalized)
	  //Average oscillation Probability
	  avgPosc_AD_217[iAD] = Posc_AD_BF_217[iAD]->GetMean();
	  //No-oscillation IDB rate (per day)
	  noOsc_IBDrate_perday_217[iAD] = IBDrate_perday_217[iAD][0]/avgPosc_AD_217[iAD];
        }
        else{
	  avgPosc_AD_217[iAD] = 0.0;
	  noOsc_IBDrate_perday_217[iAD] = 0.0;
        }
        //Printing results
        cout << "(avgPosc_AD,noOsc_IBDrate_perday_217)_" << sel << " = (" << avgPosc_AD_217[iAD]
	     << ", " << noOsc_IBDrate_perday_217[iAD] << ") " << endl;
        fprintf(file_IBDrates6AD,"%f \n", noOsc_IBDrate_perday_217[iAD]);
        //------------------------------------------------

        //------------------------------------------------
        //Filling the non-oscillated Ep spectra
        T6AD8AD->Draw(Form("Ep >> nu_nosc_spect_histo_1013_%d",iAD),Form("id==%d && per==2",sel),"");
        TotNosc_1013[iAD] =  nu_nosc_spect_histo_1013[iAD]->Integral();
        T6AD8AD->Draw(Form("Ep >> nu_nosc_spect_histo_217_%d",iAD),Form("id==%d && per==1",sel),"");
        TotNosc_217[iAD] =  nu_nosc_spect_histo_217[iAD]->Integral();
        //nu_nosc_spect_histo[iAD]->Scale(noOsc_IBDrate_perday[iAD]/TotNosc[iAD]); //normalize per day HERE? (2017-07-13)
        
        //condition to fill BF-oscillation- Ep spectra
        //cutBF = Form("(1.0 - 0.0841*((sin( 1.267 * 2.50e-3 * Ln/En ))**2) - ((cos(0.5 * asin(sqrt(0.0841))))**4) * 0.841 * (sin( 1.267 * 7.59e-5 * Ln/En ))**2)*(id==%d)",iAD);
        cutBF_6AD8AD = Form("(1.0 - 0.0841*((sin( 1.267 * 2.50e-3 * Ln/En ))**2) - ((cos(0.5 * asin(sqrt(0.0841))))**4) * 0.846 * (sin( 1.267 * 7.53e-5 * Ln/En ))**2)*(id==%d && per==2)",iAD);
        cutBF_6AD = Form("(1.0 - 0.0841*((sin( 1.267 * 2.50e-3 * Ln/En ))**2) - ((cos(0.5 * asin(sqrt(0.0841))))**4) * 0.846 * (sin( 1.267 * 7.53e-5 * Ln/En ))**2)*(id==%d && per==1)",iAD);

        //Filling and normalizing BF-oscillation Ep spectra
        T6AD8AD->Draw(Form("Ep >> BFit_spect_histo_1013_%d",sel),cutBF_6AD8AD,"");
        integ_1013 = BFit_spect_histo_1013[iAD]->Integral();
        BFit_spect_histo_1013[iAD]->Scale(1.0/integ_1013);
        T6AD8AD->Draw(Form("Ep >> BFit_spect_histo_217_%d",sel),cutBF_6AD,"");
        integ_217 = BFit_spect_histo_217[iAD]->Integral();
        BFit_spect_histo_217[iAD]->Scale(1.0/integ_217);
        //------------------------------------------------
    }
    fclose(file_IBDrates8AD);
    fclose(file_IBDrates6AD);
    

    //---------------------------------------------------
    //Definition of the grid of oscillation parameters
    double s2t_pt, dm2_pt;
    //int   sel;
    //Comment out when using the shel script!!!
    //const int     N_s2t = 5;
    //const int     N_dm2 = 5;
    //Comment out when using the shel script!!!
    //double       lo_s2t = 0.01;
    //double       hi_s2t = 0.3;
    //double DeltaLog_s2t = (log10(hi_s2t)-log10(lo_s2t))/double(N_s2t-1);
    double Delta_s2t = (hi_s2t - lo_s2t)/double(N_s2t-1);
    //Comment out when using the shel script!!!
    //double       lo_dm2 = 1e-4;
    //double       hi_dm2 = 1e-2;
    //double DeltaLog_dm2 = (log10(hi_dm2)-log10(lo_dm2))/double(N_dm2-1);
    double Delta_dm2 = (hi_dm2 - lo_dm2)/double(N_dm2-1);
    
    printf("\n");
    printf("Grid definition db_osc_spec.C\n");
    printf("---------------------------\n");
    printf("N_s2t: %d\n",N_s2t);
    printf("lo_s2t: %f\n",lo_s2t);
    printf("hi_s2t: %f\n",hi_s2t);
    printf("\n");
    printf("N_dm2: %d\n",N_dm2);
    printf("lo_dm2: %f\n",lo_dm2);
    printf("hi_dm2: %f\n",hi_dm2);
    printf("---------------------------\n");
    printf("\n");

    TCut posc_wgt1,posc_wgt2;
    
    const int dim = N_s2t*N_dm2;
    double TotWosc_1013[dim];
    double TotWosc_217[dim];
    TH1F *wosc_spect_histo_1013[dim];
    TH1F *wosc_spect_histo_217[dim];
    for (int i = 0 ; i < dim ; i++)
    {
        wosc_spect_histo_1013[i] = new TH1F(Form("wosc_spect_histo_1013_%d",i),"",NB,xbins);
        wosc_spect_histo_1013[i]->SetLineWidth(2);
        wosc_spect_histo_1013[i]->SetLineColor(i+2);
        wosc_spect_histo_1013[i]->SetLineStyle(2);

        wosc_spect_histo_217[i] = new TH1F(Form("wosc_spect_histo_217_%d",i),"",NB,xbins);
        wosc_spect_histo_217[i]->SetLineWidth(3);
        wosc_spect_histo_217[i]->SetLineColor(i+3);
        wosc_spect_histo_217[i]->SetLineStyle(3);
    }

    //File to print the oscillation paramenter, the survival prob. and oscilated spectra
    ofstream file;
    string grid_name = "files_data/db_gridOscSpectra_1230.txt";
    file.open((grid_name).c_str());
    file << fixed;
    file << setprecision(6);

    s2t_pt = 0.0;
    dm2_pt = 0.0;
    //write non-oscillated spectra for each AD to file
    for (int iAD = 0 ; iAD < nAD ; iAD++)
    {
        sel = iAD;
       //---------------------------------------------------------

        file << iAD + 1 << "\t" << s2t_pt << "\t" << dm2_pt;
        //print bin-content of non-oscilated spectra per day
        for (int ib = 0 ; ib < NB ; ib++)
        {
            double contNO_217    = nu_nosc_spect_histo_217[iAD] ->GetBinContent(ib+1);
            if (iAD==3 || iAD==7) contNO_217 = 0.0; //-- Zero content for detector 3 and 7 for the 6AD run
            file << "\t" << contNO_217;
            //fprintf(file, " %10.2f",contNO_217);
        }
        file << " \t" << TotNosc_217[iAD];
        for (int ib = 0 ; ib < NB ; ib++)
        {
            double contNO_1013   = nu_nosc_spect_histo_1013[iAD]->GetBinContent(ib+1);
                    file << "\t" << contNO_1013;
                }
                file << " \t" << TotNosc_1013[iAD] << endl;
            } // for iAD
            file << endl;

    for (int is2t = 0 ; is2t < N_s2t ; is2t++)
    {
        //s2t_pt = pow(10,log10(lo_s2t) + (double(is2t)*DeltaLog_s2t));
        s2t_pt = lo_s2t + (double(is2t)*Delta_s2t);
    
        for (int idm2 = 0 ; idm2 < N_dm2 ; idm2++)
        {
            dm2_pt = lo_dm2 + (double(idm2)*Delta_dm2);

            for (int iAD = 0 ; iAD < nAD ; iAD++)
            {
                //if (iAD < 3) sel = iAD;// this is not necessary here as there are 8 AD
                //else if (iAD >= 3) sel = iAD +1;
                sel = iAD;    
                // Condition to fill oscilated spectra for (s2t_pt,dm2_pt), i.e. wosc_spect_histo[iAD]
                posc_wgt1 = Form("(1.0 - %e*((sin( 1.267 * %e * Ln/En ))**2) - (0.25*(1 + sqrt(1 - %e))**2) * 0.846 * (sin( 1.267 * 7.53e-5 * Ln/En ))**2)*(id==%d && per==1)" ,s2t_pt,dm2_pt,s2t_pt,sel);
                posc_wgt2 = Form("(1.0 - %e*((sin( 1.267 * %e * Ln/En ))**2) - (0.25*(1 + sqrt(1 - %e))**2) * 0.846 * (sin( 1.267 * 7.53e-5 * Ln/En ))**2)*(id==%d && per==2)" ,s2t_pt,dm2_pt,s2t_pt,sel);
                // Filling oscilated spectra for (s2t_pt,dm2_pt)
                int ih = is2t*N_dm2 + idm2;
                //cout << "AD " << iAD + 1 << "\t s2t_pt =  " << s2t_pt << "\t dm2_pt =  " << dm2_pt  << "\t ih = " << ih << endl;
                //cout << "cut computed" << endl;
                //cout << "Test No. " << idm2 << endl;
                T6AD8AD->Draw(Form("Ep >> wosc_spect_histo_1013_%d",ih),posc_wgt2,"");
                TotWosc_1013[ih] =  wosc_spect_histo_1013[ih]->Integral();
                T6AD8AD->Draw(Form("Ep >> wosc_spect_histo_217_%d", ih),posc_wgt1,"");
                TotWosc_217[ih] =   wosc_spect_histo_217[ih]->Integral();

                file << iAD + 1 << " " << s2t_pt << "\t" << dm2_pt;
                //Printing bin-content for the oscilated spectra for (s2t_pt,dm2_pt)
                for (int ib = 0 ; ib < NB ; ib++)
                {
                    double cont217    = wosc_spect_histo_217[ih] ->GetBinContent(ib+1);
                    file << "\t" << cont217;
                }
                file << " \t" << TotWosc_217[ih];
                for (int ib = 0 ; ib < NB ; ib++)
                {
                    double cont1013   = wosc_spect_histo_1013[ih]->GetBinContent(ib+1);
                    file << "\t" << cont1013;
                }
                file << " \t" << TotWosc_1013[ih] << endl;

                //Printing check-points info
                if (ih%10 == 0)
                {
                    cout << ih << "  Done with detector " << iAD+1 << " for " << s2t_pt << "\t" << dm2_pt << endl;
                    //cout << " \t" << TotNosc[iAD] << endl;
                }

                //Normalizing the oscilated spectra for (s2t_pt,dm2_pt)
                integ_1013 = wosc_spect_histo_1013[ih]->Integral();
                wosc_spect_histo_1013[ih]->Scale(1.0/integ_1013);
                integ_217  = wosc_spect_histo_217[ih]->Integral();
                wosc_spect_histo_217[ih]->Scale(1.0/integ_217);
            } // for idm2
            file << endl;
        }//for is2t
        //cout << "  Done with detector " << iAD+1 << endl;
        file << endl;
    }//for iAD

    file.close();
    //break;
    //---------------------------------------------------
    //---------------------------------------------------
    
/*
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
        leg_spect->AddEntry(nu_nosc_spect_histo[i],Form("AD%d",i+1));
    }
    for (int i = 0 ; i < 4 ; i++)
    {
        leg_spect->AddEntry(wosc_spect_histo[i],Form("Spectra %d ",i));
    }
    
    
    //frame_spectra->Draw();
    BFit_spect_histo[0]->Draw("");
    nu_nosc_spect_histo[0]->Draw("same");
//    wosc_spect_histo[0]->Draw("same");
    for (int i = 0 ; i < 4 ; i++)
    {
        wosc_spect_histo[i]->Draw("same");
    }
    leg_spect->Draw();
    //break;
    
    cout << endl << "Mean value: " << BFit_spect_histo[0]->GetMean() << endl;

    //canv0->Print("files_plots/osc_test.pdf");

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
 */

} //end
