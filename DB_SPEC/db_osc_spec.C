//-- db_osc_spec.C  -  M.A.Acero O. - A.A.Alexis A. --//
//For the spectral analysis, using information from
//F.P. An et al., PRL 112 061801 (2014)
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
    // Open data file
    //TFile *fdata = new TFile("PRL112_data.root","READ");
    //---------------------------------------------------
    // define and get Data (energy) and Background histograms for the three (3) spectra (three EH)
    //const int nEH = 3;
    //TH1F *data_spect_histo[nEH];
    //TH1F *BkGd_spect_histo[nEH];
    //for (int i = 0 ; i < nEH ; i++)
    //{
    //    data_spect_histo[i] = (TH1F*) fdata->Get(Form("data_spect_histo_%d",i));//data from PRL112_data.root
    //    BkGd_spect_histo[i] = (TH1F*) fdata->Get(Form("BkGd_spect_histo_%d",i));//data from PRL112_data.root
    //}
    //---------------------------------------------------
    // Open ntuple file to read simulated data
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
    const int nAD = 6; //Number of Antineutrino detectors
    TH1F *nu_nosc_spect_histo[nAD];
    TH1F *BFit_spect_histo[nAD];
    TH1F *Posc_AD_BF[nAD];
    TH1F *Posc_AD_surv[nAD];
    for (int i = 0 ; i < nAD ; i++)
    {
        //no-oscillation Ep spectra - histograms
        nu_nosc_spect_histo[i] = new TH1F(Form("nu_nosc_spect_histo_%d",i),"",NB,xbins);// info from db-ntuple.root
        nu_nosc_spect_histo[i]->SetLineColor(i+3);
        //BF-oscillation Ep spectra - histograms
        BFit_spect_histo[i] = new TH1F(Form("BFit_spect_histo_%d",i),"",NB,xbins);// info from db-ntuple.root
        BFit_spect_histo[i]->SetLineColor(i+1);

        //Ocillation prpbability at BF - histograms
        Posc_AD_BF[i]       = new TH1F(Form("Posc_AD_BF_%d",i),"",1000,0,1);//to store <POsc(BF)>
        Posc_AD_BF[i]->SetLineColor(i+1);

        //Oscillation prpbability at (s2th,dm2) - histograms
        Posc_AD_surv[i]     = new TH1F(Form("Posc_AD_surv_%d",i),"",1000,0,1);//to store <POsc(s2th,dm2)>
    }
    //---------------------------------------------------
    //IBD rat (per day) (PRL 112 061801 (2014))
    //EH1(AD1, AD2),EH2(AD3),EH3(AD4, AD5, AD6)
    double IBDrate_perday[nAD][2] =
    {
        {653.30,2.31},{664.15,2.33},
        {581.97,2.33},
        { 73.31,0.66},{ 73.03,0.66},{72.20,0.66}
    };
  
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

    int sel;
    double TotNosc[nAD];
    double avgPosc_AD[nAD]; //<POsc(s2t_BF,dm2_31)>
    double noOsc_IBDrate_perday[nAD];
    double integ;
    for (int iAD = 0 ; iAD < nAD ; iAD++)
    {
        if (iAD < 3) sel = iAD;
        else if (iAD >= 3) sel = iAD+1;
        //------------------------------------------------
        //Filling Ocillation prpbability at BF - histogram
        T->Draw(Form("(1.0 - 0.089*((sin( 1.267 * 2.32e-3 * Ln/En ))**2) - ((cos(0.5 * asin(sqrt(0.089))))**4) * 0.861 * (sin( 1.267 * 7.59e-5 * Ln/En ))**2) >> Posc_AD_BF_%d",iAD),Form("id==%d",sel));
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
        //Filling the non-oscillated Ep spectra
        T->Draw(Form("Ep >> nu_nosc_spect_histo_%d",iAD),Form("(id==%d)",sel),"");
        TotNosc[iAD] =  nu_nosc_spect_histo[iAD]->Integral();
        //nu_nosc_spect_histo[iAD]->Scale(noOsc_IBDrate_perday[iAD]/TotNosc[iAD]); //normalize per day HERE? (2017-07-13)
        
        //condition to fill BF-oscillation- Ep spectra
        cutBF = Form("(1.0 - 0.089*((sin( 1.267 * 2.32e-3 * Ln/En ))**2) - ((cos(0.5 * asin(sqrt(0.089))))**4) * 0.861 * (sin( 1.267 * 7.59e-5 * Ln/En ))**2)*(id==%d)",iAD);

        //Filling and normalizing BF-oscillation Ep spectra
        T->Draw(Form("Ep >> BFit_spect_histo_%d",sel),cutBF,"");
        integ = BFit_spect_histo[iAD]->Integral();
        BFit_spect_histo[iAD]->Scale(1.0/integ);
        //------------------------------------------------
    }
    
    //---------------------------------------------------
    //Definition of the grid of oscillation parameters
    double s2t_pt, dm2_pt;
    //int   sel;

    const int     N_s2t = 200;
    const int     N_dm2 = 200;

    double       lo_s2t = 0.01;
    double       hi_s2t = 0.3;
    double DeltaLog_s2t = (log10(hi_s2t)-log10(lo_s2t))/double(N_s2t-1);

    double       lo_dm2 = 1e-4;
    double       hi_dm2 = 1e-2;
    double DeltaLog_dm2 = (log10(hi_dm2)-log10(lo_dm2))/double(N_dm2-1);
    
    TCut cut;
    
    const int dim = N_s2t*N_dm2;
    double TotWosc[dim];
    TH1F *wosc_spect_histo[dim];
    for (int i = 0 ; i < dim ; i++)
    {
        wosc_spect_histo[i] = new TH1F(Form("wosc_spect_histo_%d",i),"",NB,xbins);
        wosc_spect_histo[i]->SetLineWidth(2);
        wosc_spect_histo[i]->SetLineColor(i+2);
        wosc_spect_histo[i]->SetLineStyle(2);
    }


    //File to print the oscillation paramenter, the survival prob. and oscilated spectra
    //ofstream file;
    //string result = "db_SurvParams_V2.txt";
    //file.open (result.c_str());
    //file << setprecision(5);
    
    FILE *file;
    //file = fopen("files_data/db_gridOscSpectra_5M.txt","w");
    file = fopen("files_data/db_gridOscSpectra_test.txt","w");

    //write non-oscillated spectra for each AD to file
    for (int iAD = 0 ; iAD < nAD ; iAD++)
    {
       if      (iAD <  3) sel = iAD;
       else if (iAD >= 3) sel = iAD+1;

      //file << iAD + 1 << " " << s2t_pt << "\t" << dm2_pt;
       fprintf(file,"%d %8.2e %8.2e",iAD+1,s2t_pt,dm2_pt);
      //print bin-content of non-oscilated spectra per day
      for (int ib = 0 ; ib < 26 ; ib++)
      {
         double contNO   = nu_nosc_spect_histo[iAD]->GetBinContent(ib+1);
         //file << "\t" << contNO;
          fprintf(file," %10.2f",contNO);
      }
      //file << " \t" << TotNosc[iAD] /*<< "\t" << "1.00000"*/ << endl;
        fprintf(file," %10.2f\n",TotNosc[iAD]);
    } // for iAD
    //file << endl;
    fprintf(file,"\n");


    for (int is2t = 0 ; is2t < N_s2t ; is2t++)
    {
        s2t_pt = pow(10,log10(lo_s2t) + (double(is2t)*DeltaLog_s2t));
    
        for (int idm2 = 0 ; idm2 < N_dm2 ; idm2++)
        {
            dm2_pt = pow(10,log10(lo_dm2) + (double(idm2)*DeltaLog_dm2));

            for (int iAD = 0 ; iAD < nAD ; iAD++)
            {
                if (iAD < 3) sel = iAD;
                else if (iAD >= 3) sel = iAD +1;
                    
                // Condition to fill oscilated spectra for (s2t_pt,dm2_pt), i.e. wosc_spect_histo[iAD]
                cut = Form("(1.0 - %e*((sin( 1.267 * %e * Ln/En ))**2) - (0.25*(1 + sqrt(1 - %e))**2) * 0.861 * (sin( 1.267 * 7.59e-5 * Ln/En ))**2)*(id==%d)" ,s2t_pt,dm2_pt,s2t_pt,sel);
                // Filling oscilated spectra for (s2t_pt,dm2_pt)
                int ih = is2t*N_dm2 + idm2;
                //cout << "AD " << iAD + 1 << "\t s2t_pt =  " << s2t_pt << "\t dm2_pt =  " << dm2_pt  << "\t ih = " << ih << endl;
                //cout << "cut computed" << endl;
                //cout << "Test No. " << idm2 << endl;
                T->Draw(Form("Ep >> wosc_spect_histo_%d",ih),cut,"");
                TotWosc[ih] =  wosc_spect_histo[ih]->Integral();

                // Survival probability for (s2t_pt,dm2_pt) to fill histogram at each detector, i.e. Posc_AD_surv[iAD]
                //T->Draw(Form("(1.0 - %e*((sin( 1.267 * %e * Ln/En ))**2) - (0.25*(1 + sqrt(1 - %e))**2) * 0.861 * (sin( 1.267 * 7.59e-5 * Ln/En ))**2) >> Posc_AD_surv_%d",s2t_pt,dm2_pt,s2t_pt,iAD),Form("id==%d",sel));
                //Average survival probability for (s2t_pt,dm2_pt) at each detector
                //SurvP = Posc_AD_surv[iAD]->GetMean();

                //file << iAD + 1 << " " << s2t_pt << "\t" << dm2_pt;
                fprintf(file,"%d %8.2e %8.2e",iAD+1,s2t_pt,dm2_pt);
                //Printing bin-content for the oscilated spectra for (s2t_pt,dm2_pt)
                for (int ib = 0 ; ib < 26 ; ib++)
                {
                    double cont   = wosc_spect_histo[ih]->GetBinContent(ib+1);
                    //file << "\t" << cont;
                    fprintf(file," %10.2f",cont);
                }
                //file << " \t" << TotWosc[ih] << endl;
                fprintf(file," %10.2f\n",TotWosc[ih]);
                
                //Printing check-points info
                if (ih%10 == 0)
                {
                    cout << ih << "  Done with detector " << iAD+1 << " for " << s2t_pt << "\t" << dm2_pt << endl;
                    //cout << " \t" << TotNosc[iAD] << endl;
                }

                //Normalizing the oscilated spectra for (s2t_pt,dm2_pt)
                //*CHECK if this is nneded .. guess is NO 
                integ = wosc_spect_histo[ih]->Integral();
                wosc_spect_histo[ih]->Scale(1.0/integ);
            } // for idm2
            //file << endl;
            fprintf(file,"\n");
        }//for is2t
        //cout << "  Done with detector " << iAD+1 << endl;
        //file << endl;
        fprintf(file,"\n");
    }//for iAD

    //file.close();
    fclose(file);
    
    //break;
    //---------------------------------------------------
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

} //end
