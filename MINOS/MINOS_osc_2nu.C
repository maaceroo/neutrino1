//------------------------------------------------------------------//
//--       MINOS_osc_2nu.C  -  M.A.Acero O. - A.A.Alexis A.       --//
//For the MINOS (anti)numu data analysis, using information from  --//
// - MINOS Col. PRL 110, 251801 (2013) - (Anti)NuMu disappearance --//
//------------------------------------------------------------------//
#include "constants.h"
#include <iostream>
#include <fstream>
#include <string>

void MINOS_osc_2nu()
{ //begin

    TString filePath = dirName;

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
    TFile *fntuple = new TFile(filePath + "/data/minos_ntuple.root","READ");
    TTree *Tnumu  = (TTree*)fntuple->Get("Tnumu");
    TCut cutBF_numu;
    TTree *Tnumub = (TTree*)fntuple->Get("Tnumub");
    TCut cutBF_numub;
    //---------------------------------------------------
    //Neutrinos
    double xbins_numu110[NB_numu110+1];
    double delta_bins110 = (9.0 - 0.5)/17.0; // 0.5 GeV/bin
    for (int i = 0 ; i < (NB_numu110 - 5) ; i++)
        xbins_numu110[i] = lo110 + i*delta_bins110;
    xbins_numu110[18] = xbins_numu110[17] + 0.75;
    xbins_numu110[19] = xbins_numu110[18] + 0.75;
    xbins_numu110[20] = xbins_numu110[19] + 0.75;
    xbins_numu110[21] = xbins_numu110[20] + 0.75;
    xbins_numu110[22] = xbins_numu110[21] + 1.0;
    xbins_numu110[23] = hi110;
    //AntiNeutrinos
    double xbins_nub110[NB_numubar110+1];
    double delta_binsb110 = (12.0 - lob110)/(NB_numubar110-1); // 1.0 GeV/bin
    for (int i = 0 ; i < NB_numubar110 ; i++)
        xbins_nub110[i] = lob110 + delta_binsb110*i;
    xbins_nub110[12] = xbins_nub110[11] + 2.0;
    //---------------------------------------------------
    TH1F *numu_nosc_spect_histo;
    TH1F *numu_BFit_spect_histo;
    TH1F *numu_Posc_BF;
    TH1F *numu_Posc_surv;
    TH1F *numub_nosc_spect_histo;
    TH1F *numub_BFit_spect_histo;
    TH1F *numub_Posc_BF;
    TH1F *numub_Posc_surv;
    //no-oscillation RecoE spectra - histograms
    numu_nosc_spect_histo = new TH1F("numu_nosc_spect_histo","",NB_numu110,xbins_numu110);
    numu_nosc_spect_histo->SetLineColor(3);
    numub_nosc_spect_histo = new TH1F("numub_nosc_spect_histo","",NB_numubar110,xbins_nub110);
    numub_nosc_spect_histo->SetLineColor(3);
    //BF-oscillation Ep spectra - histograms
    numu_BFit_spect_histo = new TH1F("numu_BFit_spect_histo","",NB_numu110,xbins_numu110);
    numu_BFit_spect_histo->SetLineColor(1);
    numub_BFit_spect_histo = new TH1F("numub_BFit_spect_histo","",NB_numubar110,xbins_nub110);
    numub_BFit_spect_histo->SetLineColor(1);
    //Ocillation prpbability at BF - histograms
    numu_Posc_BF          = new TH1F("numu_Posc_BF","",1000,0,1);   //to store <POsc(BF)>
    numu_Posc_BF->SetLineColor(1);
    numub_Posc_BF         = new TH1F("numub_Posc_BF","",1000,0,1);   //to store <POsc(BF)>
    numub_Posc_BF->SetLineColor(1);
    //Oscillation prpbability at (s2th,dm2) - histograms
    numu_Posc_surv        = new TH1F("numu_Posc_surv","",1000,0,1); //to store <POsc(s2th,dm2)>
    numub_Posc_surv       = new TH1F("numub_Posc_surv","",1000,0,1); //to store <POsc(s2th,dm2)>

    //---------------------------------------------------
    double numu_TotNosc;
    double numu_avgPosc; //<POsc(s2t_BF,dm2_31)>
    double numu_integ;
    double numub_TotNosc;
    double numub_avgPosc; //<POsc(s2t_BF,dm2_31)>
    double numub_integ;
    //------------------------------------------------
    //Filling Ocillation prpbability at BF - neutrinos
    Tnumu->Draw("1.0 - 0.950*((sin( 1.267 * 2.41e-3 * 735.0/Etrue ))**2) >> numu_Posc_BF");
    numu_integ = numu_Posc_BF->Integral();
    numu_Posc_BF->Scale(1.0/numu_integ); //Used to plot the Survival Probabilities with BF parameters (normalized)
    //Average oscillation Probability
    numu_avgPosc = numu_Posc_BF->GetMean();
    //Filling the non-oscillated RecoE spectra
    Tnumu->Draw("Ereco >> numu_nosc_spect_histo");
    numu_TotNosc = numu_nosc_spect_histo->Integral();
    //numu_nosc_spect_histo->Scale(1.0/numu_TotNosc);
    //------------------------------------------------
    //Filling Ocillation prpbability at BF - antineutrinos
    Tnumub->Draw("1.0 - 0.97*((sin( 1.267 * 2.50e-3 * 735.0/Etrueb ))**2) >> numub_Posc_BF");
    numub_integ = numub_Posc_BF->Integral();
    numub_Posc_BF->Scale(1.0/numub_integ); //Used to plot the Survival Probabilities with BF parameters (normalized)
    //Average oscillation Probability
    numub_avgPosc = numub_Posc_BF->GetMean();
    //------------------------------------------------
    //Filling the non-oscillated RecoE spectra
    Tnumub->Draw("Erecob >> numub_nosc_spect_histo");
    numub_TotNosc = numub_nosc_spect_histo->Integral();

    //condition to fill BF-oscillation- Ereco spectra
    cutBF_numu  = "1.0 - 0.95*((sin( 1.267 * 2.41e-3 * 735.0/Etrue ))**2)";
    cutBF_numub = "1.0 - 0.97*((sin( 1.267 * 2.50e-3 * 735.0/Etrueb ))**2)";

    //Filling and normalizing BF-oscillation Ereco spectra
    Tnumu ->Draw("Ereco >>  numu_BFit_spect_histo", cutBF_numu, "");
    numu_integ  = numu_BFit_spect_histo-> Integral();
    Tnumub->Draw("Erecob >> numub_BFit_spect_histo",cutBF_numub,"");
    numub_integ = numub_BFit_spect_histo->Integral();

    //---------------------------------------------------
    //Definition of the grid of oscillation parameters
    double s2t_pt, dm2_pt;

//    const int     N_s2t = 5;
//    const int     N_dm2 = 5;
//
//    double       lo_s2t = 0.75;
//    double       hi_s2t = 1.0;
//    double DeltaLog_s2t = (log10(hi_s2t)-log10(lo_s2t))/double(N_s2t-1);
    double Delta_s2t = (hi_s2t - lo_s2t)/double(N_s2t-1);

//    double       lo_dm2 = 1e-3;
//    double       hi_dm2 = 4e-3;
//    double DeltaLog_dm2 = (log10(hi_dm2)-log10(lo_dm2))/double(N_dm2-1);
    double Delta_dm2 = (hi_dm2 - lo_dm2)/double(N_dm2-1);


    printf("\n");
    printf("Grid definition MINOS_osc_2nu.C\n");
    printf("-------------------------------\n");
    printf("N_s2t: %d\n",N_s2t);
    printf("lo_s2t: %f\n",lo_s2t);
    printf("hi_s2t: %f\n",hi_s2t);
    printf("\n");
    printf("N_dm2: %d\n",N_dm2);
    printf("lo_dm2: %f\n",lo_dm2);
    printf("hi_dm2: %f\n",hi_dm2);
    printf("-------------------------------\n");
    printf("\n");

    TCut cut_numu;
    TCut cut_numub;

    const int dim = N_s2t*N_dm2;
    double numu_TotWosc[dim];
    TH1F *numu_wosc_spect_histo[dim];
    double numub_TotWosc[dim];
    TH1F *numub_wosc_spect_histo[dim];
    for (int i = 0 ; i < dim ; i++)
    {
        numu_wosc_spect_histo[i] = new TH1F(Form("numu_wosc_spect_histo_%d",i),"",NB_numu110,xbins_numu110);
        numu_wosc_spect_histo[i]->SetLineWidth(2);
        numu_wosc_spect_histo[i]->SetLineColor(i+2);
        numu_wosc_spect_histo[i]->SetLineStyle(2);

        numub_wosc_spect_histo[i] = new TH1F(Form("numub_wosc_spect_histo_%d",i),"",NB_numubar110,xbins_nub110);
        numub_wosc_spect_histo[i]->SetLineWidth(2);
        numub_wosc_spect_histo[i]->SetLineColor(i+2);
        numub_wosc_spect_histo[i]->SetLineStyle(2);
    }

    ofstream numu_file;
    string numu_grid_name = "/data/minosNuMu_gridOscSpectra.txt";
    numu_file.open(filePath + (numu_grid_name).c_str());
    numu_file << fixed;
    numu_file << setprecision(6);

    ofstream numub_file;
    string numub_grid_name = "data/minosNuMuB_gridOscSpectra.txt";
    numub_file.open(filePath + (numub_grid_name).c_str());
    numub_file << fixed;
    numub_file << setprecision(6);

    s2t_pt = 0.0;
    dm2_pt = 0.0;
    //write non-oscillated spectra
    numu_file << s2t_pt << "\t" << dm2_pt;
    numub_file << s2t_pt << "\t" << dm2_pt;
    //print bin-content of non-oscilated spectra per day
    for (int ib = 0 ; ib < 23 ; ib++)
    {
        double contNO   = numu_nosc_spect_histo->GetBinContent(ib+1);
        numu_file << "\t" << contNO;
        if (ib < 12) {
            contNO = numub_nosc_spect_histo->GetBinContent(ib+1);
            numub_file << "\t" << contNO;
        }
    }
    numu_file  << " \t" << numu_TotNosc  << endl;
    numu_file  << endl;
    numub_file << " \t" << numub_TotNosc << endl;
    numub_file << endl;

    for (int is2t = 0 ; is2t < N_s2t ; is2t++)
    {
        //s2t_pt = pow(10,log10(lo_s2t) + (double(is2t)*DeltaLog_s2t));
        s2t_pt = lo_s2t + (double(is2t)*Delta_s2t);
        for (int idm2 = 0 ; idm2 < N_dm2 ; idm2++)
        {
            dm2_pt = lo_dm2 + (double(idm2)*Delta_dm2);
            // Condition to fill oscilated spectra for (s2t_pt,dm2_pt), i.e. wosc_spect_histo
            cut_numu  = Form("1.0 - %e*((sin( 1.267 * %e * 735.0/Etrue ))**2)" ,s2t_pt,dm2_pt);
            cut_numub = Form("1.0 - %e*((sin( 1.267 * %e * 735.0/Etrueb ))**2)" ,s2t_pt,dm2_pt);
            // Filling oscilated spectra for (s2t_pt,dm2_pt)
            int ih = is2t*N_dm2 + idm2;
            //cout "\t s2t_pt =  " << s2t_pt << "\t dm2_pt =  " << dm2_pt  << "\t ih = " << ih << endl;
            //cout << "cut computed" << endl;
            //cout << "Test No. " << idm2 << endl;
            Tnumu->Draw(Form("Ereco >> numu_wosc_spect_histo_%d",ih),cut_numu,"");
            numu_TotWosc[ih] =  numu_wosc_spect_histo[ih]->Integral();
            Tnumub->Draw(Form("Erecob >> numub_wosc_spect_histo_%d",ih),cut_numub,"");
            numub_TotWosc[ih] =  numub_wosc_spect_histo[ih]->Integral();

            numu_file << s2t_pt << "\t" << dm2_pt;
            numub_file << s2t_pt << "\t" << dm2_pt;
            //Printing bin-content for the oscilated spectra for (s2t_pt,dm2_pt)
            for (int ib = 0 ; ib < 23 ; ib++)
            {
                double cont   = numu_wosc_spect_histo[ih]->GetBinContent(ib+1);
                numu_file << "\t" << cont;
                if (ib < 12) {
                    cont = numub_wosc_spect_histo[ih]->GetBinContent(ib+1);
                    numub_file << "\t" << cont;
                }
            }
            numu_file  << " \t" << numu_TotWosc[ih]  << endl;
            numub_file << " \t" << numub_TotWosc[ih] << endl;

            //Printing check-points info
            if (ih%10 == 0)
            {
                cout << ih << "  Done with spectrum for " << s2t_pt << "\t" << dm2_pt << endl;
                //cout << " \t" << TotNosc << endl;
            }

            //Normalizing the oscilated spectra for (s2t_pt,dm2_pt)
            numu_integ = numu_wosc_spect_histo[ih]->Integral();
            numu_wosc_spect_histo[ih]->Scale(1.0/numu_integ);
            numub_integ = numub_wosc_spect_histo[ih]->Integral();
            numub_wosc_spect_histo[ih]->Scale(1.0/numub_integ);
        } // for idm2
        numu_file << endl;
        numub_file << endl;
        //fprintf(file,"\n");
    }//for is2t
    numu_file << endl;
    numub_file << endl;

    numu_file.close();
    numub_file.close();


    //---------------------------------------------------
    // Drawing section

    TH2F *frame_spectra = new TH2F("frame_spectra","",NB_numu110,lo110,hi110,10,0,200e4);
    frame_spectra->GetXaxis()->SetTitle("Reco. Energy (GeV)");
    frame_spectra->GetYaxis()->SetTitle("Events");

    TCanvas *canv0 = new TCanvas("canv0","canv0",775,500);

    TLegend *leg_spect = new TLegend(0.7,0.6,0.9,0.9);
    leg_spect->SetFillColor(0);
    leg_spect->AddEntry(numu_BFit_spect_histo,"Best Fit");
    leg_spect->AddEntry(numu_nosc_spect_histo,"No Osc.");
    for (int i = 0 ; i < 4 ; i++)
        leg_spect->AddEntry(numu_wosc_spect_histo[i],Form("Spectra %d ",i));

    numu_BFit_spect_histo->Scale(1.0/(numu_BFit_spect_histo->Integral()));
    numu_nosc_spect_histo->Scale(1.0/(numu_nosc_spect_histo->Integral()));
    
    //frame_spectra->Draw();
    numu_BFit_spect_histo->Draw("same");
    numu_nosc_spect_histo->Draw("same");
    for (int i = 0 ; i < 4 ; i++)
        numu_wosc_spect_histo[i]->Draw("same");
    leg_spect->Draw();
    //break;

    cout << endl << "Mean value: " << numu_BFit_spect_histo->GetMean() << endl;
//
    //canv0->Print("files_plots/MINOS_osc_test.pdf");
//
//    //---------------------------------------------------
//    // Drawing Survival Probabilities at the six ADs
//    TH2F *frame_POscBF = new TH2F("frame_POscBF","",1000,0.0,1.01,10,0.0,0.013);
//    frame_POscBF->GetXaxis()->SetTitleFont(ft);
//    frame_POscBF->GetXaxis()->SetTitleSize(sz);
//    frame_POscBF->GetXaxis()->SetLabelFont(ft);
//    frame_POscBF->GetXaxis()->SetLabelSize(sz);
//    frame_POscBF->GetXaxis()->SetTitle("Survival Probability");
//    frame_POscBF->GetYaxis()->SetTitleFont(ft);
//    frame_POscBF->GetYaxis()->SetTitleSize(sz);
//    frame_POscBF->GetYaxis()->SetLabelFont(ft);
//    frame_POscBF->GetYaxis()->SetLabelSize(sz);
//    //frame_POscBF->GetYaxis()->SetTitle("a.u.");
//
//    TLegend *leg_avgs = new TLegend(0.15,0.7,0.35,0.9);
//    leg_avgs->SetTextFont(ft);
//    leg_avgs->SetTextSize(0.8*sz);
//    leg_avgs->SetFillColor(0);
//    leg_avgs->SetLineColor(0);
//    leg_avgs->AddEntry(numu_Posc_BF,Form("P_{avg} = %f",numu_avgPosc));
//
//    TCanvas *canv1 = new TCanvas("canv1","canv1",775,500);
//    frame_POscBF->Draw();
//    numu_Posc_BF->Draw("same");
//    leg_avgs->Draw();
//
//    canv1->Print("files_plots/MINOS_POsc_avg.pdf");
    //---------------------------------------------------

    fntuple->Close();
} //end
