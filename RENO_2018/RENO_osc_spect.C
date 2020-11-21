//-------------------------------------------------------------------------//
//--- RENO_osc_spect.C - By M.A. Acero O., A.A. Aguilar-A. - 2020-02-07 ---//
//-------------------------------------------------------------------------//
// This macro can be executed under ROOT typing                            //
// "root[0] .x RENO_osc_spect.C"                                           //
// For the spect only analysis, using information from                     //
// G. Bak et al., RENO Coll. (2018) PRL121, 201801                         //
//-------------------------------------------------------------------------//

#include "constants.h"
#include <iostream>
#include <fstream>
#include <string>

void RENO_osc_spect()
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
  
    // Open ntuple file to read simulated data
    TString filePath = dirName;
    TFile *fntuple = new TFile(filePath + "/files_root/RENO-ntuple_BFosc.root","READ");
    TTree *T = (TTree*)fntuple->Get("T");
    TCut cutBF;
  
    double daqTime[nDet] = {1807.88,2193.04};
  
    double wrd_array_near[nRea];
    double wrd_array_far[nRea];
  
    TFile *wrd_File = new TFile(filePath + "/files_root/ldist_RENO.root","READ");
    TH1F *wrd_histo_near = ((TH1F*)(wrd_File->Get("histo_ldist_RENO_near")));;
    TH1F *wrd_histo_far = ((TH1F*)(wrd_File->Get("histo_ldist_RENO_far")));;
    
    for (int ir = 0 ; ir < nRea ; ir++)
    {
        wrd_array_near[ir] = wrd_histo_near->GetBinContent(ir+1);
        // std::cout << blid << "  " << wrd_array_near[ir]  << std::endl;
        wrd_array_far[ir] = wrd_histo_far->GetBinContent(ir+1);
        // std::cout << blid << "  " << wrd_array_far[ir]  << std::endl;
    }

    //---------------------------------------------------
    // histogram binning for the simulated data
    //const double    NB = 26;
    //const double    lo = 1.2;
    //const double    hi = 8.4;
    double xbins[NB+1];
    double delta_bins2 = (5.6 - 1.2)/22; // 0.2 MeV/bin
  
    for (int i = 0 ; i < (NB-3) ; i++){
        xbins[i] = 1.2 + delta_bins2*i;
    }
    xbins[23] = xbins[22] + 0.4;
    xbins[24] = xbins[23] + 0.4;
    xbins[25] = 8.4 - 1.4;
    xbins[26] = 8.4;
    //---------------------------------------------------
  
    TH1F *nosc_spect_histo[nDet];
    TH1F *BFit_spect_histo[nDet];
    TH1F *Posc_AD_BF[nDet];
    TH1F *Posc_AD_surv[nDet][nRea];
    TH1F *BFit_spect_histo_d[nDet];
    for (int i = 0 ; i < nDet ; i++)
        {
            //BF-oscillation Ep spectra
            //- This is taken directly from the n-tuple: 5e7 events simulated following
            //- the RENO MC spectra with the BF oscillation parameters
            BFit_spect_histo[i] = new TH1F(Form("BFit_spect_histo_%d",i),"",NB,xbins);// info from RENO-ntuple_BFosc.root
            BFit_spect_histo[i]->SetLineColor(kGreen+i);
      
            //BF-oscillation Ep spectra
            BFit_spect_histo_d[i] = new TH1F(Form("BFit_spect_histo_d_%d",i),"",NB,xbins);//
            BFit_spect_histo_d[i]->SetLineColor(kGreen+i);
      
            //no-oscillation Ep spectra - histograms
            //- This histogram is built by "un-oscillating" the BF spectra (using RENO BF parameters)
            nosc_spect_histo[i] = new TH1F(Form("nosc_spect_histo_%d",i),"",NB,xbins);// info from RENO-ntuple_BFosc.root
            nosc_spect_histo[i]->SetLineColor(kRed+i);
      
            //Oscillation probability at BF - histogram
            //to store <POsc(BF)>
            Posc_AD_BF[i]  = new TH1F(Form("Posc_AD_BF_%d",i),"",1000,0,1);
            Posc_AD_BF[i]->SetLineColor(i+1);

            for (int j = 0 ; j < nRea ; j++)
            {
                //Oscillation probability at (s2th,dm2) - histogram
                Posc_AD_surv[i][j]  = new TH1F(Form("Posc_AD_surv_%d_%d",i,j),"",1000,0,1);//to store <POsc(s2th,dm2)>
            }
        }
    //---------------------------------------------------
    //IBD rate (per day) (PRL121 2018)
    double IBDrate_perday[nDet][2];
    for (int i = 0 ; i < nDet ; i++) {
        for (int j = 0 ; j < 2; j++) {
            IBDrate_perday[i][j] = IBDrate_data[i][j]; //-- Set to observed values from Table I (PRL121 2018)
        }
    }
    double TotNosc[nDet];
    double avgPosc[nDet];              //<POsc(s2t_BF,dm2_ee)>
    double avgPosc_AD[nDet];           //<POsc(s2t_BF,dm2_ee)>
    double noOsc_IBDrate_perday[nDet]; //IBD rate per day w/o oscillations
    double integ;
  
    FILE *file_IBDrates;
    file_IBDrates = fopen(filePath + "/files/RENO_noOsc_IBDrates_perday.txt","w");
    
    //int iAD = 2;
    //std::cout << "\n Checking the oscillation Probability" << std::endl;
    //printf("(1.0 - %e*((sin( 1.267 * %e * Ln/En ))**2) - ((cos(0.5 * asin(sqrt(%e))))**4) * %e * (sin( 1.267 * %e * Ln/En ))**2) >> Posc_AD_BF_%d \n \n",ssq2th13RENO,dmsqeeRENO,ssq2th13RENO,ssq2th12RENO,dmsq21RENO,iAD);
    
    for (int iAD = 0 ; iAD < nDet ; iAD++)
    {
        //------------------------------------------------
        //Filling the un-oscillated number of events (using the BF histogram and the oscillation Probability at BF)
        //T->Draw(Form("(1.0 - %e*((sin( 1.267 * %e * Ln/En ))**2) - ((cos(0.5 * asin(sqrt(%e))))**4) * %e * (sin( 1.267 * %e * Ln/En ))**2) >> Posc_AD_BF_%d",ssq2th13RENO,dmsqeeRENO,ssq2th13RENO,ssq2th12RENO,dmsq21RENO,iAD),Form("id==%d",iAD));
        //-- 03-04-2020
        T->Draw(Form("Ep >> nosc_spect_histo_%d",iAD),Form("(1.0/(1.0 - %e*((sin( 1.267 * %e * Ln/En ))**2) - ((cos(0.5 * asin(sqrt(%e))))**4) * %e * (sin( 1.267 * %e * Ln/En ))**2))*(id==%d)",ssq2th13RENO,dmsqeeRENO,ssq2th13RENO,ssq2th12RENO,dmsq21RENO,iAD));
        TotNosc[iAD] =  nosc_spect_histo[iAD]->Integral();
        std::cout << "TEST1:\t TotNosc[" << iAD << "] = " << TotNosc[iAD] << std::endl;
        //------------------------------------------------
        //------------------------------------------------
        //(BF-oscillation probability) condition to fill BF-oscillation- Ep spectra
        //cutBF = Form("(1.0 - %e*((sin( 1.267 * %e * Ln/En ))**2) - ((cos(0.5 * asin(sqrt(%e))))**4) * %e * (sin( 1.267 * %e * Ln/En ))**2)*(id==%d)",ssq2th13RENO,dmsqeeRENO,ssq2th13RENO,ssq2th12RENO,dmsq21RENO,iAD);
        //Filling BF-oscillation Ep spectra
        cutBF = Form("(id==%d)",iAD);
        //Filling and normalizing BF-oscillation Ep spectra
        T->Draw(Form("Ep >> BFit_spect_histo_%d",iAD),cutBF,"");
        integ = BFit_spect_histo[iAD]->Integral();
        std::cout << "TEST2:\t BFit_spect_histo[" << iAD << "] = " << integ << std::endl;
        //Filling and normalizing BF-oscillation Ep spectra
        T->Draw(Form("Ep >> BFit_spect_histo_d_%d",iAD),cutBF,"");
        integ = BFit_spect_histo_d[iAD]->Integral();
        std::cout << "TEST3:\t BFit_spect_histo_d[" << iAD << "] = " << integ << std::endl;
        //std::cout << "integ = " << integ <<std::endl;
        BFit_spect_histo[iAD]->Scale(1.0/integ);
        //BFit_spect_histo_d[iAD]->Scale(integ);
        //------------------------------------------------
        //------------------------------------------------
        //Average oscillation Probability - I CHANGED THIS HERE TO CORRECTLY COMPUTE THE AVG. OSCILLATION PROBABILITY 2020.04.07
        avgPosc[iAD] = integ/TotNosc[iAD];
        //No-oscillation IDB rate (per day)
        //- We are using the measured IBD rate per day, given that RENO does not provide the IBD rate per day at the BF.
        noOsc_IBDrate_perday[iAD] = BFtoObs[iAD]*IBDrate_data[iAD][0]/avgPosc[iAD];
        //Printing results
        std::cout << "(avgPosc,noOsc_IBDrate_perday)_" << iAD << " = (" << avgPosc[iAD] << ", " << noOsc_IBDrate_perday[iAD] << ") " << std::endl;
        fprintf(file_IBDrates,"%f \n", noOsc_IBDrate_perday[iAD]);
        //------------------------------------------------
    
    }
    fclose(file_IBDrates);
    TCanvas *canv0 = new TCanvas("canv0","AD1");
    canv0->cd();
    nosc_spect_histo[0]->Draw("histo");
    BFit_spect_histo_d[0]->Draw("same histo");
    TCanvas *canv1 = new TCanvas("canv1","AD2");
    canv1->cd();
    nosc_spect_histo[1]->Draw("histo");
    BFit_spect_histo_d[1]->Draw("same histo");
//    std::cout << "\n Checking the oscillation Probability" << std::endl;
//    for (int iAD = 0 ; iAD < nDet ; iAD++)
//    {
//        printf("(1.0 - %e*((sin( 1.267 * %e * Ln/En ))**2) - ((cos(0.5 * asin(sqrt(%e))))**4) * %e * (sin( 1.267 * %e * Ln/En ))**2)*(id==%d)",ssq2th13RENO,dmsqeeRENO,ssq2th13RENO,ssq2th12RENO,dmsq21RENO,iAD);
//    }


    //---------------------------------------------------
    //Definition of the grid of oscillation parameters
    double s2t_pt, dm2_pt;
    double DeltaLin_s2t = (hi_s2t - lo_s2t)/double(N_s2t-1);
    double DeltaLin_dm2 = (hi_dm2 - lo_dm2)/double(N_dm2-1);

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


    ofstream file;
    string grid_name = "/files/RENO_gridOscSpectra_test.txt";
    file.open(filePath + (grid_name).c_str());
    file << fixed;
    file << setprecision(6);

    s2t_pt = 0.0;
    dm2_pt = 0.0;
    for (int iAD = 0 ; iAD < nDet ; iAD++)
    {
        //write non-oscillated spectra for each AD and NR to file
        file << iAD + 1 << "\t" << s2t_pt << "\t" << dm2_pt;
        //print bin-content of non-oscillated spectra per day
        for (int ib = 0 ; ib < NB ; ib++)
        {
            double contNO = nosc_spect_histo[iAD]->GetBinContent(ib+1);
            file << "\t" << contNO;
            //std::cout <<  contNO << "  ";
        }

        file << "\t" << TotNosc[iAD] << std::endl;

    }// for iAD

    file << std::endl;

    for (int is2t = 0 ; is2t < N_s2t ; is2t++)
        {
            s2t_pt = lo_s2t + double(is2t)*DeltaLin_s2t;

            for (int idm2 = 0 ; idm2 < N_dm2 ; idm2++)
            {
                dm2_pt = lo_dm2 + double(idm2)*DeltaLin_dm2;

                for (int iAD = 0 ; iAD < nDet ; iAD++)
                {
                    // Condition to fill oscillated spectra for (s2t_pt,dm2_pt), i.e. wosc_spect_histo[iAD]
                    //cut = Form("((1.0 - %e*((sin( 1.267 * %e * Ln/En ))**2) - (0.25*(1 + sqrt(1 - %e))**2) * %e * (sin( 1.267 * %e * Ln/En ))**2))*((id==%d))",s2t_pt,dm2_pt,s2t_pt,ssq2th12RENO,dmsq21RENO,iAD);

                    //- 2020.04.09 - There was an error with the denominator, as it was written exactly as the numerator, which is not correct because of the ssq2th13RENO term should be used as in the BF cut (line 127)
                    cut = Form("((1.0 - %e*((sin( 1.267 * %e * Ln/En ))**2) - (0.25*(1 + sqrt(1 - %e))**2) * %e * (sin( 1.267 * %e * Ln/En ))**2))*(1.0/(1.0 - %e*((sin( 1.267 * %e * Ln/En ))**2) - ((cos(0.5 * asin(sqrt(%e))))**4) * %e * (sin( 1.267 * %e * Ln/En ))**2))*(id==%d)",s2t_pt,dm2_pt,s2t_pt,ssq2th12RENO,dmsq21RENO,ssq2th13RENO,dmsqeeRENO,ssq2th13RENO,ssq2th12RENO,dmsq21RENO,iAD);

                    // Filling oscilated spectra for (s2t_pt,dm2_pt)
                    int ih = is2t*N_dm2 + idm2;
                    //std::cout << "AD " << iAD + 1 << "\t s2t_pt =  " << s2t_pt << "\t dm2_pt =  " << dm2_pt  << "\t ih = " << ih << std::endl;
                    //std::cout << "cut computed" << std::endl;
                    //std::cout << "Test No. " << idm2 << std::endl;
                    T->Draw(Form("Ep >> wosc_spect_histo_%d",ih),cut,"");
                    TotWosc[ih] =  wosc_spect_histo[ih]->Integral();
                    //std::cout << TotWosc1[10] <<std::endl;

                    file << iAD + 1 << "\t" << s2t_pt << "\t" << dm2_pt;
                    //fprintf(file,"%d %8.2e %8.2e",iAD+1,s2t_pt,dm2_pt);
                    //Printing bin-content for the oscilated spectra for (s2t_pt,dm2_pt)
                    for (int ib = 0 ; ib < NB ; ib++)
                    {
                        double cont   = wosc_spect_histo[ih]->GetBinContent(ib+1);
                        //std::cout << " cont = " << cont << " ih = " << ih <<std::endl;

                        file << "\t" << cont;
                        //fprintf(file," %10.5f",cont);
                        //std::cout << "cont =  " << cont1 <<std::endl;
                    }

                    file << " \t" << TotWosc[ih] << std::endl;
                    //fprintf(file," %10.2f\n",TotWosc[ih]);

                    //Printing check-points info
                    if (ih%10 == 0)
                    {
                        std::cout << ih << "  Done with detector " << iAD+1  << " for "<< s2t_pt << "\t" << dm2_pt << std::endl;
                        // std::cout << " \t" << TotNosc[iAD] << std::endl;
                    }

                    //Normalizing the oscilated spectra for (s2t_pt,dm2_pt)
                    integ = wosc_spect_histo[ih]->Integral();
                    wosc_spect_histo[ih]->Scale(1.0/integ);
                }//for iAD
                file << std::endl;
                //fprintf(file,"\n");
            }//for idm2
            //std::cout << "  Done with detector " << iAD+1 << std::endl;
            file << std::endl;
            //fprintf(file,"");
        }//for is2t

    file.close();
    //fclose(file);

} //end
