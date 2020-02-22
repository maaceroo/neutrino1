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
  TFile *fntuple = new TFile("files_root/RENO-ntuple.root","READ");
  TTree *T = (TTree*)fntuple->Get("T");
  TCut cutBF;
  //const int nDet = 2; //Number of Antineutrino detectors
  //const int nRea = 6; //Number of Nuclear Reactors
  
  double daqTime[nDet] = {1807.88,2193.04};
  
  double wrd_array_near[nRea];
  double wrd_array_far[nRea];
  
  TFile *wrd_File = new TFile("files_root/ldist_RENO.root","READ");
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
    //TH1F *Posc_AD_BF[nDet][nRea];
    TH1F *Posc_AD_BF[nDet];
    TH1F *Posc_AD_surv[nDet][nRea];
    TH1F *BFit_spect_histo_d[nDet];
    TH1F *Nevtexp[nDet];
    for (int i = 0 ; i < nDet ; i++)
        {
            Nevtexp[i] = new TH1F(Form("Nevtexp_%d",i),"",NB,xbins);// info from db-ntuple.root
            Nevtexp[i]->SetLineColor(i+3);
      
            //BF-oscillation Ep spectra
            BFit_spect_histo[i] = new TH1F(Form("BFit_spect_histo_%d",i),"",NB,xbins);// info from RENO-ntuple.root
            BFit_spect_histo[i]->SetLineColor(i+1);
      
            //BF-oscillation Ep spectra
            BFit_spect_histo_d[i] = new TH1F(Form("BFit_spect_histo_d_%d",i),"",NB,xbins);//
            BFit_spect_histo_d[i]->SetLineColor(i+1);
      
            //no-oscillation Ep spectra - histograms
            nosc_spect_histo[i] = new TH1F(Form("nosc_spect_histo_%d",i),"",NB,xbins);// info from db-ntuple.root
            nosc_spect_histo[i]->SetLineColor(i+3);
      
            //Oscillation probability at BF - histogram
            Posc_AD_BF[i]  = new TH1F(Form("Posc_AD_BF_%d",i),"",1000,0,1);//to store <POsc(BF)>
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
    file_IBDrates = fopen("files/RENO_noOsc_IBDrates_perday.txt","w");
    
    //int iAD = 2;
    //std::cout << "\n Checking the oscillation Probability" << std::endl;
    //printf("(1.0 - %e*((sin( 1.267 * %e * Ln/En ))**2) - ((cos(0.5 * asin(sqrt(%e))))**4) * %e * (sin( 1.267 * %e * Ln/En ))**2) >> Posc_AD_BF_%d \n \n",ssq2th13,dmsqee,ssq2th13,ssq2th12,dmsq21,iAD);
    
    for (int iAD = 0 ; iAD < nDet ; iAD++)
    {
        //------------------------------------------------
        //Filling Oscillation probability at BF - histogram
        //T->Draw(Form("(1.0 - 0.087*((sin( 1.267 * 2.49e-3 * Ln/En ))**2) - ((cos(0.5 * asin(sqrt(0.087))))**4) * 0.846 * (sin( 1.267 * 7.53e-5 * Ln/En ))**2) >> Posc_AD_BF_%d_%d",iAD,iNR),Form("id==%d  && ir==%d",iAD,iNR));
        T->Draw(Form("(1.0 - %e*((sin( 1.267 * %e * Ln/En ))**2) - ((cos(0.5 * asin(sqrt(%e))))**4) * %e * (sin( 1.267 * %e * Ln/En ))**2) >> Posc_AD_BF_%d",ssq2th13,dmsqee,ssq2th13,ssq2th12,dmsq21,iAD),Form("id==%d",iAD));
        //integ = Posc_AD_BF[iAD][iNR]->Integral();
        integ = Posc_AD_BF[iAD]->Integral();
        //std::cout << "integ = " << integ <<std::endl;
        Posc_AD_BF[iAD]->Scale(1.0/integ); //Used to plot the Survival Probabilities for each AD and each NR with BF parameters (normalized)
        //Average oscillation Probability
        avgPosc[iAD] = Posc_AD_BF[iAD]->GetMean();
        //No-oscillation IDB rate (per day)
        noOsc_IBDrate_perday[iAD] = IBDrate_perday[iAD][0]/avgPosc[iAD];
        //Printing results
        std::cout << "(avgPosc,noOsc_IBDrate_perday)_" << iAD << " = (" << avgPosc[iAD]
                  << ", " << noOsc_IBDrate_perday[iAD] << ") " << std::endl;
        fprintf(file_IBDrates,"%f \n", noOsc_IBDrate_perday[iAD]);
        //------------------------------------------------
    
        //Filling the non-oscillated Ep spectra
        T->Draw(Form("Ep >> nosc_spect_histo_%d",iAD),Form("(id==%d)",iAD),"");
        TotNosc[iAD] =  nosc_spect_histo[iAD]->Integral();
    
        //------------------------------------------------
        //(BF-oscillation probability) condition to fill BF-oscillation- Ep spectra
        cutBF = Form("(1.0 - %e*((sin( 1.267 * %e * Ln/En ))**2) - ((cos(0.5 * asin(sqrt(%e))))**4) * %e * (sin( 1.267 * %e * Ln/En ))**2)*(id==%d)",ssq2th13,dmsqee,ssq2th13,ssq2th12,dmsq21,iAD);
        //Filling and normalizing BF-oscillation Ep spectra
        T->Draw(Form("Ep >> BFit_spect_histo_%d",iAD),cutBF,"");
        integ = BFit_spect_histo[iAD]->Integral();
        //Filling and normalizing BF-oscillation Ep spectra
        T->Draw(Form("Ep >> BFit_spect_histo_d_%d",iAD),cutBF,"");
        integ = BFit_spect_histo_d[iAD]->Integral();
        //std::cout << "integ = " << integ <<std::endl;
        BFit_spect_histo[iAD]->Scale(1.0/integ);
        //BFit_spect_histo_d[iAD]->Scale(integ);
        //------------------------------------------------
    
    }
    fclose(file_IBDrates);
    
    //std::cout << "\n Checking the oscillation Probability" << std::endl;
    //printf("(1.0 - %e*((sin( 1.267 * %e * Ln/En ))**2) - ((cos(0.5 * asin(sqrt(%e))))**4) * %e * (sin( 1.267 * %e * Ln/En ))**2)*(id==%d)",ssq2th13,dmsqee,ssq2th13,ssq2th12,dmsq21,iAD);

  
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
    
    double nebf = 0.0;
    double neno = 0.0;
    for (int iAD = 0 ; iAD < nDet ; iAD++)
    {
        for (int ib = 0 ; ib < NB ; ib++)
        {
            nebf = BFit_spect_histo_d[iAD]->GetBinContent(ib+1);
            //neno = nebf/avgPosc_AD[iAD]/*IBDrate_perday[iAD][0]*daqTime[iAD]*/;
            neno = nebf/avgPosc[iAD]/*IBDrate_perday[iAD][0]*daqTime[iAD]*/;
            //std::cout << " neno = " << neno <<std::endl;
            Nevtexp[iAD]->SetBinContent(ib+1,neno);
        }
    }
  
    ofstream file;
    string grid_name = "files/RENO_gridOscSpectra_test.txt";
    file.open((grid_name).c_str());
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
                    cut = Form("((1.0 - %e*((sin( 1.267 * %e * Ln/En ))**2) - (0.25*(1 + sqrt(1 - %e))**2) * %e * (sin( 1.267 * %e * Ln/En ))**2))*((id==%d))",s2t_pt,dm2_pt,s2t_pt,ssq2th12,dmsq21,iAD);
	      
                    // Filling oscilated spectra for (s2t_pt,dm2_pt)
                    int ih = is2t*N_dm2 + idm2;
                    //std::cout << "AD " << iAD + 1 << "\t s2t_pt =  " << s2t_pt << "\t dm2_pt =  " << dm2_pt  << "\t ih = " << ih << std::endl;
                    //std::cout << "cut computed" << std::endl;
                    //std::cout << "Test No. " << idm2 << std::endl;
                    T->Draw(Form("Ep >> wosc_spect_histo_%d",ih),cut,"");
                    TotWosc[ih] =  wosc_spect_histo[ih]->Integral();
                    //std::cout << TotWosc1[10] <<std::endl;
                    // Survival probability for (s2t_pt,dm2_pt) to fill histogram at each detector, i.e. Posc_AD_surv[iAD]
                    //T->Draw(Form("(1.0 - %e*((sin( 1.267 * %e * Ln/En ))**2) - (0.25*(1 + sqrt(1 - %e))**2) * 0.861 * (sin( 1.267 * 7.53e-5 * Ln/En ))**2) >> Posc_AD_surv_%d",s2t_pt,dm2_pt,s2t_pt,iAD),Form("id==%d",sel));
                    //Average survival probability for (s2t_pt,dm2_pt) at each detector
                    //SurvP = Posc_AD_surv[iAD]->GetMean();
	      
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
                    //*CHECK if this is nneded .. guess is NO
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
