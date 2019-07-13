//------------------------------------------------------------------------------------------//
//--- RENO_osc_spect.C - By M.A. Acero O., A.A. Aguilar-A. and D.J. Polo T. - 2018-15-09 ---//
//------------------------------------------------------------------------------------------//
// This macro can be executed under ROOT typing   "root[0] .x RENO_osc_spect.C"             //
// For the spect only analysis, using information from F.P. An et al., PRL 112 061801 (2014)//
//------------------------------------------------------------------------------------------//
// 2017-02-03                                                                               //
// This macro performs a simple chi^2 analysis of the RENO data by comparing the IBD        //
// rate at the six AD reported in article 1610.04326 (2017)                                 //
//------------------------------------------------------------------------------------------//

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
  const int nAD = 2; //Number of Antineutrino detectors
  const int nNR = 6; //Number of Nuclear Reactors
  
  double daqTime[nAD] = {458.49,489.93};
  ////////////////////////////////////////////////////////////////*************************//////////////////////////////////////////////
  
  double wrd_array_near[nNR];
  double wrd_array_far[nNR];
  
  TFile *wrd_File = new TFile("files_root/ldist_RENO_2x6.root","READ");
  TH1F *wrd_histo_near = ((TH1F*)(wrd_File->Get("histo_ldist_RENO_near")));;
  TH1F *wrd_histo_far = ((TH1F*)(wrd_File->Get("histo_ldist_RENO_far")));;
    
  for (int ir = 0 ; ir < nNR ; ir++)
    {
      
      wrd_array_near[ir] = wrd_histo_near->GetBinContent(ir+1);
      // std::cout << blid << "  " << wrd_array_near[ir]  << std::endl;
      wrd_array_far[ir] = wrd_histo_far->GetBinContent(ir+1);
      // std::cout << blid << "  " << wrd_array_far[ir]  << std::endl;
      
    }

  //---------------------------------------------------
  // histogram binning for the simulated data
  //const double    NB = 27; 
  //const double    lo = 1.2;  
  //const double    hi = 8.4;
  double xbins[NB+1];
  double delta_bins2 = (6.0 - 1.2)/24; // 0.2 MeV/bin
  
  for (int i = 0 ; i < (NB-2) ; i++)
    {
      xbins[i] = 1.2 + delta_bins2*i;
    }
  xbins[25] = xbins[24] + 0.4;
  xbins[26] = 8.4 - 1.4;
  xbins[27] = 8.4;
  
  //---------------------------------------------------
  
  TH1F *nosc_spect_histo[nAD];
  TH1F *BFit_spect_histo[nAD];
  TH1F *Posc_AD_BF[nAD][nNR];
  TH1F *Posc_AD_surv[nAD][nNR];
  TH1F *BFit_spect_histo_d[nAD];
  TH1F *Nevtexp[nAD];
  for (int i = 0 ; i < nAD ; i++)
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
      
      for (int j = 0 ; j < nNR ; j++)
	{
	  
	  //Oscillation probability at BF - histogram
	  Posc_AD_BF[i][j]  = new TH1F(Form("Posc_AD_BF_%d_%d",i,j),"",1000,0,1);//to store <POsc(BF)>
	  Posc_AD_BF[i][j]->SetLineColor(i+1);
	  
	  //Oscillation probability at (s2th,dm2) - histogram
	  Posc_AD_surv[i][j]  = new TH1F(Form("Posc_AD_surv_%d_%d",i,j),"",1000,0,1);//to store <POsc(s2th,dm2)>
	}
    }
  //---------------------------------------------------
  //IBD rate (per day), total background and efficiencies (1610.04326v4 (2017))
  double IBDrate_perday[nAD][2] =
    {
      {616.67,1.44},{61.24,0.42}
    };
  /*
  //(IBD candidates)/(DAQ live time -days-) from 1610.04326v4 (2017)
  double IBDrate_data[nAD][2] =
  {
  {634.20,1.18},{64.38,0.36}
  };
  double totalBgd[nAD][2] =
  {
  {17.54,0.83},{3.14,0.23}
  };
  // double emuem[nAD] = {0.7644,0.7644};
  /*/
  
  double TotNosc[nAD];
  double avgPosc[nAD][nNR];              //<POsc(s2t_BF,dm2_ee)>
  double avgPosc_AD[nAD];           //<POsc(s2t_BF,dm2_ee)>
  double noOsc_IBDrate_perday[nAD][nNR]; //IBD rate per day w/o oscillations
  double noOsc_IBDrate_perday_AD[nAD]; //IBD rate per day w/o oscillations
  double integ;
  
  for (int iAD = 0 ; iAD < nAD ; iAD++)
    {
      for (int iNR = 0 ; iNR < nNR ; iNR++)
	{
	  //------------------------------------------------
	  //Filling Oscillation probability at BF - histogram
	  T->Draw(Form("(1.0 - 0.087*((sin( 1.267 * 2.49e-3 * Ln/En ))**2) - ((cos(0.5 * asin(sqrt(0.087))))**4) * 0.846 * (sin( 1.267 * 7.53e-5 * Ln/En ))**2) >> Posc_AD_BF_%d_%d",iAD,iNR),Form("id==%d  && ir==%d",iAD,iNR));
	  integ = Posc_AD_BF[iAD][iNR]->Integral();
	  //std::cout << "integ = " << integ <<std::endl;
	  Posc_AD_BF[iAD][iNR]->Scale(1.0/integ); //Used to plot the Survival Probabilities for each AD and each NR with BF parameters (normalized)
	  //Average oscillation Probability
	  avgPosc[iAD][iNR] = Posc_AD_BF[iAD][iNR]->GetMean();
	  //No-oscillation IDB rate (per day)
	  noOsc_IBDrate_perday[0][iNR] = (IBDrate_perday[0][0]/avgPosc[0][iNR])*wrd_array_near[iNR];
	  noOsc_IBDrate_perday[1][iNR] = (IBDrate_perday[1][0]/avgPosc[1][iNR])*wrd_array_far[iNR];
	  //Printing results
	  std::cout << "(avgPosc,noOsc_IBDrate_perday)_" << iAD << "_" << iNR << " = (" << avgPosc[iAD][iNR]
	       << ", " << noOsc_IBDrate_perday[iAD][iNR] << ") " << std::endl;
	  //------------------------------------------------
	  
	}
      
      //Filling the non-oscillated Ep spectra
      T->Draw(Form("Ep >> nosc_spect_histo_%d",iAD),Form("(id==%d)",iAD),"");
      TotNosc[iAD] =  nosc_spect_histo[iAD]->Integral();
      
      //------------------------------------------------
      //(BF-oscillation probability) condition to fill BF-oscillation- Ep spectra
      cutBF = Form("(1.0 - 0.087*((sin( 1.267 * 2.49e-3 * Ln/En ))**2) - ((cos(0.5 * asin(sqrt(0.087))))**4) * 0.846 * (sin( 1.267 * 7.53e-5 * Ln/En ))**2)*(id==%d)",iAD);
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
  
  //---------------------------------------------------
  //Definition of the grid of oscillation parameters
  double s2t_pt, dm2_pt;
  //const int     N_s2t = 100;
  //const int     N_dm2 = 100;
  //double       lo_s2t = 0.01;
  //double       hi_s2t = 0.2;
  double DeltaLin_s2t = (hi_s2t - lo_s2t)/double(N_s2t-1);
  //double DeltaLog_s2t = (log10(hi_s2t)-log10(lo_s2t))/double(N_s2t-1);
  //double       lo_dm2 = 1.2e-3;
  //double       hi_dm2 = 3.5e-3;
  double DeltaLin_dm2 = (hi_dm2 - lo_dm2)/double(N_dm2-1);
  //double DeltaLog_dm2 = (log10(hi_dm2)-log10(lo_dm2))/double(N_dm2-1);
  
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
  for (int iAD = 0 ; iAD < nAD ; iAD++)
    {
      for (int ib = 0 ; ib < NB ; ib++)
	{
	  nebf = BFit_spect_histo_d[iAD]->GetBinContent(ib+1);
	  neno = nebf/avgPosc_AD[iAD]/*IBDrate_perday[iAD][0]*daqTime[iAD]*/;
	  //std::cout << " neno = " << neno <<std::endl;
	  Nevtexp[iAD]->SetBinContent(ib+1,neno);
	}
    }
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
    //FILE *file;
    //file = fopen("files/RENO_gridOscSpectra_test.txt","w");
    ofstream file;
    string grid_name = "files/RENO_gridOscSpectra_test.txt";
    file.open((grid_name).c_str());
    file << fixed;
    file << setprecision(6);

    s2t_pt = 0.0;
    dm2_pt = 0.0;
    for (int iAD = 0 ; iAD < nAD ; iAD++)
    {
        //write non-oscillated spectra for each AD and NR to file
        file << iAD + 1 << "\t" << s2t_pt << "\t" << dm2_pt;
        //fprintf(file,"%d %8.2e %8.2e",iAD+1,s2t_pt,dm2_pt);
        //print bin-content of non-oscillated spectra per day
        for (int ib = 0 ; ib < NB ; ib++)
        {
            double contNO = nosc_spect_histo[iAD]->GetBinContent(ib+1);
            file << "\t" << contNO;
            //	  std::cout <<  contNO << "  ";
            //fprintf(file," %10.2f ",contNO);
        }
      
        //fprintf(file," %10.2f\n",TotNosc[iAD]);
        file << "\t" << TotNosc[iAD] << std::endl;
      
    }
    //} // for iAD
    file << std::endl;
    //fprintf(file,"\n");
  
    for (int is2t = 0 ; is2t < N_s2t ; is2t++)
        {
            //s2t_pt = 10**(log10(lo_s2t) + double(is2t)*DeltaLog_s2t);
            s2t_pt = lo_s2t + double(is2t)*DeltaLin_s2t;
      
            for (int idm2 = 0 ; idm2 < N_dm2 ; idm2++)
            {
                //dm2_pt = 10**(log10(lo_dm2) + double(idm2)*DeltaLog_dm2);
                dm2_pt = lo_dm2 + double(idm2)*DeltaLin_dm2;
	  
                for (int iAD = 0 ; iAD < nAD ; iAD++)
                {
                    // Condition to fill oscillated spectra for (s2t_pt,dm2_pt), i.e. wosc_spect_histo[iAD]
                    cut = Form("((1.0 - %e*((sin( 1.267 * %e * Ln/En ))**2) - (0.25*(1 + sqrt(1 - %e))**2) * 0.861 * (sin( 1.267 * 7.53e-5 * Ln/En ))**2))*((id==%d))",s2t_pt,dm2_pt,s2t_pt,iAD);
	      
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
  
    const int columns = 31;
    //const int rows2 = rows-3;
    const int rows = nAD*N_s2t*N_dm2 + nAD;
    const int rows2 = rows/2;
    ofstream Ratio;
    string Ratio_spectra = "files/RENO_gridOscSpectra_Ratio.txt";
    Ratio.open((Ratio_spectra).c_str());
  
    ifstream matriz("files/RENO_gridOscSpectra_test.txt");
    double ** matr;
    matr = new double*[rows];
    for(int k = 0 ; k < rows ; k++)
    {
        matr[k] = new double[columns];
    }
    
    for(int j = 0 ; j < rows ; j++)
    {
        for(int l = 0 ; l < columns ; l++)
        {
            matriz >> matr[j][l];
            //std::cout << matr[j][l] << " ";
        }
    //std::cout << std::endl;
    }

    for(int i = 0; i < rows2 ; i++ )
    {
        int k = 2*i;
        int n = 2*i + 1;
        Ratio << matr[k][1] << "  " << matr[k][2] << "  ";
        for(int l = 3 ; l < columns ; l++)
        {
            //std::cout << matr[n][l] << "  " << matr[k][l] << "  ";
            double ratio = matr[n][l]/matr[k][l];
            Ratio << ratio << "  " ;
        }
        Ratio << std::endl;
    
    }
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////7
  
  //break;
  /*
  //---------------------------------------------------
  // Drawing section
  
  TH2F *frame_spectra = new TH2F("frame_spectra","",NB,lo,hi,10,0,18);
  frame_spectra->GetXaxis()->SetTitle("Prompt energy (MeV)");
  frame_spectra->GetYaxis()->SetTitle("Events/day (bkgd-subtracted)");
  
  TCanvas *canv0 = new TCanvas("canv0","canv0",775,500);
  
  
  TLegend *leg_spect = new TLegend(0.7,0.6,0.9,0.9);
  leg_spect->SetFillColor(4);
  
  TLegend *leg_spect1 = new TLegend(0.7,0.6,0.9,0.9);
  leg_spect1->SetFillColor(5);
  for (int i = 0 ; i < 1 ; i++)
    {
      leg_spect->AddEntry(BFit_spect_histo[i],Form("AD%d",i+1));
      
    }
  
  for (int i = 0 ; i < 1 ; i++)
    {
      leg_spect1->AddEntry(nosc_spect_histo[i],Form("AD%d",i+1));
      
    }
  
  //for (int i = 0 ; i < 4 ; i++)
  //	{
  //  leg_spect->AddEntry(wosc_spect_histo[i],Form("Spectra %d ",i));
  //}
  TCanvas *c = new TCanvas("c","c",775,500);
  //frame_spectra->Draw();
  
  nosc_spect_histo[0]->Draw();
  //wosc_spect_histo[0]->Draw("same");
  BFit_spect_histo_d[0]->Draw("same");
  //Nevtexp[0]->Draw("same");
  //leg_spect->Draw("same");
  TCanvas *c3 = new TCanvas("c3","c3",775,500);
  
  nosc_spect_histo[1]->Draw();
  //wosc_spect_histo[1]->Draw("same");
  BFit_spect_histo_d[1]->Draw("same");
  //Nevtexp[1]->Draw("same");
  //      leg_spect->Draw("SAME");
  c->Print("Plots/PLOT_BF_Near.pdf");
  c3->Print("Plots/PLOT_BF_far.pdf");
  
  printf("****************************************************************************************\n");
  printf("* Declaration of arrays needed for RENO_minuit_spect.C lines 67. Copy and paste there. *\n");
  printf("****************************************************************************************\n\n");
  //-- Printing out for RENO_minuit_spect.C
  
  printf("double noOsc_IBDrate_perday[nAD] ={ ");
  for (int iAD = 0 ; iAD < nAD ; iAD++)
    {
      printf("%7.2f", noOsc_IBDrate_perday_AD[iAD]);
      if(iAD < nAD-1) 
	printf(",");
      else
	printf("};\n");
    }
  
  printf("******************************************************************************************\n");
  
  */
  //break;
  /*
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
      for (int j = 0 ; j < nNR ; j++)
	{
	  leg_avgs->AddEntry(Posc_AD_BF[i][j],Form("AD%d NR%d P_{avg} = %f",i+1, j+1,avgPosc[i][j]));
	}
      leg_spect->AddEntry(nosc_spect_histo[i],Form("AD%d ",i+1));
    }
  TCanvas *canv0 = new TCanvas("canv0","canv0",775,500);
  frame_POscBF->Draw();
  for (int j = 0 ; j < nAD ; j++) {
    for (int j = 0 ; j < nAD ; j++) {
      Posc_AD_BF[j]->Draw("same");
    }
  }
  leg_avgs->Draw();
  
  canv0->Print("POsc_avg.pdf");
  //---------------------------------------------------
  */
  
  
  
} //end
