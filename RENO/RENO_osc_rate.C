//------------------------------------------------------------------------------------------//
//--- RENO_osc_rate.C - By M.A. Acero O., A.A. Aguilar-A. and D.J. Polo T. - 2017-22-10 ----//
//------------------------------------------------------------------------------------------//
// This macro can be executed under ROOT typing   "root[0] .x RENO_osc_rate.C"              //
// For the rate-only analysis, using information from F.P. An et al., PRL 112 061801 (2014) //
//------------------------------------------------------------------------------------------//
// 2017-23-12                                                                               //
// Modifications to perform a rate-only analysis by computing and using average values for  //
// the survival probability and taking sin2(2th) as the only free parameter (to bi fitted). //
// 2017-02-03                                                                               //
// This macro performs a simple chi^2 analysis of the RENO data by comparing the IBD        //
// rate at the six AD reported in article 1610.04326 (2017)                                 //
//------------------------------------------------------------------------------------------//

#include <constants.h>
#include <iostream>
#include <fstream>
#include <string>

void RENO_osc_rate()
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
    // Open ntuple file to read simulated data (Ep, En, Ln) for the 2 AD and 6 reactors
    TFile *fntuple = new TFile("files_root/RENO-ntuple.root","READ");
    TTree *T = (TTree*)fntuple->Get("T");
    TCut cutBF;
    //---------------------------------------------------
    // histogram binning for the simulated data
   	const double    NB = 27; 
 	const double    lo = 1.2;  
 	const double    hi = 8.4;
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
    const int nAD = 2; //Number of Antineutrino detectors
    TH1F *BFit_spect_histo[nAD];
    TH1F *Posc_AD_BF[nAD];
    TH1F *Posc_AD_surv[nAD];
    TH1F *SinDelta21[nAD];
    TH1F *SinDeltaee[nAD];
    for (int i = 0 ; i < nAD ; i++)
    {
        //BF-oscillation Ep spectra
        BFit_spect_histo[i] = new TH1F(Form("BFit_spect_histo_%d",i),"",NB,xbins);// info from RENO-ntuple.root
        BFit_spect_histo[i]->SetLineColor(i+1);

        //Ocillation probability at BF - histogram
        Posc_AD_BF[i]       = new TH1F(Form("Posc_AD_BF_%d",i),"",1000,0,1);//to store <POsc(BF)>
        Posc_AD_BF[i]->SetLineColor(i+1);

        //Ocillation probability at BF - histogram
        SinDelta21[i]       = new TH1F(Form("SinDelta21_%d",i),"",1000,0,1);//to store <sin^2(1.267 dm2_21 L/E)>
        SinDelta21[i]->SetLineColor(i+1);

        //Ocillation probability at BF - histogram
        SinDeltaee[i]       = new TH1F(Form("SinDeltaee_%d",i),"",1000,0,1);//to store <sin^2(1.267 dm2_31 L/E)>
        SinDeltaee[i]->SetLineColor(i+2);

        //Ocillation probability at (s2th,dm2) - histogram
        Posc_AD_surv[i]     = new TH1F(Form("Posc_AD_surv_%d",i),"",1000,0,1);//to store <POsc(s2th,dm2)>
    }
    //---------------------------------------------------
    //IBD rate (per day), total background and efficiencies (1610.04326v4 (2017))
    double IBDrate_perday[nAD][2] =
    {
        {616.67,1.44},{61.24,0.42}
    };
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

    double avgPosc_AD[nAD];           //<POsc(s2t_BF,dm2_ee)>
    double noOsc_IBDrate_perday[nAD]; //IBD rate per day w/o oscillations
    double integ;
    for (int iAD = 0 ; iAD < nAD ; iAD++)
    {
        //------------------------------------------------
        //Filling Oscillation probability at BF - histogram
        T->Draw(Form("(1.0 - 0.087*((sin( 1.267 * 2.49e-3 * Ln/En ))**2) - ((cos(0.5 * asin(sqrt(0.087))))**4) * 0.846 * (sin( 1.267 * 7.53e-5 * Ln/En ))**2) >> Posc_AD_BF_%d",iAD),Form("id==%d",iAD));
        integ = Posc_AD_BF[iAD]->Integral();
        Posc_AD_BF[iAD]->Scale(1.0/integ); //Used to plot the Survival Probabilities for each AD with BF parameters (normalized)
        
        //Average oscillation Probability
        avgPosc_AD[iAD] = Posc_AD_BF[iAD]->GetMean();
        //No-oscillation IDB rate (per day)
        noOsc_IBDrate_perday[iAD] = IBDrate_perday[iAD][0]/avgPosc_AD[iAD];
        //Printing results
        cout << "(avgPosc_AD,noOsc_IBDrate_perday)_" << iAD << " = (" << avgPosc_AD[iAD]
             << ", " << noOsc_IBDrate_perday[iAD] << ") " << endl;
        //------------------------------------------------

        //------------------------------------------------
        //(BF-oscillation probability) condition to fill BF-oscillation- Ep spectra
        cutBF = Form("(1.0 - 0.087*((sin( 1.267 * 2.49e-3 * Ln/En ))**2) - ((cos(0.5 * asin(sqrt(0.087))))**4) * 0.846 * (sin( 1.267 * 7.53e-5 * Ln/En ))**2)*(id==%d)",iAD);
        //Filling and normalizing BF-oscillation Ep spectra
        T->Draw(Form("Ep >> BFit_spect_histo_%d",iAD),cutBF,"");
        integ = BFit_spect_histo[iAD]->Integral();
        BFit_spect_histo[iAD]->Scale(1.0/integ);
        //------------------------------------------------
    }

    //---------------------------------------------------
    // computing <sin^2(1.267 dm2_21 L/E)> and <sin^2(1.267 dm2_ee L/E)>
    
    double avgSinDelta21[nAD]; //<sin^2(1.267 dm2_21 L/E)> for each AD
    double avgSinDeltaee[nAD]; //<sin^2(1.267 dm2_ee L/E)> for each AD

    double dm2_21 = 7.53e-5; //eV^2,                  1610.04326v4 (2017)
    double dm2_ee = 2.49e-3; //eV^2,                  1610.04326v4 (2017)
    int   sel;

    for (int iAD = 0 ; iAD < nAD ; iAD++)
    {
         T->Draw(Form("((sin(1.267 * %e * Ln/En))**2) >> SinDelta21_%d",dm2_21,iAD),Form("id==%d",iAD));
        integ = SinDelta21[iAD]->Integral();
        SinDelta21[iAD]->Scale(1.0/integ); //Used to plot
        avgSinDelta21[iAD] = SinDelta21[iAD]->GetMean();
        T->Draw(Form("((sin(1.267 * %e * Ln/En))**2) >> SinDeltaee_%d",dm2_ee,iAD),Form("id==%d",iAD));
        integ = SinDeltaee[iAD]->Integral();
        SinDeltaee[iAD]->Scale(1.0/integ); //Used to plot
        avgSinDeltaee[iAD] = SinDeltaee[iAD]->GetMean();
        cout << "avgSinDelta21_" << iAD << " = " << avgSinDelta21[iAD] << "\t"
             << "avgSinDeltaee_" << iAD << " = " << avgSinDeltaee[iAD] << endl;
    }
    //---------------------------------------------------

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
        SurvProb->SetParameters(avgSinDelta21[i],avgSinDeltaee[i]);
        double test = SurvProb->Eval(0.087);
        cout << i << "test = " << test << endl;
    }
    
    //break;
    
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
    double sig_r = 0.009;
    double sig_d = 0.002;
    
    //File to print the chi-square results
    ofstream file;
    string result = "files/RENO_chi2_rate.txt";
    file.open (result.c_str());
    file << setprecision(5);

    for (int is2t = 0 ; is2t < N_s2t ; is2t++)
    {
        for (int iAD = 0 ; iAD < nAD ; iAD++)
        {
            SurvProb->SetParameter(0,avgSinDelta21[iAD]);
            SurvProb->SetParameter(1,avgSinDeltaee[iAD]);
            
            s2t_pt = 10**(log10(lo_s2t) + double(is2t)*DeltaLog_s2t);
            //s2t_pt = lo_s2t + double(is2t)*DeltaLin_s2t;
            
            SurvP = SurvProb->Eval(s2t_pt);
            IBDrate_osc[is2t] = (SurvP * noOsc_IBDrate_perday[iAD] + totalBgd[iAD][0]);
            
            double sqrerror = (IBDrate_data[iAD][1])**2 + (totalBgd[iAD][1])**2;
            
            chi2[is2t] += ((IBDrate_osc[is2t] - IBDrate_data[iAD][0])**2)/sqrerror;
            
        }//for iAD
        file << s2t_pt << "\t" << SurvP << "\t" << chi2[is2t] << endl;
    }//for is2t
    file << endl;

    file.close();

    double chi2_min = 5e+5;
    int sel;
    for (int jj = 0 ; jj < N_s2t ; jj++)
    {
        if (chi2[jj] < chi2_min)
        {
            sel = jj;
            chi2_min = chi2[sel];
        }
    }
    cout << "chi2_min = " << chi2_min << ", for jj = " << sel << endl;
    
    //---------------------------------------------------
    //---------------------------------------------------
    // Drawing section

    TH2F *frame_spectra = new TH2F("frame_spectra","",NB,lo,hi,10,0,18);
    frame_spectra->GetXaxis()->SetTitle("Prompt energy (MeV)");
    frame_spectra->GetYaxis()->SetTitle("Events/day (bkgd-subtracted)");

    TCanvas *canv1 = new TCanvas("canv1","canv0",775,500);
    
    TLegend *leg_spect = new TLegend(0.7,0.6,0.9,0.9);
    leg_spect->SetFillColor(0);
    for (int i = 0 ; i < 1 ; i++)
    {
        leg_spect->AddEntry(BFit_spect_histo[i],Form("AD%d",i+1));
        //leg_spect->AddEntry(wosc_spect_histo[i],Form("Spectra %d ",i));
    }
    
    //frame_spectra->Draw();
    BFit_spect_histo[0]->Draw("");
    leg_spect->Draw();

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
    
    TCanvas *canv0 = new TCanvas("canv0","canv0",775,500);
    frame_POscBF->Draw();
    for (int j = 0 ; j < nAD ; j++) {
        Posc_AD_BF[j]->Draw("same");
    }
    leg_avgs->Draw();
    
    canv0->Print("Plots/POsc_avg.pdf");
    //---------------------------------------------------

    //---------------------------------------------------
    //---------------------------------------------------

    printf("******************************************************************************************\n");
    printf("* Declaration of arrays needed for RENO_minuit1.C lines 65, 67 and 69. Copy and paste there. *\n");
    printf("******************************************************************************************\n\n");
    //-- Printing out for RENO_minuit1.C
    printf("double noOsc_IBDrate_perday[nAD] = {");
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
    
    printf("double avgSinDeltaee[nAD] = {");
    for (int iAD = 0 ; iAD < nAD ; iAD++)
    {
        printf("%12.9f",avgSinDeltaee[iAD]);
        if(iAD < nAD-1)
            printf(",");
        else
            printf("};\n\n");
    }
    printf("******************************************************************************************\n");
    

    
} //end
