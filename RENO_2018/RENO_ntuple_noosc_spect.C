//------------------------------------------------------------------------------------//
//--   RENO_ntuple_noosc_spect.C - By M.A. Acero O. & A.A Aguilar-A. - 2019-12-19   --//
//------------------------------------------------------------------------------------//
// This macro can be executed under ROOT typing                        		          //
// "root[0] .x RENO_ntuple_noosc_spect.C"                                             //
// or, from a Terminal, typing
// "root -n -l -b -q RENO_ntuple_noosc_spect.C"                                       //
// It should be executed independently of the oter macros related to the analysis.    //
//------------------------------------------------------------------------------------//
// For the spectral analysis, using information from:                  		          //
// - F.P. An et al.,  RENO Coll. (2017) - arXiv:1610.04326                            //
// - F.P. An et al., arXiv:1003.1391" 6 Mar 2010. (Table 1.2)                         //
// - G. Bak et al., RENO Coll. (2018) PRL121, 201801 (data)                           //
//------------------------------------------------------------------------------------//
// With this macro we build a spectrum for each RENO antineutrino detector, filled    //
// simulated no-oscillated neutrinos events. To do so, we use the MC oscillated       //
// expected spectra reported by RENO in PRL121,201801(2018) to sample the simulated   //
// events. We then compute the average neutrino oscillation probability and           //
// "un-oscilate" the MC spectra using the Collaboration's best fit parameters.        //
//------------------------------------------------------------------------------------//
// The resulting no-osccilated spectra are used in our RENO data analysis.            //
// No-oscillated espectra (events and events / 0.2 Mev) are saved to the file         //
// "files_root/RENOplots_noosc.root"                                                  //
//------------------------------------------------------------------------------------//

#include "constants.h"

void RENO_ntuple_noosc_spect()
{ // begin
    //------------- Style --------------
    gROOT->SetStyle("Plain");
    gStyle->SetOptStat(0);
    //------------- Style --------------
    //----------  Text Style  ---------
    //ft = 10 * fontID + precision
    Int_t ft = 10 * 4 + 2;
    Double_t sz = 0.04;
    //---------------------------------

    //---------------------------------------------------
    // histogram binning for the simulated data
    //const int  NB = 26;
    //const double lo = 1.2;
    const double hi = 8.4;
    double xbins[NB+1];
    double delta_bins2 = (5.6 - 1.2)/22; // 0.2 MeV/bin
    
    for (int i = 0 ; i < (NB-3) ; i++){
        xbins[i] = 1.2 + delta_bins2*i;
    }
    xbins[23] = xbins[22] + 0.4;
    xbins[24] = xbins[23] + 0.4;
    xbins[25] = 8.4 - 1.4;
    xbins[26] = 8.4;
    double osc_prob_nd[NB] = {0.0};
    double osc_prob_fd[NB] = {0.0};
    double N_events_perbin_nd[NB] = {0.0};
    double N_events_perbin_fd[NB] = {0.0};
    
    const int nDet = 2;   // number of AD at Reno
    const int nRea = 6;   // number of reactors
  
    //---------------------------------------------------
    //-------------------
    // Energy Histograms
    //-------------------
    TString filePath = dirName;
    TFile *fenergy = new TFile(filePath + "/files_root/RENOplots.root","read");
  
    //The ND and FD MC spectra
    //From PRL116 211801 (2016). Figure 2  (Blue lines).
    TH1F *MC_spect_histo[nDet];
    TH1F *MC_spect_histo_perMeV[nDet];
    TH1F *noosc_spect_histo[nDet];
    TH1F *noosc_spect_histo_perMeV[nDet];
    for(int n = 0 ; n < nDet ; n++){
        MC_spect_histo_perMeV[n] = (TH1F*) fenergy->Get(Form("mc_histo_%d",n));  // Get data (events per MeV)
        MC_spect_histo_perMeV[n]->SetLineColor(kBlue);
        MC_spect_histo[n] = new TH1F(Form("MC_spect_histo_%d",n),"",NB,xbins);
        MC_spect_histo[n]->SetLineColor(kBlue);
        noosc_spect_histo[n] = new TH1F(Form("noosc_spect_histo_%d",n),"",NB,xbins);
        noosc_spect_histo[n]->SetLineColor(kRed);
        noosc_spect_histo_perMeV[n] = new TH1F(Form("noosc_spect_histo_perMeV_%d",n),"",NB,xbins);
        noosc_spect_histo_perMeV[n]->SetLineColor(kRed);
    }
    //- This is ti buid the number of events histogram -//
    int NumB = MC_spect_histo_perMeV[1]->GetNbinsX();
    for(int n = 0 ; n < nDet ; n++){
        for(int ib = 0; ib < NB; ib++){
            double binW = MC_spect_histo[1]->GetBinWidth(ib+1);
            double cont =  MC_spect_histo_perMeV[n]->GetBinContent(ib+1);
            MC_spect_histo[n]->SetBinContent(ib+1,cont*binW*0.2);  // Get data (events)
        }
    }
  
    //-------------------
    // Distance Histogram
    //-------------------
    TFile *fpathl = new TFile(filePath + "/files_root/ldist_RENO.root","read");
    // Get histogram - histo_ldist_RENO_2x6
    TH1F *histo_ldist_RENO_2x6 = (TH1F*) fpathl->Get("histo_ldist_RENO_2x6");
  
    //Baseline Distances (m)
    const char  *detNames[nDet] = {"near-AD1", "far-AD2"};
  
    //The 12 baselines in RENO AD*NR = 12
    double baselines[nDet][nRea] =
    {
        { 667.9, 451.8, 304.8, 336.1, 513.9, 739.1},
        {1556.5,1456.2,1395.9,1381.3,1413.8,1490.1}
    };
  //-- 18-09-2020 Trying to add different uncertainties for each baseline --
  TFile *fBLsdev = new TFile("files_root/reno_BL_deviations.root","read");
  TH1F *bldist[nDet][nRea];
  for (int ii = 0 ; ii < nDet ; ii++) {
      for (int jj = 0 ; jj < nRea ; jj++) {
          bldist[ii][jj] = (TH1F*) fBLsdev->Get(Form("bldist_%d_%d",ii,jj));
          //std::cout << "\t Mean Distance" << "(" << ii << ", " << jj << ") = " << bldist[ii][jj]->GetMean()
          //<< "\t Deviation" << "(" << ii << ", " << jj << ") = " << bldist[ii][jj]->GetStdDev() << std::endl;
      }
  }
  //-- 18-09-2020 Independent Gaussian fluctuation for each baseline with different StdDev
  TF1 *gauL[nDet][nRea];
  for (int iD = 0 ; iD < nDet ; iD++)
      for (int iR = 0 ; iR < nRea ; iR++){
          gauL[iD][iR] = new TF1(Form("gauL_%d_%d",iD,iR),"exp(-0.5*(x/[0])^2)",-10.0,10.0);
          gauL[iD][iR]->SetParameter(0,bldist[iD][iR]->GetStdDev());
      }

    //make ntuple
    //TFile *fout = new TFile("files_root/RENO-ntuple_noosc.root","RECREATE");
    TFile *fout = new TFile(filePath + "/files_root/RENO-ntuple_BFosc.root","RECREATE");
    TTree *T = new TTree("T","Monte Carlo neutrino events");
  
    float Ep, Ee, En, Ln; //Promt and neutrino energies; Baseline
    int   blid,ir,id,ad;
    float Prob_surv;

    T->Branch("Ep"  ,&Ep  ,"Ep/F");      //prompt energy
    T->Branch("Ee"  ,&Ee  ,"Ee/F");      //positron energy
    T->Branch("En"  ,&En  ,"En/F");	  //neutrino energy
    T->Branch("Ln"  ,&Ln  ,"Ln/F");	  //neutrino baseline
    T->Branch("blid",&blid,"blid/s"); //Baseline id (0-11)
  
    T->Branch("ir", &ir, "ir/s"); //reactor
    T->Branch("id", &id, "id/s"); //detector
    T->Branch("Prob_surv", &Prob_surv, "Prob_surv/F"); //Survival Probability

    double osc_prob_nd_allBins = 0.0;
    double N_events_allBins_nd = 0.0;
    double osc_prob_fd_allBins = 0.0;
    double N_events_allBins_fd = 0.0;

    //-- Energy resolution function (1610.04326) - 21.02.2020
    TF1 *gauE = new TF1("gauE","exp(-0.5*(x/[0])^2)",-2.0,2.0);
    double sigEp   = 0.0;
    double deltaEp = 0.0;
    
    int Nevents = atoi(getenv("NTUPLE_EVENTS"));
    for (int i = 0 ; i < Nevents ; i++){
        // generate a baseline (blid uniquely identifies the baseline)
        blid = histo_ldist_RENO_2x6->GetRandom();
        id =   (blid/nRea);
        ir =   (blid - id*nRea);
        Ln = baselines[id][ir] + gauL[id][ir]->GetRandom();
      
        if(id==0)  ad=0;
        else if (id==1) ad=1;
//        Ep = MC_spect_histo[ad]->GetRandom();
//        sigEp = Ep*(0.079/sqrt(Ep + 0.3));
//        gauE->SetParameter(0,sigEp);
//        deltaEp = gauE->GetRandom();
//        Ep = Ep+ deltaEp;
//        Ep = fFac1*Ep;
        do {
            Ep = MC_spect_histo[ad]->GetRandom();
            Ep = fFac1*Ep;
            sigEp = resFac*Ep*(0.079/sqrt(Ep + 0.3));
            gauE->SetParameter(0,sigEp);
            Ep = Ep + gauE->GetRandom();
        //            std::cout << "Cicle " << cicle << "\t 6AD: Ep After = " << Ep << std::endl;
                    //} while (Ep < 0.7);
        } while (Ep < 2*mel);

        if (i%10000 == 0) {
            std::cout << "NoOsc Event " << i << "   Ep = " << Ep << std::endl;
        }
        //-- We apply a incremental factor to the energy aiming to account
        //-- for an additional uncertainty on the neutrino energy and improve
        //-- our fit compared to the Collaboration's one
        //En = fFac*Ep + avg_nRecoilE + avg_constE;
        //En = fFac2*Ep + avg_nRecoilE + avg_constE;
        En = Ep + avg_nRecoilE + avg_constE;
        //En = Ep + avg_nRecoilE + avg_constE;
        
        //Using values from PRL121,201801 (2018) by RENO Coll.
        //-- ssq2th13RENO = 0.0896 +/- 0.0048(stat) +/- 0.0047(syst)
        //-- |dmsqeeRENO| =[2.68   +/- 0.12(stat)   +/- 0.07(syst)] x 10^{-3} eV^2
        //-- ssq2th12RENO = 0.307  +/- 0.013
        //-- dmsq21RENO   =[7.53   +/- 0.18] x 10^{-5} eV^2
        Prob_surv = 1.0 - ssq2th13RENO*pow(sin( 1.267 * dmsqeeRENO * Ln/En ),2) - pow(cos(0.5 * asin(sqrt(ssq2th13RENO))),4) * ssq2th12RENO * pow(sin( 1.267 * dmsq21RENO * Ln/En ),2);
        
        if( (xbins[0] <= Ep) && (Ep < xbins[NB]) ){
            if (id == 0) {
                osc_prob_nd_allBins += Prob_surv;
                N_events_allBins_nd++;
            }
            else {
                osc_prob_fd_allBins += Prob_surv;
                N_events_allBins_fd++;
            }
        }
        // Computing the average of the oscillation probability per bin
        for (int j = 0 ; j < NB ; j++){
            if( (xbins[j] <= Ep) && (Ep < xbins[j+1]) ){
                if(id==0){
                    osc_prob_nd[j] += Prob_surv;
                    N_events_perbin_nd[j]++;
                }
	      
                //if(id==1){
                else {
                    osc_prob_fd[j] += Prob_surv;
                    N_events_perbin_fd[j]++;
                }
            }
	  
        }
      T->Fill();
    }
    cout << "Done!" << endl;
    fout->Write();
    fout->Close();
    /*
     for (int j = 0 ; j < NB ; j++)
     {
         cout << " Posc_nd = " << osc_prob_nd[j] << "\tN_events = " << N_events_perbin_nd[j] << endl;
         cout << " Posc_fd = " << osc_prob_fd[j] << "\tN_events = " << N_events_perbin_fd[j] << endl;
     }
    */
    for (int j = 0 ; j < NB ; j++){
        osc_prob_nd[j] = osc_prob_nd[j]/N_events_perbin_nd[j];
        osc_prob_fd[j] = osc_prob_fd[j]/N_events_perbin_fd[j];
    }
    
    osc_prob_nd_allBins = osc_prob_nd_allBins/N_events_allBins_nd;
    osc_prob_fd_allBins = osc_prob_fd_allBins/N_events_allBins_fd;
    
    std::cout << "Osc Prob. ND = " << osc_prob_nd_allBins << std::endl;
    std::cout << "Osc Prob. FD = " << osc_prob_fd_allBins << std::endl;

    //-----------------------------
    for(int ib = 0; ib < NB; ib++){
        double binW2 = MC_spect_histo[1]->GetBinWidth(ib+1);
        double cont_nd = MC_spect_histo[0]->GetBinContent(ib+1);
        noosc_spect_histo[0]->SetBinContent(ib+1,cont_nd/osc_prob_nd[ib]);
        noosc_spect_histo_perMeV[0]->SetBinContent(ib+1,cont_nd/(0.2*binW2*osc_prob_nd[ib]));
      
        double cont_fd = MC_spect_histo[1]->GetBinContent(ib+1);
        noosc_spect_histo[1]->SetBinContent(ib+1,cont_fd/osc_prob_fd[ib]);
        noosc_spect_histo_perMeV[1]->SetBinContent(ib+1,cont_fd/(0.2*binW2*osc_prob_fd[ib]));
    }
  
    TH2F *frameND_perMeV = new TH2F("frameND_perMeV","",NB,1,hi,11,0,52500);
    frameND_perMeV->GetXaxis()->SetLabelSize(1.4*sz);
    frameND_perMeV->GetXaxis()->SetLabelFont(ft);
    frameND_perMeV->GetYaxis()->SetTitle("Events/0.2 MeV");
    frameND_perMeV->GetYaxis()->SetTitleFont(ft);
    frameND_perMeV->GetYaxis()->SetTitleOffset(0.7);
    frameND_perMeV->GetYaxis()->SetTitleSize(1.4*sz);
    frameND_perMeV->GetYaxis()->SetLabelSize(1.4*sz);
    frameND_perMeV->GetYaxis()->SetLabelFont(ft);
    TH2F *frameFD_perMeV = new TH2F("frameFD_perMeV","",NB,1,hi,11,0,6500);
    frameFD_perMeV->GetXaxis()->SetTitle("Prompt Reconstructed Energy (MeV)");
    frameFD_perMeV->GetXaxis()->SetTitleFont(ft);
    frameFD_perMeV->GetXaxis()->SetTitleOffset(0.9);
    frameFD_perMeV->GetXaxis()->SetTitleSize(1.4*sz);
    frameFD_perMeV->GetXaxis()->SetLabelSize(1.4*sz);
    frameFD_perMeV->GetXaxis()->SetLabelFont(ft);
    frameFD_perMeV->GetYaxis()->SetTitle("Events/0.2 MeV");
    frameFD_perMeV->GetYaxis()->SetTitleFont(ft);
    frameFD_perMeV->GetYaxis()->SetTitleOffset(0.7);
    frameFD_perMeV->GetYaxis()->SetTitleSize(1.4*sz);
    frameFD_perMeV->GetYaxis()->SetLabelSize(1.4*sz);
    frameFD_perMeV->GetYaxis()->SetLabelFont(ft);

    TH2F *frameND_events = new TH2F("frameND_events","",NB,1,hi,11,0,2100);
    frameND_events->GetXaxis()->SetLabelSize(1.4*sz);
    frameND_events->GetXaxis()->SetLabelFont(ft);
    frameND_events->GetYaxis()->SetTitle("Events");
    frameND_events->GetYaxis()->SetTitleFont(ft);
    frameND_events->GetYaxis()->SetTitleOffset(0.7);
    frameND_events->GetYaxis()->SetTitleSize(1.4*sz);
    frameND_events->GetYaxis()->SetLabelSize(1.4*sz);
    frameND_events->GetYaxis()->SetLabelFont(ft);
    TH2F *frameFD_events = new TH2F("frameFD_events","",NB,1,hi,11,0,260);
    frameFD_events->GetXaxis()->SetTitle("Prompt Reconstructed Energy (MeV)");
    frameFD_events->GetXaxis()->SetTitleFont(ft);
    frameFD_events->GetXaxis()->SetTitleOffset(0.9);
    frameFD_events->GetXaxis()->SetTitleSize(1.4*sz);
    frameFD_events->GetXaxis()->SetLabelSize(1.4*sz);
    frameFD_events->GetXaxis()->SetLabelFont(ft);
    frameFD_events->GetYaxis()->SetTitle("Events");
    frameFD_events->GetYaxis()->SetTitleFont(ft);
    frameFD_events->GetYaxis()->SetTitleOffset(0.7);
    frameFD_events->GetYaxis()->SetTitleSize(1.4*sz);
    frameFD_events->GetYaxis()->SetLabelSize(1.4*sz);
    frameFD_events->GetYaxis()->SetLabelFont(ft);


    TLatex *lat = new TLatex();
    lat->SetNDC();
    lat->SetTextFont(ft);
    lat->SetTextSize(1.8*sz);
    // Define the legend of the plots
    TLegend *leg1 = new TLegend(0.6,0.6,0.8,0.8);
    leg1->SetTextFont(ft);
    leg1->SetTextSize(1.4*sz);
    leg1->SetFillColor(0);
    leg1->SetLineColor(0);
    leg1->AddEntry(noosc_spect_histo_perMeV[0],"No oscillations");
    leg1->AddEntry(MC_spect_histo_perMeV[0],"RENO MC");
    
    TCanvas *canv0 = new TCanvas("canv0","Events/0.2 MeV",700,600);
    TGaxis::SetMaxDigits(3);
  
    TPad *pad1 = new TPad("pad1", "pad1", 0, 1./2., 1, 1.0);
    pad1->SetBottomMargin(0); // Upper and lower plot are joined
    pad1->Draw();             // Draw the upper pad: pad1
    pad1->cd();               // pad1 becomes the current pad
    frameND_perMeV->Draw();
    noosc_spect_histo_perMeV[0]->Draw("same h");
    MC_spect_histo_perMeV[0]->Draw("same h");
    lat->DrawLatex(0.6,0.4,"Near Detector");
    leg1->Draw();
    gPad->SetTicks(1,1);
  
    canv0->cd();          // Go back to the main canvas before defining pad2
    TPad *pad2 = new TPad("pad2", "pad2", 0, 0, 1, 1./2.);
    pad2->SetTopMargin(0);
    pad2->SetBottomMargin(0.11);
    pad2->Draw();
    pad2->cd();       // pad2 becomes the current pad
    frameFD_perMeV->Draw();
    noosc_spect_histo_perMeV[1]->Draw("same h");
    MC_spect_histo_perMeV[1]->Draw("same h");
    lat->DrawLatex(0.6,0.4,"Far Detector");
    gPad->RedrawAxis();
    gPad->SetTicks(1,1);
    
    canv0->Print(filePath + "/Plots/RENO_noosc-perMeV_spectra.pdf");
  
    TCanvas *canv1 = new TCanvas("canv1","Events",700,600);
    TGaxis::SetMaxDigits(3);
  
    TPad *pad10 = new TPad("pad10", "pad10", 0, 1./2., 1, 1.0);
    pad10->SetBottomMargin(0); // Upper and lower plot are joined
    pad10->Draw();             // Draw the upper pad: pad1
    pad10->cd();               // pad1 becomes the current pad
    frameND_events->Draw();
    noosc_spect_histo[0]->Draw("same h");
    MC_spect_histo[0]->Draw("same h");
    lat->DrawLatex(0.6,0.4,"Near Detector");
    leg1->Draw();
    gPad->SetTicks(1,1);
  
    canv1->cd();          // Go back to the main canvas before defining pad2
    TPad *pad20 = new TPad("pad20", "pad20", 0, 0, 1, 1./2.);
    pad20->SetTopMargin(0);
    pad20->SetBottomMargin(0.11);
    pad20->Draw();
    pad20->cd();       // pad2 becomes the current pad
    frameFD_events->Draw();
    noosc_spect_histo[1]->Draw("same h");
    MC_spect_histo[1]->Draw("same h");
    lat->DrawLatex(0.6,0.4,"Far Detector");
    gPad->RedrawAxis();
    gPad->SetTicks(1,1);
    
    canv1->Print(filePath + "/Plots/RENO_noosc-events_spectra.pdf");

     TFile *fout2 = new TFile(filePath + "/files_root/RENOplots_noosc.root","RECREATE");
     for (int j = 0 ; j < nDet ; j++) {
         noosc_spect_histo_perMeV[j]->Write();
         noosc_spect_histo[j]->Write();
     }
     fout2->Close();
    
    /*
     cout << endl;
     for (int j = 0 ; j < NB ; j++){
        cout << " Posc_nd = " << osc_prob_nd[j] << "\tPosc_fd = " << osc_prob_fd[j] << endl;
     }
     */

} // end
