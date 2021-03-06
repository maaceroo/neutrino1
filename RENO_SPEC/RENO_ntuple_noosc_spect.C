//------------------------------------------------------------------------------------//
//--     RENO_ntuple_noosc_spect.C - By M.A. Acero O. & D.J. Polo - 2019-07-10      --//
//------------------------------------------------------------------------------------//
// This macro can be executed under ROOT typing                        		          //
// "root[0] .x RENO_ntuple_noosc_spect.C"                                             //
// or, from a Terminal, typing
// "root -n -l -b -q RENO_ntuple_noosc_spect.C"                                       //
// It should be executed independently of the oter macros related to the analysis.    //
//------------------------------------------------------------------------------------//
// For the spectral analysis, using information from:                  		          //
// - F.P. An et al.,  RENO Coll. (2017) - arXiv:1610.04326                            //
// - F.P. An et al., arXiv:1003.1391" 6 Mar 2010. (Table 1.2)              	          //
//------------------------------------------------------------------------------------//
// With this macro we build a spectrum for each RENO antineutrino detector, filled    //
// simulated no-oscillated neutrinos events. To do so, we use the MC oscillated       //
// expected spectra reported by RENO in 1610.04326 to sample the simulated events.    //
// We then compute the average neutrino oscillation probability and "un-oscilate" the //
// MC spectra using the Collaboration's best fit parameters.                          //
//------------------------------------------------------------------------------------//
// The resulting no-osccilated spectra are used in our RENO data analysis.            //
// No-oscillated espectra (events and events / 0.2 Mev) are saved to the file         //
// "files_root/RENOplots_noosc.root"                                                  //
//------------------------------------------------------------------------------------//

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
    const int       NB = 27;
    const double    lo = 1.2;
    const double    hi = 8.4;
    double xbins[NB+1];
    double delta_bins2 = (6.0 - 1.2)/24; // 0.2 MeV/bin
  
    for (int i = 0 ; i < (NB-2) ; i++){
        xbins[i] = 1.2 + delta_bins2*i;
    }
    xbins[25] = xbins[24] + 0.4;
    xbins[26] = 8.4 - 1.4;
    xbins[27] = 8.4;
    double osc_prob_nd[NB] = {0.0};
    double osc_prob_fd[NB] = {0.0};
    double N_events_perbin_nd[NB] = {0.0};
    double N_events_perbin_fd[NB] = {0.0};
  
    //---------------------------------------------------
    //-------------------
    // Energy Histograms
    //-------------------
    TFile *fenergy = new TFile("files_root/RENOplots.root","read");
  
    //The ND and FD MC spectra
    //From PRL116 211801 (2016). Figure 2  (Blue lines).
    const int nd = 2; // Number of detectors
    TH1F *MC_spect_histo[nd];
    TH1F *MC_spect_histo_perMeV[nd];
    TH1F *noosc_spect_histo[nd];
    TH1F *noosc_spect_histo_perMeV[nd];
    for(int n = 0 ; n < nd ; n++){
        MC_spect_histo_perMeV[n] = (TH1F*) fenergy->Get(Form("spect_histo_%d",n));  // Get data
        MC_spect_histo_perMeV[n]->SetLineColor(kBlue);
        MC_spect_histo[n] = new TH1F(Form("MC_spect_histo_%d",n),"",NB,xbins);
        MC_spect_histo[n]->SetLineColor(kBlue);
        noosc_spect_histo[n] = new TH1F(Form("noosc_spect_histo_%d",n),"",NB,xbins);
        noosc_spect_histo[n]->SetLineColor(kRed);
        noosc_spect_histo_perMeV[n] = new TH1F(Form("noosc_spect_histo_perMeV_%d",n),"",NB,xbins);
        noosc_spect_histo_perMeV[n]->SetLineColor(kRed);
    }
    
    int NumB = MC_spect_histo_perMeV[1]->GetNbinsX();
    for(int n = 0 ; n < nd ; n++){
        for(int ib = 0; ib < NB; ib++){
            double binW = MC_spect_histo[1]->GetBinWidth(ib+1);
            double cont =  MC_spect_histo_perMeV[n]->GetBinContent(ib+1);
            MC_spect_histo[n]->SetBinContent(ib+1,cont*binW*0.2);  // Get data
        }
    }
  
    //-------------------
    // Distance Histogram
    //-------------------
    TFile *fpathl = new TFile("files_root/ldist_RENO_2x6.root","read");
    // Get histogram - histo_ldist_RENO_2x6
    TH1F *histo_ldist_RENO_2x6 = (TH1F*) fpathl->Get("histo_ldist_RENO_2x6");
  
    const int nDet = 2;   // number of AD at Reno
    const int nRea = 6;   // number of reactors
    //Baseline Distances (m)
    const char  *detNames[nDet] = {"near-AD1", "far-AD2"};
  
    //The 12 baselines in RENO AD*NR = 12
    double baselines[nDet][nRea] =
    {
        { 667.9, 451.8, 304.8, 336.1, 513.9, 739.1},
        {1556.5,1456.2,1395.9,1381.3,1413.8,1490.1}
    };
  
    //make ntuple
    TFile *fout = new TFile("files_root/RENO-ntuple_noosc.root","RECREATE");
    TTree *T = new TTree("T","Monte Carlo neutrino events");
  
    //Dfining some important quantities and variables
    double avg_nRecoilE = 10.0e-3; //MeV
    double avg_constE = 0.78; //MeV
    float Ep, En, Ln; //Promt and neutrino energies; Baseline
    int   blid,ir,id,ad;
    float Prob_surv;

    T->Branch("Ep"  ,&Ep  ,"Ep/F");	  //prompt energy
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

    int Nevents = 5000000;
    for (int i = 0 ; i < Nevents ; i++){
        // generate a baseline (blid uniquely identifies the baseline)
        blid = histo_ldist_RENO_2x6->GetRandom();
        id =   (blid/nRea);
        ir =   (blid - id*nRea);
        Ln = baselines[id][ir];
      
        if(id==0)  ad=0;
        else if (id==1) ad=1;
        Ep = MC_spect_histo[ad]->GetRandom();
        En = Ep + avg_nRecoilE + avg_constE;
        //En = Mn + Ep - Mp ; // Neutrino energy. Where Mp and Mn are the proton and neutron masses
        
        Prob_surv = 1.0 - 0.087*pow(sin( 1.267 * 2.49e-3 * Ln/En ),2) - pow(cos(0.5 * asin(sqrt(0.087))),4) * 0.846 * pow(sin( 1.267 * 7.53e-5 * Ln/En ),2);
        
        if (id == 0) {
            osc_prob_nd_allBins += Prob_surv;
            N_events_allBins_nd++;
        }
        else {
            osc_prob_fd_allBins += Prob_surv;
            N_events_allBins_fd++;
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

    //Adding artificial noice to the no-osccilated spectra
    TF1 *gauNoise = new TF1("gauNoise","exp(-0.5*((x)/[0])^2)",-5.0,5.0);
    gauNoise->SetParameter(0,1);
    //-----------------------------
    for(int ib = 0; ib < NB; ib++){
        //double noise = gauNoise->GetRandom();
        
        double binW2 = MC_spect_histo[1]->GetBinWidth(ib+1);
        double cont_nd = MC_spect_histo[0]->GetBinContent(ib+1);
        //double cont_nd = MC_spect_histo[0]->GetBinContent(ib+1) + noise;
        noosc_spect_histo[0]->SetBinContent(ib+1,cont_nd/osc_prob_nd[ib]);
        noosc_spect_histo_perMeV[0]->SetBinContent(ib+1,cont_nd/(0.2*binW2*osc_prob_nd[ib]));
      
        double cont_fd = MC_spect_histo[1]->GetBinContent(ib+1);
        //double cont_fd = MC_spect_histo[1]->GetBinContent(ib+1) + noise;
        noosc_spect_histo[1]->SetBinContent(ib+1,cont_fd/osc_prob_fd[ib]);
        noosc_spect_histo_perMeV[1]->SetBinContent(ib+1,cont_fd/(0.2*binW2*osc_prob_fd[ib]));
    }
  
    TH2F *frameND_perMeV = new TH2F("frameND_perMeV","",NB,1,hi,11,0,18300);
    frameND_perMeV->GetXaxis()->SetLabelSize(1.4*sz);
    frameND_perMeV->GetXaxis()->SetLabelFont(ft);
    frameND_perMeV->GetYaxis()->SetTitle("Events/0.2 MeV");
    frameND_perMeV->GetYaxis()->SetTitleFont(ft);
    frameND_perMeV->GetYaxis()->SetTitleOffset(0.7);
    frameND_perMeV->GetYaxis()->SetTitleSize(1.4*sz);
    frameND_perMeV->GetYaxis()->SetLabelSize(1.4*sz);
    frameND_perMeV->GetYaxis()->SetLabelFont(ft);
    TH2F *frameFD_perMeV = new TH2F("frameFD_perMeV","",NB,1,hi,11,0,2100);
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

    TH2F *frameND_events = new TH2F("frameND_events","",NB,1,hi,11,0,730);
    frameND_events->GetXaxis()->SetLabelSize(1.4*sz);
    frameND_events->GetXaxis()->SetLabelFont(ft);
    frameND_events->GetYaxis()->SetTitle("Events");
    frameND_events->GetYaxis()->SetTitleFont(ft);
    frameND_events->GetYaxis()->SetTitleOffset(0.7);
    frameND_events->GetYaxis()->SetTitleSize(1.4*sz);
    frameND_events->GetYaxis()->SetLabelSize(1.4*sz);
    frameND_events->GetYaxis()->SetLabelFont(ft);
    TH2F *frameFD_events = new TH2F("frameFD_events","",NB,1,hi,11,0,90);
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
    
    canv0->Print("Plots/RENO_noosc-perMeV_spectra.pdf");
  
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
    
    canv1->Print("Plots/RENO_noosc-events_spectra.pdf");

     TFile *fout2 = new TFile("files_root/RENOplots_noosc.root","RECREATE");
     for (int j = 0 ; j < nd ; j++) {
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
