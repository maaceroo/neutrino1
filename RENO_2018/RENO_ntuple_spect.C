//---------------------------------------------------------------------------//
//- RENO_ntuple_spec.C - By M.A. Acero O. & A.A. Aguilar-A. -- 2020-02-07  --//
//---------------------------------------------------------------------------//
// This macro can be executed under ROOT typing                        		 //
// "root[0] .x RENO_ntuple.C"                                            	 //
// For the spectral analysis, using information from:                  		 //
// - F.P. An et al.,  RENO Coll. (2017) - arXiv:1610.04326                   //
// - F.P. An et al., arXiv:1003.1391" 6 Mar 2010. (Table 1.2)              	 //
//---------------------------------------------------------------------------//
// Simulated events, following the RENO energy spectra from the  		     //
// G. Bak et al., RENO Coll. (2018) PRL121, 201801                           //
//---------------------------------------------------------------------------//

#include "constants.h"

void RENO_ntuple_spect()
{ // begin
    std::cout << "Begin..." << std::endl;
  
    //-------------------
    // Energy Histograms
    //-------------------
    //TFile *fenergy = new TFile("files_root/RENOplots.root","read");
    TFile *fenergy = new TFile("files_root/RENOplots_noosc.root","read");

    //The histogram of near and far data spectra
    const int nd = 2; // Number of detectors
    TH1F *noosc_spect_histo[nd];
    for(int n = 0 ; n < nd ; n++){
        noosc_spect_histo[n] = (TH1F*) fenergy->Get(Form("noosc_spect_histo_%d",n)); // Get data
    }
    //TH1F *noosc_spect_histo[nd];
    //noosc_spect_histo[0] = (TH1F*) fenergy->Get("data_spect_histo_0");  // Set to observed Near spectrum
    //noosc_spect_histo[1] = (TH1F*) fenergy->Get("reno_noosc_histo");  // Set to predicted noOsc in Far spectrum
    
    //std::cout << "Histograms..." << std::endl;
    //-- Multilping by the bin width - 01.02.2019
    //for(int n = 1 ; n < nd ; n++){
    //int NumB = noosc_spect_histo[1]->GetNbinsX();
    //for(int i = 0 ; i < NumB ; i++){
        //double binW = noosc_spect_histo[1]->GetBinWidth(i+1);
        //double cont = noosc_spect_histo[1]->GetBinContent(i+1);
        //noosc_spect_histo[1]->SetBinContent(i+1,0.2*binW*cont);
    //}

    //-------------------
    // Distance Histogram
    //-------------------
    TFile *fpathl = new TFile("files_root/ldist_RENO.root","read");
    // Get histogram - histo_ldist_RENO_2x6
    TH1F *histo_ldist_RENO_2x6 = (TH1F*) fpathl->Get("histo_ldist_RENO_2x6");
  
    //const int nDet = 2;   // number of ad at Reno
    //const int nRea = 6;   // number of reactors
    //Baseline Distances (m)
    const char  *detNames[nDet] = {"near-AD1", "far-AD2"};
  
    //The 12 baselines in RENO ad*nr = 12 ()
    double baselines[nDet][nRea] =
        {
            {667.9,451.8,304.8,336.1,513.9,739.1},
            {1556.5,1456.2,1395.9,1381.3,1413.8,1490.1}
        };
  
    //make ntuple
    std::cout << "Ntuple creation..." << std::endl;
    TFile *fout = new TFile("files_root/RENO-ntuple.root","RECREATE");
    TTree *T = new TTree("T","Monte Carlo neutrino events");
  
    float Ep, En, Ln;
    int   blid,ir,id,ad;
    T->Branch("Ep"  ,&Ep  ,"Ep/F");	  //prompt energy
    T->Branch("En"  ,&En  ,"En/F");	  //neutrino energy
    T->Branch("Ln"  ,&Ln  ,"Ln/F");	  //neutrino baseline
    T->Branch("blid",&blid,"blid/s"); //Baseline id (0-11)
  
    T->Branch("ir", &ir, "ir/s"); //reactor
    T->Branch("id", &id, "id/s"); //detector
  
    //-- A normal distribution to add a random quantity to the baseline Ln
    //-- in order to account for the reactor and detector sizes
    TF1 *gau = new TF1("gau","exp(-0.5*(x/[0])^2)",-5.0,5.0);
    gau->SetParameter(0,1);
    //-- Energy resolution function (1610.04326) - 21.02.2020
    TF1 *gauE = new TF1("gauE","exp(-0.5*(x/[0])^2)",-2.0,2.0);
    double sigEp   = 0.0;
    double deltaEp = 0.0;
    //int Nevents = 1000000; // CAUTION!! This must be commented out when using the script
    int Nevents = atoi(getenv("NTUPLE_EVENTS")); // This must be uncommented when using the script
    for (int i = 0 ; i < Nevents ; i++)
        {
            // generate a baseline (blid uniquely identifies the baseline)
            blid = histo_ldist_RENO_2x6->GetRandom();
            id =   (blid/nRea);
            ir =   (blid - id*nRea);
            Ln = baselines[id][ir] + gau->GetRandom();

            if(id==0)  ad=0;
            else if (id==1) ad=1;
            //Ep = data_spect_histo[ad]->GetRandom();
            Ep = noosc_spect_histo[ad]->GetRandom();
            sigEp = Ep*(0.079/sqrt(Ep + 0.3));
            gauE->SetParameter(0,sigEp);
            deltaEp = gauE->GetRandom();
            Ep = Ep + deltaEp;
            //-- We apply a incremental factor to the energy aiming to account
            //-- for an additional uncertainty on the neutrino energy and improve
            //-- our fit compared to the Collaboration's one
            En = fFac*Ep + avg_nRecoilE + avg_constE;
            if (i%10000 == 0)
                std::cout << " Oscillated Event " << i << "   Ep = " << Ep << std::endl;
            //En = Mn + Ep - Mp ; // Neutrino energy. Where Mp and Mn are the proton and neutron masses
            
            T->Fill();
        }
    
    std::cout << "Ntuple creation... DONE" << std::endl;
    fout->Write();
  
} // end
