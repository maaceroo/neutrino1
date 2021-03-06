//---------------------------------------------------------------------//
//--  db_ntuple.C - By M.A. Acero O. & A.A. Aguilar-A. - 2016-01-02  --//
//---------------------------------------------------------------------//
// This macro can be executed under ROOT typing                        //
// "root[0] .x db_ntuple.C"                                            //
// For the spectral analysis, using information from:                  //
// - F.P. An et al., PRL 112 061801 (2014)                             //
// - F.P. An et al., Chin.Phys.C37 011001 (2013)                       //
//---------------------------------------------------------------------//
// Simulated events, following the Daya Bay energy spectra from the    //
// paper PRL112.                                                       //
//---------------------------------------------------------------------//

#include "constants.h"

void db_ntuple()
{ // begin

    //-------------------
    // Energy Histograms
    //-------------------
    TFile *fenergy = new TFile("PRL112_217days_data.root","read");
    //Three sets of histograms one for each Experimental Hall
    const int nEH = 3;
    TH1F *nosc_spect_histo[nEH];
    TH1F *bkgd_spect_histo[nEH];
    TH1F *nu_nosc_spect_histo[nEH];
    for (int i = 0 ; i < nEH ; i++)
    {
        nosc_spect_histo[i] = (TH1F*) fenergy->Get(Form("nosc_spect_histo_%d",i));
        bkgd_spect_histo[i] = (TH1F*) fenergy->Get(Form("bkgd_spect_histo_%d",i));
        nu_nosc_spect_histo[i] = (TH1F*) nosc_spect_histo[i]->Clone();
        nu_nosc_spect_histo[i]->Add(bkgd_spect_histo[i],-1.0);
    }
    //-------------------
    // Distance Histogram
    //-------------------
    TFile *fpathl = new TFile("files_data/daya-bay-ldist6AD.root","read");
    //Only 6 detectors (PRL112 061801, 2014)
    TH1F *histo_ldist_6Det = (TH1F*) fpathl->Get("histo_ldist_6Det");

    const int nDet = 8;
    const int nRea = 6;
    //Baseline Distances (cm)
    const char  *detNames[nDet] = {"EH1-AD1", "EH1-AD2", "EH2-AD1", "EH2-AD2",
          	                   "EH3-AD1", "EH3-AD2", "EH3-AD3", "EH3-AD4"};
    //The 48 baselines in Daya-Bay ()
    double baselines[nDet][nRea] =
    {
        {362.380, 371.763,  903.466, 817.158,1353.618,1265.315},
        {357.940, 368.414,  903.347, 816.896,1354.229,1265.886},
        {1332.479,1358.148, 467.574, 489.577, 557.579, 499.207},
        {1337.429,1362.876, 472.971, 495.346, 558.707, 501.071},
        {1919.632,1894.337,1533.180,1533.628,1551.384,1524.940},
        {1917.519,1891.977,1534.919,1535.032,1554.767,1528.079},
        {1925.255,1899.861,1538.930,1539.468,1556.344,1530.078},
        {1923.149,1897.507,1540.667,1540.872,1559.721,1533.179}
    };

    //make ntuple
    TFile *fout = new TFile("files_data/db6AD-ntuple.root","RECREATE");
    TTree *T = new TTree("T","Monte Carlo neutrino events");

    float Ep, En, Ln;
    int   blid,ir,id,ad;
    T->Branch("Ep"  ,&Ep  ,"Ep/F");	  //prompt reconstructed energy
    T->Branch("En"  ,&En  ,"En/F");	  //neutrino energy
    T->Branch("Ln"  ,&Ln  ,"Ln/F");	  //neutrino baseline
    T->Branch("blid",&blid,"blid/s"); //Baseline id (0-47)
    
    T->Branch("ir", &ir, "ir/s"); //reactor
    T->Branch("id", &id, "id/s"); //detector

    //int Nevents = 1000000;
    int Nevents = atoi(getenv("NTUPLE_EVENTS"));
    printf("Ntuple Events: %d \n",Nevents);
    
    //-- 2018.12.21 -
    //-- Gaussian distribution to include the effect of the detectors and reactors dimensions
    TF1 *gau = new TF1("gau","exp(-0.5*(x/[0])^2)",-30.0,30.0);
    gau->SetParameter(0,5);
    for (int i = 0 ; i < Nevents ; i++)
        {
            // generate a baseline (blid uniquely identifies the baseline)
            blid = histo_ldist_6Det->GetRandom();
            id =   (blid/nRea);
            ir =   (blid - id*nRea);
            Ln = baselines[id][ir] + gau->GetRandom();

            // generate a neutrino energy
            //if (id is 0 or 1, espectrum from EH1; id is 2, espectrum form EH2; id is 4 to 6,espectrum from EH3) this is for Ep
            if (id < 2)       ad = 0;
            else if (id == 2) ad = 1;
            else if (id > 3 && id < 7) ad = 2;
            
            Ep = nu_nosc_spect_histo[ad]->GetRandom();
            En = Ep*1.03 + avg_nRecoilE + avg_constE;
            //En = Ep*1.0 + avg_nRecoilE + avg_constE;

            T->Fill();
            
            if(i%1000000 == 0)
                cout << "Number of events " << i << " done!" << endl;
        }
    
    fout->Write();

} // end
