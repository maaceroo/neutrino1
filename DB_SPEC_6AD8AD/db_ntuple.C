//---------------------------------------------------------------------//
//--  db_ntuple.C - By M.A. Acero O. & A.A. Aguilar-A. - 2019-10-25  --//
//---------------------------------------------------------------------//
// This macro can be executed under ROOT typing                        //
// "root[0] .x db_ntuple.C"                                            //
// For the spectral analysis, using information from:                  //
// - F.P. An et al., PRD 95 072006 (2017)                             //
// - F.P. An et al., Chin.Phys.C37 011001 (2013)                       //
//---------------------------------------------------------------------//
// Simulated events, following the Daya Bay energy spectra from the    //
// paper PRD95.                                                        //
//---------------------------------------------------------------------//

#include "constants.h"

void db_ntuple()
{ // begin

    //-------------------
    // Energy Histograms
    //-------------------
    //TFile *fenergy8AD = new TFile("PRD95_1230days_data.root","read");
    //TFile *fenergy6AD = new TFile("PRL112_217days_data.root","read");
    std::cout << "Calling the root file" << std::endl;
    TFile *fenergy = new TFile("histos_6AD217Days_8AD1013Days_data.root","read");
    //Three sets of histograms one for each Experimental Hall
    const int nEH = 3;
    //1230 days
    std::cout << "Creating histograms" << std::endl;
    TH1F *nosc_spect_histo[nEH];
    TH1F *bkgd_spect_histo[nEH];
    TH1F *nu_nosc_spect_histo[nEH];
    //217 days - 6AD
    TH1F *nosc_spect_histo6AD[nEH];
    TH1F *bkgd_spect_histo6AD[nEH];
    TH1F *nu_nosc_spect_histo6AD[nEH];
    //1013 days - 8AD
    TH1F *nosc_spect_histo8AD[nEH];
    TH1F *bkgd_spect_histo8AD[nEH];
    TH1F *nu_nosc_spect_histo8AD[nEH];
    std::cout << "Filling histograms" << std::endl;
    for (int i = 0 ; i < nEH ; i++)
    {
        std::cout << "EH " << i << std::endl;
        nosc_spect_histo[i]    = (TH1F*) fenergy->Get(Form("nosc_spect_6AD8AD_histo_%d",i));
        bkgd_spect_histo[i]    = (TH1F*) fenergy->Get(Form("bkgd_spect_6AD8AD_histo_%d",i));
        nu_nosc_spect_histo[i] = (TH1F*) nosc_spect_histo[i]->Clone();
        nu_nosc_spect_histo[i]->Add(bkgd_spect_histo[i],-1.0);
        
        nosc_spect_histo6AD[i]    = (TH1F*) fenergy->Get(Form("nosc_spect_6AD35B_histo_%d",i));
        bkgd_spect_histo6AD[i]    = (TH1F*) fenergy->Get(Form("bkgd_spect_6AD35B_histo_%d",i));
        nu_nosc_spect_histo6AD[i] = (TH1F*) nosc_spect_histo6AD[i]->Clone();
        nu_nosc_spect_histo6AD[i]->Add(bkgd_spect_histo6AD[i],-1.0);

        nosc_spect_histo8AD[i]    = (TH1F*) fenergy->Get(Form("nosc_spect_8AD35B_histo_%d",i));
        bkgd_spect_histo8AD[i]    = (TH1F*) fenergy->Get(Form("bkgd_spect_8AD35B_histo_%d",i));
        nu_nosc_spect_histo8AD[i] = (TH1F*) nosc_spect_histo8AD[i]->Clone();
        nu_nosc_spect_histo8AD[i]->Add(bkgd_spect_histo8AD[i],-1.0);
    }

    //-------------------
    // Distance Histogram
    //-------------------
    TFile *fpathl = new TFile("files_data/daya-bay-ldist.root","read");
    //8 detectors (PRD95 072006, 2017)
    TH1F *histo_ldist = (TH1F*) fpathl->Get("histo_ldist");
    //6 detectors (PRL112 )
    TH1F *histo_ldist_6Det = (TH1F*) fpathl->Get("histo_ldist_6Det");

    //const int nDet = 8;
    //const int nRea = 6;
    const int nDet = nAD;
    const int nRea = nNR;
    //Baseline Distances (cm)
    const char  *detNames[nDet] = {"EH1-AD1", "EH1-AD2", "EH2-AD1", "EH2-AD2", "EH3-AD1", "EH3-AD2", "EH3-AD3", "EH3-AD4"};
    //The 48 baselines in Daya-Bay ()
    double baselines[nDet][nRea] =
    {
        { 362.380, 371.763, 903.466, 817.158,1353.618,1265.315},
        { 357.940, 368.414, 903.347, 816.896,1354.229,1265.886},
        {1332.479,1358.148, 467.574, 489.577, 557.579, 499.207},
        {1337.429,1362.876, 472.971, 495.346, 558.707, 501.071},
        {1919.632,1894.337,1533.180,1533.628,1551.384,1524.940},
        {1917.519,1891.977,1534.919,1535.032,1554.767,1528.079},
        {1925.255,1899.861,1538.930,1539.468,1556.344,1530.078},
        {1923.149,1897.507,1540.667,1540.872,1559.721,1533.179}
    };
    //-- 01-08-2020 Trying to add different uncertainties for each baseline --
    TFile *fBLsdev = new TFile("files_data/db_BL_deviations.root","read");
    TH1F *bldist[nDet][nRea];
    for (int ii = 0 ; ii < nDet ; ii++) {
        for (int jj = 0 ; jj < nRea ; jj++) {
            bldist[ii][jj] = (TH1F*) fBLsdev->Get(Form("bldist_%d_%d",ii,jj));
            //std::cout << "\t Mean Distance" << "(" << ii << ", " << jj << ") = " << bldist[ii][jj]->GetMean()
            //<< "\t Deviation" << "(" << ii << ", " << jj << ") = " << bldist[ii][jj]->GetStdDev() << std::endl;
        }
    }
    //-- 07-08-2020 Independent Gaussian fluctuation for each baseline with different StdDev
    TF1 *gauL[nDet][nRea];
    for (int iD = 0 ; iD < nDet ; iD++)
        for (int iR = 0 ; iR < nRea ; iR++){
            gauL[iD][iR] = new TF1(Form("gauL_%d_%d",iD,iR),"exp(-0.5*(x/[0])^2)",-10.0,10.0);
            gauL[iD][iR]->SetParameter(0,bldist[iD][iR]->GetStdDev());
        }

    //make ntuple
    TFile *fout = new TFile("files_data/db-ntuple_noER_1M.root","RECREATE");
    TTree *T6AD8AD = new TTree("T6AD8AD","Monte Carlo neutrino events 6AD+8AD");

    double frac6AD = 217.0/1230.; // fraction of the total corresponding to 6AD only

    float Ep, Ee, En, Ln;
    int   blid,ir,id,ad;
    int   per;
    T6AD8AD->Branch("Ep"  ,&Ep  ,"Ep/F");      //prompt reconstructed energy
    T6AD8AD->Branch("En"  ,&En  ,"En/F");	  //neutrino energy
    T6AD8AD->Branch("Ln"  ,&Ln  ,"Ln/F");	  //neutrino baseline
    T6AD8AD->Branch("blid",&blid,"blid/s"); //Baseline id (0-47)

    T6AD8AD->Branch("ir", &ir, "ir/s"); //reactor
    T6AD8AD->Branch("id", &id, "id/s"); //detector
    T6AD8AD->Branch("per",&per,"per/s"); //period

    //int Nevents = 5000000; // CAUTION!! This must be commented out when using the script
    int Nevents = atoi(getenv("NTUPLE_EVENTS")); // This must be uncommented when using the script
    printf("Ntuple Events: %d \n",Nevents);
    
    //-- Prompt Energy resolution function
    TF1 *gauE = new TF1("gauE","exp(-0.5*(x/[0])^2)",-2.0,2.0);
    double sigEp   = 0.0;
    double deltaEp = 0.0;
    
    //TF1 *Unif = new TF1("Unif","1",0.7,1.3);
    TF1 *Unif = new TF1("Unif","1",2*mel,1.3);
    
    std::cout << "\n";
    std::cout << "Filling ntuple for 6AD+8AD combined analysis\n";
    //- Fill Ntuple for 6AD+8AD combined analysis
    for (int i = 0 ; i < Nevents ; i++)
    {
        // generate a baseline (blid uniquely identifies the baseline)
        if (i < frac6AD*Nevents) {
            blid = histo_ldist_6Det->GetRandom();
            per = 1;
        } else {
            blid = histo_ldist->GetRandom();
            per = 2;
        }
        id =   (blid/nRea);
        ir =   (blid - id*nRea);
        Ln = baselines[id][ir] + gauL[id][ir]->GetRandom();

        if ( i<frac6AD*Nevents ) {
            // generate a neutrino energy
            //if (id is 0 or 1, EH1 spectrum; id is 2, EH2 spectrum; id is 4 to 6, EH3 spectrum) this is for Ep
            if (id < 2)       ad = 0;
            else if (id == 2) ad = 1;
            else if (id > 3 && id < 7) ad = 2;
     
            Ep = nu_nosc_spect_histo6AD[ad]->GetRandom();
//            sigEp = Ep*(sqrt(pow(0.016,2) + pow(0.081,2)/Ep + pow(0.026/Ep,2))); //PRD 95,072006,2017
//            gauE->SetParameter(0,sigEp);
//            Ep = Ep + gauE->GetRandom();
            while (Ep < 2*mel){
                Ep = Unif->GetRandom();
//                sigEp = Ep*(sqrt(pow(0.016,2) + pow(0.081,2)/Ep + pow(0.026/Ep,2))); //PRD 95,072006,2017
//                gauE->SetParameter(0,sigEp);
//                Ep = Ep + gauE->GetRandom();
            }
            sigEp = Ep*(sqrt(pow(0.016,2) + pow(0.081,2)/Ep + pow(0.026/Ep,2))); //PRD 95,072006,2017
            gauE->SetParameter(0,sigEp);
            Ep = Ep + gauE->GetRandom();
//            while (Ep < 2*mel)
//                Ep = Ep + gauE->GetRandom();

            En = Ep + avg_nRecoilE + avg_constE;

        } //if 6AD period events
        else {
            // generate a neutrino energy
            //if (id is 0 or 1, EH1 spectrum; id is 2 or 3, EH2 spectrum; id is 4 to 7, EH3 spectrum) this is for Ep
            if (id < 2)       ad = 0;
            else if (id == 2 || id == 3) ad = 1;
            else if (id > 3) ad = 2;
    
            Ep = nu_nosc_spect_histo8AD[ad]->GetRandom();
            //Here: Apply fFac
//            sigEp = Ep*(sqrt(pow(0.016,2) + pow(0.081,2)/Ep + pow(0.026/Ep,2))); //PRD 95,072006,2017
//            gauE->SetParameter(0,sigEp);
//            Ep = Ep + gauE->GetRandom();
//            while (Ep < 2*mel){
//                Ep = Unif->GetRandom();
//                sigEp = Ep*(sqrt(pow(0.016,2) + pow(0.081,2)/Ep + pow(0.026/Ep,2))); //PRD 95,072006,2017
//                gauE->SetParameter(0,sigEp);
//                Ep = Ep + gauE->GetRandom();
//            }
//
            En = Ep + avg_nRecoilE + avg_constE;

        }//if 8AD period events

        T6AD8AD->Fill();
            
        if(i%(Nevents/20) == 0)
            cout << "Number of events 6AD+8AD ntuple " << i << " done!" << endl;
    }//for Nevents
    
    fout->Write();

} // end
