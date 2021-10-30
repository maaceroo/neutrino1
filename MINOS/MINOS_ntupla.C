//-----------------------------------------------------------------------//
//--  MINOS_ntuple.C - By M.A. Acero O. & A.A. Aguilar-A. - 2021-9-03  --//
//-----------------------------------------------------------------------//
// This macro can be executed under ROOT typing                          //
// "root[0] .x MINOS_ntuple.C"                                           //
// Muon and Electron (Anti)Neutrinos                                     //
// - MINOS Col. PRL 110, 171801 (2013) - NuE appearance                  //
// - MINOS Col. PRL 110, 251801 (2013) - NuMu disappearance              //
//-----------------------------------------------------------------------//

#include "constants.h"
#include <math.h>
#include <iostream>
#include <string>

void MINOS_ntupla()
{ // begin

    //------------- Style --------------
    gROOT->SetStyle("Plain");
    gStyle->SetOptStat(0);
    //------------- Style --------------

    //-------------------
    // Energy Histograms
    //-------------------
    std::cout << "Calling the root file" << std::endl;
    TFile *fenergy = new TFile("MINOS_spectra_PRL108-PRL110.root","read");

    // define and fill the histograms for Muon (Anti)neutrinos PRL110 ------------------------
    std::cout << "Creating and Filling histograms" << std::endl;
    //- Muon neutrinos
    TH1F *numu110_noosc_histo;
    TH1F *numu110_bkgd_histo;
    numu110_noosc_histo = (TH1F*) fenergy->Get("numu110_noosc_histo");
    numu110_bkgd_histo  = (TH1F*) fenergy->Get("numu110_bkgd_histo");
    //- Muon antineutrinos
    TH1F *numub110_noosc_histo;
    TH1F *numub110_bkgd_histo;
    numub110_noosc_histo = (TH1F*) fenergy->Get("numub110_noosc_histo");
    numub110_bkgd_histo  = (TH1F*) fenergy->Get("numub110_bkgd_histo");
    

    //Root file for the ntuple
    TString filePath = dirName;
    TFile *fout   = new TFile(filePath + "/data/minos_ntuple.root","RECREATE");
    TTree *Tnumu  = new TTree("Tnumu" ,"MC neutrino events");
    TTree *Tnumub = new TTree("Tnumub","MC neutrino events");
    //-------------------
    // True Vs. Reco Energy Matrix
    //-------------------
    std::cout << "Calling the Matrix root file" << std::endl;
    TFile *fmatrix = new TFile("./data/MINOS_Matrix_TvsR_Energy.root","read");
    TH2F *TrueReco_Matrix;
    TrueReco_Matrix = (TH2F*) fmatrix->Get("TRmat_histo");
    TH1D *TrueReco_px;
    TH1D *TrueReco_py;
    TrueReco_px = TrueReco_Matrix->ProjectionX("TrueReco_px");
    TrueReco_py = TrueReco_Matrix->ProjectionY("TrueReco_py");
    //fout->cd();
    //TrueReco_px->Write();
    //TrueReco_py->Write();
    int Nbins   = 80;
    int n_histo = Nbins;
    TH1F *True_array[80];
    for (int j = 0 ; j < n_histo ; j++) {
        True_array[j] = new TH1F(Form("True_array_%d",j),"",Nbins,0,20);
        for (int i = 0 ; i < Nbins ; i++) {
            double value = TrueReco_Matrix->GetBinContent(i+1,j+1);
            True_array[j]->SetBinContent(i+1,value);
        }
//        fout->cd();
//        True_array[j]->Write();
    }

    float Ereco, Etrue;
    float Erecob, Etrueb;
    float BL = 735.0; // km
    //int Nevents = 1e7;
    int Nevents = atoi(getenv("NTUPLE_EVENTS")); // This must be uncommented when using the script


    Tnumu-> Branch("Ereco", &Ereco, "Ereco/F" );
    Tnumu-> Branch("Etrue", &Etrue, "Etrue/F" );
    Tnumub->Branch("Erecob",&Erecob,"Erecob/F");
    Tnumub->Branch("Etrueb",&Etrueb,"Etrueb/F");
    
    double e10   = 10.0;
    double sig10 = 1.2;
    TF1 *gauE = new TF1("gauE","exp(-0.5*(x/[0])^2)",-20.0,20.0);
    double sigEr   = 0.0;
    double deltaEr = 0.0;

    double ehi = 14.0;
    TH2F *histoNu  = new TH2F("histoNu", "", 56,0,ehi,56,0,ehi);
    TH2F *histoNub = new TH2F("histoNub", "", 56,0,ehi,56,0,ehi);

    int recoBin;
    for (int i = 0 ; i < Nevents ; i++) {
        Ereco  = numu110_noosc_histo->GetRandom();
        recoBin = int(Ereco/0.25);
        //sigEr  = sig10*sqrt(Ereco/e10);
//        sigEr  = (0.043 + 0.213*Ereco - 0.0051*pow(Ereco,2))/(2*sqrt(2*log(2)));
//        gauE->SetParameter(0,sigEr);
//        deltaEr = gauE->GetRandom();
//        Etrue   = Ereco + deltaEr;
        Etrue = True_array[recoBin]->GetRandom();
        //std::cout << "Ereco = " << Ereco << "   Etrue = " << Etrue << endl;
        
        histoNu->Fill(Etrue,Ereco);

        Erecob = numub110_noosc_histo->GetRandom();
        //sigEr  = sig10*sqrt(Erecob/e10);
//        sigEr  = (0.043 + 0.213*Erecob - 0.0051*pow(Erecob,2))/(2*sqrt(2*log(2)));
//        gauE->SetParameter(0,sigEr);
//        deltaEr = gauE->GetRandom();
//        Etrueb  = Erecob + deltaEr;
        recoBin = int(Erecob/0.25);
        Etrueb = True_array[recoBin]->GetRandom();
        //std::cout << "Ereco = " << Erecob << "   Etrue = " << Etrueb << endl;

        histoNub->Fill(Etrueb,Erecob);

        Tnumu ->Fill();
        Tnumub->Fill();
    }
    
    fout->cd();
    histoNu->Write();
    histoNub->Write();

    fout->Write();
    
    fenergy->Close();
    fout->Close();
} // end
