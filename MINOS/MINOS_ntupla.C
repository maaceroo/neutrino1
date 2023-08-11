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
    //TRandom Seed
    gRandom = new TRandom3(108);

    //-------------------
    // Energy Histograms
    //-------------------
    std::cout << "Calling the root file" << std::endl;
    TFile *fenergy = new TFile("MINOS_spectra_PRL108-PRL110.root","read");
    TFile *feneFit = new TFile("fitMINOS.root","read");

    // define and fill the histograms for Muon (Anti)neutrinos PRL110 ------------------------
    std::cout << "Creating and Filling histograms" << std::endl;
    //- Muon neutrinos
    TH1F *numu110_noosc_histo;
    TH1F *numu110_nooscbs_histo;
    TH1F *numu110_bkgd_histo;
    //numu110_noosc_histo = (TH1F*) fenergy->Get("numu110_noosc_histo");
    numu110_noosc_histo   = (TH1F*) feneFit->Get("numu_spect_histo_fine"); //2023.07.04 - This histogram was built with the bkgd substracted
    numu110_bkgd_histo    = (TH1F*) fenergy->Get("numu110_bkgd_histo"); //Not used now (11.02.22). Should be substracted from the noosc spectrum.
    numu110_nooscbs_histo = (TH1F*)(numu110_noosc_histo->Clone("numu110_nooscbs_histo"));
    //numu110_nooscbs_histo ->Add(numu110_bkgd_histo,-1); //NOTE: 30.06.2023 - there might be an error here!
    //- Muon antineutrinos
    TH1F *numub110_noosc_histo;
    TH1F *numub110_nooscbs_histo;
    TH1F *numub110_bkgd_histo;
    //numub110_noosc_histo   = (TH1F*) fenergy->Get("numub110_noosc_histo");
    numub110_noosc_histo   = (TH1F*) feneFit->Get("numub_spect_histo_fine"); //2023.07.04 - This histogram was built with the bkgd substracted
    numub110_bkgd_histo    = (TH1F*) fenergy->Get("numub110_bkgd_histo");
    numub110_nooscbs_histo = (TH1F*)(numub110_noosc_histo->Clone("numub110_nooscbs_histo"));
    //numub110_nooscbs_histo ->Add(numub110_bkgd_histo,-1);
    //- Muon antineutrinos WS
    TH1F *numubWS110_noosc_histo;
    TH1F *numubWS110_nooscbs_histo;
    TH1F *numubWS110_bkgd_histo;
    //numubWS110_noosc_histo   = (TH1F*) fenergy->Get("numubWS110_noosc_histo");
    numubWS110_noosc_histo   = (TH1F*) feneFit->Get("numubWS_spect_histo_fine"); //2023.07.04 - This histogram was built with the bkgd substracted
    numubWS110_bkgd_histo    = (TH1F*) fenergy->Get("numubWS110_bkgd_histo");
    numubWS110_nooscbs_histo = (TH1F*)(numubWS110_noosc_histo->Clone("numubWS110_nooscbs_histo"));
    //numubWS110_nooscbs_histo ->Add(numubWS110_bkgd_histo,-1);
    

    //Root file for the ntuple
    TString filePath = dirName;
    TFile *fout     = new TFile(filePath + "/data/minos_ntuple.root","RECREATE");
    TTree *Tnumu    = new TTree("Tnumu" ,  "MC neutrino events");
    TTree *Tnumub   = new TTree("Tnumub",  "MC neutrino events");
    TTree *TnumubWS = new TTree("TnumubWS","MC neutrino events");
    //-------------------
    // True Vs. Reco Energy Matrix
    //--------------------
    std::cout << "Calling the Matrix root file - Neutrinos" << std::endl;
    TFile *fmatrix = new TFile("./data/MINOS_Matrix_TvsR_Energy.root","read");
    TH2F *TrueReco_Matrix;
    TrueReco_Matrix = (TH2F*) fmatrix->Get("TRmat_histo");
    //dirty fix - Neutrinos
    //double tt1 = TrueReco_Matrix->GetBinContent(12,12);
    //TrueReco_Matrix->SetBinContent(12,12,tt1*0.90);
    //double tt2 = TrueReco_Matrix->GetBinContent(13,13);
    //TrueReco_Matrix->SetBinContent(13,13,tt2*1.05);
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
    //----------------------------------------
    std::cout << "Calling the Matrix root file - Anti-Neutrinos" << std::endl;
    //    TFile *fmatrix_nubar = new TFile("./data/MINOS_Matrix_TvsR_Energy_antiNuMu.root","read");
    TH2F *TrueReco_Matrix_nubar;
    TrueReco_Matrix_nubar = (TH2F*) fmatrix->Get("TRmat_b_histo");
    Nbins   = 60;
    n_histo = Nbins;
    TH1F *True_array_nubar[60];
    for (int j = 0 ; j < n_histo ; j++) {
        True_array_nubar[j] = new TH1F(Form("True_array_nubar_%d",j),"",Nbins,0,30);
        for (int i = 0 ; i < Nbins ; i++) {
            double value = TrueReco_Matrix_nubar->GetBinContent(i+1,j+1);
            True_array_nubar[j]->SetBinContent(i+1,value);
        }
    }

    std::cout << "Calling the Matrix root file - Anti-Neutrinos - WS sample" << std::endl;
    TH2F *TrueReco_Matrix_nubarWS;
    TrueReco_Matrix_nubarWS = (TH2F*) fmatrix->Get("TRmat_b_histo_Mod");
    //dirty fix - Neutrinos
    //double tt1b = TrueReco_Matrix_nubarWS->GetBinContent(3,10);
    //TrueReco_Matrix_nubarWS->SetBinContent(10,3,1.0);
    //double tt2b = TrueReco_Matrix_nubarWS->GetBinContent(3,11);
    //TrueReco_Matrix_nubarWS->SetBinContent(11,3,1.1);

    Nbins   = 60;
    n_histo = 15;
    TH1F *True_array_nubarWS[15];
    for (int j = 0 ; j < n_histo ; j++) {
      True_array_nubarWS[j] = new TH1F(Form("True_array_nubarWS_%d",j),"",Nbins,0,30);
      for (int i = 0 ; i < Nbins ; i++) {
	double value = TrueReco_Matrix_nubarWS->GetBinContent(i+1,j+1);
	True_array_nubarWS[j]->SetBinContent(i+1,value);
      }
    }
    
    //----------------------------------------
    //----------------------------------------

    float Ereco, Etrue;
    float Erecob, Etrueb;
    float ErecobWS, EtruebWS;
    float BL = 735.0; // km
    //int Nevents = 1e7;
    int Nevents = atoi(getenv("NTUPLE_EVENTS")); // This must be uncommented when using the script

    Tnumu   ->Branch("Ereco",   &Ereco,   "Ereco/F"   );
    Tnumu   ->Branch("Etrue",   &Etrue,   "Etrue/F"   );
    Tnumub  ->Branch("Erecob",  &Erecob,  "Erecob/F"  );
    Tnumub  ->Branch("Etrueb",  &Etrueb,  "Etrueb/F"  );
    TnumubWS->Branch("ErecobWS",&ErecobWS,"ErecobWS/F");
    TnumubWS->Branch("EtruebWS",&EtruebWS,"EtruebWS/F");
    
    double ehi = 14.0;
    TH2F *histoNu    = new TH2F("histoNu",    "", 56,0,ehi,56,0,ehi);
    TH2F *histoNub   = new TH2F("histoNub",   "", 28,0,ehi,28,0,ehi);
    double ehiWS = 25.0;
    TH2F *histoNubWS = new TH2F("histoNubWS", "", 28,0,ehi,28,0,ehi);

    int recoBin;
    TF1 *fFixing = new TF1("fFixing","gaus",-5,5);
    for (int i = 0 ; i < Nevents ; i++) {
      //neutrinos
      //Ereco  = numu110_noosc_histo->GetRandom();
      Ereco  = numu110_nooscbs_histo->GetRandom();
      recoBin = int(Ereco/0.25);
      Etrue = fudge1*True_array[recoBin]->GetRandom();
      //std::cout << "Ereco = " << Ereco << "   Etrue = " << Etrue << endl;
        
      histoNu->Fill(Etrue,Ereco);

      //antineutrinos
      //Erecob = numub110_noosc_histo->GetRandom();
      Erecob = numub110_nooscbs_histo->GetRandom();
      recoBin = int(Erecob/0.5);
      Etrueb = fudgeb1*True_array_nubar[recoBin]->GetRandom();
      //std::cout << "Ereco = " << Erecob << "   Etrue = " << Etrueb << endl;
      
      histoNub->Fill(Etrueb,Erecob);
      
      //antineutrinos
      ErecobWS = numubWS110_nooscbs_histo->GetRandom();
      //recoBin = int(ErecobWS/0.5);
      recoBin = int(ErecobWS/2.0);
      EtruebWS = fudgebWS1*True_array_nubarWS[recoBin]->GetRandom();
      /*      
      if(recoBin == 1) {
	fFixing->SetParameters(1,0,1.0);
	EtruebWS += fFixing->GetRandom();
      }
      
      if(recoBin == 2) {
	fFixing->SetParameters(1,0,1.0);    
	EtruebWS += fFixing->GetRandom();
      }
      
      if(recoBin == 3) {
	fFixing->SetParameters(1,0,1.0);    
	EtruebWS += fFixing->GetRandom();
      }
      */
      //std::cout << "Ereco = " << Erecob << "   Etrue = " << Etrueb << endl;
      
      histoNubWS->Fill(EtruebWS,ErecobWS);
      
      Tnumu   ->Fill(); //- Numu         from numu-dominated beam
      Tnumub  ->Fill(); //- NumuBar      from numuBar-enhanced beam
      TnumubWS->Fill(); //- NumuBar (WS) from numu-dominated beam
    }
    
    fout->cd();
    histoNu->Write();
    histoNub->Write();
    histoNubWS->Write();

    fout->Write();
    
    fenergy->Close();
    fout->Close();
    fmatrix->Close();
    //fmatrix_nubar->Close();

    cout << "Exiting ntuple execution." << endl;
} // end
