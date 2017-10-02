// csv files are digitalizations of data plots in the article, by using Engauge software.
// plot from Daya Bay Coll. PRL 112, 061801 (2014) - arXiv:13106732

void db_PRL112_plots()
{// begin

    //------------- Style --------------
    gROOT->SetStyle("Plain");
    gStyle->SetOptStat(0);
    //------------- Style --------------
    //----------  Text Style  ---------
    //ft = 10 * fontID + precision
    Int_t ft = 10 * 4 + 2;
    Double_t sz = 0.04;
    //---------------------------------

    //---------------------------------------------------------------------------------------------------
    //get data from digitized files - EH1
    TGraph *DB_specEH1_Data  = new TGraph("data_files/EH1/DB_specEH1_Data.csv","%lg,%lg","");
    TGraph *DB_specEH1_NoOsc = new TGraph("data_files/EH1/DB_specEH1_NoOsc.csv","%lg,%lg","");
    TGraph *DB_specEH1_BF    = new TGraph("data_files/EH1/DB_specEH1_BF.csv","%lg,%lg","");
    TGraph *DB_specEH1_BG    = new TGraph("data_files/EH1/DB_specEH1_BG.csv","%lg,%lg","");
    
    TGraph *DB_ratioEH1_BF   = new TGraph("data_files/EH1/DB_ratioEH1_BF.csv","%lg,%lg","");
    TGraph *DB_ratioEH1_Data = new TGraph("data_files/EH1/DB_ratioEH1_Data.csv","%lg,%lg","");
    TGraph *DB_ratioEH1_Erro = new TGraph("data_files/EH1/DB_ratioEH1_Erro.csv","%lg,%lg","");
    
    //Note: "DB_ratioEH1_NoOsc" is not digitized. Trivial horizontal line
    //---------------------------------------------------------------------------------------------------
    //get data from digitized files - EH2
    TGraph *DB_specEH2_Data  = new TGraph("data_files/EH2/DB_specEH2_Data.csv","%lg,%lg","");
    TGraph *DB_specEH2_NoOsc = new TGraph("data_files/EH2/DB_specEH2_NoOsc.csv","%lg,%lg","");
    TGraph *DB_specEH2_BF    = new TGraph("data_files/EH2/DB_specEH2_BF.csv","%lg,%lg","");
    TGraph *DB_specEH2_BG    = new TGraph("data_files/EH2/DB_specEH2_BG.csv","%lg,%lg","");
    
    TGraph *DB_ratioEH2_BF   = new TGraph("data_files/EH2/DB_ratioEH2_BF.csv","%lg,%lg","");
    TGraph *DB_ratioEH2_Data = new TGraph("data_files/EH2/DB_ratioEH2_Data.csv","%lg,%lg","");
    TGraph *DB_ratioEH2_Erro = new TGraph("data_files/EH2/DB_ratioEH2_Erro.csv","%lg,%lg","");
    
    //Note: "DB_ratioEH1_NoOsc" is not digitized. Trivial horizontal line
    //---------------------------------------------------------------------------------------------------
    //get data from digitized files - EH3
    TGraph *DB_specEH3_Data  = new TGraph("data_files/EH3/DB_specEH3_Data.csv","%lg,%lg","");
    TGraph *DB_specEH3_NoOsc = new TGraph("data_files/EH3/DB_specEH3_NoOsc.csv","%lg,%lg","");
    TGraph *DB_specEH3_BF    = new TGraph("data_files/EH3/DB_specEH3_BF.csv","%lg,%lg","");
    TGraph *DB_specEH3_BG    = new TGraph("data_files/EH3/DB_specEH3_BG.csv","%lg,%lg","");
    
    TGraph *DB_ratioEH3_BF   = new TGraph("data_files/EH3/DB_ratioEH3_BF.csv","%lg,%lg","");
    TGraph *DB_ratioEH3_Data = new TGraph("data_files/EH3/DB_ratioEH3_Data.csv","%lg,%lg","");
    TGraph *DB_ratioEH3_Erro = new TGraph("data_files/EH3/DB_ratioEH3_Erro.csv","%lg,%lg","");
    
    //Note: "DB_ratioEH1_NoOsc" is not digitized. Trivial horizontal line
    //---------------------------------------------------------------------------------------------------

    //define histograms
    double    NB = 26;
    double    lo = 0.7;
    double    hi = 12.0;

    double xbins[27];
    xbins[0] = 0.7;
    double delta_bins2 = (7.3 - 1.3)/24; // 0.25 MeV/bin
    for (int i = 0 ; i < (NB-1) ; i++)
        {
            xbins[i+1] = 1.3 + delta_bins2*i;
        }
    xbins[26] = hi;

    //-----------
    const int nEH = 3;
    //-- Events / MeV
    TH1F *data_spect_histoPerMeV[nEH];
    TH1F *nosc_spect_histoPerMeV[nEH];
    TH1F *BFit_spect_histoPerMeV[nEH];
    TH1F *BkGd_spect_histoPerMeV[nEH];
    //-- Events
    TH1F *data_spect_histo[nEH];
    TH1F *nosc_spect_histo[nEH];
    TH1F *BFit_spect_histo[nEH];
    TH1F *BkGd_spect_histo[nEH];
    //-- Ratio
    TH1F *data_ratio_histo[nEH];
    TH1F *nosc_ratio_histo[nEH];
    TH1F *BFit_ratio_histo[nEH];
    
    for (int i = 0 ; i < nEH ; i++)
    {
        // spectra (Events / Mev)
        //-----------
        data_spect_histoPerMeV[i] = new TH1F(Form("data_spect_histoPerMeV_%d",i),"",NB,xbins);
        data_spect_histoPerMeV[i]->SetLineWidth(2);
        data_spect_histoPerMeV[i]->SetMarkerStyle(34);
        data_spect_histoPerMeV[i]->SetMarkerSize(1.1);
        
        nosc_spect_histoPerMeV[i] = new TH1F(Form("nosc_spect_histoPerMeV_%d",i),"",NB,xbins);
        nosc_spect_histoPerMeV[i]->SetLineWidth(1.5);
        nosc_spect_histoPerMeV[i]->SetLineColor(4);
        
        BFit_spect_histoPerMeV[i] = new TH1F(Form("BFit_spect_histoPerMeV_%d",i),"",NB,xbins);
        BFit_spect_histoPerMeV[i]->SetLineWidth(3);
        BFit_spect_histoPerMeV[i]->SetLineColor(2);
        
        BkGd_spect_histoPerMeV[i] = new TH1F(Form("BkGd_spect_histoPerMeV_%d",i),"",NB,xbins);
        BkGd_spect_histoPerMeV[i]->SetLineWidth(2);
        BkGd_spect_histoPerMeV[i]->SetLineColor(6);

        // spectra (Events)
        //-----------
        data_spect_histo[i] = new TH1F(Form("data_spect_histo_%d",i),"",NB,xbins);
        data_spect_histo[i]->SetLineWidth(2);
        data_spect_histo[i]->SetMarkerStyle(34);
        data_spect_histo[i]->SetMarkerSize(1.1);
        
        nosc_spect_histo[i] = new TH1F(Form("nosc_spect_histo_%d",i),"",NB,xbins);
        nosc_spect_histo[i]->SetLineWidth(1.5);
        nosc_spect_histo[i]->SetLineColor(4);
        
        BFit_spect_histo[i] = new TH1F(Form("BFit_spect_histo_%d",i),"",NB,xbins);
        BFit_spect_histo[i]->SetLineWidth(3);
        BFit_spect_histo[i]->SetLineColor(2);
        
        BkGd_spect_histo[i] = new TH1F(Form("BkGd_spect_histo_%d",i),"",NB,xbins);
        BkGd_spect_histo[i]->SetLineWidth(2);
        BkGd_spect_histo[i]->SetLineColor(6);
        
        //Data - Background ratios
        //------------------
        data_ratio_histo[i] = new TH1F(Form("data_ratio_histo_%d",i),"",NB,xbins);
        data_ratio_histo[i]->SetLineWidth(2);
        data_ratio_histo[i]->SetMarkerStyle(8);
        data_ratio_histo[i]->SetMarkerSize(0.8);
        
        nosc_ratio_histo[i] = new TH1F(Form("nosc_ratio_histo_%d",i),"",NB,xbins);
        nosc_ratio_histo[i]->SetLineWidth(2);
        nosc_ratio_histo[i]->SetLineColor(4);
        
        BFit_ratio_histo[i] = new TH1F(Form("BFit_ratio_histo_%d",i),"",NB,xbins);
        BFit_ratio_histo[i]->SetLineWidth(2);
        BFit_ratio_histo[i]->SetLineColor(2);
    }
    
    // set histogram contents
    double ctnt = 0;
    double erro = 0;
    
    double binW = 0.0;

    for (int i = 0 ; i < NB ; i++)
        {
            // spectra EH1
            ctnt = DB_specEH1_Data->GetY()[i];
            data_spect_histoPerMeV[0]->SetBinContent(i+1,ctnt);
            binW = data_spect_histoPerMeV[0]->GetBinWidth(i+1);
            data_spect_histo[0]->SetBinContent(i+1,ctnt*binW);
            
            ctnt = DB_specEH1_NoOsc->GetY()[i];
            nosc_spect_histoPerMeV[0]->SetBinContent(i+1,ctnt);
            binW = nosc_spect_histoPerMeV[0]->GetBinWidth(i+1);
            nosc_spect_histo[0]->SetBinContent(i+1,ctnt*binW);
            
            ctnt = DB_specEH1_BF->GetY()[i];
            BFit_spect_histoPerMeV[0]->SetBinContent(i+1,ctnt);
            binW = BFit_spect_histoPerMeV[0]->GetBinWidth(i+1);
            BFit_spect_histo[0]->SetBinContent(i+1,ctnt*binW);
            
            ctnt = DB_specEH1_BG->GetY()[i];
            BkGd_spect_histoPerMeV[0]->SetBinContent(i+1,ctnt);
            binW = BkGd_spect_histoPerMeV[0]->GetBinWidth(i+1);
            BkGd_spect_histo[0]->SetBinContent(i+1,ctnt*binW);
            
            // spectra EH2
            ctnt = DB_specEH1_Data->GetY()[i];
            data_spect_histoPerMeV[1]->SetBinContent(i+1,ctnt);
            binW = data_spect_histoPerMeV[1]->GetBinWidth(i+1);
            data_spect_histo[1]->SetBinContent(i+1,ctnt*binW);
            
            ctnt = DB_specEH1_NoOsc->GetY()[i];
            nosc_spect_histoPerMeV[1]->SetBinContent(i+1,ctnt);
            binW = nosc_spect_histoPerMeV[1]->GetBinWidth(i+1);
            nosc_spect_histo[1]->SetBinContent(i+1,ctnt*binW);
            
            ctnt = DB_specEH1_BF->GetY()[i];
            BFit_spect_histoPerMeV[1]->SetBinContent(i+1,ctnt);
            binW = BFit_spect_histoPerMeV[1]->GetBinWidth(i+1);
            BFit_spect_histo[1]->SetBinContent(i+1,ctnt*binW);
            
            ctnt = DB_specEH1_BG->GetY()[i];
            BkGd_spect_histoPerMeV[1]->SetBinContent(i+1,ctnt);
            binW = BkGd_spect_histoPerMeV[1]->GetBinWidth(i+1);
            BkGd_spect_histo[1]->SetBinContent(i+1,ctnt*binW);
            
            // spectra EH3
            ctnt = DB_specEH1_Data->GetY()[i];
            data_spect_histoPerMeV[2]->SetBinContent(i+1,ctnt);
            binW = data_spect_histoPerMeV[2]->GetBinWidth(i+1);
            data_spect_histo[2]->SetBinContent(i+1,ctnt*binW);
            
            ctnt = DB_specEH1_NoOsc->GetY()[i];
            nosc_spect_histoPerMeV[2]->SetBinContent(i+1,ctnt);
            binW = nosc_spect_histoPerMeV[2]->GetBinWidth(i+1);
            nosc_spect_histo[2]->SetBinContent(i+1,ctnt*binW);
            
            ctnt = DB_specEH1_BF->GetY()[i];
            BFit_spect_histoPerMeV[2]->SetBinContent(i+1,ctnt);
            binW = BFit_spect_histoPerMeV[2]->GetBinWidth(i+1);
            BFit_spect_histo[2]->SetBinContent(i+1,ctnt*binW);
            
            ctnt = DB_specEH1_BG->GetY()[i];
            BkGd_spect_histoPerMeV[2]->SetBinContent(i+1,ctnt);
            binW = BkGd_spect_histoPerMeV[2]->GetBinWidth(i+1);
            BkGd_spect_histo[2]->SetBinContent(i+1,ctnt*binW);

            //------------------------------------------------------
            // ratios EH1
            ctnt = DB_ratioEH1_Data->GetY()[i];
            erro = DB_ratioEH1_Erro->GetY()[i];
            data_ratio_histo[0]->SetBinContent(i+1,ctnt);
            data_ratio_histo[0]->SetBinError(i+1,erro);
            
            ctnt = 1.0;
            nosc_ratio_histo[0]->SetBinContent(i+1,ctnt);
            
            ctnt = DB_ratioEH1_BF->GetY()[i];
            BFit_ratio_histo[0]->SetBinContent(i+1,ctnt);

            // ratios EH2
            ctnt = DB_ratioEH2_Data->GetY()[i];
            erro = DB_ratioEH2_Erro->GetY()[i];
            data_ratio_histo[1]->SetBinContent(i+1,ctnt);
            data_ratio_histo[1]->SetBinError(i+1,erro);
            
            ctnt = 1.0;
            nosc_ratio_histo[1]->SetBinContent(i+1,ctnt);
            
            ctnt = DB_ratioEH2_BF->GetY()[i];
            BFit_ratio_histo[1]->SetBinContent(i+1,ctnt);
            
            // ratios EH3
            ctnt = DB_ratioEH3_Data->GetY()[i];
            erro = DB_ratioEH3_Erro->GetY()[i];
            data_ratio_histo[2]->SetBinContent(i+1,ctnt);
            data_ratio_histo[2]->SetBinError(i+1,erro);
            
            ctnt = 1.0;
            nosc_ratio_histo[2]->SetBinContent(i+1,ctnt);
            
            ctnt = DB_ratioEH3_BF->GetY()[i];
            BFit_ratio_histo[2]->SetBinContent(i+1,ctnt);
            
/*            // tests that errors look OK
            ctnt = wtd_near_site_ratio_bfit_lo->GetY()[i];
            ns_wtd_ratio_bfit_lo_histo->SetBinContent(i+1,ctnt);
            //
            ctnt = wtd_near_site_ratio_bfit_hi->GetY()[i];
            ns_wtd_ratio_bfit_hi_histo->SetBinContent(i+1,ctnt);
*/
        }

    // Drawng section
    //-------------------
    
    TH2F *frame_spectra1 = new TH2F("frame_spectra1","",NB,lo,hi,10,0,64.5e3);
    frame_spectra1->GetXaxis()->SetTitle("Prompt Reconstructed Energy (MeV)");
    frame_spectra1->GetYaxis()->SetTitle("Events/MeV");
    
    TH2F *frame_spectra2 = new TH2F("frame_spectra2","",NB,lo,hi,10,0,29.5e3);
    frame_spectra2->GetXaxis()->SetTitle("Prompt Reconstructed Energy (MeV)");
    frame_spectra2->GetYaxis()->SetTitle("Events/MeV");
    
    TH2F *frame_spectra3 = new TH2F("frame_spectra3","",NB,lo,hi,10,0,13.55e3);
    frame_spectra3->GetXaxis()->SetTitle("Prompt Reconstructed Energy (MeV)");
    frame_spectra3->GetYaxis()->SetTitle("Events/MeV");
    
    TH2F *frame_ratios = new TH2F("frame_ratios","",NB,lo,hi,10,0.83,1.12);
    frame_ratios->GetXaxis()->SetTitle("Prompt Reconstructed Energy (MeV)");
    frame_ratios->GetYaxis()->SetTitle("Data - Background / Prediction");

    TCanvas *canv0 = new TCanvas("canv0","",700,450);
    canv0->Divide(1,2);

    canv0->cd(1);
    gPad->SetPad(0.005,0.300,0.995,0.995);
    frame_spectra1->Draw();
    BkGd_spect_histoPerMeV[0]->Draw("same hist");
    BFit_spect_histoPerMeV[0]->Draw("same hist");
    data_spect_histoPerMeV[0]->Draw("P same hist");
    nosc_spect_histoPerMeV[0]->Draw("same");
    gPad->SetTicks(1,1);
    
    canv0->cd(2);
    gPad->SetPad(0.005,0.006,0.995,0.300);
    frame_ratios->Draw();
    BFit_ratio_histo[0]->Draw("same hist");
    data_ratio_histo[0]->Draw("P same");
    nosc_ratio_histo[0]->Draw("same");
    gPad->RedrawAxis();
    gPad->SetTicks(1,1);

    canv0->Print("canv_DB-EH1.pdf");

    canv0->Clear("D");
    
    canv0->cd(1);
    gPad->SetPad(0.005,0.300,0.995,0.995);
    frame_spectra2->Draw();
    BkGd_spect_histoPerMeV[1]->Draw("same hist");
    BFit_spect_histoPerMeV[1]->Draw("same hist");
    data_spect_histoPerMeV[1]->Draw("P same hist");
    nosc_spect_histoPerMeV[1]->Draw("same");
    gPad->SetTicks(1,1);
    
    canv0->cd(2);
    gPad->SetPad(0.005,0.005,0.995,0.300);
    frame_ratios->Draw();
    BFit_ratio_histo[1]->Draw("same hist");
    data_ratio_histo[1]->Draw("P same");
    nosc_ratio_histo[1]->Draw("same");
    gPad->SetTicks(1,1);
    gPad->RedrawAxis();
    
    canv0->Print("canv_DB-EH2.pdf");
    
    canv0->Clear("D");
    
    canv0->cd(1);
    gPad->SetPad(0.005,0.300,0.995,0.995);
    frame_spectra3->Draw();
    BkGd_spect_histoPerMeV[2]->Draw("same hist");
    BFit_spect_histoPerMeV[2]->Draw("same hist");
    data_spect_histoPerMeV[2]->Draw("P same hist");
    nosc_spect_histoPerMeV[2]->Draw("same");
    gPad->SetTicks(1,1);
    
    canv0->cd(2);
    gPad->SetPad(0.005,0.005,0.995,0.300);
    frame_ratios->Draw();
    BFit_ratio_histo[2]->Draw("same hist");
    data_ratio_histo[2]->Draw("P same");
    nosc_ratio_histo[2]->Draw("same");
    gPad->SetTicks(1,1);
    gPad->RedrawAxis();
    
    canv0->Print("canv_DB-EH3.pdf");
    
    //---------------------------------------------------------
    // write to output file
    TFile *fout = new TFile("PRL112_data.root","recreate");
    fout->cd();

    for (int i = 0 ; i < nEH ; i++)
        {
            //spectra (Events / MeV)
            data_spect_histoPerMeV[i]->Write();
            nosc_spect_histoPerMeV[i]->Write();
            BFit_spect_histoPerMeV[i]->Write();
            BkGd_spect_histoPerMeV[i]->Write();
            //spectra (Events)
            data_spect_histo[i]->Write();
            nosc_spect_histo[i]->Write();
            BFit_spect_histo[i]->Write();
            BkGd_spect_histo[i]->Write();
            
            //ratio
            data_ratio_histo[i]->Write();
            nosc_ratio_histo[i]->Write();
            BFit_ratio_histo[i]->Write();
        }

    fout->Close();


}// end

