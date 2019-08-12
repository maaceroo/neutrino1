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
    TGraph *DB_specEH1_DataperMev  = new TGraph("files_data/EH1/DB_specEH1_Data.csv","%lg,%lg","");
    TGraph *DB_specEH1_NoOscperMev = new TGraph("files_data/EH1/DB_specEH1_NoOsc.csv","%lg,%lg","");
    TGraph *DB_specEH1_BFperMev    = new TGraph("files_data/EH1/DB_specEH1_BF.csv","%lg,%lg","");
    TGraph *DB_specEH1_BGperMev    = new TGraph("files_data/EH1/DB_specEH1_BG.csv","%lg,%lg","");
    
    TGraph *DB_ratioEH1_BF   = new TGraph("files_data/EH1/DB_ratioEH1_BF.csv","%lg,%lg","");
    TGraph *DB_ratioEH1_Data = new TGraph("files_data/EH1/DB_ratioEH1_Data.csv","%lg,%lg","");
    TGraph *DB_ratioEH1_Erro = new TGraph("files_data/EH1/DB_ratioEH1_Erro.csv","%lg,%lg","");
    
    //Note: "DB_ratioEH1_NoOsc" is not digitized. Trivial horizontal line
    //---------------------------------------------------------------------------------------------------
    //get data from digitized files - EH2
    TGraph *DB_specEH2_DataperMev  = new TGraph("files_data/EH2/DB_specEH2_DataNew.csv","%lg,%lg","");
    TGraph *DB_specEH2_NoOscperMev = new TGraph("files_data/EH2/DB_specEH2_NoOscNew.csv","%lg,%lg","");
    TGraph *DB_specEH2_BFperMev    = new TGraph("files_data/EH2/DB_specEH2_BF.csv","%lg,%lg","");
    TGraph *DB_specEH2_BGperMev    = new TGraph("files_data/EH2/DB_specEH2_BG.csv","%lg,%lg","");
    
    TGraph *DB_ratioEH2_BF   = new TGraph("files_data/EH2/DB_ratioEH2_BF.csv","%lg,%lg","");
    TGraph *DB_ratioEH2_Data = new TGraph("files_data/EH2/DB_ratioEH2_Data.csv","%lg,%lg","");
    TGraph *DB_ratioEH2_Erro = new TGraph("files_data/EH2/DB_ratioEH2_Erro.csv","%lg,%lg","");
    
    //Note: "DB_ratioEH1_NoOsc" is not digitized. Trivial horizontal line
    //---------------------------------------------------------------------------------------------------
    //get data from digitized files - EH3
    TGraph *DB_specEH3_DataperMev  = new TGraph("files_data/EH3/DB_specEH3_Data.csv","%lg,%lg","");
    TGraph *DB_specEH3_NoOscperMev = new TGraph("files_data/EH3/DB_specEH3_NoOsc.csv","%lg,%lg","");
    TGraph *DB_specEH3_BFperMev    = new TGraph("files_data/EH3/DB_specEH3_BF.csv","%lg,%lg","");
    TGraph *DB_specEH3_BGperMev    = new TGraph("files_data/EH3/DB_specEH3_BG.csv","%lg,%lg","");
    
    TGraph *DB_ratioEH3_BF   = new TGraph("files_data/EH3/DB_ratioEH3_BF.csv","%lg,%lg","");
    TGraph *DB_ratioEH3_Data = new TGraph("files_data/EH3/DB_ratioEH3_Data.csv","%lg,%lg","");
    TGraph *DB_ratioEH3_Erro = new TGraph("files_data/EH3/DB_ratioEH3_Erro.csv","%lg,%lg","");
    
    //Note: "DB_ratioEH1_NoOsc" is not digitized. Trivial horizontal line
    //------------------------------------------------------------------------------------------

    //define histograms
    const int       NB = 26;
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
    TH1F *bkgd_spect_histoPerMeV[nEH];
    //-- Events
    TH1F *data_spect_histo[nEH];
    TH1F *nosc_spect_histo[nEH];
    TH1F *BFit_spect_histo[nEH];
    TH1F *bkgd_spect_histo[nEH];
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
        nosc_spect_histoPerMeV[i]->SetLineWidth(1);
        nosc_spect_histoPerMeV[i]->SetLineColor(4);
        
        BFit_spect_histoPerMeV[i] = new TH1F(Form("BFit_spect_histoPerMeV_%d",i),"",NB,xbins);
        BFit_spect_histoPerMeV[i]->SetLineWidth(3);
        BFit_spect_histoPerMeV[i]->SetLineColor(2);
        
        bkgd_spect_histoPerMeV[i] = new TH1F(Form("bkgd_spect_histoPerMeV_%d",i),"",NB,xbins);
        bkgd_spect_histoPerMeV[i]->SetLineWidth(2);
        bkgd_spect_histoPerMeV[i]->SetLineColor(6);

        // spectra (Events)
        //-----------
        data_spect_histo[i] = new TH1F(Form("data_spect_histo_%d",i),"",NB,xbins);
        data_spect_histo[i]->SetLineWidth(2);
        data_spect_histo[i]->SetMarkerStyle(34);
        data_spect_histo[i]->SetMarkerSize(1.1);
        
        nosc_spect_histo[i] = new TH1F(Form("nosc_spect_histo_%d",i),"",NB,xbins);
        nosc_spect_histo[i]->SetLineWidth(1);
        nosc_spect_histo[i]->SetLineColor(4);
        
        BFit_spect_histo[i] = new TH1F(Form("BFit_spect_histo_%d",i),"",NB,xbins);
        BFit_spect_histo[i]->SetLineWidth(3);
        BFit_spect_histo[i]->SetLineColor(2);
        
        bkgd_spect_histo[i] = new TH1F(Form("bkgd_spect_histo_%d",i),"",NB,xbins);
        bkgd_spect_histo[i]->SetLineWidth(2);
        bkgd_spect_histo[i]->SetLineColor(6);
        
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
            ctnt = DB_specEH1_DataperMev->GetY()[i];
            data_spect_histoPerMeV[0]->SetBinContent(i+1,ctnt);
            binW = data_spect_histoPerMeV[0]->GetBinWidth(i+1);
            data_spect_histo[0]->SetBinContent(i+1,ctnt*binW);
            
            ctnt = DB_specEH1_NoOscperMev->GetY()[i];
            nosc_spect_histoPerMeV[0]->SetBinContent(i+1,ctnt);
            binW = nosc_spect_histoPerMeV[0]->GetBinWidth(i+1);
            nosc_spect_histo[0]->SetBinContent(i+1,ctnt*binW);
            
            ctnt = DB_specEH1_BFperMev->GetY()[i];
            BFit_spect_histoPerMeV[0]->SetBinContent(i+1,ctnt);
            binW = BFit_spect_histoPerMeV[0]->GetBinWidth(i+1);
            BFit_spect_histo[0]->SetBinContent(i+1,ctnt*binW);
            
            ctnt = DB_specEH1_BGperMev->GetY()[i];
            bkgd_spect_histoPerMeV[0]->SetBinContent(i+1,ctnt);
            binW = bkgd_spect_histoPerMeV[0]->GetBinWidth(i+1);
            bkgd_spect_histo[0]->SetBinContent(i+1,ctnt*binW);
            
            // spectra EH2
            ctnt = DB_specEH2_DataperMev->GetY()[i];
            data_spect_histoPerMeV[1]->SetBinContent(i+1,ctnt);
            binW = data_spect_histoPerMeV[1]->GetBinWidth(i+1);
            data_spect_histo[1]->SetBinContent(i+1,ctnt*binW);
            
            ctnt = DB_specEH2_NoOscperMev->GetY()[i];
            nosc_spect_histoPerMeV[1]->SetBinContent(i+1,ctnt);
            binW = nosc_spect_histoPerMeV[1]->GetBinWidth(i+1);
            nosc_spect_histo[1]->SetBinContent(i+1,ctnt*binW);
            
            ctnt = DB_specEH2_BFperMev->GetY()[i];
            BFit_spect_histoPerMeV[1]->SetBinContent(i+1,ctnt);
            binW = BFit_spect_histoPerMeV[1]->GetBinWidth(i+1);
            BFit_spect_histo[1]->SetBinContent(i+1,ctnt*binW);
            
            ctnt = DB_specEH2_BGperMev->GetY()[i];
            bkgd_spect_histoPerMeV[1]->SetBinContent(i+1,ctnt);
            binW = bkgd_spect_histoPerMeV[1]->GetBinWidth(i+1);
            bkgd_spect_histo[1]->SetBinContent(i+1,ctnt*binW);
            
            // spectra EH3
            ctnt = DB_specEH3_DataperMev->GetY()[i];
            data_spect_histoPerMeV[2]->SetBinContent(i+1,ctnt);
            binW = data_spect_histoPerMeV[2]->GetBinWidth(i+1);
            data_spect_histo[2]->SetBinContent(i+1,ctnt*binW);
            
            ctnt = DB_specEH3_NoOscperMev->GetY()[i];
            nosc_spect_histoPerMeV[2]->SetBinContent(i+1,ctnt);
            binW = nosc_spect_histoPerMeV[2]->GetBinWidth(i+1);
            nosc_spect_histo[2]->SetBinContent(i+1,ctnt*binW);
            
            ctnt = DB_specEH3_BFperMev->GetY()[i];
            BFit_spect_histoPerMeV[2]->SetBinContent(i+1,ctnt);
            binW = BFit_spect_histoPerMeV[2]->GetBinWidth(i+1);
            BFit_spect_histo[2]->SetBinContent(i+1,ctnt*binW);
            
            ctnt = DB_specEH3_BGperMev->GetY()[i];
            bkgd_spect_histoPerMeV[2]->SetBinContent(i+1,ctnt);
            binW = bkgd_spect_histoPerMeV[2]->GetBinWidth(i+1);
            bkgd_spect_histo[2]->SetBinContent(i+1,ctnt*binW);

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
    frame_spectra1->GetYaxis()->SetTitle("Events/MeV");
    frame_spectra1->GetYaxis()->SetTitleFont(ft);
    frame_spectra1->GetYaxis()->SetTitleOffset(0.7);
    frame_spectra1->GetYaxis()->SetTitleSize(1.4*sz);
    frame_spectra1->GetYaxis()->SetLabelSize(1.4*sz);
    frame_spectra1->GetYaxis()->SetLabelFont(ft);

    TH2F *frame_spectra2 = new TH2F("frame_spectra2","",NB,lo,hi,10,0,29.5e3);
    frame_spectra2->GetYaxis()->SetTitle("Events/MeV");
    frame_spectra2->GetYaxis()->SetTitleFont(ft);
    frame_spectra2->GetYaxis()->SetTitleOffset(0.7);
    frame_spectra2->GetYaxis()->SetTitleSize(1.4*sz);
    frame_spectra2->GetYaxis()->SetLabelSize(1.4*sz);
    frame_spectra2->GetYaxis()->SetLabelFont(ft);

    TH2F *frame_spectra3 = new TH2F("frame_spectra3","",NB,lo,hi,10,0,13.55e3);
    frame_spectra3->GetXaxis()->SetTitle("Prompt Reconstructed Energy (MeV)");
    frame_spectra3->GetXaxis()->SetTitleFont(ft);
    frame_spectra3->GetXaxis()->SetTitleOffset(0.9);
    frame_spectra3->GetXaxis()->SetTitleSize(1.4*sz);
    frame_spectra3->GetXaxis()->SetLabelSize(1.4*sz);
    frame_spectra3->GetXaxis()->SetLabelFont(ft);
    frame_spectra3->GetYaxis()->SetTitle("Events/MeV");
    frame_spectra3->GetYaxis()->SetTitleFont(ft);
    frame_spectra3->GetYaxis()->SetTitleOffset(0.7);
    frame_spectra3->GetYaxis()->SetTitleSize(1.4*sz);
    frame_spectra3->GetYaxis()->SetLabelSize(1.4*sz);
    frame_spectra3->GetYaxis()->SetLabelFont(ft);

    TH2F *frame_ratios = new TH2F("frame_ratios","",NB,lo,hi,10,0.83,1.12);
    frame_ratios->GetXaxis()->SetTitle("Prompt Reconstructed Energy (MeV)");
    frame_ratios->GetYaxis()->SetTitle("Data - Background / Prediction");

    TCanvas *canv0 = new TCanvas("canv0","",700,450);
    canv0->Divide(1,2);

    canv0->cd(1);
    gPad->SetPad(0.005,0.300,0.995,0.995);
    frame_spectra1->Draw();
    bkgd_spect_histoPerMeV[0]->Draw("same hist");
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

    canv0->Print("files_plots/canv_DB.pdf");

    canv0->Clear("D");
    
    canv0->cd(1);
    gPad->SetPad(0.005,0.300,0.995,0.995);
    frame_spectra2->Draw();
    bkgd_spect_histoPerMeV[1]->Draw("same hist");
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
    
    canv0->Print("files_plots/canv_DB-EH2.pdf");
    
    canv0->Clear("D");
    
    canv0->cd(1);
    gPad->SetPad(0.005,0.300,0.995,0.995);
    frame_spectra3->Draw();
    bkgd_spect_histoPerMeV[2]->Draw("same hist");
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
    
    canv0->Print("files_plots/canv_DB-EH3.pdf");
    
    //---------------------------------------------------------
    // write to output file
    TFile *fout = new TFile("PRL112_217days_data.root","recreate");
    fout->cd();

    for (int i = 0 ; i < nEH ; i++)
        {
            //spectra (Events / MeV)
            data_spect_histoPerMeV[i]->Write();
            nosc_spect_histoPerMeV[i]->Write();
            BFit_spect_histoPerMeV[i]->Write();
            bkgd_spect_histoPerMeV[i]->Write();
            //spectra (Events)
            data_spect_histo[i]->Write();
            nosc_spect_histo[i]->Write();
            BFit_spect_histo[i]->Write();
            bkgd_spect_histo[i]->Write();
            
            //ratio
            data_ratio_histo[i]->Write();
            nosc_ratio_histo[i]->Write();
            BFit_ratio_histo[i]->Write();
        }

    fout->Close();

    //-- Our Bes fit espectra - BEGIN------------------------------------------------------------//
    double bfOscAD1_array[NB] = {44968.589844, 48651.742188, 65845.304688, 80212.937500, 91929.820312, 100853.164062, 104951.203125, 106869.843750, 105027.710938, 100476.445312, 91899.554688, 84035.843750, 76399.820312, 70017.156250, 62114.343750, 54348.109375, 47021.695312, 39262.187500, 32153.701172, 25413.869141, 19849.662109, 14527.382812, 10495.794922, 7417.966797, 4629.693848, 6937.430176};
    double bfOscAD2_array[NB] = {45798.835938, 49802.015625, 66966.273438, 81310.820312, 92438.156250, 102606.867188, 106976.484375, 108740.359375, 106891.429688, 101904.625000, 94120.523438, 85332.226562, 78009.054688, 71388.085938, 63426.437500, 55402.273438, 47999.664062, 39843.132812, 32604.730469, 25825.503906, 20385.707031, 14810.106445, 10360.466797, 7593.097656, 4732.743164, 7040.079102};
    double bfOscAD3_array[NB] = {39646.667969, 44560.691406, 59413.000000, 72553.234375, 82326.945312, 91207.617188, 95359.109375, 97101.359375, 95388.695312, 90927.539062, 83537.578125, 76337.468750, 70078.000000, 63335.417969, 56824.496094, 50187.023438, 42872.617188, 35609.285156, 28885.322266, 23285.650391, 18354.080078, 13330.003906, 9481.607422, 6389.695801, 4432.358398, 4804.482422};
    double bfOscAD4_array[NB] = {5419.828613, 5621.157715, 7572.794434, 9034.805664, 10458.491211, 11191.432617, 11929.741211, 11927.525391, 11798.802734, 11094.764648, 10271.388672, 9269.329102, 8726.391602, 7939.876953, 7042.828125, 6356.278320, 5429.204590, 4457.473633, 3580.233887, 2939.905029, 2242.790527, 1692.819092, 1199.760254, 766.403015, 520.860657, 594.856934};
    double bfOscAD5_array[NB] = {5279.687988, 5689.091309, 7536.551758, 9056.147461, 10293.578125, 11267.240234, 11892.850586, 12118.606445, 11968.355469, 11414.507812, 10515.269531, 9511.712891, 8720.475586, 7901.094238, 7109.237305, 6377.724609, 5326.480957, 4569.768555, 3612.557861, 2987.107178, 2365.126221, 1716.851562, 1139.369263, 816.917358, 546.164246, 637.580933};
    double bfOscAD6_array[NB] = {5388.755859, 5746.058105, 7512.610352, 9199.747070, 10270.195312, 11032.949219, 11644.550781, 11945.448242, 11672.162109, 11147.319336, 10234.995117, 9234.738281, 8582.391602, 7745.772461, 7123.770020, 6118.626953, 5371.090332, 4378.521973, 3680.504150, 2957.345215, 2350.883789, 1613.286621, 1181.495850, 766.518250, 545.671204, 608.436157};
    
    const int nAD = 6;
    double totalBgd[nAD][2]     = { {13.20,0.98},{13.01,0.98},{ 9.57,0.71},{ 3.52,0.14},{ 3.48,0.14},{3.43,0.14} };
    double NoscTot[nAD] = {1496310.973633, 1522309.700195, 1356229.946777, 169079.744629, 170370.055481, 168053.844360};
    double noNoscTot[nAD] = {1522533.0, 1548558.0, 1383643.0, 181663.0, 183053.0, 180550.0};
    double SurvPavg[nAD] = {0.0};
    for (int iAD = 0 ; iAD < nAD ; iAD++) {
        SurvPavg[iAD] = NoscTot[iAD]/noNoscTot[iAD];
    }
    double noOsc_IBDrate_perday[nAD];
    //-- File to get noOsc normalizations
    ifstream IBDrates_file("files_data/db6AD_noOsc_IBDrates_perday.txt");
    for (int i=0; i< nAD; i ++)
    {
        IBDrates_file >> noOsc_IBDrate_perday[i];
        //cout << "noOscIBD_rates_perday " << i << ": " << noOsc_IBDrate_perday[i] << endl;
    }//for
    double emuem[nAD]           = {0.7957,0.7927,0.8282,0.9577,0.9568,0.9566};
    double daqTime[nAD]         = {191.001,191.001,189.645,189.779,189.779,189.779};
    
    double Td = 0.0;
    TH1F *bfOscAD1_histo = new TH1F("bfOscAD1_histo","",NB,xbins);
    TH1F *bfOscAD2_histo = new TH1F("bfOscAD2_histo","",NB,xbins);
    TH1F *bfOscAD3_histo = new TH1F("bfOscAD3_histo","",NB,xbins);
    bfOscAD3_histo->SetLineColor(2);
    bfOscAD3_histo->SetLineWidth(2);
    TH1F *bfOscAD4_histo = new TH1F("bfOscAD4_histo","",NB,xbins);
    TH1F *bfOscAD5_histo = new TH1F("bfOscAD5_histo","",NB,xbins);
    TH1F *bfOscAD6_histo = new TH1F("bfOscAD6_histo","",NB,xbins);
    for (int iBIN = 0 ; iBIN < NB ; iBIN++) {
        Td = bfOscAD1_array[iBIN]*(SurvPavg[0]*noOsc_IBDrate_perday[0]/NoscTot[0])*emuem[0]*daqTime[0];
        bfOscAD1_histo->SetBinContent(iBIN+1,Td);
        
        Td = bfOscAD2_array[iBIN]*(SurvPavg[1]*noOsc_IBDrate_perday[1]/NoscTot[1])*emuem[1]*daqTime[1];
        bfOscAD2_histo->SetBinContent(iBIN+1,Td);
        
        Td = bfOscAD3_array[iBIN]*(SurvPavg[2]*noOsc_IBDrate_perday[2]/NoscTot[2])*emuem[2]*daqTime[2];
        bfOscAD3_histo->SetBinContent(iBIN+1,Td);
        
        Td = bfOscAD4_array[iBIN]*(SurvPavg[3]*noOsc_IBDrate_perday[3]/NoscTot[3])*emuem[3]*daqTime[3];
        bfOscAD4_histo->SetBinContent(iBIN+1,Td);
        
        Td = bfOscAD5_array[iBIN]*(SurvPavg[4]*noOsc_IBDrate_perday[4]/NoscTot[4])*emuem[4]*daqTime[4];
        bfOscAD5_histo->SetBinContent(iBIN+1,Td);
        
        Td = bfOscAD6_array[iBIN]*(SurvPavg[5]*noOsc_IBDrate_perday[5]/NoscTot[5])*emuem[5]*daqTime[5];
        bfOscAD6_histo->SetBinContent(iBIN+1,Td);
    }
    
    TH1F *our_BFit_histo[nEH];
    for (int i = 0 ; i < nEH ; i++){
        our_BFit_histo[i] = new TH1F(Form("BFit_spect_histo_%d",i),"",NB,xbins);
        our_BFit_histo[i]->SetLineWidth(2);
        our_BFit_histo[i]->SetLineColor(2);
        
    }

    TH1F *temp_histo = new TH1F("temp_histo","",NB,xbins);;
    
    our_BFit_histo[0]->Add(bfOscAD1_histo,bfOscAD2_histo);
    our_BFit_histo[1] = (TH1F*)bfOscAD3_histo->Clone("our_BFit_histo_1");
    temp_histo->Add(bfOscAD4_histo,bfOscAD5_histo);
    our_BFit_histo[2]->Add(temp_histo,bfOscAD6_histo);
    
    for (int j = 0 ; j < nEH ; j++) {
        for (int i = 0 ; i < NB ; i++) {
            double bgnd = bkgd_spect_histoPerMeV[j]->GetBinContent(i+1);
            double wid  = our_BFit_histo[j]->GetBinWidth(i+1);
            double cont = our_BFit_histo[j]->GetBinContent(i+1);
            our_BFit_histo[j]->SetBinContent(i+1,(cont/wid) + bgnd);
        }
    }
    
    //-- Our Bes fit espectra - END--//

    /////////////////////////
    TLatex *lat = new TLatex();
    lat->SetNDC();
    lat->SetTextFont(ft);
    lat->SetTextSize(2.6*sz);

    TLegend *leg1 = new TLegend(0.6,0.5,0.8,0.8);
    leg1->SetTextFont(ft);
    leg1->SetTextSize(1.7*sz);
    leg1->SetFillColor(0);
    leg1->SetLineColor(0);
    
    leg1->AddEntry(data_spect_histo[0],"Daya Bay Data","pl");
    leg1->AddEntry(nosc_spect_histo[0],"No Oscillations","l");
    leg1->AddEntry(our_BFit_histo[0],"Best fit","l");
    leg1->AddEntry(bkgd_spect_histo[0],"Total Background","l");

    TCanvas *ca = new TCanvas("ca", "canvas", 700, 900);
    TGaxis::SetMaxDigits(3);

    TPad *pad1 = new TPad("pad1", "pad1", 0, 2./3., 1, 1.0);
    pad1->SetBottomMargin(0); // Upper and lower plot are joined
    pad1->Draw();             // Draw the upper pad: pad1
    pad1->cd();               // pad1 becomes the current pad
    frame_spectra1->Draw();
    //BFit_spect_histoPerMeV[0]->Draw("same hist");
    data_spect_histoPerMeV[0]->Draw("P same");
    our_BFit_histo[0]->Draw("same hist");
    nosc_spect_histoPerMeV[0]->Draw("same");
    bkgd_spect_histoPerMeV[0]->Draw("same");
    gPad->SetTicks(1,1);
    leg1->Draw();
    lat->DrawLatex(0.7,0.3,"EH1");
    
    // lower plot will be in pad
    ca->cd();          // Go back to the main canvas before defining pad2
    TPad *pad2 = new TPad("pad2", "pad2", 0, 1./3., 1, 2./3.);
    pad2->SetTopMargin(0);
    pad2->SetBottomMargin(0);
    pad2->Draw();
    pad2->cd();       // pad2 becomes the current pad
    frame_spectra2->Draw();
    //BFit_spect_histoPerMeV[1]->Draw("same hist");
    data_spect_histoPerMeV[1]->Draw("P same");
    //our_BFit_histo[1]->Draw("same hist");
    nosc_spect_histoPerMeV[1]->Draw("same");
    bkgd_spect_histoPerMeV[1]->Draw("same");
    gPad->SetTicks(1,1);
    lat->DrawLatex(0.7,0.3,"EH2");
    
    // lower plot will be in pad
    ca->cd();          // Go back to the main canvas before defining pad2
    TPad *pad3 = new TPad("pad3", "pad3", 0, 0.0, 1, 1./3.);
    pad3->SetTopMargin(0);
    pad3->SetBottomMargin(0.11);
    pad3->Draw();
    pad3->cd();       // pad2 becomes the current pad
    frame_spectra3->Draw();
    //BFit_spect_histoPerMeV[2]->Draw("same hist");
    data_spect_histoPerMeV[2]->Draw("P same");
    our_BFit_histo[2]->Draw("same hist");
    nosc_spect_histoPerMeV[2]->Draw("same");
    bkgd_spect_histoPerMeV[2]->Draw("same");
    gPad->SetTicks(1,1);
    lat->DrawLatex(0.7,0.3,"EH3");
    
    ca->Print("files_plots/DB_spect-PRL112.pdf");
    ca->Print("files_plots/DB_spect-PRL112.eps");

    TH2F *frame_eventsEH1 = new TH2F("frame_eventsEH1","",NB,lo,hi,10,0,64000*0.25);
    frame_eventsEH1->GetXaxis()->SetTitle("Energ#'ia Reconstruida (MeV)");
    frame_eventsEH1->GetYaxis()->SetTitle("Eventos");
    
    TH2F *frame_eventsEH2 = new TH2F("frame_eventsEH2","",NB,lo,hi,10,0,29500*0.25);
    frame_eventsEH2->GetXaxis()->SetTitle("Energ#'ia Reconstruida (MeV)");
    frame_eventsEH2->GetYaxis()->SetTitle("Eventos");
    
    TH2F *frame_eventsEH3 = new TH2F("frame_eventsEH3","",NB,lo,hi,10,0,13750*0.25);
    frame_eventsEH3->GetXaxis()->SetTitle("Energ#'ia Reconstruida (MeV)");
    frame_eventsEH3->GetYaxis()->SetTitle("Eventos");

    TCanvas *ce = new TCanvas("ce", "canvas", 700, 900);
    
    TPad *pad4 = new TPad("pad4", "pad4", 0, 2./3., 1, 1.0);
    pad4->SetBottomMargin(0); // Upper and lower plot are joined
    pad4->Draw();             // Draw the upper pad: pad1
    pad4->cd();               // pad1 becomes the current pad
    frame_eventsEH1->Draw();
    BFit_spect_histo[0]->Draw("same hist");
    data_spect_histo[0]->Draw("P same hist");
    nosc_spect_histo[0]->Draw("same");
    bkgd_spect_histo[0]->Draw("same");
    gPad->SetTicks(1,1);
    leg1->Draw();
    lat->DrawLatex(0.7,0.3,"EH1");
    
    // lower plot will be in pad
    ce->cd();          // Go back to the main canvas before defining pad2
    TPad *pad5 = new TPad("pad5", "pad5", 0, 1./3., 1, 2./3.);
    pad5->SetTopMargin(0);
    pad5->SetBottomMargin(0);
    pad5->Draw();
    pad5->cd();       // pad2 becomes the current pad
    frame_eventsEH2->Draw();
    BFit_spect_histo[1]->Draw("same hist");
    data_spect_histo[1]->Draw("P same hist");
    nosc_spect_histo[1]->Draw("same");
    bkgd_spect_histo[1]->Draw("same");
    gPad->SetTicks(1,1);
    lat->DrawLatex(0.7,0.3,"EH2");
    
    // lower plot will be in pad
    ce->cd();          // Go back to the main canvas before defining pad2
    TPad *pad6 = new TPad("pad6", "pad6", 0, 0.0, 1, 1./3.);
    pad6->SetTopMargin(0);
    pad6->SetBottomMargin(0.1);
    pad6->Draw();
    pad6->cd();       // pad2 becomes the current pad
    frame_eventsEH3->Draw();
    BFit_spect_histo[2]->Draw("same hist");
    data_spect_histo[2]->Draw("P same hist");
    nosc_spect_histo[2]->Draw("same");
    bkgd_spect_histo[2]->Draw("same");
    gPad->SetTicks(1,1);
    lat->DrawLatex(0.7,0.3,"EH3");
    
    ce->Print("files_plots/DB_Events-PRL112.pdf");

}// end

