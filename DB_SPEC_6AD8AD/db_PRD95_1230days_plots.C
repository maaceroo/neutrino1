// plot from Daya Bay Coll. PRD 95, 072006 (2017) - arXiv:1610.04802

void db_PRD95_1230days_plots()
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
    TGraph *DB_specEH1_Data  = new TGraph("files_data/EH1/DB1230_specEH1_Data.txt","%lg %lg","");
    TGraph *DB_specEH1_NoOsc = new TGraph("files_data/EH1/DB1230_specEH1_NoOsc.txt","%lg %lg","");
    TGraph *DB_specEH1_BF    = new TGraph("files_data/EH1/DB1230_specEH1_BF.txt","%lg %lg","");
    TGraph *DB_specEH1_BG    = new TGraph("files_data/EH1/DB1230_specEH1_BG.txt","%lg %lg","");
    
    //---------------------------------------------------------------------------------------------------
    //get data from digitized files - EH2
    TGraph *DB_specEH2_Data  = new TGraph("files_data/EH2/DB1230_specEH2_Data.txt","%lg %lg","");
    TGraph *DB_specEH2_NoOsc = new TGraph("files_data/EH2/DB1230_specEH2_NoOsc.txt","%lg %lg","");
    TGraph *DB_specEH2_BF    = new TGraph("files_data/EH2/DB1230_specEH2_BF.txt","%lg %lg","");
    TGraph *DB_specEH2_BG    = new TGraph("files_data/EH2/DB1230_specEH2_BG.txt","%lg %lg","");
    
    //---------------------------------------------------------------------------------------------------
    //get data from digitized files - EH3
    TGraph *DB_specEH3_Data  = new TGraph("files_data/EH3/DB1230_specEH3_Data.txt","%lg %lg","");
    TGraph *DB_specEH3_NoOsc = new TGraph("files_data/EH3/DB1230_specEH3_NoOsc.txt","%lg %lg","");
    TGraph *DB_specEH3_BF    = new TGraph("files_data/EH3/DB1230_specEH3_BF.txt","%lg %lg","");
    TGraph *DB_specEH3_BG    = new TGraph("files_data/EH3/DB1230_specEH3_BG.txt","%lg %lg","");
    
    //---------------------------------------------------------------------------------------------------

    cout << "1. Load files... Done" << endl;
    
    //define histograms
    double    NB = 35;
    double    lo = 0.7;
    double    hi = 12.0;
    
    double xbins[36];
    xbins[0] = 0.7;
    double delta_bins2 = (7.9 - 1.3)/33; // 0.20 MeV/bin
    for (int i = 0 ; i < (NB-1) ; i++)
    {
        xbins[i+1] = 1.3 + delta_bins2*i;
    }
    xbins[35] = hi;

    cout << "2. Bins definition... Done" << endl;

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
    }
    
    cout << "3. Histogram definition... Done" << endl;

    // set histogram contents
    double ctnt = 0;
    double binW = 0.0;

    for (int i = 0 ; i < NB ; i++)
        {
            // spectra EH1
            ctnt = DB_specEH1_Data->GetY()[i];
            data_spect_histo[0]->SetBinContent(i+1,ctnt);
            binW = data_spect_histo[0]->GetBinWidth(i+1);
            data_spect_histoPerMeV[0]->SetBinContent(i+1,ctnt/(1.0e5*binW));
            
            ctnt = DB_specEH1_NoOsc->GetY()[i];
            nosc_spect_histo[0]->SetBinContent(i+1,ctnt);
            binW = nosc_spect_histo[0]->GetBinWidth(i+1);
            nosc_spect_histoPerMeV[0]->SetBinContent(i+1,ctnt/(1.0e5*binW));
            
            ctnt = DB_specEH1_BF->GetY()[i];
            BFit_spect_histo[0]->SetBinContent(i+1,ctnt);
            binW = BFit_spect_histo[0]->GetBinWidth(i+1);
            BFit_spect_histoPerMeV[0]->SetBinContent(i+1,ctnt/(1.0e5*binW));
            
            ctnt = DB_specEH1_BG->GetY()[i];
            bkgd_spect_histo[0]->SetBinContent(i+1,ctnt);
            binW = bkgd_spect_histo[0]->GetBinWidth(i+1);
            bkgd_spect_histoPerMeV[0]->SetBinContent(i+1,ctnt/(1.0e5*binW));
            
            // spectra EH2
            ctnt = DB_specEH2_Data->GetY()[i];
            data_spect_histo[1]->SetBinContent(i+1,ctnt);
            binW = data_spect_histo[1]->GetBinWidth(i+1);
            data_spect_histoPerMeV[1]->SetBinContent(i+1,ctnt/(1.0e5*binW));
            
            ctnt = DB_specEH2_NoOsc->GetY()[i];
            nosc_spect_histo[1]->SetBinContent(i+1,ctnt);
            binW = nosc_spect_histo[1]->GetBinWidth(i+1);
            nosc_spect_histoPerMeV[1]->SetBinContent(i+1,ctnt/(1.0e5*binW));
            
            ctnt = DB_specEH2_BF->GetY()[i];
            BFit_spect_histo[1]->SetBinContent(i+1,ctnt);
            binW = BFit_spect_histo[1]->GetBinWidth(i+1);
            BFit_spect_histoPerMeV[1]->SetBinContent(i+1,ctnt/(1.0e5*binW));
            
            ctnt = DB_specEH2_BG->GetY()[i];
            bkgd_spect_histo[1]->SetBinContent(i+1,ctnt);
            binW = bkgd_spect_histo[1]->GetBinWidth(i+1);
            bkgd_spect_histoPerMeV[1]->SetBinContent(i+1,ctnt/(1.0e5*binW));
            
            // spectra EH3
            ctnt = DB_specEH3_Data->GetY()[i];
            data_spect_histo[2]->SetBinContent(i+1,ctnt);
            binW = data_spect_histo[2]->GetBinWidth(i+1);
            data_spect_histoPerMeV[2]->SetBinContent(i+1,ctnt/(1.0e5*binW));
            
            ctnt = DB_specEH3_NoOsc->GetY()[i];
            nosc_spect_histo[2]->SetBinContent(i+1,ctnt);
            binW = nosc_spect_histo[2]->GetBinWidth(i+1);
            nosc_spect_histoPerMeV[2]->SetBinContent(i+1,ctnt/(1.0e5*binW));
            
            ctnt = DB_specEH3_BF->GetY()[i];
            BFit_spect_histo[2]->SetBinContent(i+1,ctnt);
            binW = BFit_spect_histo[2]->GetBinWidth(i+1);
            BFit_spect_histoPerMeV[2]->SetBinContent(i+1,ctnt/(1.0e5*binW));
            
            ctnt = DB_specEH3_BG->GetY()[i];
            bkgd_spect_histo[2]->SetBinContent(i+1,ctnt);
            binW = bkgd_spect_histo[2]->GetBinWidth(i+1);
            bkgd_spect_histoPerMeV[2]->SetBinContent(i+1,ctnt/(1.0e5*binW));
        }

    cout << "4. Fill Histograms... Done" << endl;

    for (int i = 0 ; i < nEH ; i++)
    {
        //Data/NoOsc ratios
        //------------------
        data_ratio_histo[i] = (TH1F*)data_spect_histo[i]->Clone(Form("data_ratio_histo_%d",i));
        //data_ratio_histo[i] = new TH1F(Form("data_ratio_histo_%d",i),"",NB,xbins);
        data_ratio_histo[i]->Divide(nosc_spect_histo[i]);
        data_ratio_histo[i]->SetLineWidth(2);
        data_ratio_histo[i]->SetMarkerStyle(8);
        data_ratio_histo[i]->SetMarkerSize(0.8);
        
        nosc_ratio_histo[i] = (TH1F*)nosc_spect_histo[i]->Clone(Form("nosc_ratio_histo_%d",i));
        //nosc_ratio_histo[i] = new TH1F(Form("nosc_ratio_histo_%d",i),"",NB,xbins);
        nosc_ratio_histo[i]->Divide(nosc_spect_histo[i]);
        nosc_ratio_histo[i]->SetLineWidth(2);
        nosc_ratio_histo[i]->SetLineColor(4);
        
        BFit_ratio_histo[i] = (TH1F*)BFit_spect_histo[i]->Clone(Form("BFit_ratio_histo_%d",i));
        //BFit_ratio_histo[i] = new TH1F(Form("BFit_ratio_histo_%d",i),"",NB,xbins);
        BFit_ratio_histo[i]->Divide(nosc_spect_histo[i]);
        BFit_ratio_histo[i]->SetLineWidth(2);
        BFit_ratio_histo[i]->SetLineColor(2);
    }

    cout << "5. Ratio Histograms... Done" << endl;

    // Drawng section
    //-------------------
    
    TH2F *frame_spectra1 = new TH2F("frame_spectra1","",NB,lo,hi,10,0,3.5);
    frame_spectra1->GetXaxis()->SetTitle("Prompt Energy (MeV)");
    frame_spectra1->GetYaxis()->SetTitle("Events/(MeV #times 10^{5})");
    
    TH2F *frame_spectra2 = new TH2F("frame_spectra2","",NB,lo,hi,10,0,3.15);
    frame_spectra2->GetXaxis()->SetTitle("Prompt Energy (MeV)");
    frame_spectra2->GetYaxis()->SetTitle("Events/(MeV #times 10^{5})");
    
    TH2F *frame_spectra3 = new TH2F("frame_spectra3","",NB,lo,hi,10,0,1.0);
    frame_spectra3->GetXaxis()->SetTitle("Prompt Energy (MeV)");
    frame_spectra3->GetYaxis()->SetTitle("Events/(MeV #times 10^{5})");
    
    TH2F *frame_ratios12 = new TH2F("frame_ratios12","",NB,lo,hi,10,0.93,1.05);
    frame_ratios12->GetXaxis()->SetTitle("Prompt Energy (MeV)");
    frame_ratios12->GetYaxis()->SetTitle("R^{obs} / R^{pred}_{no-osc}");
    
    TH2F *frame_ratios3 = new TH2F("frame_ratios3","",NB,lo,hi,10,0.90,1.02);
    frame_ratios3->GetXaxis()->SetTitle("Prompt Energy (MeV)");
    frame_ratios3->GetYaxis()->SetTitle("R^{obs} / R^{pred}_{no-osc}");
    
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
    frame_ratios12->Draw();
    BFit_ratio_histo[0]->Draw("same");
    data_ratio_histo[0]->Draw("P same");
    nosc_ratio_histo[0]->Draw("same");
    gPad->RedrawAxis();
    gPad->SetTicks(1,1);

    canv0->Print("files_plots/canv_DB1230-EH1.pdf");

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
    frame_ratios12->Draw();
    BFit_ratio_histo[1]->Draw("same");
    data_ratio_histo[1]->Draw("P same");
    nosc_ratio_histo[1]->Draw("same");
    gPad->SetTicks(1,1);
    gPad->RedrawAxis();
    
    canv0->Print("files_plots/canv_DB1230-EH2.pdf");
    
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
    frame_ratios3->Draw();
    BFit_ratio_histo[2]->Draw("same");
    data_ratio_histo[2]->Draw("P same");
    nosc_ratio_histo[2]->Draw("same");
    gPad->SetTicks(1,1);
    gPad->RedrawAxis();
    
    canv0->Print("files_plots/canv_DB1230-EH3.pdf");
    
    //---------------------------------------------------------
    // write to output file
    TFile *fout = new TFile("PRD95_1230days_data.root","recreate");
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

/*
    /////////////////////////
    TLatex *lat = new TLatex();
    lat->SetNDC();
    lat->SetTextFont(ft);
    lat->SetTextSize(1.6*sz);

    TLegend *leg1 = new TLegend(0.6,0.5,0.8,0.8);
    leg1->SetTextFont(ft);
    leg1->SetTextSize(1.1*sz);
    leg1->SetFillColor(0);
    leg1->SetLineColor(0);
    
    leg1->AddEntry(data_spect_histo[0],"Datos Daya Bay","pl");
    leg1->AddEntry(nosc_spect_histo[0],"Sin Oscilaci#acute{o}n","l");
    leg1->AddEntry(BFit_spect_histo[0],"Mejor Ajuste","l");
    leg1->AddEntry(bkgd_spect_histo[0],"Background Total","f");

    TCanvas *ca = new TCanvas("ca", "canvas", 700, 900);
    
    TPad *pad1 = new TPad("pad1", "pad1", 0, 2./3., 1, 1.0);
    pad1->SetBottomMargin(0); // Upper and lower plot are joined
    pad1->Draw();             // Draw the upper pad: pad1
    pad1->cd();               // pad1 becomes the current pad
    frame_spectra1->Draw();
    BFit_spect_histoPerMeV[0]->Draw("same hist");
    data_spect_histoPerMeV[0]->Draw("P same hist");
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
    BFit_spect_histoPerMeV[1]->Draw("same hist");
    data_spect_histoPerMeV[1]->Draw("P same hist");
    nosc_spect_histoPerMeV[1]->Draw("same");
    bkgd_spect_histoPerMeV[1]->Draw("same");
    gPad->SetTicks(1,1);
    lat->DrawLatex(0.7,0.3,"EH2");
    
    // lower plot will be in pad
    ca->cd();          // Go back to the main canvas before defining pad2
    TPad *pad3 = new TPad("pad3", "pad3", 0, 0.0., 1, 1./3.);
    pad3->SetTopMargin(0);
    pad3->SetBottomMargin(0.1);
    pad3->Draw();
    pad3->cd();       // pad2 becomes the current pad
    frame_spectra3->Draw();
    BFit_spect_histoPerMeV[2]->Draw("same hist");
    data_spect_histoPerMeV[2]->Draw("P same hist");
    nosc_spect_histoPerMeV[2]->Draw("same");
    bkgd_spect_histoPerMeV[2]->Draw("same");
    gPad->SetTicks(1,1);
    lat->DrawLatex(0.7,0.3,"EH3");
    
    ca->Print("files_plots/DB_spect-PRL112.pdf");

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
    
    TPad *pad1 = new TPad("pad1", "pad1", 0, 2./3., 1, 1.0);
    pad1->SetBottomMargin(0); // Upper and lower plot are joined
    pad1->Draw();             // Draw the upper pad: pad1
    pad1->cd();               // pad1 becomes the current pad
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
    TPad *pad2 = new TPad("pad2", "pad2", 0, 1./3., 1, 2./3.);
    pad2->SetTopMargin(0);
    pad2->SetBottomMargin(0);
    pad2->Draw();
    pad2->cd();       // pad2 becomes the current pad
    frame_eventsEH2->Draw();
    BFit_spect_histo[1]->Draw("same hist");
    data_spect_histo[1]->Draw("P same hist");
    nosc_spect_histo[1]->Draw("same");
    bkgd_spect_histo[1]->Draw("same");
    gPad->SetTicks(1,1);
    lat->DrawLatex(0.7,0.3,"EH2");
    
    // lower plot will be in pad
    ce->cd();          // Go back to the main canvas before defining pad2
    TPad *pad3 = new TPad("pad3", "pad3", 0, 0.0., 1, 1./3.);
    pad3->SetTopMargin(0);
    pad3->SetBottomMargin(0.1);
    pad3->Draw();
    pad3->cd();       // pad2 becomes the current pad
    frame_eventsEH3->Draw();
    BFit_spect_histo[2]->Draw("same hist");
    data_spect_histo[2]->Draw("P same hist");
    nosc_spect_histo[2]->Draw("same");
    bkgd_spect_histo[2]->Draw("same");
    gPad->SetTicks(1,1);
    lat->DrawLatex(0.7,0.3,"EH3");
    
    ce->Print("files_plots/DB_Events-PRL112.pdf");
*/
}// end

