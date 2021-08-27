// csv files are digitalizations of data plots in the articles, by using Engauge software.
// Spectra from
// - MINOS Coll. PRL108, 191801 (2012) - anti_numu        CC events
// - MINOS Coll. PRL110, 171801 (2013) - nue and anti_nue CC events
// This macro can be executed using in root [0] .x MINOS_spectra.C
// The output from code is a Root file "MINOS_plot.root"
void MINOS_spectra(){

    //------------- Style --------------
    gROOT->SetStyle("Plain");
    gStyle->SetOptStat(0);
    //------------- Style --------------
    //----------  Text Style  ---------
    //ft = 10 * fontID + precision
    Int_t ft = 10 * 4 + 2;
    Double_t sz = 0.04;
    //---------------------------------

	// READING FILES - Muon Antineutrino PRL108 ------------------
	//Muon Antineutrino Data (PRL108) and No Osc. Prediction.
	TGraph *minos_numubar_data     = new TGraph("data/PRL108_numuBar_FD_Data.csv","%lg,%lg","");
	TGraph *minos_numubar_noosc    = new TGraph("data/PRL108_numuBar_FD_NoOscPred.csv","%lg,%lg","");
	//Muon Antineutrino Background (PRL108)
	TGraph *minos_numubar_bkgd     = new TGraph("data/PRL108_numuBar_FD_Bkgd.csv","%lg,%lg","");
    //Muon Antineutrino Best Fint (PRL108)
    TGraph *minos_numubar_bfit     = new TGraph("data/PRL108_numuBar_FD_BF.csv","%lg,%lg","");

    // READING FILES - Electron (Anti)neutrino PRL110 ------------------
    std::string nue_name[3]  = {"TopLeft" ,"MidLeft" ,"LowLeft"};
    std::string nueb_name[3] = {"TopRight","MidRight","LowRight"};
    //Electron (Anti)Neutrino Data (PRL110).
    TGraph *minos_nue_data[3];
    TGraph *minos_nueb_data[3];
    //Electron (Anti)Neutrino Background (PRL110)
    TGraph *minos_nue_bkgd[3];
    TGraph *minos_nueb_bkgd[3];
    //Electron (Anti)Neutrino nue CC (PRL110)
    TGraph *minos_nue_nueCC[3];
    TGraph *minos_nueb_nueCC[3];
    //Electron (Anti)Neutrino nuebar CC (PRL110)
    TGraph *minos_nue_nuebCC[3];
    TGraph *minos_nueb_nuebCC[3];
    for (int i = 0 ; i < 3 ; i++) {
        //-- Electron neutrinos ---------------------------------------------------------
        minos_nue_data[i]   = new TGraph(("data/PRL110_nue_FD_Data_"+nue_name[i]+".csv").c_str(),"%lg,%lg","");
        minos_nue_bkgd[i]   = new TGraph(("data/PRL110_nue_FD_Bkgd_"+nue_name[i]+".csv").c_str(),"%lg,%lg","");
        minos_nue_nueCC[i]  = new TGraph(("data/PRL110_nue_FD_nueCCSignal_"+nue_name[i]+".csv").c_str(),"%lg,%lg","");
        minos_nue_nuebCC[i] = new TGraph(("data/PRL110_nue_FD_antinueCCSignal_"+nue_name[i]+".csv").c_str(),"%lg,%lg","");
        //-- Electron antineutrinos ------------------------------------------------------
        minos_nueb_data[i]   = new TGraph(("data/PRL110_nue_FD_Data_"+nueb_name[i]+".csv").c_str(),"%lg,%lg","");
        minos_nueb_bkgd[i]   = new TGraph(("data/PRL110_nue_FD_Bkgd_"+nueb_name[i]+".csv").c_str(),"%lg,%lg","");
        minos_nueb_nueCC[i]  = new TGraph(("data/PRL110_nue_FD_nueCCSignal_"+nueb_name[i]+".csv").c_str(),"%lg,%lg","");
        minos_nueb_nuebCC[i] = new TGraph(("data/PRL110_nue_FD_antinueCCSignal_"+nueb_name[i]+".csv").c_str(),"%lg,%lg","");
    }

    // define binning for Muon Antineutrinos PRL108
    const int  NB_numubar = 14;
    const double lo = 0.0;
    const double hi = 50.0;
    double xbins_numubar[NB_numubar+1];
    double delta_bins = (10.0 - 0.0)/10.0; // 1.0 GeV/bin
    
    for (int i = 0 ; i < (NB_numubar-3) ; i++)
        xbins_numubar[i] = lo + delta_bins*i;
    
    xbins_numubar[11] = xbins_numubar[10] + 5.0;
    xbins_numubar[12] = xbins_numubar[11] + 5.0;
    xbins_numubar[13] = xbins_numubar[12] + 10.0;
    xbins_numubar[14] = hi;

    // define binning for Electron (Anti)neutrinos PRL110
    const int  NB_nue = 7;
    const double loe = 1.0;
    const double hie = 8.0;
    double xbins_nue[NB_nue+1];
    double delta_binse = (hie - loe)/NB_nue; // 1.0 GeV/bin
    
    for (int i = 0 ; i <= (NB_nue) ; i++)
        xbins_nue[i] = loe + delta_binse*i;
    //------------------------------------------------------------------------------------------//

    // define and fill the histograms for Muon Antineutrinos PRL108 ------------------------
    TH1F *numubar_data_histo  = new TH1F("numubar_data_histo" ,"",NB_numubar,xbins_numubar);
    numubar_data_histo->SetLineWidth(2);
    numubar_data_histo->SetMarkerStyle(8);
    numubar_data_histo->SetMarkerSize(0.8);
    numubar_data_histo->SetBinErrorOption(TH1::kPoisson);
    TH1F *numubar_noosc_histo = new TH1F("numubar_noosc_histo","",NB_numubar,xbins_numubar);
    numubar_noosc_histo->SetLineWidth(2);
    numubar_noosc_histo->SetLineColor(kRed);
    TH1F *numubar_bkgd_histo  = new TH1F("numubar_bkgd_histo" ,"",NB_numubar,xbins_numubar);
    numubar_bkgd_histo->SetLineWidth(1);
    numubar_bkgd_histo->SetLineColor(kBlack);
    numubar_bkgd_histo->SetFillColor(kGray);
    TH1F *numubar_bfit_histo  = new TH1F("numubar_bfit_histo" ,"",NB_numubar,xbins_numubar);;
    numubar_bfit_histo->SetLineWidth(2);
    numubar_bfit_histo->SetLineColor(kBlue);

    double ctnt = 0;
    double cont1;
    double sum = 0;
    for (int j = 0 ; j < NB_numubar ; j++)
    {
        //Data -------------------//
        ctnt = minos_numubar_data->GetY()[j];
        numubar_data_histo->SetBinContent(j+1,ctnt);
        //NoOsc ------------------//
        ctnt = minos_numubar_noosc->GetY()[j];
        numubar_noosc_histo->SetBinContent(j+1,ctnt);
        //Bkgd _------------------//
        ctnt = minos_numubar_bkgd->GetY()[j];
        numubar_bkgd_histo->SetBinContent(j+1,ctnt);
        //BFit _------------------//
        ctnt = minos_numubar_bfit->GetY()[j];
        numubar_bfit_histo->SetBinContent(j+1,ctnt);
    }

    // define and fill the histograms for Electron (Anti)neutrinos PRL110 ------------------
    TH1F *nue_data_histo[3];
    TH1F *nue_bkgd_histo[3];
    TH1F *nue_nueCC_histo[3];
    TH1F *nue_nuebCC_histo[3];
    TH1F *nueb_data_histo[3];
    TH1F *nueb_bkgd_histo[3];
    TH1F *nueb_nueCC_histo[3];
    TH1F *nueb_nuebCC_histo[3];
    for (int i = 0 ; i < 3 ; i++) {
        nue_data_histo[i] = new TH1F(Form("nue_data_histo_%d",i),"",NB_nue,xbins_nue);
        nue_data_histo[i]->SetLineWidth(2);
        nue_data_histo[i]->SetMarkerStyle(8);
        nue_data_histo[i]->SetMarkerSize(0.8);
        nue_data_histo[i]->SetBinErrorOption(TH1::kPoisson);
        nueb_data_histo[i] = new TH1F(Form("nueb_data_histo_%d",i),"",NB_nue,xbins_nue);
        nueb_data_histo[i]->SetLineWidth(2);
        nueb_data_histo[i]->SetMarkerStyle(8);
        nueb_data_histo[i]->SetMarkerSize(0.8);
        nueb_data_histo[i]->SetBinErrorOption(TH1::kPoisson);

        nue_bkgd_histo[i] = new TH1F(Form("nue_bkgd_histo_%d",i),"",NB_nue,xbins_nue);
        nue_bkgd_histo[i]->SetLineWidth(2);
        nue_bkgd_histo[i]->SetLineColor(kRed);
        nue_bkgd_histo[i]->SetFillColor(10);
        nueb_bkgd_histo[i] = new TH1F(Form("nueb_bkgd_histo_%d",i),"",NB_nue,xbins_nue);
        nueb_bkgd_histo[i]->SetLineWidth(2);
        nueb_bkgd_histo[i]->SetLineColor(kRed);
        nueb_bkgd_histo[i]->SetFillColor(10);

        nue_nueCC_histo[i] = new TH1F(Form("nue_nueCC_histo_%d",i),"",NB_nue,xbins_nue);
        nue_nueCC_histo[i]->SetLineColor(kViolet);
        nue_nueCC_histo[i]->SetFillColor(kViolet-1);
        nueb_nueCC_histo[i] = new TH1F(Form("nueb_nueCC_histo_%d",i),"",NB_nue,xbins_nue);
        nueb_nueCC_histo[i]->SetLineColor(kViolet-9);
        nueb_nueCC_histo[i]->SetFillColor(kViolet-8);

        nue_nuebCC_histo[i] = new TH1F(Form("nue_nuebCC_histo_%d",i),"",NB_nue,xbins_nue);
        nue_nuebCC_histo[i]->SetLineColor(kViolet);
        nue_nuebCC_histo[i]->SetFillColor(kViolet-1);
        nueb_nuebCC_histo[i] = new TH1F(Form("nueb_nuebCC_histo_%d",i),"",NB_nue,xbins_nue);
        nueb_nuebCC_histo[i]->SetLineColor(kViolet-9);
        nueb_nuebCC_histo[i]->SetFillColor(kViolet-8);

        for (int j = 0 ; j < NB_nue ; j++)
        {
            //Data -------------------//
            ctnt = minos_nue_data[i]->GetY()[j];
            nue_data_histo[i]->SetBinContent(j+1,ctnt);
            ctnt = minos_nueb_data[i]->GetY()[j];
            nueb_data_histo[i]->SetBinContent(j+1,ctnt);
            //Bkgd ------------------//
            ctnt = minos_nue_bkgd[i]->GetY()[j];
            nue_bkgd_histo[i]->SetBinContent(j+1,ctnt);
            ctnt = minos_nueb_bkgd[i]->GetY()[j];
            nueb_bkgd_histo[i]->SetBinContent(j+1,ctnt);
            //NueCC -----------------//
            ctnt = minos_nue_nueCC[i]->GetY()[j];
            nue_nueCC_histo[i]->SetBinContent(j+1,ctnt);
            ctnt = minos_nueb_nueCC[i]->GetY()[j];
            nueb_nueCC_histo[i]->SetBinContent(j+1,ctnt);
            //NuebCC ----------------//
            ctnt = minos_nue_nuebCC[i]->GetY()[j];
            nue_nuebCC_histo[i]->SetBinContent(j+1,ctnt);
            ctnt = minos_nueb_nuebCC[i]->GetY()[j];
            nueb_nuebCC_histo[i]->SetBinContent(j+1,ctnt);
        }
    }
    //---------------------------------
	
    //- Frame for Muon Antineutrinos PRL108
	TH2F *frame_spectra = new TH2F("frame_spectra","",NB_numubar,lo,hi,10,0,65.0);
    frame_spectra->GetXaxis()->SetTitle("Reconstructed Antineutrino Energy (GeV)");
    frame_spectra->GetXaxis()->SetTitleFont(ft);
    frame_spectra->GetXaxis()->SetTitleOffset(1.15);
    frame_spectra->GetXaxis()->CenterTitle();
    frame_spectra->GetXaxis()->SetTitleSize(1.0*sz);
    frame_spectra->GetXaxis()->SetLabelSize(0.9*sz);
    frame_spectra->GetXaxis()->SetLabelFont(ft);
    frame_spectra->GetYaxis()->SetTitle("Events / GeV");
    frame_spectra->GetYaxis()->SetTitleFont(ft);
    frame_spectra->GetYaxis()->SetTitleOffset(1.1);
    frame_spectra->GetYaxis()->CenterTitle();
    frame_spectra->GetYaxis()->SetTitleSize(1.0*sz);
    frame_spectra->GetYaxis()->SetLabelSize(0.9*sz);
    frame_spectra->GetYaxis()->SetLabelFont(ft);
	
    //- Frame for Electron neutrinos PRL110
    TH2F *frame_spectra_nue = new TH2F("frame_spectra_nue","",NB_nue,loe,hie,10,0,36.0);
    frame_spectra_nue->GetXaxis()->SetTitle("Reconstructed Energy (GeV)");
    frame_spectra_nue->GetXaxis()->SetTitleFont(ft);
    frame_spectra_nue->GetXaxis()->SetTitleOffset(1.15);
    frame_spectra_nue->GetXaxis()->CenterTitle();
    frame_spectra_nue->GetXaxis()->SetTitleSize(1.0*sz);
    frame_spectra_nue->GetXaxis()->SetLabelSize(0.9*sz);
    frame_spectra_nue->GetXaxis()->SetLabelFont(ft);
    frame_spectra_nue->GetYaxis()->SetTitle("Events");
    frame_spectra_nue->GetYaxis()->SetTitleFont(ft);
    frame_spectra_nue->GetYaxis()->SetTitleOffset(1.1);
    frame_spectra_nue->GetYaxis()->CenterTitle();
    frame_spectra_nue->GetYaxis()->SetTitleSize(1.0*sz);
    frame_spectra_nue->GetYaxis()->SetLabelSize(0.9*sz);
    frame_spectra_nue->GetYaxis()->SetLabelFont(ft);
    //- Frame for Electron antineutrinos PRL110
    TH2F *frame_spectra_nueb = new TH2F("frame_spectra_nueb","",NB_nue,loe,hie,10,0,8.5);
    frame_spectra_nueb->GetXaxis()->SetTitle("Reconstructed Energy (GeV)");
    frame_spectra_nueb->GetXaxis()->SetTitleFont(ft);
    frame_spectra_nueb->GetXaxis()->SetTitleOffset(1.15);
    frame_spectra_nueb->GetXaxis()->CenterTitle();
    frame_spectra_nueb->GetXaxis()->SetTitleSize(1.0*sz);
    frame_spectra_nueb->GetXaxis()->SetLabelSize(0.9*sz);
    frame_spectra_nueb->GetXaxis()->SetLabelFont(ft);
    frame_spectra_nueb->GetYaxis()->SetTitle("Events");
    frame_spectra_nueb->GetYaxis()->SetTitleFont(ft);
    frame_spectra_nueb->GetYaxis()->SetTitleOffset(1.1);
    frame_spectra_nueb->GetYaxis()->CenterTitle();
    frame_spectra_nueb->GetYaxis()->SetTitleSize(1.0*sz);
    frame_spectra_nueb->GetYaxis()->SetLabelSize(0.9*sz);
    frame_spectra_nueb->GetYaxis()->SetLabelFont(ft);

    // Drawing section
    //------------------------------------
    TLatex *lat = new TLatex();
    lat->SetNDC();
    lat->SetTextFont(ft);
    lat->SetTextSize(1.5*sz);
    
    TLegend *leg11 = new TLegend(0.35,0.6,0.8,0.8);
    leg11->SetTextFont(ft);
    leg11->SetTextSize(0.8*sz);
    leg11->SetFillColor(0);
    leg11->SetLineColor(0);
    
    leg11->AddEntry(numubar_data_histo,"MINOS FD Data","pe0l");
    leg11->AddEntry(numubar_noosc_histo,"Prediction w/o oscillations","l");
    leg11->AddEntry(numubar_bfit_histo,"Best fit","l");
    leg11->AddEntry(numubar_bkgd_histo,"Background w/ oscillations","f");

    TCanvas *canv0 = new TCanvas("canv0","",700,600);
    TGaxis::SetMaxDigits(3);
    
    frame_spectra->Draw();
    numubar_data_histo->Draw("pe1 same");
    numubar_noosc_histo->Draw("h same");
    numubar_bfit_histo->Draw("h same");
    numubar_bkgd_histo->Draw("h same");
    leg11->Draw();
    gPad->RedrawAxis();
    gPad->SetTicks(1,1);

    //--
    TLegend *legnue = new TLegend(0.5,0.5,0.85,0.8);
    legnue->SetTextFont(ft);
    legnue->SetTextSize(1.2*sz);
    legnue->SetFillColor(0);
    legnue->SetLineColor(0);
    
    legnue->AddEntry(nue_data_histo[0],  "MINOS FD Data","pe0l");
    legnue->AddEntry(nue_bkgd_histo[0],  "Background","l");
    legnue->AddEntry(nue_nueCC_histo[0], "#nu_{e}-CC Signal","f");
    legnue->AddEntry(nueb_nuebCC_histo[0],"#bar{#nu}_{e}-CC Signal","f");
    //------------------------------------

    TCanvas *canvnue = new TCanvas("canvnue","",2*450,3*400);
    canvnue->Divide(2,3);
    TGaxis::SetMaxDigits(3);
    
    canvnue->cd(1);
    frame_spectra_nue->Draw();
    nue_nueCC_histo[0]->Draw("h same");
    nue_nuebCC_histo[0]->Draw("h same");
    nue_bkgd_histo[0]->Draw("h same");
    nue_data_histo[0]->Draw("pe0e1 same");
    lat->DrawLatex(0.15,0.8,"0.6 < #alpha_{LEM} < 0.7");
    lat->DrawLatex(0.15,0.7,"#nu Mode");
    gPad->RedrawAxis();
    gPad->SetTicks(1,1);

    canvnue->cd(3);
    frame_spectra_nue->Draw();
    nue_nueCC_histo[1]->Draw("h same");
    nue_nuebCC_histo[1]->Draw("h same");
    nue_bkgd_histo[1]->Draw("hf same");
    nue_data_histo[1]->Draw("pe0e1 same");
    lat->DrawLatex(0.15,0.8,"0.7 < #alpha_{LEM} < 0.8");
    lat->DrawLatex(0.15,0.7,"#nu Mode");
    gPad->RedrawAxis();
    gPad->SetTicks(1,1);

    canvnue->cd(5);
    frame_spectra_nue->Draw();
    nue_nueCC_histo[2]->Draw("h same");
    nue_nuebCC_histo[2]->Draw("h same");
    nue_bkgd_histo[2]->Draw("hf same");
    nue_data_histo[2]->Draw("pe0e1 same");
    legnue->Draw();
    lat->DrawLatex(0.15,0.8,"#alpha_{LEM} > 0.8");
    lat->DrawLatex(0.15,0.7,"#nu Mode");
    gPad->RedrawAxis();
    gPad->SetTicks(1,1);

    canvnue->cd(2);
    frame_spectra_nueb->Draw();
    nueb_nueCC_histo[0]->Draw("h same");
    nueb_nuebCC_histo[0]->Draw("h same");
    nueb_bkgd_histo[0]->Draw("hf same");
    nueb_data_histo[0]->Draw("pe0e1 same");
    lat->DrawLatex(0.15,0.8,"0.6 < #alpha_{LEM} < 0.7");
    lat->DrawLatex(0.15,0.7,"#bar{#nu} Mode");
    gPad->RedrawAxis();
    gPad->SetTicks(1,1);

    canvnue->cd(4);
    frame_spectra_nueb->Draw();
    nueb_nueCC_histo[1]->Draw("h same");
    nueb_nuebCC_histo[1]->Draw("h same");
    nueb_bkgd_histo[1]->Draw("hf same");
    nueb_data_histo[1]->Draw("pe0e1 same");
    lat->DrawLatex(0.15,0.8,"0.7 < #alpha_{LEM} < 0.8");
    lat->DrawLatex(0.15,0.7,"#bar{#nu} Mode");
    gPad->RedrawAxis();
    gPad->SetTicks(1,1);

    canvnue->cd(6);
    frame_spectra_nueb->Draw();
    nueb_nueCC_histo[2]->Draw("h same");
    nueb_nuebCC_histo[2]->Draw("h same");
    nueb_bkgd_histo[2]->Draw("hf same");
    nueb_data_histo[2]->Draw("pe0e1 same");
    lat->DrawLatex(0.15,0.8,"#alpha_{LEM} > 0.8");
    lat->DrawLatex(0.15,0.7,"#bar{#nu} Mode");
    gPad->RedrawAxis();
    gPad->SetTicks(1,1);

/*
    canv1->Print("Plots/RENO_bump.eps");
    canv1->Print("Plots/RENO_bump.pdf");
*/
}
