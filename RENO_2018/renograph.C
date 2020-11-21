// txt files are digitalizations of data plots in the article, by using Engauge software.
// plot from RENO Coll. (2017) - arXiv:1610.04326
// This macro can be executed using in root [0] .x renograph.C
// The output from code is a Root file "RENOplots.root" 
#include "constants.h"

void renograph()
{
    
    //------------- Style --------------
    gROOT->SetStyle("Plain");
    gStyle->SetOptStat(0);
    //------------- Style --------------
    //----------  Text Style  ---------
    //ft = 10 * fontID + precision
    Int_t ft = 10 * 4 + 2;
    Double_t sz = 0.04;
    //---------------------------------

	// READING FILES.  //

    //Antineutrino detector #1- near detector. Data (point) and Monte Carlo (histogram).
    TGraph *reno_data_neardet     = new TGraph("data/ND_Data.csv","%lg,%lg","");
    TGraph *reno_mc_neardet       = new TGraph("data/ND_MC.csv","%lg,%lg","");

    //Antineutrino detector #2- far detector. Data (point) and Monte Carlo (histogram).
    TGraph *reno_data_fardet      = new TGraph("data/FD_Data.csv","%lg,%lg","");
    TGraph *reno_mc_fardet        = new TGraph("data/FD_MC.csv","%lg,%lg","");

    // Background near detector
    TGraph *reno_bg_near_fast     = new TGraph("data/ND_FstN.csv","%lg,%lg","");
    TGraph *reno_bg_near_lihe     = new TGraph("data/ND_LiHe.csv","%lg,%lg","");
    TGraph *reno_bg_near_accident = new TGraph("data/ND_Acc.csv","%lg,%lg","");
    TGraph *reno_bg_near_cf       = new TGraph("data/ND_Cf.csv","%lg,%lg","");

    // Background far detector
    TGraph *reno_bg_far_fast     = new TGraph("data/FD_FstN.csv","%lg,%lg","");
    TGraph *reno_bg_far_lihe     = new TGraph("data/FD_LiHe.csv","%lg,%lg","");
    TGraph *reno_bg_far_accident = new TGraph("data/FD_Acc.csv","%lg,%lg","");
    TGraph *reno_bg_far_cf       = new TGraph("data/FD_Cf.csv","%lg,%lg","");

    // define number of bins
    //const int  NB = 26;
    //const double lo = 1.2;
    //const double hi = 8.4;
    //const int nDet = 2;
    double xbins[NB+1];
    double delta_bins2 = (5.6 - 1.2)/22; // 0.2 MeV/bin
    
    for (int i = 0 ; i < (NB-3) ; i++){
        xbins[i] = 1.2 + delta_bins2*i;
    }
    xbins[23] = xbins[22] + 0.4;
    xbins[24] = xbins[23] + 0.4;
    xbins[25] = 8.4 - 1.4;
    xbins[26] = 8.4;

	// define the histograms
	// nDet = 0 -> Near detector; nDet = 1 -> Far detector
	TH1F *data_spect_histo[nDet];
    TH1F *mc_histo[nDet];
    TH1F *bkgd_histo[nDet];
	for (int n = 0 ; n < nDet ; n++){
        data_spect_histo[n] = new TH1F(Form("data_spect_histo_%d",n),"",NB,xbins);
        data_spect_histo[n]->SetLineWidth(2);
        data_spect_histo[n]->SetMarkerStyle(34);
        data_spect_histo[n]->SetMarkerSize(1.1);
	
        mc_histo[n] = new TH1F(Form("mc_histo_%d",n),"",NB,xbins);
        mc_histo[n]->SetLineWidth(2);
        mc_histo[n]->SetLineColor(kBlue);

        bkgd_histo[n] = new TH1F(Form("bkgd_histo_%d",n),"",NB,xbins);
        bkgd_histo[n]->SetLineWidth(2);
        bkgd_histo[n]->SetLineColor(6);
	}	
	
    // nDet = 0 -> Near detector; nDet = 1 -> Far detector
    /*
	TH1F *ratio_histo[nDet];
	for (int n = 0 ; n < nDet ; n++){
        ratio_histo[n] = new TH1F(Form("ratio_histo_%d",n),"",NB,xbins);
        ratio_histo[n]->SetLineWidth(2);
        ratio_histo[n]->SetMarkerStyle(34);
        ratio_histo[n]->SetMarkerSize(1.1);
	}
	*/
	// Define the legend of the plots
	TLegend *leg1 = new TLegend(.75,.80,.95,.95,"Near Detector");
	leg1->AddEntry(mc_histo[0],"MC");
	leg1->AddEntry(data_spect_histo[0],"Data");
		
	// Near detector backgrounds
	TH1F *reno_bg_near_fast_histo      = new TH1F("reno_bg_near_fast_histo","",NB,xbins);
	TH1F *reno_bg_near_accident_histo  = new TH1F("reno_bg_near_accident_histo", "",NB,xbins);
	TH1F *reno_bg_near_lihe_histo      = new TH1F("reno_bg_near_lihe_histo", "",NB,xbins);
	TH1F *reno_bg_near_cf_histo        = new TH1F("reno_bg_near_cf_histo", "",NB,xbins);
	
	TLegend *leg2 = new TLegend(.75,.80,.95,.95,"Near Detector");
	leg2->AddEntry(reno_bg_near_cf_histo,"Cf");
	leg2->AddEntry(reno_bg_near_fast_histo," Fast Neutron");
	leg2->AddEntry(reno_bg_near_accident_histo,"Accidental");
	leg2->AddEntry(reno_bg_near_lihe_histo,"Li He");
	
	THStack *hst = new THStack("hst","Backgrounds Near Detector");
	reno_bg_near_fast_histo->SetLineColor(7);
	reno_bg_near_fast_histo->SetLineWidth(2);
	reno_bg_near_accident_histo->SetLineColor(kBlue);
	reno_bg_near_accident_histo->SetLineWidth(2);
	reno_bg_near_lihe_histo->SetLineColor(kGreen);
	reno_bg_near_lihe_histo->SetLineWidth(2);
	reno_bg_near_cf_histo->SetLineColor(kRed);
	reno_bg_near_cf_histo->SetLineWidth(2);
	
	TLegend *leg4 = new TLegend(.75,.80,.95,.95,"Far Detector");
	leg4->AddEntry(mc_histo[1],"MC");
	leg4->AddEntry(data_spect_histo[1],"data");
	
    //Far detector backgrounds
	TH1F *reno_bg_far_fast_histo = new TH1F("reno_bg_far_fast_histo","",NB,xbins);
	TH1F *reno_bg_far_acci_histo  = new TH1F("reno_bg_far_acci_histo", "",NB,xbins);
	TH1F *reno_bg_far_lihe_histo  = new TH1F("reno_bg_far_lihe_histo", "",NB,xbins);		
	TH1F *reno_bg_far_cf_histo  = new TH1F("reno_bg_far_cf_histo", "",NB,xbins);
	
	TLegend *leg3 = new TLegend(.75,.80,.95,.95,"Far Detector");
	leg3->AddEntry(reno_bg_far_fast_histo," Fast Neutron");
	leg3->AddEntry(reno_bg_far_acci_histo,"Accidental");
	leg3->AddEntry(reno_bg_far_lihe_histo,"Li He");
	leg3->AddEntry(reno_bg_far_cf_histo,"Cf");

	THStack *hst1 = new THStack("hst1","Backgrounds far Detector");
	reno_bg_far_fast_histo->SetLineColor(7);
	reno_bg_far_fast_histo->SetLineWidth(2);
	reno_bg_far_acci_histo->SetLineColor(kBlue);
	reno_bg_far_acci_histo->SetLineWidth(2);
	reno_bg_far_lihe_histo->SetLineColor(kGreen);
	reno_bg_far_lihe_histo->SetLineWidth(2);
	reno_bg_far_cf_histo->SetLineColor(kRed);
	reno_bg_far_cf_histo->SetLineWidth(2);

	// Far detector  -  Predictions
    /*
	TH1F *reno_bestfit_histo = new TH1F("reno_bestfit_histo","",NB,xbins);
	reno_bestfit_histo->SetLineColor(kRed);
	TH1F *reno_noosc_histo = new TH1F("reno_noosc_histo", "",NB,xbins);
	reno_noosc_histo->SetLineColor(kBlue);
	TLegend *leg5 = new TLegend(.75,.80,.95,.95,"Far Detector");
	leg5->AddEntry(reno_bestfit_histo,"Prediction (Best fit)");
	leg5->AddEntry(reno_noosc_histo,"Prediction (No oscillation)");
    
	ofstream ratio_obs,ratio_spect;
	string ratios_obs = "files/ratio_obs.txt";
	string ratios_spect = "files/ratio_spect.txt";
	ratio_obs.open((ratios_obs).c_str());
	ratio_spect.open((ratios_spect).c_str());
	*/
	double ctnt = 0;
	double ctnt1 = 0;
	double cont1,cont2;
	double sum = 0;
	for (int j = 0 ; j < NB ; j++)
	  {// ratio histo
          //ctnt = reno_data_neardet->GetY()[j];
          //ctnt1 = reno_data_fardet->GetY()[j];
          //cont1 = ctnt1/ctnt;
          //ratio_histo[0]->SetBinContent(j+1,cont1);
          //cout << "j + 1 = " << j + 1 << "   ctnt1 = " << ctnt1 << "ctnt = " << ctnt << "cont1 = " << cont1 <<endl;
          //ratio_obs  << cont1 << "\t"  << endl;
	    
          //ctnt = reno_espect_neardet->GetY()[j];
          //ctnt1 = reno_espect_fardet->GetY()[j];
          //cont2 = ctnt1/ctnt;
          //ratio_histo[1]->SetBinContent(j+1,cont2);
          //cout << "j + 1 = " << j + 1 << "   ctnt1 = " << ctnt1 << "ctnt = " << ctnt << "cont2 = " << cont2 <<endl;
	    
          //ratio_spect << cont2 << "\t" << endl;
          //	ratio <<endl;
	    
          //Near detector
          // data and Monte carlo
          
          ctnt = reno_data_neardet->GetY()[j];
          data_spect_histo[0]->SetBinContent(j+1,ctnt);
          // cout << "j + 1 = " << j + 1 << "   ctnt = " << ctnt << "ctnt = " << ctnt1 <<endl;
	    
          ctnt = reno_mc_neardet->GetY()[j];
          mc_histo[0]->SetBinContent(j+1,ctnt);
          
          //////////////////////////////////////////////////////////////////
          
          // Backgrounds
          
          // lihe
          ctnt = reno_bg_near_lihe->GetY()[j];
          reno_bg_near_lihe_histo->SetBinContent(j+1,ctnt);
          bkgd_histo[0]->SetBinContent(j+1,0.0);
          
          // Accidental
          ctnt = reno_bg_near_accident->GetY()[j];
          reno_bg_near_accident_histo->SetBinContent(j+1,ctnt);
          
          // Fast
          ctnt = reno_bg_near_fast->GetY()[j];
          reno_bg_near_fast_histo->SetBinContent(j+1,ctnt);
          
          // cf
          ctnt = reno_bg_near_cf->GetY()[j];
          reno_bg_near_cf_histo->SetBinContent(j+1,ctnt);
          
          // Far detector
          ctnt = reno_mc_fardet->GetY()[j];
          //cout << "j + 1 = " << j + 1 << "   ctnt = " << ctnt << endl;
          mc_histo[1]->SetBinContent(j+1,ctnt);
          
          ctnt = reno_data_fardet->GetY()[j];
          data_spect_histo[1]->SetBinContent(j+1,ctnt);
          
          // Backgrounds
          // lihe
          ctnt = reno_bg_far_lihe->GetY()[j];
          reno_bg_far_lihe_histo->SetBinContent(j+1,ctnt);
          bkgd_histo[1]->SetBinContent(j+1,0.0);

          // Accidental
          ctnt = reno_bg_far_accident->GetY()[j];
          reno_bg_far_acci_histo->SetBinContent(j+1,ctnt);
          
          // cf
          ctnt = reno_bg_far_cf->GetY()[j];
          reno_bg_far_cf_histo->SetBinContent(j+1,ctnt);
          
          // Fast
          ctnt = reno_bg_far_fast->GetY()[j];
          reno_bg_far_fast_histo->SetBinContent(j+1,ctnt);
          //cout << "j + 1 = " << j + 1 << "   ctnt = " << ctnt << endl;
          
          //Prediction best fit and Prediction no oscillation
          //ctnt = reno_bestfit->GetY()[j];
          //reno_bestfit_histo->SetBinContent(j+1,ctnt);
          //ctnt = reno_noosc->GetY()[j];
          //reno_noosc_histo->SetBinContent(j+1,ctnt);
          
	  }

    reno_bg_near_cf_histo->Add(reno_bg_near_lihe_histo,-1);
    reno_bg_near_lihe_histo->Add(reno_bg_near_accident_histo,-1);
    reno_bg_near_accident_histo->Add(reno_bg_near_fast_histo,-1);
    hst->Add(reno_bg_near_fast_histo);
    hst->Add(reno_bg_near_accident_histo);
    hst->Add(reno_bg_near_lihe_histo);
    hst->Add(reno_bg_near_cf_histo);

    reno_bg_far_cf_histo->Add(reno_bg_far_lihe_histo,-1);
    reno_bg_far_lihe_histo->Add(reno_bg_far_acci_histo,-1);
    reno_bg_far_acci_histo->Add(reno_bg_far_fast_histo,-1);
    hst1->Add(reno_bg_far_fast_histo);
    hst1->Add(reno_bg_far_acci_histo);
    hst1->Add(reno_bg_far_lihe_histo);
    hst1->Add(reno_bg_far_cf_histo);

    std::cout << "Background" << std::endl;
    //-- Total background histograms --
    bkgd_histo[0]->Add(reno_bg_near_fast_histo);
    bkgd_histo[0]->Add(reno_bg_near_accident_histo);
    bkgd_histo[0]->Add(reno_bg_near_lihe_histo);
    bkgd_histo[0]->Add(reno_bg_near_cf_histo);

    bkgd_histo[1]->Add(reno_bg_far_fast_histo);
    bkgd_histo[1]->Add(reno_bg_far_acci_histo);
    bkgd_histo[1]->Add(reno_bg_far_lihe_histo);
    bkgd_histo[1]->Add(reno_bg_far_cf_histo);
    //---------------------------------
	
	TH2F *frame_spect_histo = new TH2F("frame_spect_histo","",NB,1,hi,10,0,52500);
	frame_spect_histo->GetXaxis()->SetTitle("Prompt Energy (MeV)");
	frame_spect_histo->GetYaxis()->SetTitle("Events/0.2 MeV");	 
	
	TH2F *frame_spectrafd = new TH2F("frame_spectrafd","",NB,1,hi,10,0,6100);
    frame_spectrafd->GetXaxis()->SetTitle("Prompt Reconstructed Energy (MeV)");
    frame_spectrafd->GetXaxis()->SetTitleFont(ft);
    frame_spectrafd->GetXaxis()->SetTitleOffset(0.9);
    frame_spectrafd->GetXaxis()->SetTitleSize(1.4*sz);
    frame_spectrafd->GetXaxis()->SetLabelSize(1.4*sz);
    frame_spectrafd->GetXaxis()->SetLabelFont(ft);
    frame_spectrafd->GetYaxis()->SetTitle("Events/0.2 MeV");
    frame_spectrafd->GetYaxis()->SetTitleFont(ft);
    frame_spectrafd->GetYaxis()->SetTitleOffset(0.7);
    frame_spectrafd->GetYaxis()->SetTitleSize(1.4*sz);
    frame_spectrafd->GetYaxis()->SetLabelSize(1.4*sz);
    frame_spectrafd->GetYaxis()->SetLabelFont(ft);
	
	TH2F *frame_spectrand = new TH2F("frame_spectrand","",NB,1,hi,10,0,52500);
    frame_spectrand->GetYaxis()->SetTitle("Events/0.2 MeV");
    frame_spectrand->GetYaxis()->SetTitleFont(ft);
    frame_spectrand->GetYaxis()->SetTitleOffset(0.7);
    frame_spectrand->GetYaxis()->SetTitleSize(1.4*sz);
    frame_spectrand->GetYaxis()->SetLabelSize(1.4*sz);
    frame_spectrand->GetYaxis()->SetLabelFont(ft);
	
	TH2F *frame_backnd = new TH2F("frame_backnd","",NB,1,hi,10,0,1250);
	frame_backnd->GetXaxis()->SetTitle("Prompt Energy (MeV)");
	frame_backnd->GetYaxis()->SetTitle("Events/0.2 MeV");
	
	TH2F *frame_backfd = new TH2F("frame_backfd","",NB,1,hi,10,0,525);
	frame_backfd->GetXaxis()->SetTitle("Prompt Energy (MeV)");
	frame_backfd->GetYaxis()->SetTitle("Events/0.2 MeV");
	
	//TH2F *frame_pred = new TH2F("frame_pred","",NB,1,hi,10,0,2010);
	//frame_pred->GetXaxis()->SetTitle("Prompt Energy (MeV)");
	//frame_pred->GetYaxis()->SetTitle("Events/0.2 MeV");

	// Drawing section
	
	// near Drawing seccion
	// Graph of near detector
	TCanvas* c = new TCanvas("c","NEAR DETECTOR",700,500);
	frame_spectrand->Draw();
	
	data_spect_histo[0]->Draw("PE same");
	mc_histo[0]->Draw("same");
	leg1->Draw();
	//c->Print("Plots/ND_Spectra.pdf");
	
	TCanvas *c2 = new TCanvas("c2","Near Detector",700,500);
	frame_backnd->Draw();
    hst->Draw();
	//reno_bg_near_cf_histo->Draw("same");
	//reno_bg_near_lihe_histo->Draw("same");
	//reno_bg_near_accident_histo->Draw("same");
	//reno_bg_near_fast_histo->Draw("same");
	leg2->Draw();
	//c2->Print("Plots/ND_BG.pdf");

	// far Drawing section
	
	TCanvas *c3 = new TCanvas("c3","",700,500);
	frame_backfd->Draw();
    hst1->Draw();
	//reno_bg_far_lihe_histo->Draw("same");
	//reno_bg_far_acci_histo->Draw("same");
	//reno_bg_far_cf_histo->Draw("same");
	//reno_bg_far_fast_histo->Draw("same");
	leg3->Draw();
	//c3->Print("Plots/FD_BG.pdf");
	
	// Graph of far detector
	TCanvas* c1 = new TCanvas("c1","",700,500);
	//h->Draw("");
	frame_spectrafd->Draw();
	data_spect_histo[1]->Draw("PE same");
	mc_histo[1]->Draw("same");
	leg4->Draw();
	//c1->Print("Plots/FD_Spectra.pdf");

    /*
	// Graph of far detector
	TCanvas* c4 = new TCanvas("c4","",700,500);
	frame_pred->Draw();
	reno_bestfit_histo->Draw("same");
	reno_noosc_histo->Draw("same");
	leg5->Draw();
	c4->Print("Plots/Predictions.pdf");
	 
	TCanvas* c7 = new TCanvas("c7","Ratio_data_detector",700,500);
	ratio_histo[0]->Draw();
	c7->Print("Plots/ratio_obs.pdf");
	TCanvas* c8 = new TCanvas("c8","Ratio_spect_detector",700,500);
	ratio_histo[1]->Draw();
	c8->Print("Plots/ratio_spect.pdf");
     */
	
	// write to output file
	//TFile *fout = new TFile("files_root/RENOplots.root","recreate");
	//const char name = getenv("NS2T");
	TString filePath = dirName;
	TString fileName = "/files_root/RENOplots.root";
	TFile *fout = new TFile(filePath + fileName,"recreate");
	fout->cd();
	
    c1->Write();
    c2->Write();
    
	//reno_bestfit_histo->Write();
	//reno_noosc_histo->Write();
    
    //------------------------------------
    TLatex *lat = new TLatex();
    lat->SetNDC();
    lat->SetTextFont(ft);
    lat->SetTextSize(2.6*sz);
    
    TLegend *leg11 = new TLegend(0.6,0.5,0.8,0.8);
    leg11->SetTextFont(ft);
    leg11->SetTextSize(1.7*sz);
    leg11->SetFillColor(0);
    leg11->SetLineColor(0);
    
    leg11->AddEntry(data_spect_histo[0],"RENO Data","p");
    leg11->AddEntry(mc_histo[0],"RENO MC","l");
    leg11->AddEntry(bkgd_histo[0],"Total Background","l");

    TCanvas *canv0 = new TCanvas("canv0","",700,600);
    TGaxis::SetMaxDigits(3);
    
    TPad *pad1 = new TPad("pad1", "pad1", 0, 1./2., 1, 1.0);
    pad1->SetBottomMargin(0); // Upper and lower plot are joined
    pad1->Draw();             // Draw the upper pad: pad1
    pad1->cd();               // pad1 becomes the current pad
    frame_spectrand->Draw();
    data_spect_histo[0]->Draw("P same");
    bkgd_histo[0]->Draw("same");
    //hst->Draw("same");
    mc_histo[0]->Draw("same");
    leg11->Draw();
    gPad->SetTicks(1,1);
    
    canv0->cd();          // Go back to the main canvas before defining pad2
    TPad *pad2 = new TPad("pad2", "pad2", 0, 0, 1, 1./2.);
    pad2->SetTopMargin(0);
    pad2->SetBottomMargin(0.11);
    pad2->Draw();
    pad2->cd();       // pad2 becomes the current pad
    frame_spectrafd->Draw();
    data_spect_histo[1]->Draw("P same");
    bkgd_histo[1]->Draw("same");
    //hst1->Draw("same");
    mc_histo[1]->Draw("same");
    gPad->RedrawAxis();
    gPad->SetTicks(1,1);
    
    //canv0->Print("Plots/RENO_Spectra.pdf");
    canv0->Print(filePath + "/Plots/RENO_Spectra.pdf");
    //------------------------------------
	
    //-- 2019.02.04 - Begin
    //-- Normalizing data histograms: it would have unity area and each bin would have units of
    //-- number of events
    double area = 0.0;
	for (int i = 0 ; i < nDet ; i++)
	  {
          for (int nbin = 0 ; nbin < NB; nbin++) {
              double widt = data_spect_histo[i]->GetBinWidth(nbin+1);
              double cont = data_spect_histo[i]->GetBinContent(nbin+1);
              data_spect_histo[i]->SetBinContent(nbin+1,0.2*widt*cont);
          }
          //area = data_spect_histo[i]->Integral();
          //data_spect_histo[i]->Scale(1.0/area);
          data_spect_histo[i]->Write();
          mc_histo[i]->Write();
          //ratio_histo[i]->Write();
	  }
    //-- 2019.02.04 - End

	reno_bg_near_lihe_histo->Write();
	reno_bg_near_accident_histo->Write();
	reno_bg_near_cf_histo->Write();
	reno_bg_near_fast_histo->Write();
	reno_bg_far_lihe_histo->Write();
	reno_bg_far_acci_histo->Write();
	reno_bg_far_cf_histo->Write();
	reno_bg_far_fast_histo->Write();
        bkgd_histo[0]->Write();
        bkgd_histo[1]->Write();
	
	fout->Close();
	
	
}
