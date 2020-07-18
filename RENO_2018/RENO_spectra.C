// txt files are digitalizations of data plots in the article, by using Engauge software.
// plot from RENO Coll. (2018) - arXiv:1806.00248
// This macro can be executed using in root [0] .x renograph.C
// The output from code is a Root file "RENOplots.root" 
void RENO_spectra(){

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
	const int  NB = 26;
 	const double lo = 1.2;  
 	const double hi = 8.4;
	const int nd = 2;
	double xbins[NB+1];
	double delta_bins2 = (5.6 - 1.2)/22; // 0.2 MeV/bin
	
	for (int i = 0 ; i < (NB-3) ; i++){
        xbins[i] = 1.2 + delta_bins2*i;
    }
    xbins[23] = xbins[22] + 0.4;
    xbins[24] = xbins[23] + 0.4;
    xbins[25] = 8.4 - 1.4;
    xbins[26] = 8.4;

    //Information from files/RENO_gridOscSpectra_test.txt -----//
    /*
    double noOscND_array[NB] = {
        98494.0, 138037.0, 174801.0, 208353.0, 237357.0, 258357.0, 275161.0, 281847.0, 282413.0,
        279508.0, 267610.0, 252612.0, 237130.0, 217169.0, 197348.0, 178047.0, 157754.0, 140650.0,
        123423.0, 108318.0, 94332.0, 80303.0, 67598.0, 56924.0, 83381.0, 66776.0, 28800.0
    };
    double noOscFD_array[NB] = {
        8854.0, 12333.0, 15416.0, 18469.0, 21167.0, 22764.0, 24473.0, 24924.0, 25133.0,
        24580.0, 23765.0, 22473.0, 20687.0, 18918.0, 17575.0, 15664.0, 13926.0, 12518.0,
        10918.0, 9450.0, 8381.0, 7153.0, 6126.0, 5220.0, 7494.0, 5758.0, 3358.0
    };
    //---------------------------------------------------------//
    //Information from files/RENO_gridOscSpectra_test.txt -----//
    double bfOscND_array[NB] = {
        95485.476562, 134415.687500, 170686.687500, 204011.062500, 232884.437500, 253991.312500, 271485.500000, 278526.250000, 279644.843750, 277010.906250, 265352.968750, 250575.890625, 235393.968750, 215778.984375, 196249.109375, 177153.250000, 157001.812500, 139997.750000, 122859.703125, 107847.601562, 93923.070312, 79967.867188, 67341.226562, 56726.988281, 83138.679688, 66607.570312, 28731.505859
    };
    double bfOscFD_array[NB] = {
        8395.631836, 11514.088867, 14243.035156, 16961.671875, 19391.025391, 20852.914062, 22454.701172, 22929.558594, 23198.109375, 22771.792969, 22100.732422, 20977.773438, 19383.070312, 17789.500000, 16581.134766, 14825.713867, 13221.453125, 11917.731445, 10421.828125, 9043.328125, 8038.103027, 6875.809082, 5900.748047, 5036.811035, 7249.881836, 5591.371582, 3279.748535
    };
     */
    //------------------------------------------------------------------------------------------//

	// define the histograms
	// nd = 0 -> Near detector; nd = 1 -> Far detector
    TH1F *data_spect_histo[nd];
    TH1F *mc_histo[nd];
    TH1F *bkgd_histo[nd];
    //TH1F *noosc_histo[nd];
    //TH1F *bfit_histo[nd];
	for (int n = 0 ; n < nd ; n++){
        data_spect_histo[n] = new TH1F(Form("data_spect_histo_%d",n),"",NB,xbins);
        data_spect_histo[n]->SetLineWidth(2);
        data_spect_histo[n]->SetMarkerStyle(34);
        data_spect_histo[n]->SetMarkerSize(1.1);
	
        mc_histo[n] = new TH1F(Form("mc_histo_%d",n),"",NB,xbins);
        mc_histo[n]->SetLineWidth(2);
        mc_histo[n]->SetLineColor(kBlue);

        //noosc_histo[n] = new TH1F(Form("noosc_histo_%d",n),"",NB,xbins);
        //noosc_histo[n]->SetLineColor(kBlue);
        
        //bfit_histo[n] = new TH1F(Form("bfit_histo_%d",n),"",NB,xbins);
        //bfit_histo[n]->SetLineColor(kRed);
        
        bkgd_histo[n] = new TH1F(Form("bkgd_histo_%d",n),"",NB,xbins);
        bkgd_histo[n]->SetLineWidth(2);
        bkgd_histo[n]->SetLineColor(6);
	}	
	
    // nd = 0 -> Near detector; nd = 1 -> Far detector
	// Near detector backgrounds
	TH1F *reno_bg_near_fast_histo      = new TH1F("reno_bg_near_fast_histo","",NB,xbins);
	TH1F *reno_bg_near_accident_histo  = new TH1F("reno_bg_near_accident_histo", "",NB,xbins);
	TH1F *reno_bg_near_lihe_histo      = new TH1F("reno_bg_near_lihe_histo", "",NB,xbins);
	TH1F *reno_bg_near_cf_histo        = new TH1F("reno_bg_near_cf_histo", "",NB,xbins);
	
	THStack *hst = new THStack("hst","Backgrounds Near Detector");
	reno_bg_near_fast_histo->SetLineColor(7);
	reno_bg_near_fast_histo->SetLineWidth(2);
	hst->Add(reno_bg_near_fast_histo);
	reno_bg_near_accident_histo->SetLineColor(kBlue);
	reno_bg_near_accident_histo->SetLineWidth(2);
	hst->Add(reno_bg_near_accident_histo);
	reno_bg_near_lihe_histo->SetLineColor(kGreen);
	reno_bg_near_lihe_histo->SetLineWidth(2);
	hst->Add(reno_bg_near_lihe_histo);
	reno_bg_near_cf_histo->SetLineColor(kRed);
	reno_bg_near_cf_histo->SetLineWidth(2);
	hst->Add(reno_bg_near_cf_histo);		
	
    //Far detector backgrounds
	TH1F *reno_bg_far_fast_histo = new TH1F("reno_bg_far_fast_histo","",NB,xbins);
	TH1F *reno_bg_far_acci_histo  = new TH1F("reno_bg_far_acci_histo", "",NB,xbins);
	TH1F *reno_bg_far_lihe_histo  = new TH1F("reno_bg_far_lihe_histo", "",NB,xbins);		
	TH1F *reno_bg_far_cf_histo  = new TH1F("reno_bg_far_cf_histo", "",NB,xbins);
	
	THStack *hst1 = new THStack("hst1","Backgrounds far Detector");
	reno_bg_far_fast_histo->SetLineColor(7);
	reno_bg_far_fast_histo->SetLineWidth(2);
	hst1->Add(reno_bg_far_fast_histo);
	reno_bg_far_acci_histo->SetLineColor(kBlue);
	reno_bg_far_acci_histo->SetLineWidth(2);
	hst1->Add(reno_bg_far_acci_histo);
	reno_bg_far_lihe_histo->SetLineColor(kGreen);
	reno_bg_far_lihe_histo->SetLineWidth(2);
	hst1->Add(reno_bg_far_lihe_histo);
	reno_bg_far_cf_histo->SetLineColor(kRed);
	reno_bg_far_cf_histo->SetLineWidth(2);
	hst1->Add(reno_bg_far_cf_histo);		

	double ctnt = 0;
	double ctnt1 = 0;
	double cont1,cont2;
	double sum = 0;
	for (int j = 0 ; j < NB ; j++)
    {
        //Near detector -------------------//
        ctnt = reno_data_neardet->GetY()[j];
        data_spect_histo[0]->SetBinContent(j+1,ctnt);
        // cout << "j + 1 = " << j + 1 << "   ctnt = " << ctnt << "ctnt = " << ctnt1 <<endl;
        ctnt = reno_mc_neardet->GetY()[j];
        mc_histo[0]->SetBinContent(j+1,ctnt);

        // Backgrounds
        // lihe
        ctnt = reno_bg_near_lihe->GetY()[j];
        reno_bg_near_lihe_histo->SetBinContent(j+1,ctnt);
        bkgd_histo[0]->SetBinContent(j+1,ctnt);
        
        // Accidental
        ctnt = reno_bg_near_accident->GetY()[j];
        reno_bg_near_accident_histo->SetBinContent(j+1,ctnt);
        
        // Fast
        ctnt = reno_bg_near_fast->GetY()[j];
        reno_bg_near_fast_histo->SetBinContent(j+1,ctnt);
        
        // cf
        ctnt = reno_bg_near_cf->GetY()[j];
        reno_bg_near_cf_histo->SetBinContent(j+1,ctnt);
        
        // Far detector -------------------//
        ctnt = reno_data_fardet->GetY()[j];
        data_spect_histo[1]->SetBinContent(j+1,ctnt);
        ctnt = reno_mc_fardet->GetY()[j];
        mc_histo[1]->SetBinContent(j+1,ctnt);

        // Backgrounds
        // lihe
        ctnt = reno_bg_far_lihe->GetY()[j];
        reno_bg_far_lihe_histo->SetBinContent(j+1,ctnt);
        bkgd_histo[1]->SetBinContent(j+1,ctnt);

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
        
        //BF ND
        //ctnt = bfOscND_array[j];
        //bfit_histo[0]->SetBinContent(j+1,ctnt);
        //BF FD
        //ctnt = bfOscFD_array[j];
        //bfit_histo[1]->SetBinContent(j+1,ctnt);

        //non oscillation ND
        //ctnt = noOscND_array[j];
        //noosc_histo[0]->SetBinContent(j+1,ctnt);
        //non oscillation FD
        //ctnt = noOscFD_array[j];
        //noosc_histo[1]->SetBinContent(j+1,ctnt);
    }
/*
	for (int i = 0 ; i < nd ; i++)
    {
        for (int nbin = 0 ; nbin < NB; nbin++) {
            double widt = noosc_histo[i]->GetBinWidth(nbin+1);
            double cont = noosc_histo[i]->GetBinContent(nbin+1);
            noosc_histo[i]->SetBinContent(nbin+1,cont/(0.2*widt));

            widt = bfit_histo[i]->GetBinWidth(nbin+1);
            cont = bfit_histo[i]->GetBinContent(nbin+1);
            bfit_histo[i]->SetBinContent(nbin+1,cont/(0.2*widt));
        }
    }
 */
	/*
	double areaNDdata = 623.26*458.49;
	double areaFDdata = 65.57*489.93;
    double areaNDdata_histo = areaNDdata*(4542790.111328/4592503);
    double areaFDdata_histo = areaFDdata*(380947.268066/407497);
    std::cout << data_spect_histo[0]->Integral() << "\t" << data_spect_histo[1]->Integral() << "\n"
              << areaNDdata << "\t" << areaFDdata << std::endl;
     */
    /*
    double areaND = noosc_histo[0]->Integral();
    double areaFD = noosc_histo[1]->Integral();
    noosc_histo[0]->Scale(areaNDdata/areaND);
    noosc_histo[1]->Scale(areaFDdata/areaFD);

    double areaNDbf = bfit_histo[0]->Integral();
    double areaFDbf = bfit_histo[1]->Integral();
    bfit_histo[0]->Scale(areaNDdata_histo/areaNDbf);
    bfit_histo[1]->Scale(areaFDdata_histo/areaFDbf);
    
    TH1F *bfit_histoFD_bump;
    bfit_histoFD_bump = (TH1F*) bfit_histo[1]->Clone("bfit_histoFD_bump");
    double areaFDbf_bump = bfit_histoFD_bump->Integral();
    double areaFDdt_bump = data_spect_histo[1]->Integral();
    double ratioFDdt_bf_bump = areaFDdt_bump/areaFDbf_bump;
    bfit_histoFD_bump->Scale(ratioFDdt_bf_bump);
    
    double areaNDbf_frac = bfit_histo[0]->Integral(1,12) + bfit_histo[0]->Integral(26,27);
    double areaFDbf_frac = bfit_histo[1]->Integral(1,12) + bfit_histo[1]->Integral(26,27);
    double areaNDdt_frac = data_spect_histo[0]->Integral(1,12) + data_spect_histo[0]->Integral(26,27);
    double areaFDdt_frac = data_spect_histo[1]->Integral(1,12) + data_spect_histo[1]->Integral(26,27);
    
    double ratioNDdt_bf = areaNDdt_frac/areaNDbf_frac;
    double ratioFDdt_bf = areaFDdt_frac/areaFDbf_frac;

    areaNDbf = bfit_histo[0]->Integral();
    areaFDbf = bfit_histo[1]->Integral();
    bfit_histo[0]->Scale(ratioNDdt_bf);
    bfit_histo[1]->Scale(ratioFDdt_bf);
    noosc_histo[0]->Scale(ratioNDdt_bf);
    noosc_histo[1]->Scale(ratioFDdt_bf);
     */

    std::cout << "Background" << std::endl;
    //-- Total background histograms --
    bkgd_histo[0]->Add(reno_bg_near_accident_histo);
    bkgd_histo[0]->Add(reno_bg_near_cf_histo);
    bkgd_histo[0]->Add(reno_bg_near_fast_histo);

    bkgd_histo[1]->Add(reno_bg_far_acci_histo);
    bkgd_histo[1]->Add(reno_bg_far_cf_histo);
    bkgd_histo[1]->Add(reno_bg_far_fast_histo);
    //---------------------------------
	
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
	
    // Drawing section
	
	// near Drawing seccion
	// Graph of near detector
    
    //------------------------------------
    /////////////////////////
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
    leg11->AddEntry(mc_histo[0],"MC (osc.)","l");
    //leg11->AddEntry(noosc_histo[0],"No Oscillations","l");
    //leg11->AddEntry(bfit_histo[0],"Best fit","l");
    leg11->AddEntry(bkgd_histo[0],"Total Background","l");

    TCanvas *canv0 = new TCanvas("canv0","",700,600);
    TGaxis::SetMaxDigits(3);
    
    TPad *pad1 = new TPad("pad1", "pad1", 0, 1./2., 1, 1.0);
    pad1->SetBottomMargin(0); // Upper and lower plot are joined
    pad1->Draw();             // Draw the upper pad: pad1
    pad1->cd();               // pad1 becomes the current pad
    frame_spectrand->Draw();
    data_spect_histo[0]->Draw("P same");
    mc_histo[0]->Draw("h same");
    bkgd_histo[0]->Draw("h same");
    //noosc_histo[0]->Draw("h same");
    //bfit_histo[0]->Draw("h same");
    leg11->Draw();
    lat->DrawLatex(0.7,0.3,"ND");
    gPad->SetTicks(1,1);
    
    canv0->cd();          // Go back to the main canvas before defining pad2
    TPad *pad2 = new TPad("pad2", "pad2", 0, 0, 1, 1./2.);
    pad2->SetTopMargin(0);
    pad2->SetBottomMargin(0.11);
    pad2->Draw();
    pad2->cd();       // pad2 becomes the current pad
    frame_spectrafd->Draw();
    data_spect_histo[1]->Draw("P same");
    mc_histo[1]->Draw("h same");
    bkgd_histo[1]->Draw("h same");
    //noosc_histo[1]->Draw("h same");
    //bfit_histo[1]->Draw("h same");
    lat->DrawLatex(0.7,0.3,"FD");
    gPad->RedrawAxis();
    gPad->SetTicks(1,1);
    
    //canv0->Print("Plots/canv_RENO.pdf");
    //canv0->Print("Plots/canv_RENO.eps");
    //------------------------------------

/*
    TCanvas *canv1 = new TCanvas("canv1","",700,300);
    TGaxis::SetMaxDigits(3);
    
    frame_spectrafd->Draw();
    data_spect_histo[1]->Draw("P same");
    //noosc_histo[1]->Draw("h same");
    bfit_histoFD_bump->Draw("h same");
    lat->DrawLatex(0.7,0.3,"FD");
    gPad->RedrawAxis();
    gPad->SetTicks(1,1);
    
    canv1->Print("Plots/RENO_bump.eps");
    canv1->Print("Plots/RENO_bump.pdf");
*/
}
