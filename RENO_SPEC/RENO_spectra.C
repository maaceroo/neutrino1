// txt files are digitalizations of data plots in the article, by using Engauge software.
// plot from RENO Coll. (2017) - arXiv:1610.04326
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
	TGraph *reno_data_neardet  = new TGraph("data/dataneardet.txt","%lg,%lg","");
	TGraph *reno_espect_neardet  = new TGraph("data/espneardet.txt","%lg,%lg","");

	//Antineutrino detector #2- far detector. Data (point) and Monte Carlo (histogram).  
	TGraph *reno_data_fardet  = new TGraph("data/datafardet.txt","%lg,%lg","");
	TGraph *reno_espect_fardet  = new TGraph("data/espfardet.txt","%lg,%lg","");

	// Background near detector
	TGraph *reno_bg_near_fast = new TGraph("data/fast_near.txt","%lg,%lg","");
	TGraph *reno_bg_near_lihe  = new TGraph("data/lihe_near.txt","%lg,%lg","");
	TGraph *reno_bg_near_accident  = new TGraph("data/acci_near.txt","%lg,%lg","");
	TGraph *reno_bg_near_cf  = new TGraph("data/cf_near.txt","%lg,%lg","");

	// Background far detector
	TGraph *reno_bg_far_fast  = new TGraph("data/fast_far.txt","%lg,%lg","");
	TGraph *reno_bg_far_lihe  = new TGraph("data/lihe_far.txt","%lg,%lg","");
	TGraph *reno_bg_far_accident  = new TGraph("data/acc_far.txt","%lg,%lg","");
	TGraph *reno_bg_far_cf  = new TGraph("data/cf_far.txt","%lg,%lg","");

	// define number of bins
	const int  NB = 27; 
 	const double lo = 1.2;  
 	const double hi = 8.4;
	const int nd = 2;
	double xbins[NB+1];
	double delta_bins2 = (6.0 - 1.2)/24; // 0.2 MeV/bin
	
	for (int i = 0 ; i < (NB-2) ; i++){
        xbins[i] = 1.2 + delta_bins2*i;
    }
    xbins[25] = xbins[24] + 0.4;
    xbins[26] = 8.4 - 1.4;
    xbins[27] = 8.4;

    //Information from files/RENO_gridOscSpectra_test.txt -------------------------------------//
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
    //------------------------------------------------------------------------------------------//
    //Information from files/RENO_gridOscSpectra_test.txt -------------------------------------//
    double bfOscND_array[NB] = {
        95942.718750, 134933.000000, 171330.046875, 204706.843750, 233750.796875, 255043.078125, 272334.843750,
        279281.437500, 280082.843750, 277358.375000, 265712.218750, 251099.000000, 235903.718750, 216166.359375,
        196509.937500, 177345.312500, 157145.078125, 140087.828125, 122944.578125, 107938.984375, 94030.210938,
        80062.023438, 67405.906250, 56773.175781, 83197.023438, 66631.164062, 28746.216797
    };
    double bfOscFD_array[NB] = {
        8405.608398, 11561.667969, 14338.573242, 17109.140625, 19587.277344, 21081.503906, 22710.755859, 23195.269531,
        23465.232422, 23029.289062, 22343.611328, 21200.369141, 19580.912109, 17962.941406, 16736.125000, 14957.951172,
        13333.698242, 12014.272461, 10502.206055, 9109.700195, 8094.432617, 6921.635742, 5938.180176, 5067.385742,
        7290.916504, 5619.614258, 3293.185303
    };
    //------------------------------------------------------------------------------------------//

	// define the histograms
	// nd = 0 -> Near detector; nd = 1 -> Far detector
	TH1F *data_spect_histo[nd];
    TH1F *bkgd_histo[nd];
    TH1F *noosc_histo[nd];
    TH1F *bfit_histo[nd];
	for (int n = 0 ; n < nd ; n++){
        data_spect_histo[n] = new TH1F(Form("data_spect_histo_%d",n),"",NB,xbins);
        data_spect_histo[n]->SetLineWidth(2);
        data_spect_histo[n]->SetMarkerStyle(34);
        data_spect_histo[n]->SetMarkerSize(1.1);
	
        noosc_histo[n] = new TH1F(Form("noosc_histo_%d",n),"",NB,xbins);
        noosc_histo[n]->SetLineColor(kBlue);
        
        bfit_histo[n] = new TH1F(Form("bfit_histo_%d",n),"",NB,xbins);
        bfit_histo[n]->SetLineColor(kRed);
        
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
        //Near detector
        ctnt = reno_data_neardet->GetY()[j];
        data_spect_histo[0]->SetBinContent(j+1,ctnt);
        // cout << "j + 1 = " << j + 1 << "   ctnt = " << ctnt << "ctnt = " << ctnt1 <<endl;
    
        //////////////////////////////////////////////////////////////////
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
        
        // Far detector
        
        ctnt = reno_data_fardet->GetY()[j];
        data_spect_histo[1]->SetBinContent(j+1,ctnt);
        
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
        
        //non oscillation ND
        ctnt = bfOscND_array[j];
        bfit_histo[0]->SetBinContent(j+1,ctnt);
        //non oscillation FD
        ctnt = bfOscFD_array[j];
        bfit_histo[1]->SetBinContent(j+1,ctnt);

        //non oscillation ND
        ctnt = noOscND_array[j];
        noosc_histo[0]->SetBinContent(j+1,ctnt);
        //non oscillation FD
        ctnt = noOscFD_array[j];
        noosc_histo[1]->SetBinContent(j+1,ctnt);
    }

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
	//** PROBLEM HERE!! **//
	double areaNDdata = 623.26*458.49;
	double areaFDdata = 65.57*489.93;
    double areaNDdata_histo = areaNDdata*0.989029;
    double areaFDdata_histo = areaFDdata*0.93325;
    std::cout << data_spect_histo[0]->Integral() << "\t" << data_spect_histo[1]->Integral() << "\n"
              << areaNDdata << "\t" << areaFDdata << std::endl;

    double areaND = noosc_histo[0]->Integral();
    double areaFD = noosc_histo[1]->Integral();
    noosc_histo[0]->Scale(areaNDdata/areaND);
    noosc_histo[1]->Scale(areaFDdata/areaFD);

    double areaNDbf = bfit_histo[0]->Integral();
    double areaFDbf = bfit_histo[1]->Integral();
    bfit_histo[0]->Scale(areaNDdata_histo/areaNDbf);
    bfit_histo[1]->Scale(areaFDdata_histo/areaFDbf);
    
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

    std::cout << "Background" << std::endl;
    //-- Total background histograms --
    bkgd_histo[0]->Add(reno_bg_near_accident_histo);
    bkgd_histo[0]->Add(reno_bg_near_cf_histo);
    bkgd_histo[0]->Add(reno_bg_near_fast_histo);

    bkgd_histo[1]->Add(reno_bg_far_acci_histo);
    bkgd_histo[1]->Add(reno_bg_far_cf_histo);
    bkgd_histo[1]->Add(reno_bg_far_fast_histo);
    //---------------------------------
	
	TH2F *frame_spectrafd = new TH2F("frame_spectrafd","",NB,1,hi,10,0,2000);
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
	
	TH2F *frame_spectrand = new TH2F("frame_spectrand","",NB,1,hi,10,0,18500);
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
    leg11->AddEntry(noosc_histo[0],"No Oscillations","l");
    leg11->AddEntry(bfit_histo[0],"Best fit","l");
    leg11->AddEntry(bkgd_histo[0],"Total Background","l");

    TCanvas *canv0 = new TCanvas("canv0","",700,600);
    TGaxis::SetMaxDigits(3);
    
    TPad *pad1 = new TPad("pad1", "pad1", 0, 1./2., 1, 1.0);
    pad1->SetBottomMargin(0); // Upper and lower plot are joined
    pad1->Draw();             // Draw the upper pad: pad1
    pad1->cd();               // pad1 becomes the current pad
    frame_spectrand->Draw();
    data_spect_histo[0]->Draw("P same");
    bkgd_histo[0]->Draw("h same");
    noosc_histo[0]->Draw("h same");
    bfit_histo[0]->Draw("h same");
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
    bkgd_histo[1]->Draw("h same");
    noosc_histo[1]->Draw("h same");
    bfit_histo[1]->Draw("h same");
    lat->DrawLatex(0.7,0.3,"FD");
    gPad->RedrawAxis();
    gPad->SetTicks(1,1);
    
    canv0->Print("Plots/canv_RENO.pdf");
    canv0->Print("Plots/canv_RENO.eps");
    //------------------------------------
	
}
