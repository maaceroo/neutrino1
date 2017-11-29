// txt files are digitalizations of data plots in the article, by using Engauge software.
// plot from RENO Coll. (2017) - arXiv:1610.04326
// This macro can be executed using in root [0] .x renograph.C
// The output from code is a Root file "RENOplots.root" 
void renograph(){

   //------------- Style --------------
    gROOT->SetStyle("Plain");
    gStyle->SetOptStat(0);
    //------------- Style --------------
    //----------  Text Style  ---------
    //ft = 10 * fontID + precision
    Int_t ft = 10 * 4 + 2;
    Double_t sz = 0.04;
    //---------------------------------

	// READING FILES.  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////

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

	// Predictions 
	TGraph *reno_bestfit  = new TGraph("data/bestfit.txt","%lg,%lg","");
	TGraph *reno_noosc  = new TGraph("data/noosc.txt","%lg,%lg","");
		

	// define number of bins //////////////////////////////////////////////////////////////////////////////////////////////77
	const int  NB = 27; 
 	const double lo = 1.2;  
 	const double hi = 8.4;
	const int nd = 2;
	double xbins[NB+1];
	double delta_bins2 = (6.0 - 1.2)/24; // 0.2 MeV/bin
	
	for (int i = 0 ; i < (NB-2) ; i++)
        {
            xbins[i+1] = 1.2 + delta_bins2*i;
	
	}
	xbins[25] = xbins[24] + 0.4;
    	xbins[26] = 8.4 - 1.4;
    	xbins[27] = 8.4;

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 

	// define the histograms
	
	//If nd = 0 then we are talking about data of the Near detector; if nd =1 then is data of the far detector
	TH1F *data_spect_histo[nd];
	for (int n = 0 ; n < nd ; n++){

	data_spect_histo[n] = new TH1F(Form("data_spect_histo_%d",n),"",NB,xbins);
	data_spect_histo[n]->SetLineWidth(2);
	data_spect_histo[n]->SetMarkerStyle(34);
	data_spect_histo[n]->SetMarkerSize(1.1);
	}	

	//  Monte Carlo
	TH1F *near_spect_histo = new TH1F("near_spect_histo","",NB,xbins); //// to fill near detector Monte Carlo simulation
	near_spect_histo->SetLineColor(kBlue);
	
	TH1F *far_spect_histo = new TH1F("far_spect_histo","",NB,xbins);//to Fill Far detector Monte Carlo
	far_spect_histo->SetLineColor(kBlue);
		
	
	//If nd = 0 then we are talking about the ratio between far- Near data; if nd =1 then is spect ratio far-near detector
	 TH1F *ratios_histo[nd];
	for (int n = 0 ; n < nd ; n++){

	  ratios_histo[n] = new TH1F(Form("ratios_histo_%d",n),"",NB,xbins);
	  ratios_histo[n]->SetLineWidth(2);
	  ratios_histo[n]->SetMarkerStyle(34);
	  ratios_histo[n]->SetMarkerSize(1.1);
	}

	// Define the legend of the plots
	TLegend *leg1 = new TLegend(.75,.80,.95,.95,"Near Detector");
	leg1->AddEntry(near_spect_histo,"MC");
	leg1->AddEntry(data_spect_histo[0],"data");
	

	// Near detector backgrounds
	TH1F *reno_bg_near_fast_histo = new TH1F("reno_bg_near_fast_histo","",NB,xbins);
	TH1F *reno_bg_near_accident_histo  = new TH1F("reno_bg_near_accident_histo", "",NB,xbins);
	TH1F *reno_bg_near_lihe_histo  = new TH1F("reno_bg_near_lihe_histo", "",NB,xbins);
	TH1F *reno_bg_near_cf_histo  = new TH1F("reno_bg_near_cf_histo", "",NB,xbins);

	TLegend *leg2 = new TLegend(.75,.80,.95,.95,"Near Detector");
	leg2->AddEntry(reno_bg_near_cf_histo,"Cf");
	leg2->AddEntry(reno_bg_near_fast_histo," Fast Neutron");
	leg2->AddEntry(reno_bg_near_accident_histo,"Accidental");
	leg2->AddEntry(reno_bg_near_lihe_histo,"Li He");
	
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

	TLegend *leg4 = new TLegend(.75,.80,.95,.95,"Far Detector");
	leg4->AddEntry(far_spect_histo,"MC");
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

	// Far detector  -  Predictions
	TH1F *reno_bestfit_histo = new TH1F("reno_bestfit_histo","",NB,xbins);
	reno_bestfit_histo->SetLineColor(kRed);
	TH1F *reno_noosc_histo = new TH1F("reno_noosc_histo", "",NB,xbins);
	reno_noosc_histo->SetLineColor(kBlue);
	TLegend *leg5 = new TLegend(.75,.80,.95,.95,"Far Detector");
	leg5->AddEntry(reno_bestfit_histo,"Prediction (Best fit)");
	leg5->AddEntry(reno_noosc_histo,"Prediction (No oscillation)");

	ofstream ratiodata,ratiospect;
	string ratiosdata = "files/ratiodata.txt";
	string ratiosspect = "files/ratiospect.txt";
	ratiodata.open((ratiosdata).c_str());
	ratiospect.open((ratiosspect).c_str());
	
	double ctnt = 0;
	double ctnt1 = 0;
	double cont1,cont2;
	double sum = 0;
	for (int j = 0 ; j < NB ; j++)
	  {

	    ////////////////////////////////////////////7/// ratio histo
	    
	    ctnt = reno_data_neardet->GetY()[j];
	    ctnt1 = reno_data_fardet->GetY()[j];
	    cont1 = ctnt1/ctnt;
	    ratios_histo[0]->SetBinContent(j+1,cont1);
	    //cout << "j + 1 = " << j + 1 << "   ctnt1 = " << ctnt1 << "ctnt = " << ctnt << "cont1 = " << cont1 <<endl;	
	    ratiodata  << cont1 << "\t"  << endl;
	    
	    ctnt = reno_espect_neardet->GetY()[j];
	    ctnt1 = reno_espect_fardet->GetY()[j];
	    cont2 = ctnt1/ctnt;
	    ratios_histo[1]->SetBinContent(j+1,cont2);
	    //cout << "j + 1 = " << j + 1 << "   ctnt1 = " << ctnt1 << "ctnt = " << ctnt << "cont2 = " << cont2 <<endl;	
	    
	    ratiospect << cont2 << "\t" << endl;
	    //	ratios <<endl;
	    
	    /////////////////////////////////////////////////////////////////////////////////////////
	    
	    //Near detector
	    // data and Monte carlo
	    
	    ctnt = reno_data_neardet->GetY()[j];
	    //	ctnt1 =  reno_data_neardet->GetY()[12];
	    data_spect_histo[0]->SetBinContent(j+1,ctnt);
	    // cout << "j + 1 = " << j + 1 << "   ctnt = " << ctnt << "ctnt = " << ctnt1 <<endl;
	    
	    ctnt = reno_espect_neardet->GetY()[j];
	    near_spect_histo->SetBinContent(j+1,ctnt);

	    //////////////////////////////////////////////////////////////////
	   		
	    // Backgrounds

	    // lihe
	    ctnt = reno_bg_near_lihe->GetY()[j];
	    reno_bg_near_lihe_histo->SetBinContent(j+1,ctnt);
	    
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
	    
	    ctnt = reno_espect_fardet->GetY()[j];
	    //cout << "j + 1 = " << j + 1 << "   ctnt = " << ctnt << endl;
	    far_spect_histo->SetBinContent(j+1,ctnt);	
	    
	    ctnt = reno_data_fardet->GetY()[j];
	    data_spect_histo[1]->SetBinContent(j+1,ctnt);
	    
	    // Backgrounds
	    // lihe
	    ctnt = reno_bg_far_lihe->GetY()[j];
	    reno_bg_far_lihe_histo->SetBinContent(j+1,ctnt);
	    
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
	    ctnt = reno_bestfit->GetY()[j];
	    reno_bestfit_histo->SetBinContent(j+1,ctnt);
	    ctnt = reno_noosc->GetY()[j];
	    reno_noosc_histo->SetBinContent(j+1,ctnt);
	    
	  }	
	
	
	TH2F *frame_spect_histo = new TH2F("frame_spect_histo","",NB,lo,hi,10,0,18300);
	frame_spect_histo->GetXaxis()->SetTitle("Prompt Energy (MeV)");
	frame_spect_histo->GetYaxis()->SetTitle("Events/0.2 MeV");	 
	
	TH2F *frame_spectrafd = new TH2F("frame_spectrafd","",NB,lo,hi,10,0,1890);
	frame_spectrafd->GetXaxis()->SetTitle("Prompt Energy (MeV)");
	frame_spectrafd->GetYaxis()->SetTitle("Events/0.2 MeV");
   
	TH2F *frame_spectrand = new TH2F("frame_spectrand","",NB,lo,hi,10,0,18300);
	frame_spectrand->GetXaxis()->SetTitle("Prompt Energy (MeV)");
	frame_spectrand->GetYaxis()->SetTitle("Events/0.2 MeV");
	
	TH2F *frame_backnd = new TH2F("frame_backnd","",NB,lo,hi,10,0,800);
	frame_backnd->GetXaxis()->SetTitle("Prompt Energy (MeV)");
	frame_backnd->GetYaxis()->SetTitle("Events/0.2 MeV");
	
	TH2F *frame_backfd = new TH2F("frame_backfd","",NB,lo,hi,10,0,130);
	frame_backfd->GetXaxis()->SetTitle("Prompt Energy (MeV)");
	frame_backfd->GetYaxis()->SetTitle("Events/0.2 MeV");
	
	TH2F *frame_pred = new TH2F("frame_pred","",NB,lo,hi,10,0,2010);
	frame_pred->GetXaxis()->SetTitle("Prompt Energy (MeV)");
	frame_pred->GetYaxis()->SetTitle("Events/0.2 MeV");

	// Drawing section
	
	// near Drawing seccion
	// Graph of near detector
	TCanvas* c = new TCanvas("c","NEAR DETECTOR",700,500);
	frame_spectrand->Draw();
	
	data_spect_histo[0]->Draw("P same");
	near_spect_histo->Draw("same");
	leg1->Draw();
	c->Print("Plots/Neardetector.pdf");	
	
	TCanvas *c2 = new TCanvas("c2","Near Detector",700,500);
	//gPad->SetLogy();
	//hst->Draw();
	frame_backnd->Draw();
	reno_bg_near_cf_histo->Draw("same");
	reno_bg_near_lihe_histo->Draw("same");
	reno_bg_near_accident_histo->Draw("same");
	reno_bg_near_fast_histo->Draw("same");
	leg2->Draw();
	// far Drawing section
	c2->Print("Plots/Neardetectorbg.pdf");
	 
	TCanvas *c3 = new TCanvas("c3","",700,500);
	//gPad->SetLogy();
	//hst1->Draw();
	frame_backfd->Draw();
	reno_bg_far_lihe_histo->Draw("same");
	reno_bg_far_acci_histo->Draw("same");
	reno_bg_far_cf_histo->Draw("same");
	reno_bg_far_fast_histo->Draw("same");
	leg3->Draw();
	c3->Print("Plots/fardetector.pdf");
	
	// Graph of far detector
	TCanvas* c1 = new TCanvas("c1","",700,500);
	//h->Draw("");
	frame_spectrafd->Draw();
	data_spect_histo[1]->Draw("P same");
	far_spect_histo->Draw("same");
	leg4->Draw();
	c1->Print("Plots/fardetectorbg.pdf");
	
	// Graph of far detector
	TCanvas* c4 = new TCanvas("c4","",700,500);
	//h->Draw("");
	frame_pred->Draw();
	reno_bestfit_histo->Draw("same");
	reno_noosc_histo->Draw("same");
	leg5->Draw();
	c4->Print("Plots/Predictions.pdf");
	 
	TCanvas* c7 = new TCanvas("c7","Ratio_data_detector",700,500);
	ratios_histo[0]->Draw();
	c7->Print("Plots/ratio1.pdf");
	TCanvas* c8 = new TCanvas("c8","Ratio_spect_detector",700,500);
	ratios_histo[1]->Draw();
	c8->Print("Plots/ratio2.pdf");
	// write to output file
	TFile *fout = new TFile("files_root/RENOplots.root","recreate");
	fout->cd();
	
	reno_bestfit_histo->Write();
	reno_noosc_histo->Write();
	
	for (int i = 0 ; i < nd ; i++)
	  {
	    //spectra (Events / MeV)
	    data_spect_histo[i]->Write();
	    ratios_histo[i]->Write();
	  }
	// near_data_histo->Write();
	near_spect_histo->Write();
	
	//far_data_histo->Write();
	far_spect_histo->Write();
	
	reno_bg_far_lihe_histo->Write();
	reno_bg_far_acci_histo->Write();
	reno_bg_far_cf_histo->Write();
	reno_bg_far_fast_histo->Write();
	reno_bg_far_lihe_histo->Write();
	reno_bg_far_acci_histo->Write();
	reno_bg_far_cf_histo->Write();
	reno_bg_far_fast_histo->Write();
	
	fout->Close();
	
	
}
