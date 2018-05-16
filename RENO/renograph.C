//------------------------------------------------------------------------------------------//
//---------------------  renograph.C - By  D.J. Polo T. - 2017-26-09 -----------------------//
//------------------------------------------------------------------------------------------//
// txt files are digitalizations of data plots in the article, by using Engauge software.   //
//------------------- plot from RENO Coll. (2017) - arXiv:1610.04326v4 ---------------------//                                         
//---------------- This macro can be executed using in root [0] .x renograph.C -------------//
//--------------------- The output from code is a Root file "RENOplots.root" ---------------// 
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
	TGraph *reno_data_neardet  = new TGraph("files/data/dataneardet.txt","%lg,%lg","");
	TGraph *reno_espect_neardet  = new TGraph("files/data/espneardet.txt","%lg,%lg","");

	//Antineutrino detector #2- far detector. Data (point) and Monte Carlo (histogram).  
	TGraph *reno_data_fardet  = new TGraph("files/data/datafardet.txt","%lg,%lg","");
	TGraph *reno_espect_fardet  = new TGraph("files/data/espfardet.txt","%lg,%lg","");

	// Background near detector
	TGraph *reno_bg_near_fast = new TGraph("files/data/fast_near.txt","%lg,%lg","");
	TGraph *reno_bg_near_lihe  = new TGraph("files/data/lihe_near.txt","%lg,%lg","");
	TGraph *reno_bg_near_accident  = new TGraph("files/data/acci_near.txt","%lg,%lg","");
	TGraph *reno_bg_near_cf  = new TGraph("files/data/cf_near.txt","%lg,%lg","");

	// Background far detector
	TGraph *reno_bg_far_fast  = new TGraph("files/data/fast_far.txt","%lg,%lg","");
	TGraph *reno_bg_far_lihe  = new TGraph("files/data/lihe_far.txt","%lg,%lg","");
	TGraph *reno_bg_far_accident  = new TGraph("files/data/acc_far.txt","%lg,%lg","");
	TGraph *reno_bg_far_cf  = new TGraph("files/data/cf_far.txt","%lg,%lg","");

	// Predictions 
	TGraph *reno_bestfit  = new TGraph("files/data/bestfit.txt","%lg,%lg","");
	TGraph *reno_noosc  = new TGraph("files/data/noosc.txt","%lg,%lg","");
		

	// define number of bins //////////////////////////////////////////////////////////////////////////////////////////////
	
	const int  NB = 27; 
 	const double lo = 1.2;  
 	const double hi = 8.4;
	double xbins[NB+1];
	xbins[0] = 1.2;
	double delta_bins2 = (6.0 - 1.2)/24; // 0.2 MeV/bin
	
	for (int i = 0 ; i < (NB-2) ; i++)
        {
	  xbins[i+1] = 1.2 + delta_bins2*(i+1);
	
	}
	xbins[25] = xbins[24] + 0.4;
    	xbins[26] = 8.4 - 1.4;
    	xbins[27] = 8.4;

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 

	// define the histograms

	const int nd = 2;
	//If nd = 0 then we are talking about data of the Near detector; if nd =1 then is data of the far detector
	TH1F *data_spect_histo[nd];
	TH1F *spect_histo[nd];
	for (int n = 0 ; n < nd ; n++){

	data_spect_histo[n] = new TH1F(Form("data_spect_histo_%d",n),"",NB,xbins);
	data_spect_histo[n]->SetLineWidth(2);
	data_spect_histo[n]->SetMarkerStyle(34);
	data_spect_histo[n]->SetMarkerSize(1.1);

	spect_histo[n] = new TH1F(Form("spect_histo_%d",n),"",NB,xbins);
	spect_histo[n]->SetLineColor(kBlue);
	}	
	
	//If nd = 0 then we are talking about the ratio between far- Near data; if nd =1 then is MC ratio far-near detector
	TH1F *ratio_histo[nd];
	for (int n = 0 ; n < nd ; n++){
	  
	  ratio_histo[n] = new TH1F(Form("ratio_histo_%d",n),"",NB,xbins);
	  ratio_histo[n]->SetLineWidth(2);
	  ratio_histo[n]->SetMarkerStyle(34);
	  ratio_histo[n]->SetMarkerSize(1.1);
	}
	
	// Define the legend of the plots
	TLegend *leg1 = new TLegend(.2,.5,0.3,.6," Detector cercano");
	leg1->AddEntry(spect_histo[0],"MC");
	leg1->AddEntry(data_spect_histo[0],"datos");
	

	// Near detector backgrounds
	TH1F *reno_bg_near_fast_histo = new TH1F("reno_bg_near_fast_histo","",NB,xbins);
	TH1F *reno_bg_near_accident_histo  = new TH1F("reno_bg_near_accident_histo", "",NB,xbins);
	TH1F *reno_bg_near_lihe_histo  = new TH1F("reno_bg_near_lihe_histo", "",NB,xbins);
	TH1F *reno_bg_near_cf_histo  = new TH1F("reno_bg_near_cf_histo", "",NB,xbins);

	TLegend *leg2 = new TLegend(.40,.60,0.89,.89,"Detector Cercano");
	//leg2->AddEntry(reno_bg_near_cf_histo,"Cf");
	leg2->AddEntry(reno_bg_near_fast_histo,"  Neutrones R#acute{a}pidos ");
	leg2->AddEntry(reno_bg_near_accident_histo,"Accidental");
	leg2->AddEntry(reno_bg_near_lihe_histo,"Li/He");
	
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

	TLegend *leg4 = new TLegend(.2,.5,0.3,.6, "Detector Lejano");
	leg4->AddEntry(spect_histo[1],"MC");
	leg4->AddEntry(data_spect_histo[1],"datos");


	//Far detector backgrounds
	TH1F *reno_bg_far_fast_histo = new TH1F("reno_bg_far_fast_histo","",NB,xbins);
	TH1F *reno_bg_far_acci_histo  = new TH1F("reno_bg_far_acci_histo", "",NB,xbins);
	TH1F *reno_bg_far_lihe_histo  = new TH1F("reno_bg_far_lihe_histo", "",NB,xbins);		
	TH1F *reno_bg_far_cf_histo  = new TH1F("reno_bg_far_cf_histo", "",NB,xbins);
	
	TLegend *leg3 = new TLegend(.4,.6,0.89,.89," Detector Lejano");
	leg3->AddEntry(reno_bg_far_fast_histo," Neutrones R#acute{a}pidos");
	leg3->AddEntry(reno_bg_far_acci_histo,"Accidental");
	leg3->AddEntry(reno_bg_far_lihe_histo,"Li/He");
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
	    ratio_histo[0]->SetBinContent(j+1,cont1);
	    //cout << "j + 1 = " << j + 1 << "   ctnt1 = " << ctnt1 << "ctnt = " << ctnt << "cont1 = " << cont1 <<endl;	
		    
	    ctnt = reno_espect_neardet->GetY()[j];
	    ctnt1 = reno_espect_fardet->GetY()[j];
	    cont2 = ctnt1/ctnt;
	    ratio_histo[1]->SetBinContent(j+1,cont2);
	    //cout << "j + 1 = " << j + 1 << "   ctnt1 = " << ctnt1 << "ctnt = " << ctnt << "cont2 = " << cont2 <<endl;	
	    
	    
	    /////////////////////////////////////////////////////////////////////////////////////////
	    
	    //Near detector
	    // data and Monte carlo
	    
	    ctnt = reno_data_neardet->GetY()[j];
	    //	ctnt1 =  reno_data_neardet->GetY()[12];
	    data_spect_histo[0]->SetBinContent(j+1,ctnt);
	    // cout << "j + 1 = " << j + 1 << "   ctnt = " << ctnt << "ctnt = " << ctnt1 <<endl;
	    
	    ctnt = reno_espect_neardet->GetY()[j];
	    spect_histo[0]->SetBinContent(j+1,ctnt);

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
	    spect_histo[1]->SetBinContent(j+1,ctnt);	
	    
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
	
	TH2F *frame_spectrafd = new TH2F("frame_spectrafd","",NB,1,hi,10,0,1890);
	frame_spectrafd->GetXaxis()->SetTitle("                 Energ#acute{i}a de la se#tilde{n}al r#acute{a}pida (MeV)                 ");
	frame_spectrafd->GetYaxis()->SetTitle("                                  Eventos/0.2 MeV                              ");
		
	TH2F *frame_spectrand = new TH2F("frame_spectrand","",NB,1,hi,10,0,18300);
	frame_spectrand->GetXaxis()->SetTitle("                 Energ#acute{i}a de la se#tilde{n}al r#acute{a}pida (MeV)                 ");
	frame_spectrand->GetYaxis()->SetTitle("                                                       Eventos/0.2 MeV                         ");
	
	TH2F *frame_backnd = new TH2F("frame_backnd","",NB,1,hi,10,0,800);
	frame_backnd->GetXaxis()->SetTitle("                 Energ#acute{i}a de la se#tilde{n}al r#acute{a}pida (MeV)                    ");
	frame_backnd->GetYaxis()->SetTitle("                                ");
	
	TH2F *frame_backfd = new TH2F("frame_backfd","",NB,1,hi,10,0,130);
	frame_backfd->GetXaxis()->SetTitle("                 Energ#acute{i}a de la se#tilde{n}al r#acute{a}pida (MeV)                    ");
	frame_backfd->GetYaxis()->SetTitle("                                     ");
	
	TH2F *frame_pred = new TH2F("frame_pred","",NB,1,hi,10,0,2010);
	frame_pred->GetXaxis()->SetTitle("                 Energ#acute{i}a de la se#tilde{n}al r#acute{a}pida (MeV)                      ");
	frame_pred->GetYaxis()->SetTitle("                                       Eventos/0.2 MeV                                            ");

	// Drawing section
	
	// near Drawing seccion
	// Graph of near detector
	TCanvas* c = new TCanvas("c","NEAR DETECTOR",1400,900);
	frame_spectrand->Draw();
	data_spect_histo[0]->Draw("PE same");
	spect_histo[0]->Draw("same");
	leg1->Draw();
		
	TPad *subpad = new TPad("subpad","",0.54,0.54,0.89,0.89); 
	subpad->Draw(); 
	subpad->cd(); 
	//	TCanvas *c2 = new TCanvas("c2","Near Detector",700,500);
	frame_backnd->Draw();
	reno_bg_near_cf_histo->Draw("same");
	reno_bg_near_lihe_histo->Draw("same");
	reno_bg_near_accident_histo->Draw("same");
	reno_bg_near_fast_histo->Draw("same");
	leg2->Draw();
	//c2->Print("Plots/Neardetectorbg.pdf");
	c->Print("Plots/Neardetector.pdf");	

	// Graph of far detector
	TCanvas* c1 = new TCanvas("c1","",1400,900);
	frame_spectrafd->Draw();
	data_spect_histo[1]->Draw("PE same");
	spect_histo[1]->Draw("same");
	leg4->Draw();
	
	TPad *subpad2 = new TPad("subpad2","",0.54,0.54,0.89,0.89); 
	subpad2->Draw(); 
	subpad2->cd(); 
	//TCanvas *c3 = new TCanvas("c3","",700,500);
	frame_backfd->Draw();
	reno_bg_far_lihe_histo->Draw("same");
	reno_bg_far_acci_histo->Draw("same");
	reno_bg_far_cf_histo->Draw("same");
	reno_bg_far_fast_histo->Draw("same");
	leg3->Draw();
	//c3->Print("Plots/fardetectorbg.pdf");
	c1->Print("Plots/fardetector.pdf");
	
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
	
	// write to output file
	TFile *fout = new TFile("files/RENOplots.root","recreate");
	fout->cd();
	
	reno_bestfit_histo->Write();
	reno_noosc_histo->Write();
	
	for (int i = 0 ; i < nd ; i++)
	  {
	    //spectra (Events / MeV)
	    data_spect_histo[i]->Write();
	    spect_histo[i]->Write();
	    ratio_histo[i]->Write();
	  }
	
	reno_bg_near_lihe_histo->Write();
	reno_bg_near_accident_histo->Write();
	reno_bg_near_cf_histo->Write();
	reno_bg_near_fast_histo->Write();
	reno_bg_far_lihe_histo->Write();
	reno_bg_far_acci_histo->Write();
	reno_bg_far_cf_histo->Write();
	reno_bg_far_fast_histo->Write();
	
	fout->Close();
	
	
}
