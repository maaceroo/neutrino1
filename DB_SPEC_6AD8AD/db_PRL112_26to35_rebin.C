// --- By M.A.Acero & A.A.Aguilar-Arévalo --- 04-Jun.2019 --- //
// data points are digitalizations of histograms in the article, by using Engauge software.
// plot from Daya Bay Coll. PRL 112, 061801 (2014) - arXiv:13106732
// Here we translate the 217Days-6AD data spectra from 26-bins to a 35-bins

void db_PRL112_26to35_rebin()
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
    TFile *fdata6AD = new TFile("PRL112_217days_data.root");
    TFile *fdata8AD = new TFile("PRD95_1230days_data.root");

    const int nEH = 3;
    //-- Events
    TH1F *data_spect_6AD26B_histoPerMeV[nEH];
    TH1F *data_spect_6AD8AD_histoPerMeV[nEH];
    for (int i = 0 ; i < nEH ; i++) {
        data_spect_6AD26B_histoPerMeV[i] = (TH1F*)fdata6AD->Get(Form("data_spect_histoPerMeV_%d",i));
        data_spect_6AD8AD_histoPerMeV[i] = (TH1F*)fdata8AD->Get(Form("data_spect_histoPerMeV_%d",i));
        data_spect_6AD8AD_histoPerMeV[i]->SetLineColor(kRed);
    }
    //----------------------------------------------------------------------------------------
    TLegend *leg0 = new TLegend(0.6,0.7,0.8,0.87);
    leg0->SetTextFont(ft);
    leg0->SetTextSize(1.1*sz);
    leg0->SetFillColor(0);
    leg0->SetLineColor(0);
    
    leg0->AddEntry(data_spect_6AD26B_histoPerMeV[0],"DB 217 Days - 6AD - 26 bins","l");

    //-- 26-bins data spectra
    TCanvas *canv0 = new TCanvas("canv0","6AD data - 26 Bins",0,0,700,3*450);
    canv0->Divide(1,3);
    
    canv0->cd(1);
    data_spect_6AD26B_histoPerMeV[0]->Draw("hist");
    leg0->Draw();
    gPad->SetTicks(1,1);
    
    canv0->cd(2);
    data_spect_6AD26B_histoPerMeV[1]->Draw("hist");
    gPad->SetTicks(1,1);

    canv0->cd(3);
    data_spect_6AD26B_histoPerMeV[2]->Draw("hist");
    gPad->SetTicks(1,1);
    
    //----------------------------------------------------------------------------------------
    //define 35-bins histogram
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

    TH1F *data_spect_6AD35B_histoPerMeV[nEH];
    TH1F *data_spect_8AD35B_histoPerMeV[nEH];
    for (int i = 0 ; i < nEH ; i++) {
        data_spect_6AD35B_histoPerMeV[i] = new TH1F(Form("data_spect_6AD35B_histoPerMeV_%d",i),"",NB,xbins);
        data_spect_6AD35B_histoPerMeV[i]->SetLineColor(kBlue);
        data_spect_8AD35B_histoPerMeV[i] = new TH1F(Form("data_spect_8AD35B_histoPerMeV_%d",i),"",NB,xbins);
        data_spect_8AD35B_histoPerMeV[i]->SetLineColor(kGreen);
    }
    double cont26 = 0.0;
    for (int j = 0 ; j < nEH ; j++) {
        //- First bin has the same content for both histograms
        cont26 = data_spect_6AD26B_histoPerMeV[j]->GetBinContent(1);
        //std::cout << "Hist26Bins 1: " << cont26 << std::endl;
        data_spect_6AD35B_histoPerMeV[j]->SetBinContent(1,cont26);
        //std::cout << "Hist35Bins 1: " << data_spect_6AD35B_histoPerMeV[j]->GetBinContent(1) << "\n" << std::endl;
        
        for (int k = 0 ; k < 6 ; k++)
        {
            cont26 = (4./5.)*data_spect_6AD26B_histoPerMeV[j]->GetBinContent(2+4*k);
            data_spect_6AD35B_histoPerMeV[j] -> SetBinContent(2+5*k,cont26);
            
            cont26 = (1./5.)*data_spect_6AD26B_histoPerMeV[j]->GetBinContent(2+4*k) +                       (3./5.)*data_spect_6AD26B_histoPerMeV[j]->GetBinContent(3+4*k);
            data_spect_6AD35B_histoPerMeV[j] -> SetBinContent(3+5*k,cont26);
            
            cont26 = (2./5.)*data_spect_6AD26B_histoPerMeV[j]->GetBinContent(3+4*k) +                       (2./5.)*data_spect_6AD26B_histoPerMeV[j]->GetBinContent(4+4*k);
            data_spect_6AD35B_histoPerMeV[j] -> SetBinContent(4+5*k,cont26);
            
            cont26 = (3./5.)*data_spect_6AD26B_histoPerMeV[j]->GetBinContent(4+4*k) +                       (1./5.)*data_spect_6AD26B_histoPerMeV[j]->GetBinContent(5+4*k);
            data_spect_6AD35B_histoPerMeV[j] -> SetBinContent(5+5*k,cont26);
            
            cont26 = (4./5.)*data_spect_6AD26B_histoPerMeV[j]->GetBinContent(5+4*k);
            data_spect_6AD35B_histoPerMeV[j] -> SetBinContent(6+5*k,cont26);
        }
        //- Last four bins of the 35-bins-histograms contain data from the last bin of the
        //- 26-bins-histograms. Last bin is the biggest (a factor of 86/98) and the previous
        //- 3 are equally filled with a factor of (4/98) each
        cont26 = data_spect_6AD26B_histoPerMeV[j]->GetBinContent(26);
        //std::cout << "Hist26Bins 26: " << cont26 << std::endl;
        data_spect_6AD35B_histoPerMeV[j]->SetBinContent(32, 4*cont26/98);
        //std::cout << "Hist35Bins 32: " << data_spect_6AD35B_histoPerMeV[j]->GetBinContent(32) << "\n" << std::endl;
        data_spect_6AD35B_histoPerMeV[j]->SetBinContent(33, 4*cont26/98);
        //std::cout << "Hist35Bins 33: " << data_spect_6AD35B_histoPerMeV[j]->GetBinContent(33) << "\n" << std::endl;
        data_spect_6AD35B_histoPerMeV[j]->SetBinContent(34, 4*cont26/98);
        //std::cout << "Hist35Bins 34: " << data_spect_6AD35B_histoPerMeV[j]->GetBinContent(34) << "\n" << std::endl;
        data_spect_6AD35B_histoPerMeV[j]->SetBinContent(35,86*cont26/98);
        //std::cout << "Hist35Bins 35: " << data_spect_6AD35B_histoPerMeV[j]->GetBinContent(35) << "\n" << std::endl;
    }
    
    double sclFac = 1e5;
    for (int i = 0 ; i < nEH ; i++) {
        data_spect_6AD8AD_histoPerMeV[i]->Scale(sclFac);
        data_spect_8AD35B_histoPerMeV[i]
        ->Add(data_spect_6AD8AD_histoPerMeV[i],data_spect_6AD35B_histoPerMeV[i],1,-1);
    }
    
    //-- comparing the total number of events in each spectrum
    double Nevents26, Nevents35;
    for (int l = 0 ; l < nEH ; l++){
	Nevents26 = data_spect_6AD26B_histoPerMeV[l]->Integral();
	Nevents35 = data_spect_6AD35B_histoPerMeV[l]->Integral();
	std::cout << "26Bins Spectrum EH " << l << ":\t" << Nevents26 << std::endl;
	std::cout << "35Bins Spectrum EH " << l << ":\t" << Nevents35 << std::endl;
	std::cout << std::endl;
	Nevents26 = 0.0;
	Nevents35 = 0.0;
    }
    
    //----------------------------------------------------------------------------------------
    TLegend *leg1 = new TLegend(0.6,0.7,0.8,0.87);
    leg1->SetTextFont(ft);
    leg1->SetTextSize(1.1*sz);
    leg1->SetFillColor(0);
    leg1->SetLineColor(0);
    
    leg1->AddEntry(data_spect_6AD8AD_histoPerMeV[0],"DB 2130 Days - 8AD 35 - Bins","l");
    leg1->AddEntry(data_spect_6AD35B_histoPerMeV[0],"DB  217 Days - 6AD 35 - Bins","l");
    leg1->AddEntry(data_spect_8AD35B_histoPerMeV[0],"DB 1913 Days - 8AD 35 - Bins","l");

    //-- 35-bins data spectra
    TCanvas *canv1 = new TCanvas("canv1","35 Bins",700,0,700,3*450);
    canv1->Divide(1,3);
    
    //sclFac = data_spect_6AD8AD_histoPerMeV[0]->Integral()/data_spect_8AD35B_histoPerMeV[0]->Integral();
    //data_spect_8AD35B_histoPerMeV[0]->Scale(sclFac);
    canv1->cd(1);
    data_spect_6AD8AD_histoPerMeV[0]->Draw("hist");
    data_spect_8AD35B_histoPerMeV[0]->Draw("same hist");
    data_spect_6AD35B_histoPerMeV[0]->Draw("same hist");
    leg1->Draw();
    gPad->SetTicks(1,1);
    
    //sclFac = data_spect_6AD8AD_histoPerMeV[1]->Integral()/data_spect_8AD35B_histoPerMeV[1]->Integral();
    //data_spect_8AD35B_histoPerMeV[1]->Scale(sclFac);
    canv1->cd(2);
    data_spect_6AD8AD_histoPerMeV[1]->Draw("hist");
    data_spect_8AD35B_histoPerMeV[1]->Draw("same hist");
    data_spect_6AD35B_histoPerMeV[1]->Draw("same hist");
    gPad->SetTicks(1,1);
    
    //sclFac = data_spect_6AD8AD_histoPerMeV[2]->Integral()/data_spect_8AD35B_histoPerMeV[2]->Integral();
    //data_spect_8AD35B_histoPerMeV[2]->Scale(sclFac);
    canv1->cd(3);
    data_spect_6AD8AD_histoPerMeV[2]->Draw("hist");
    data_spect_8AD35B_histoPerMeV[2]->Draw("same hist");
    data_spect_6AD35B_histoPerMeV[2]->Draw("same hist");
    gPad->SetTicks(1,1);

    //---------------------------------------------------------
    // write to output file
    TFile *fout = new TFile("histos_6AD217Days_8AD1913Days_data.root","recreate");
    fout->cd();
    for (int j = 0 ; j < nEH ; j++) {
        data_spect_6AD8AD_histoPerMeV[j]->Write();
        data_spect_6AD35B_histoPerMeV[j]->Write();
        data_spect_8AD35B_histoPerMeV[j]->Write();
    }
    fout->Close();

}// end
