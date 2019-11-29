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
    TH1F *data_spect_6AD26B_histo[nEH];
    TH1F *data_spect_6AD8AD_histo[nEH];
    TH1F *nosc_spect_6AD26B_histo[nEH];
    TH1F *nosc_spect_6AD8AD_histo[nEH];
    TH1F *bkgd_spect_6AD26B_histo[nEH];
    TH1F *bkgd_spect_6AD8AD_histo[nEH];
    for (int i = 0 ; i < nEH ; i++) {
        data_spect_6AD26B_histo[i] = (TH1F*)fdata6AD->Get(Form("data_spect_histo_%d",i));
        data_spect_6AD8AD_histo[i] = (TH1F*)fdata8AD->Get(Form("data_spect_histo_%d",i));
        data_spect_6AD26B_histo[i]->SetName(Form("data_spect_6AD26B_histo_%d",i));
        data_spect_6AD8AD_histo[i]->SetName(Form("data_spect_6AD8AD_histo_%d",i));
        data_spect_6AD8AD_histo[i]->SetLineColor(kRed);
        data_spect_6AD8AD_histo[i]->Sumw2();

        nosc_spect_6AD26B_histo[i] = (TH1F*)fdata6AD->Get(Form("nosc_spect_histo_%d",i));
        nosc_spect_6AD8AD_histo[i] = (TH1F*)fdata8AD->Get(Form("nosc_spect_histo_%d",i));
        nosc_spect_6AD26B_histo[i]->SetName(Form("nosc_spect_6AD26B_histo_%d",i));
        nosc_spect_6AD8AD_histo[i]->SetName(Form("nosc_spect_6AD8AD_histo_%d",i));
        nosc_spect_6AD8AD_histo[i]->SetLineColor(kRed+2);

        bkgd_spect_6AD26B_histo[i] = (TH1F*)fdata6AD->Get(Form("bkgd_spect_histo_%d",i));
        bkgd_spect_6AD8AD_histo[i] = (TH1F*)fdata8AD->Get(Form("bkgd_spect_histo_%d",i));
        bkgd_spect_6AD26B_histo[i]->SetName(Form("bkgd_spect_6AD26B_histo_%d",i));
        bkgd_spect_6AD8AD_histo[i]->SetName(Form("bkgd_spect_6AD8AD_histo_%d",i));
        bkgd_spect_6AD8AD_histo[i]->SetLineColor(kRed+1);
    }
    //----------------------------------------------------------------------------------------
    TLegend *leg0 = new TLegend(0.6,0.7,0.8,0.87);
    leg0->SetTextFont(ft);
    leg0->SetTextSize(1.1*sz);
    leg0->SetFillColor(0);
    leg0->SetLineColor(0);
    
    leg0->AddEntry(data_spect_6AD26B_histo[0],"DB 217 Days - 6AD - 26 bins","l");

    //-- 26-bins data spectra
    TCanvas *canv0 = new TCanvas("canv0","6AD data - 26 Bins",0,0,700,3*450);
    canv0->Divide(1,3);
    
    canv0->cd(1);
    data_spect_6AD26B_histo[0]->Draw("hist");
    nosc_spect_6AD26B_histo[0]->Draw("hist same");
    bkgd_spect_6AD26B_histo[0]->Draw("same hist");
    leg0->Draw();
    gPad->SetTicks(1,1);
    
    canv0->cd(2);
    data_spect_6AD26B_histo[1]->Draw("hist");
    nosc_spect_6AD26B_histo[1]->Draw("hist same");
    bkgd_spect_6AD26B_histo[1]->Draw("same hist");
    gPad->SetTicks(1,1);

    canv0->cd(3);
    data_spect_6AD26B_histo[2]->Draw("hist");
    nosc_spect_6AD26B_histo[2]->Draw("hist same");
    bkgd_spect_6AD26B_histo[2]->Draw("same hist");
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

    TH1F *data_spect_6AD35B_histo[nEH];
    TH1F *data_spect_8AD35B_histo[nEH];

    TH1F *nosc_spect_6AD35B_histo[nEH];
    TH1F *nosc_spect_8AD35B_histo[nEH];

    TH1F *bkgd_spect_6AD35B_histo[nEH];
    TH1F *bkgd_spect_8AD35B_histo[nEH];
    
    TH1F *data_spect_6AD35B_histoPerMeV[nEH];
    TH1F *data_spect_8AD35B_histoPerMeV[nEH];
    TH1F *data_spect_6AD8AD_histoPerMeV[nEH];
    
    TH1F *nosc_spect_6AD35B_histoPerMeV[nEH];
    TH1F *nosc_spect_8AD35B_histoPerMeV[nEH];
    TH1F *nosc_spect_6AD8AD_histoPerMeV[nEH];
    
    TH1F *bkgd_spect_6AD35B_histoPerMeV[nEH];
    TH1F *bkgd_spect_8AD35B_histoPerMeV[nEH];
    TH1F *bkgd_spect_6AD8AD_histoPerMeV[nEH];
    for (int i = 0 ; i < nEH ; i++) {
        data_spect_6AD35B_histo[i] = new TH1F(Form("data_spect_6AD35B_histo_%d",i),"",NB,xbins);
        data_spect_6AD35B_histo[i]->SetLineColor(kBlue);
        data_spect_8AD35B_histo[i] = new TH1F(Form("data_spect_8AD35B_histo_%d",i),"",NB,xbins);
        data_spect_8AD35B_histo[i]->SetLineColor(kGreen);

        nosc_spect_6AD35B_histo[i] = new TH1F(Form("nosc_spect_6AD35B_histo_%d",i),"",NB,xbins);
        nosc_spect_6AD35B_histo[i]->SetLineColor(kBlue+2);
        nosc_spect_8AD35B_histo[i] = new TH1F(Form("nosc_spect_8AD35B_histo_%d",i),"",NB,xbins);
        nosc_spect_8AD35B_histo[i]->SetLineColor(kGreen+2);

        bkgd_spect_6AD35B_histo[i] = new TH1F(Form("bkgd_spect_6AD35B_histo_%d",i),"",NB,xbins);
        bkgd_spect_6AD35B_histo[i]->SetLineColor(kBlue+1);
        bkgd_spect_8AD35B_histo[i] = new TH1F(Form("bkgd_spect_8AD35B_histo_%d",i),"",NB,xbins);
        bkgd_spect_8AD35B_histo[i]->SetLineColor(kGreen+1);

        data_spect_6AD35B_histoPerMeV[i] = new TH1F(Form("data_spect_6AD35B_histoPerMeV_%d",i),"",NB,xbins);
        data_spect_6AD35B_histoPerMeV[i]->SetLineColor(kBlue);
        data_spect_8AD35B_histoPerMeV[i] = new TH1F(Form("data_spect_8AD35B_histoPerMeV_%d",i),"",NB,xbins);
        data_spect_8AD35B_histoPerMeV[i]->SetLineColor(kGreen);
        
        nosc_spect_6AD35B_histoPerMeV[i] = new TH1F(Form("nosc_spect_6AD35B_histoPerMeV_%d",i),"",NB,xbins);
        nosc_spect_6AD35B_histoPerMeV[i]->SetLineColor(kBlue+2);
        nosc_spect_8AD35B_histoPerMeV[i] = new TH1F(Form("nosc_spect_8AD35B_histoPerMeV_%d",i),"",NB,xbins);
        nosc_spect_8AD35B_histoPerMeV[i]->SetLineColor(kGreen+2);
        
        bkgd_spect_6AD35B_histoPerMeV[i] = new TH1F(Form("bkgd_spect_6AD35B_histoPerMeV_%d",i),"",NB,xbins);
        bkgd_spect_6AD35B_histoPerMeV[i]->SetLineColor(kBlue+1);
        bkgd_spect_8AD35B_histoPerMeV[i] = new TH1F(Form("bkgd_spect_8AD35B_histoPerMeV_%d",i),"",NB,xbins);
        bkgd_spect_8AD35B_histoPerMeV[i]->SetLineColor(kGreen+1);

        data_spect_6AD8AD_histoPerMeV[i] = new TH1F(Form("data_spect_6AD8AD_histoPerMeV_%d",i),"",NB,xbins);
        data_spect_6AD8AD_histoPerMeV[i]->SetLineColor(kRed);
        nosc_spect_6AD8AD_histoPerMeV[i] = new TH1F(Form("nosc_spect_6AD8AD_histoPerMeV_%d",i),"",NB,xbins);
        nosc_spect_6AD8AD_histoPerMeV[i]->SetLineColor(kRed+2);
        bkgd_spect_6AD8AD_histoPerMeV[i] = new TH1F(Form("bkgd_spect_6AD8AD_histoPerMeV_%d",i),"",NB,xbins);
        bkgd_spect_6AD8AD_histoPerMeV[i]->SetLineColor(kRed+1);
        //PerMev Versions of these spectra are defined above (lines 23-33)
    }
    
    double cont26 = 0.0;
    double bkgd26 = 0.0;
    double nosc26 = 0.0;
    for (int j = 0 ; j < nEH ; j++) {
        //- First bin has the same content for both histograms
        cont26 = data_spect_6AD26B_histo[j]->GetBinContent(1);
        bkgd26 = bkgd_spect_6AD26B_histo[j]->GetBinContent(1);
        nosc26 = nosc_spect_6AD26B_histo[j]->GetBinContent(1);
        //std::cout << "Hist26Bins 1: " << cont26 << std::endl;
        data_spect_6AD35B_histo[j]->SetBinContent(1,cont26);
        bkgd_spect_6AD35B_histo[j]->SetBinContent(1,bkgd26);
        nosc_spect_6AD35B_histo[j]->SetBinContent(1,nosc26);
        //std::cout << "Hist35Bins 1: " << data_spect_6AD35B_histoPerMeV[j]->GetBinContent(1) << "\n" << std::endl;
        
        for (int k = 0 ; k < 6 ; k++)
        {
            cont26 = (4./5.)*data_spect_6AD26B_histo[j]->GetBinContent(2+4*k);
            data_spect_6AD35B_histo[j] -> SetBinContent(2+5*k,cont26);
            bkgd26 = (4./5.)*bkgd_spect_6AD26B_histo[j]->GetBinContent(2+4*k);
            bkgd_spect_6AD35B_histo[j] -> SetBinContent(2+5*k,bkgd26);
            nosc26 = (4./5.)*nosc_spect_6AD26B_histo[j]->GetBinContent(2+4*k);
            nosc_spect_6AD35B_histo[j] -> SetBinContent(2+5*k,nosc26);

            cont26 = (1./5.)*data_spect_6AD26B_histo[j]->GetBinContent(2+4*k) + (3./5.)*data_spect_6AD26B_histo[j]->GetBinContent(3+4*k);
            data_spect_6AD35B_histo[j] -> SetBinContent(3+5*k,cont26);
            bkgd26 = (1./5.)*bkgd_spect_6AD26B_histo[j]->GetBinContent(2+4*k) + (3./5.)*bkgd_spect_6AD26B_histo[j]->GetBinContent(3+4*k);
            bkgd_spect_6AD35B_histo[j] -> SetBinContent(3+5*k,bkgd26);
            nosc26 = (1./5.)*nosc_spect_6AD26B_histo[j]->GetBinContent(2+4*k) + (3./5.)*nosc_spect_6AD26B_histo[j]->GetBinContent(3+4*k);
            nosc_spect_6AD35B_histo[j] -> SetBinContent(3+5*k,nosc26);

            cont26 = (2./5.)*data_spect_6AD26B_histo[j]->GetBinContent(3+4*k) + (2./5.)*data_spect_6AD26B_histo[j]->GetBinContent(4+4*k);
            data_spect_6AD35B_histo[j] -> SetBinContent(4+5*k,cont26);
            bkgd26 = (2./5.)*bkgd_spect_6AD26B_histo[j]->GetBinContent(3+4*k) + (2./5.)*bkgd_spect_6AD26B_histo[j]->GetBinContent(4+4*k);
            bkgd_spect_6AD35B_histo[j] -> SetBinContent(4+5*k,bkgd26);
            nosc26 = (2./5.)*nosc_spect_6AD26B_histo[j]->GetBinContent(3+4*k) + (2./5.)*nosc_spect_6AD26B_histo[j]->GetBinContent(4+4*k);
            nosc_spect_6AD35B_histo[j] -> SetBinContent(4+5*k,nosc26);

            cont26 = (3./5.)*data_spect_6AD26B_histo[j]->GetBinContent(4+4*k) + (1./5.)*data_spect_6AD26B_histo[j]->GetBinContent(5+4*k);
            data_spect_6AD35B_histo[j] -> SetBinContent(5+5*k,cont26);
            bkgd26 = (3./5.)*bkgd_spect_6AD26B_histo[j]->GetBinContent(4+4*k) + (1./5.)*bkgd_spect_6AD26B_histo[j]->GetBinContent(5+4*k);
            bkgd_spect_6AD35B_histo[j] -> SetBinContent(5+5*k,bkgd26);
            nosc26 = (3./5.)*nosc_spect_6AD26B_histo[j]->GetBinContent(4+4*k) + (1./5.)*nosc_spect_6AD26B_histo[j]->GetBinContent(5+4*k);
            nosc_spect_6AD35B_histo[j] -> SetBinContent(5+5*k,nosc26);

            cont26 = (4./5.)*data_spect_6AD26B_histo[j]->GetBinContent(5+4*k);
            data_spect_6AD35B_histo[j] -> SetBinContent(6+5*k,cont26);
            bkgd26 = (4./5.)*bkgd_spect_6AD26B_histo[j]->GetBinContent(5+4*k);
            bkgd_spect_6AD35B_histo[j] -> SetBinContent(6+5*k,bkgd26);
            nosc26 = (4./5.)*nosc_spect_6AD26B_histo[j]->GetBinContent(5+4*k);
            nosc_spect_6AD35B_histo[j] -> SetBinContent(6+5*k,nosc26);
        }
        //- Last four bins of the 35-bins-histograms contain data from the last bin of the
        //- 26-bins-histograms. Last bin is the biggest (a factor of 86/94) and the previous
        //- 3 are equally filled with a factor of (4/98) each
        cont26 = data_spect_6AD26B_histo[j]->GetBinContent(26);
        bkgd26 = bkgd_spect_6AD26B_histo[j]->GetBinContent(26);
        nosc26 = nosc_spect_6AD26B_histo[j]->GetBinContent(26);
        //std::cout << "Hist26Bins 26: " << cont26 << std::endl;
        data_spect_6AD35B_histo[j]->SetBinContent(32, 4*cont26/94);
        bkgd_spect_6AD35B_histo[j]->SetBinContent(32, 4*bkgd26/94);
        nosc_spect_6AD35B_histo[j]->SetBinContent(32, 4*nosc26/94);
        //std::cout << "Hist35Bins 32: " << data_spect_6AD35B_histo[j]->GetBinContent(32) << "\n" << std::endl;
        data_spect_6AD35B_histo[j]->SetBinContent(33, 4*cont26/94);
        bkgd_spect_6AD35B_histo[j]->SetBinContent(33, 4*bkgd26/94);
        nosc_spect_6AD35B_histo[j]->SetBinContent(33, 4*nosc26/94);
        //std::cout << "Hist35Bins 33: " << data_spect_6AD35B_histo[j]->GetBinContent(33) << "\n" << std::endl;
        data_spect_6AD35B_histo[j]->SetBinContent(34, 4*cont26/94);
        bkgd_spect_6AD35B_histo[j]->SetBinContent(34, 4*bkgd26/94);
        nosc_spect_6AD35B_histo[j]->SetBinContent(34, 4*nosc26/94);
        //std::cout << "Hist35Bins 34: " << data_spect_6AD35B_histo[j]->GetBinContent(34) << "\n" << std::endl;
        data_spect_6AD35B_histo[j]->SetBinContent(35,82*cont26/94);
        bkgd_spect_6AD35B_histo[j]->SetBinContent(35,82*bkgd26/94);
        nosc_spect_6AD35B_histo[j]->SetBinContent(35,82*nosc26/94);
        //std::cout << "Hist35Bins 35: " << data_spect_6AD35B_histo[j]->GetBinContent(35) << "\n" << std::endl;
    }
    
    //double sclFac = 1e5;
    double sclFac = 1.0;
    for (int i = 0 ; i < nEH ; i++) {
        data_spect_6AD8AD_histo[i]->Scale(sclFac);
        data_spect_8AD35B_histo[i]->Add(data_spect_6AD8AD_histo[i],data_spect_6AD35B_histo[i],1,-1);
        
        nosc_spect_6AD8AD_histo[i]->Scale(sclFac);
        nosc_spect_8AD35B_histo[i]->Add(nosc_spect_6AD8AD_histo[i],nosc_spect_6AD35B_histo[i],1,-1);
        
        bkgd_spect_6AD8AD_histo[i]->Scale(sclFac);
        bkgd_spect_8AD35B_histo[i]->Add(bkgd_spect_6AD8AD_histo[i],bkgd_spect_6AD35B_histo[i],1,-1);
    }
    
    //-- comparing the total number of events in each spectrum
    double Nevents26, Nevents35;
    double Nbkgd26, Nbkgd35;
    double Nnosc26, Nnosc35;
    for (int l = 0 ; l < nEH ; l++){
        Nevents26 = data_spect_6AD26B_histo[l]->Integral();
        Nevents35 = data_spect_6AD35B_histo[l]->Integral();
        Nbkgd26   = bkgd_spect_6AD26B_histo[l]->Integral();
        Nbkgd35   = bkgd_spect_6AD35B_histo[l]->Integral();
        Nnosc26   = nosc_spect_6AD26B_histo[l]->Integral();
        Nnosc35   = nosc_spect_6AD35B_histo[l]->Integral();
        std::cout << "26Bins Spectrum EH " << l << ":\t" << Nevents26 << std::endl;
        std::cout << "35Bins Spectrum EH " << l << ":\t" << Nevents35 << std::endl;
        std::cout << std::endl;
        std::cout << "26Bins Nosc EH " << l << ":\t" << Nnosc26 << std::endl;
        std::cout << "35Bins Nosc EH " << l << ":\t" << Nnosc35 << std::endl;
        std::cout << std::endl;
        std::cout << "26Bins Bkgd EH " << l << ":\t" << Nbkgd26 << std::endl;
        std::cout << "35Bins Bkgd EH " << l << ":\t" << Nbkgd35 << std::endl;
        std::cout << std::endl;
        Nevents26 = 0.0;
        Nevents35 = 0.0;
        Nnosc26 = 0.0;
        Nnosc35 = 0.0;
        Nbkgd26 = 0.0;
        Nbkgd35 = 0.0;
    }

    for (int i = 0 ; i < nEH ; i++) {
        for (int j = 0 ; j < NB ; j++) {
            double wid = data_spect_6AD35B_histo[i]->GetBinWidth(j+1);
            double con = data_spect_6AD35B_histo[i]->GetBinContent(j+1);
            data_spect_6AD35B_histoPerMeV[i]->SetBinContent(j+1,con/wid);

            wid = nosc_spect_6AD35B_histo[i]->GetBinWidth(j+1);
            con = nosc_spect_6AD35B_histo[i]->GetBinContent(j+1);
            nosc_spect_6AD35B_histoPerMeV[i]->SetBinContent(j+1,con/wid);

            wid = bkgd_spect_6AD35B_histo[i]->GetBinWidth(j+1);
            con = bkgd_spect_6AD35B_histo[i]->GetBinContent(j+1);
            bkgd_spect_6AD35B_histoPerMeV[i]->SetBinContent(j+1,con/wid);

            wid = data_spect_8AD35B_histo[i]->GetBinWidth(j+1);
            con = data_spect_8AD35B_histo[i]->GetBinContent(j+1);
            data_spect_8AD35B_histoPerMeV[i]->SetBinContent(j+1,con/wid);
  
            wid = nosc_spect_8AD35B_histo[i]->GetBinWidth(j+1);
            con = nosc_spect_8AD35B_histo[i]->GetBinContent(j+1);
            nosc_spect_8AD35B_histoPerMeV[i]->SetBinContent(j+1,con/wid);
    
            wid = bkgd_spect_8AD35B_histo[i]->GetBinWidth(j+1);
            con = bkgd_spect_8AD35B_histo[i]->GetBinContent(j+1);
            bkgd_spect_8AD35B_histoPerMeV[i]->SetBinContent(j+1,con/wid);

            wid = data_spect_6AD8AD_histo[i]->GetBinWidth(j+1);
            con = data_spect_6AD8AD_histo[i]->GetBinContent(j+1);
            data_spect_6AD8AD_histoPerMeV[i]->SetBinContent(j+1,con/wid);

            wid = nosc_spect_6AD8AD_histo[i]->GetBinWidth(j+1);
            con = nosc_spect_6AD8AD_histo[i]->GetBinContent(j+1);
            nosc_spect_6AD8AD_histoPerMeV[i]->SetBinContent(j+1,con/wid);

            wid = bkgd_spect_6AD8AD_histo[i]->GetBinWidth(j+1);
            con = bkgd_spect_6AD8AD_histo[i]->GetBinContent(j+1);
            bkgd_spect_6AD8AD_histoPerMeV[i]->SetBinContent(j+1,con/wid);
        }
    }

    //---------------------------------------------------------
    // write to output file
    TFile *fout = new TFile("histos_6AD217Days_8AD1013Days_data.root","recreate");
    fout->cd();
    for (int j = 0 ; j < nEH ; j++) {
        data_spect_6AD8AD_histoPerMeV[j]->Write();
        data_spect_6AD35B_histoPerMeV[j]->Write();
        data_spect_8AD35B_histoPerMeV[j]->Write();
        nosc_spect_6AD8AD_histoPerMeV[j]->Write();
        nosc_spect_6AD35B_histoPerMeV[j]->Write();
        nosc_spect_8AD35B_histoPerMeV[j]->Write();
        bkgd_spect_6AD8AD_histoPerMeV[j]->Write();
        bkgd_spect_6AD35B_histoPerMeV[j]->Write();
        bkgd_spect_8AD35B_histoPerMeV[j]->Write();

        data_spect_6AD8AD_histo[j]->Write();
        data_spect_6AD35B_histo[j]->Write();
        data_spect_8AD35B_histo[j]->Write();
        nosc_spect_6AD8AD_histo[j]->Write();
        nosc_spect_6AD35B_histo[j]->Write();
        nosc_spect_8AD35B_histo[j]->Write();
        bkgd_spect_6AD8AD_histo[j]->Write();
        bkgd_spect_6AD35B_histo[j]->Write();
        bkgd_spect_8AD35B_histo[j]->Write();
    }
    fout->Close();

    //----------------------------------------------------------------------------------------
    TLegend *leg1 = new TLegend(0.6,0.7,0.8,0.87);
    leg1->SetTextFont(ft);
    leg1->SetTextSize(1.1*sz);
    leg1->SetFillColor(0);
    leg1->SetLineColor(0);
    leg1->SetHeader("Data Spectra","C");
    
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

    //-- 35-bins nosc spectra
    TLegend *leg2 = new TLegend(0.6,0.7,0.8,0.87);
    leg2->SetTextFont(ft);
    leg2->SetTextSize(1.1*sz);
    leg2->SetFillColor(0);
    leg2->SetLineColor(0);
    leg2->SetHeader("Nosc Spectra","C");
    
    leg2->AddEntry(nosc_spect_6AD8AD_histoPerMeV[0],"DB 2130 Days - 8AD 35 - Bins","l");
    leg2->AddEntry(nosc_spect_6AD35B_histoPerMeV[0],"DB  217 Days - 6AD 35 - Bins","l");
    leg2->AddEntry(nosc_spect_8AD35B_histoPerMeV[0],"DB 1913 Days - 8AD 35 - Bins","l");

    TCanvas *canv2 = new TCanvas("canv2","35 Bins - Nosc",700,0,700,3*450);
    canv2->Divide(1,3);
    
    canv2->cd(1);
    nosc_spect_6AD8AD_histoPerMeV[0]->Draw("hist");
    nosc_spect_8AD35B_histoPerMeV[0]->Draw("same hist");
    nosc_spect_6AD35B_histoPerMeV[0]->Draw("same hist");
    leg2->Draw();
    gPad->SetTicks(1,1);
    
    canv2->cd(2);
    nosc_spect_6AD8AD_histoPerMeV[1]->Draw("hist");
    nosc_spect_8AD35B_histoPerMeV[1]->Draw("same hist");
    nosc_spect_6AD35B_histoPerMeV[1]->Draw("same hist");
    gPad->SetTicks(1,1);
    
    canv2->cd(3);
    nosc_spect_6AD8AD_histoPerMeV[2]->Draw("hist");
    nosc_spect_8AD35B_histoPerMeV[2]->Draw("same hist");
    nosc_spect_6AD35B_histoPerMeV[2]->Draw("same hist");
    gPad->SetTicks(1,1);

}// end
