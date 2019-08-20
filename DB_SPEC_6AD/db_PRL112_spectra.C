// csv files are digitalizations of data plots in the article, by using Engauge software.
// plot from Daya Bay Coll. PRL 112, 061801 (2014) - arXiv:13106732

void db_PRL112_spectra()
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
    //----------------------------------------------------------------------------------
    //get data from digitized files - EH1
    TGraph *DB_specEH1_DataperMev  = new TGraph("files_data/EH1/DB_specEH1_Data.csv","%lg,%lg","");
    TGraph *DB_specEH1_NoOscperMev = new TGraph("files_data/EH1/DB_specEH1_NoOsc.csv","%lg,%lg","");
    TGraph *DB_specEH1_BFperMev    = new TGraph("files_data/EH1/DB_specEH1_BF.csv","%lg,%lg","");
    TGraph *DB_specEH1_BGperMev    = new TGraph("files_data/EH1/DB_specEH1_BG.csv","%lg,%lg","");
    //----------------------------------------------------------------------------------

    //get data from digitized files - EH2
    TGraph *DB_specEH2_DataperMev  = new TGraph("files_data/EH2/DB_specEH2_DataNew.csv","%lg,%lg","");
    TGraph *DB_specEH2_NoOscperMev = new TGraph("files_data/EH2/DB_specEH2_NoOscNew.csv","%lg,%lg","");
    TGraph *DB_specEH2_BFperMev    = new TGraph("files_data/EH2/DB_specEH2_BFNew.csv","%lg,%lg","");
    TGraph *DB_specEH2_BGperMev    = new TGraph("files_data/EH2/DB_specEH2_BGNew.csv","%lg,%lg","");
    //---------------------------------------------------------------------------------------------------
    //get data from digitized files - EH3
    TGraph *DB_specEH3_DataperMev  = new TGraph("files_data/EH3/DB_specEH3_Data.csv","%lg,%lg","");
    TGraph *DB_specEH3_NoOscperMev = new TGraph("files_data/EH3/DB_specEH3_NoOsc.csv","%lg,%lg","");
    TGraph *DB_specEH3_BFperMev    = new TGraph("files_data/EH3/DB_specEH3_BF.csv","%lg,%lg","");
    TGraph *DB_specEH3_BGperMev    = new TGraph("files_data/EH3/DB_specEH3_BGNew.csv","%lg,%lg","");
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
    TH1F *nosc_spect_histoPerMeV_v2[nEH];
    TH1F *BFit_spect_histoPerMeV[nEH];
    TH1F *bkgd_spect_histoPerMeV[nEH];
    //-- Events
    TH1F *data_spect_histo[nEH];
    TH1F *nosc_spect_histo[nEH];
    TH1F *BFit_spect_histo[nEH];
    TH1F *bkgd_spect_histo[nEH];
    
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
        BFit_spect_histoPerMeV[i]->SetLineWidth(1);
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
    
    //---------------------------------------------------------

    //-- Our Best fit espectra - BEGIN------------------------------------------------------------//
    double bfOscAD1_array[NB] = {44923.179688, 48640.695312, 65846.484375, 80217.257812, 91964.742188, 100931.882812, 105007.835938, 106902.523438, 105053.367188, 100511.992188, 91955.726562, 84075.023438, 76425.867188, 70040.695312, 62126.417969, 54363.859375, 47031.996094, 39275.808594, 32170.667969, 25420.900391, 19855.759766, 14531.860352, 10498.394531, 7419.506836, 4630.621094, 6938.362305};
    double bfOscAD2_array[NB] = {45753.960938, 49790.562500, 66968.796875, 81314.257812, 92476.765625, 102681.539062, 107033.632812, 108769.773438, 106918.703125, 101940.265625, 94181.234375, 85365.640625, 78034.429688, 71411.421875, 63438.160156, 55418.308594, 48009.636719, 39859.222656, 32620.236328, 25832.160156, 20392.111328, 14814.250000, 10362.870117, 7594.559082, 4733.692383, 7041.067383};
    double bfOscAD3_array[NB] = {39544.410156, 43747.273438, 58613.476562, 71693.867188, 81308.460938, 90179.468750, 94306.164062, 96184.609375, 94410.101562, 90107.164062, 82660.015625, 75791.046875, 69316.937500, 62899.125000, 56124.015625, 50032.566406, 42616.984375, 35477.511719, 28607.634766, 23233.873047, 18383.054688, 13853.395508, 9755.124023, 6859.923828, 4729.486328, 16705.773438};
    double bfOscAD4_array[NB] = {5512.877441, 5539.559082, 7469.577148, 8896.383789, 10338.937500, 11110.057617, 11870.690430, 11852.689453, 11739.267578, 11064.242188, 10250.256836, 9268.994141, 8700.360352, 7938.599609, 7031.878418, 6363.790527, 5421.415527, 4454.815918, 3582.439941, 2939.592773, 2245.497559, 1688.112305, 1203.147827, 766.343445, 522.436096, 601.256531};
    double bfOscAD5_array[NB] = {5368.732422, 5622.143066, 7425.176270, 8925.012695, 10180.092773, 11159.664062, 11829.381836, 12068.737305, 11908.782227, 11373.059570, 10503.091797, 9504.885742, 8703.018555, 7884.785645, 7097.388672, 6380.542969, 5326.126465, 4563.375977, 3614.833984, 2989.694824, 2362.248047, 1720.823364, 1136.899780, 821.734741, 543.913086, 644.981018};
    double bfOscAD6_array[NB] = {5482.606934, 5661.974609, 7419.580078, 9081.479492, 10129.814453, 10943.884766, 11571.001953, 11875.160156, 11635.333008, 11106.860352, 10196.017578, 9240.456055, 8562.715820, 7738.964355, 7124.865234, 6104.131348, 5375.571777, 4376.656738, 3680.919434, 2957.970215, 2351.785156, 1610.425659, 1183.924683, 763.532043, 550.173157, 615.812134};
    //-- Our No-Osc. espectra - //
    double NoOscAD1_array[NB] = {46994.0, 50441.0, 67978.0, 82506.0, 94304.0, 103122.0, 107005.0, 108732.0, 106722.0, 101997.0, 93167.0, 85049.0, 77234.0, 70704.0, 62676.0, 54816.0, 47401.0, 39562.0, 32372.0, 25570.0, 19963.0, 14608.0, 10549.0, 7453.0, 4651.0, 6957.0};
    double NoOscAD2_array[NB] = {47834.0, 51600.0, 69122.0, 83595.0, 94786.0, 104867.0, 109016.0, 110610.0, 108594.0, 103420.0, 95404.0, 86346.0, 78853.0, 72082.0, 63993.0, 55871.0, 48376.0, 40144.0, 32818.0, 25981.0, 20502.0, 14889.0, 10412.0, 7629.0, 4754.0, 7060.0};
    double NoOscAD3_array[NB] = {41865.0, 45682.0, 60841.0, 74055.0, 83655.0, 92428.0, 96394.0, 98111.0, 96097.0, 91527.0, 83816.0, 76738.0, 70120.0, 63581.0, 56686.0, 50490.0, 42971.0, 35752.0, 28811.0, 23384.0, 18497.0, 13933.0, 9808.0, 6895.0, 4752.0, 16754.0};
    double NoOscAD4_array[NB] = {5664.0, 5890.0, 8093.0, 9752.0, 11397.0, 12262.0, 13080.0, 13016.0, 12836.0, 12039.0, 11099.0, 9986.0, 9331.0, 8475.0, 7476.0, 6739.0, 5721.0, 4684.0, 3756.0, 3073.0, 2341.0, 1756.0, 1249.0, 794.0, 540.0, 614.0};
    double NoOscAD5_array[NB] = {5512.0, 5977.0, 8042.0, 9783.0, 11220.0, 12318.0, 13036.0, 13255.0, 13021.0, 12377.0, 11373.0, 10241.0, 9334.0, 8420.0, 7545.0, 6758.0, 5620.0, 4800.0, 3790.0, 3126.0, 2463.0, 1790.0, 1180.0, 851.0, 562.0, 659.0};
    double NoOscAD6_array[NB] = {5628.0, 6016.0, 8033.0, 9949.0, 11165.0, 12079.0, 12752.0, 13043.0, 12725.0, 12089.0, 11043.0, 9959.0, 9185.0, 8264.0, 7577.0, 6466.0, 5674.0, 4603.0, 3861.0, 3093.0, 2453.0, 1675.0, 1229.0, 791.0, 569.0, 629.0};

    double noOsc_area = 0.0;
    //Clon no-oscillated spectra
    for (int i = 0 ; i < nEH ; i++) {
        nosc_spect_histoPerMeV_v2[i] = (TH1F*)nosc_spect_histoPerMeV[i]->Clone(Form("nosc_spect_histoPerMeV_v2_%d",i));
        noOsc_area = nosc_spect_histoPerMeV_v2[i]->Integral();
        nosc_spect_histoPerMeV_v2[i]->Scale(1.0/noOsc_area);
        nosc_spect_histoPerMeV_v2[i]->SetLineStyle(2);
    }
    
    std::cout << "area = " << nosc_spect_histoPerMeV_v2[0]->Integral() << std::endl;
    std::cout << "area = " << nosc_spect_histoPerMeV_v2[1]->Integral() << std::endl;
    std::cout << "area = " << nosc_spect_histoPerMeV_v2[2]->Integral() << std::endl;

    const int nAD = 6;
    double totalBgd[nAD][2]     = { {13.20,0.98},{13.01,0.98},{ 9.57,0.71},{ 3.52,0.14},{ 3.48,0.14},{3.43,0.14} };
    double NoscTot[nAD] = {1496761.428711, 1522757.259277, 1357141.464844, 168373.220032, 169659.126892, 167341.617188};
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
        cout << "noOscIBD_rates_perday " << i << ": " << noOsc_IBDrate_perday[i] << endl;
    }//for
    double emuem[nAD]           = {0.7957,0.7927,0.8282,0.9577,0.9568,0.9566};
    double daqTime[nAD]         = {191.001,191.001,189.645,189.779,189.779,189.779};
    
    double Td = 0.0;
    TH1F *bfOscAD1_histo = new TH1F("bfOscAD1_histo","",NB,xbins);
    TH1F *bfOscAD2_histo = new TH1F("bfOscAD2_histo","",NB,xbins);
    TH1F *bfOscAD3_histo = new TH1F("bfOscAD3_histo","",NB,xbins);
    TH1F *bfOscAD4_histo = new TH1F("bfOscAD4_histo","",NB,xbins);
    TH1F *bfOscAD5_histo = new TH1F("bfOscAD5_histo","",NB,xbins);
    TH1F *bfOscAD6_histo = new TH1F("bfOscAD6_histo","",NB,xbins);
    double NoO = 0.0;
    TH1F *NoOscAD1_histo = new TH1F("NoOscAD1_histo","",NB,xbins);
    TH1F *NoOscAD2_histo = new TH1F("NoOscAD2_histo","",NB,xbins);
    TH1F *NoOscAD3_histo = new TH1F("NoOscAD3_histo","",NB,xbins);
    TH1F *NoOscAD4_histo = new TH1F("NoOscAD4_histo","",NB,xbins);
    TH1F *NoOscAD5_histo = new TH1F("NoOscAD5_histo","",NB,xbins);
    TH1F *NoOscAD6_histo = new TH1F("NoOscAD6_histo","",NB,xbins);
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

        //Noosc
        Td = NoOscAD1_array[iBIN]*(1.0*noOsc_IBDrate_perday[0]/noNoscTot[0])*emuem[0]*daqTime[0];
        NoOscAD1_histo->SetBinContent(iBIN+1,Td);
        
        Td = NoOscAD2_array[iBIN]*(1.0*noOsc_IBDrate_perday[1]/noNoscTot[1])*emuem[1]*daqTime[1];
        NoOscAD2_histo->SetBinContent(iBIN+1,Td);
        
        Td = NoOscAD3_array[iBIN]*(1.0*noOsc_IBDrate_perday[2]/noNoscTot[2])*emuem[2]*daqTime[2];
        NoOscAD3_histo->SetBinContent(iBIN+1,Td);
        
        Td = NoOscAD4_array[iBIN]*(1.0*noOsc_IBDrate_perday[3]/noNoscTot[3])*emuem[3]*daqTime[3];
        NoOscAD4_histo->SetBinContent(iBIN+1,Td);
        
        Td = NoOscAD5_array[iBIN]*(1.0*noOsc_IBDrate_perday[4]/noNoscTot[4])*emuem[4]*daqTime[4];
        NoOscAD5_histo->SetBinContent(iBIN+1,Td);
        
        Td = NoOscAD6_array[iBIN]*(1.0*noOsc_IBDrate_perday[5]/noNoscTot[5])*emuem[5]*daqTime[5];
        NoOscAD6_histo->SetBinContent(iBIN+1,Td);

    }
    
    TH1F *our_BFit_histo[nEH];
    TH1F *our_noOsc_histo[nEH];
    for (int i = 0 ; i < nEH ; i++){
        our_BFit_histo[i] = new TH1F(Form("our_BFit_histo_%d",i),"",NB,xbins);
        our_BFit_histo[i]->SetLineWidth(2);
        our_BFit_histo[i]->SetLineColor(2);
        our_BFit_histo[i]->SetLineStyle(2);
        our_noOsc_histo[i] = new TH1F(Form("our_noOsc_histo_%d",i),"",NB,xbins);
        our_noOsc_histo[i]->SetLineWidth(2);
        our_noOsc_histo[i]->SetLineColor(4);
        our_noOsc_histo[i]->SetLineStyle(2);
    }

    TH1F *temp_histo = new TH1F("temp_histo","",NB,xbins);
    TH1F *zero_histo = new TH1F("zero_histo","",NB,xbins);
    for (int k = 0 ; k < NB ; k++) {
        zero_histo->SetBinContent(k+1,0.0);
    }

    our_BFit_histo[0]->Add(bfOscAD1_histo,bfOscAD2_histo);
    our_BFit_histo[1]->Add(bfOscAD3_histo,zero_histo);
    temp_histo->Add(bfOscAD4_histo,bfOscAD5_histo);
    our_BFit_histo[2]->Add(temp_histo,bfOscAD6_histo);

    our_noOsc_histo[0]->Add(NoOscAD1_histo,NoOscAD2_histo);
    our_noOsc_histo[1]->Add(NoOscAD3_histo,zero_histo);
    temp_histo->Add(NoOscAD4_histo,NoOscAD5_histo);
    our_noOsc_histo[2]->Add(temp_histo,NoOscAD6_histo);

    for (int j = 0 ; j < nEH ; j++) {
        for (int i = 0 ; i < NB ; i++) {
            double bgnd = bkgd_spect_histoPerMeV[j]->GetBinContent(i+1);
            double wid  = our_BFit_histo[j]->GetBinWidth(i+1);
            double cont = our_BFit_histo[j]->GetBinContent(i+1);
            our_BFit_histo[j]->SetBinContent(i+1,(cont/wid) + bgnd);

            cont = our_noOsc_histo[j]->GetBinContent(i+1);
            our_noOsc_histo[j]->SetBinContent(i+1,(cont/wid) + bgnd);
        }
    }
    
    TH1F *ratio_histo[nEH];
    TH1F *ratioBF_histo[nEH];
    for (int j = 0 ; j < nEH ; j++) {
        ratio_histo[j] = (TH1F*)nosc_spect_histoPerMeV[j]->Clone(Form("ratio_histo_%d",j));
        ratio_histo[j]->Divide(our_noOsc_histo[j]);
        ratio_histo[j]->SetLineColor(4);
        ratioBF_histo[j] = (TH1F*)BFit_spect_histoPerMeV[j]->Clone(Form("ratioBF_histo_%d",j));
        ratioBF_histo[j]->Divide(our_BFit_histo[j]);
        ratioBF_histo[j]->SetLineColor(2);
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
    BFit_spect_histoPerMeV[0]->Draw("same hist");
    data_spect_histoPerMeV[0]->Draw("P same");
    our_BFit_histo[0]->Draw("same hist");
    nosc_spect_histoPerMeV[0]->Draw("same");
    our_noOsc_histo[0]->Draw("same");//check
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
    data_spect_histoPerMeV[1]->Draw("P same");
    our_BFit_histo[1]->Draw("same hist");
    nosc_spect_histoPerMeV[1]->Draw("same");
    our_noOsc_histo[1]->Draw("same");//check
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
    BFit_spect_histoPerMeV[2]->Draw("same hist");
    data_spect_histoPerMeV[2]->Draw("P same");
    our_BFit_histo[2]->Draw("same hist");
    nosc_spect_histoPerMeV[2]->Draw("same");
    our_noOsc_histo[2]->Draw("same");//check
    bkgd_spect_histoPerMeV[2]->Draw("same");
    gPad->SetTicks(1,1);
    lat->DrawLatex(0.7,0.3,"EH3");
    
    ca->Print("files_plots/DB_spect-PRL112.pdf");
    ca->Print("files_plots/DB_spect-PRL112.eps");

    //-----
    TH2F *frame_ratio1 = new TH2F("frame_ratio1","",NB,lo,hi,10,0.90,1.10);
    frame_ratio1->GetYaxis()->SetTitle("DB/OUR");
    frame_ratio1->GetYaxis()->SetTitleFont(ft);
    frame_ratio1->GetYaxis()->SetTitleOffset(0.7);
    frame_ratio1->GetYaxis()->SetTitleSize(1.4*sz);
    frame_ratio1->GetYaxis()->SetLabelSize(1.4*sz);
    frame_ratio1->GetYaxis()->SetLabelFont(ft);
    TH2F *frame_ratio2 = new TH2F("frame_ratio2","",NB,lo,hi,10,0.90,1.030);
    frame_ratio2->GetYaxis()->SetTitle("DB/OUR");
    frame_ratio2->GetYaxis()->SetTitleFont(ft);
    frame_ratio2->GetYaxis()->SetTitleOffset(0.7);
    frame_ratio2->GetYaxis()->SetTitleSize(1.4*sz);
    frame_ratio2->GetYaxis()->SetLabelSize(1.4*sz);
    frame_ratio2->GetYaxis()->SetLabelFont(ft);
    TH2F *frame_ratio3 = new TH2F("frame_ratio3","",NB,lo,hi,10,0.960,1.040);
    frame_ratio3->GetXaxis()->SetTitle("Prompt Reconstructed Energy (MeV)");
    frame_ratio3->GetXaxis()->SetTitleFont(ft);
    frame_ratio3->GetXaxis()->SetTitleOffset(0.9);
    frame_ratio3->GetXaxis()->SetTitleSize(1.4*sz);
    frame_ratio3->GetXaxis()->SetLabelSize(1.4*sz);
    frame_ratio3->GetXaxis()->SetLabelFont(ft);
    frame_ratio3->GetYaxis()->SetTitle("DB/OUR");
    frame_ratio3->GetYaxis()->SetTitleFont(ft);
    frame_ratio3->GetYaxis()->SetTitleOffset(0.7);
    frame_ratio3->GetYaxis()->SetTitleSize(1.4*sz);
    frame_ratio3->GetYaxis()->SetLabelSize(1.4*sz);
    frame_ratio3->GetYaxis()->SetLabelFont(ft);

    TCanvas *ra = new TCanvas("ra", "ratios", 700, 600);
    TGaxis::SetMaxDigits(3);
    TLine* line1 = new TLine(0.7,1,12,1);
    line1->SetLineColor(kGray);

    TPad *padr1 = new TPad("padr1", "padr1", 0, 2./3., 1, 1.0);
    padr1->SetBottomMargin(0); // Upper and lower plot are joined
    padr1->Draw();             // Draw the upper pad: pad1
    padr1->cd();               // pad1 becomes the current pad
    frame_ratio1->Draw();
    ratio_histo[0]->Draw("same h");
    ratioBF_histo[0]->Draw("same h");
    gPad->SetTicks(1,1);
    line1->Draw();
    lat->DrawLatex(0.7,0.3,"EH1");
    
    // lower plot will be in pad
    ra->cd();          // Go back to the main canvas before defining pad2
    TPad *padr2 = new TPad("padr2", "padr2", 0, 1./3., 1, 2./3.);
    padr2->SetTopMargin(0);
    padr2->SetBottomMargin(0);
    padr2->Draw();
    padr2->cd();       // pad2 becomes the current pad
    frame_ratio2->Draw();
    ratio_histo[1]->Draw("same h");
    ratioBF_histo[1]->Draw("same h");
    gPad->SetTicks(1,1);
    line1->Draw();
    lat->DrawLatex(0.7,0.3,"EH2");
    
    // lower plot will be in pad
    ra->cd();          // Go back to the main canvas before defining pad2
    TPad *padr3 = new TPad("padr3", "padr3", 0, 0.0, 1, 1./3.);
    padr3->SetTopMargin(0);
    padr3->SetBottomMargin(0.11);
    padr3->Draw();
    padr3->cd();       // pad2 becomes the current pad
    frame_ratio3->Draw();
    ratio_histo[2]->Draw("same h");
    ratioBF_histo[2]->Draw("same h");
    gPad->SetTicks(1,1);
    line1->Draw();
    lat->DrawLatex(0.7,0.3,"EH3");

    ra->Print("files_plots/DB_noOsc_Ratios-PRL112.pdf");

    // write to output file
    TFile *fout = new TFile("PRL112_217days_spectra.root","recreate");
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
        //our BF spectr (Events / MeV)
        our_BFit_histo[i]->Write();
    }
    
    fout->Close();
    
}// end

