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
    TGraph *DB_specEH2_BFperMev    = new TGraph("files_data/EH2/DB_specEH2_BF.csv","%lg,%lg","");
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
        BFit_spect_histoPerMeV[i]->SetLineColor(kOrange);
        
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
    double bfOscAD1_array[NB] = {44830.445312, 48525.503906, 65690.703125, 80047.203125, 91762.281250, 100622.945312, 104772.742188, 106734.234375, 104915.367188, 100403.921875, 91802.539062, 83923.070312, 76338.742188, 69961.921875, 62079.875000, 54320.734375, 46995.878906, 39243.304688, 32131.652344, 25395.611328, 19840.394531, 14519.090820, 10490.928711, 7414.633789, 4627.903320, 6935.844238};
    double bfOscAD2_array[NB] = {45659.839844, 49673.339844, 66808.632812, 81154.804688, 92266.351562, 102364.796875, 106803.710938, 108614.546875, 106782.851562, 101830.375000, 94013.593750, 85225.710938, 77949.195312, 71334.664062, 63393.441406, 55375.582031, 47973.539062, 39824.195312, 32579.154297, 25808.324219, 20375.724609, 14801.668945, 10355.750000, 7589.788086, 4731.012695, 7038.440430};
    double bfOscAD3_array[NB] = {39366.269531, 43607.937500, 58450.476562, 71456.437500, 81102.109375, 89890.742188, 94097.867188, 96008.195312, 94196.140625, 89930.945312, 82487.898438, 75641.898438, 69216.898438, 62831.027344, 56060.109375, 49975.859375, 42563.449219, 35437.519531, 28585.085938, 23200.431641, 18366.640625, 13834.255859, 9748.004883, 6854.117188, 4725.411133, 16698.916016};
    double bfOscAD4_array[NB] = {5417.008301, 5620.425781, 7560.541504, 9006.519531, 10413.240234, 11133.336914, 11862.096680, 11856.900391, 11727.947266, 11029.192383, 10212.175781, 9217.815430, 8679.847656, 7899.536621, 7008.798340, 6327.159668, 5405.603027, 4439.207520, 3566.301514, 2929.092529, 2234.992676, 1687.221802, 1195.978027, 764.104980, 519.395874, 593.776367};
    double bfOscAD5_array[NB] = {5277.126465, 5688.483887, 7524.675293, 9027.935547, 10249.531250, 11208.653320, 11825.549805, 12046.866211, 11896.350586, 11347.019531, 10454.595703, 9458.767578, 8673.900391, 7860.806641, 7074.909180, 6348.411133, 5303.349609, 4550.911133, 3598.490234, 2976.067383, 2356.879639, 1711.177246, 1135.795166, 814.496704, 544.642151, 636.395935};
    double bfOscAD6_array[NB] = {5385.978027, 5745.570801, 7500.989258, 9171.907227, 10226.217773, 10975.889648, 11578.887695, 11874.760742, 11602.110352, 11081.333008, 10175.903320, 9183.279297, 8536.461914, 7706.292480, 7089.224609, 6090.473633, 5347.640137, 4360.496582, 3666.051270, 2946.402344, 2342.623535, 1607.964844, 1177.772827, 764.228638, 544.112671, 607.320618};
    //-- Our No-Osc. espectra - //
    double NoOscAD1_array[NB] = {46994.0, 50441.0, 67978.0, 82506.0, 94304.0, 103122.0, 107005.0, 108732.0, 106722.0, 101997.0, 93167.0, 85049.0, 77234.0, 70704.0, 62676.0, 54816.0, 47401.0, 39562.0, 32372.0, 25570.0, 19963.0, 14608.0, 10549.0, 7453.0, 4651.0, 6957.0};
    double NoOscAD2_array[NB] = {47834.0, 51600.0, 69122.0, 83595.0, 94786.0, 104867.0, 109016.0, 110610.0, 108594.0, 103420.0, 95404.0, 86346.0, 78853.0, 72082.0, 63993.0, 55871.0, 48376.0, 40144.0, 32818.0, 25981.0, 20502.0, 14889.0, 10412.0, 7629.0, 4754.0, 7060.0};
    double NoOscAD3_array[NB] = {41843.0, 45699.0, 60877.0, 74047.0, 83659.0, 92423.0, 96405.0, 98112.0, 96088.0, 91539.0, 83808.0, 76736.0, 70117.0, 63577.0, 56685.0, 50496.0, 42967.0, 35743.0, 28819.0, 23375.0, 18496.0, 13925.0, 9808.0, 6894.0, 4751.0, 16754.0};
    double NoOscAD4_array[NB] = {5537.0, 5880.0, 8076.0, 9771.0, 11405.0, 12250.0, 13067.0, 13043.0, 12864.0, 12051.0, 11111.0, 9983.0, 9360.0, 8480.0, 7493.0, 6737.0, 5735.0, 4692.0, 3758.0, 3077.0, 2341.0, 1763.0, 1247.0, 795.0, 539.0, 608.0};
    double NoOscAD5_array[NB] = {5390.0, 5950.0, 8034.0, 9793.0, 11222.0, 12334.0, 13027.0, 13253.0, 13049.0, 12400.0, 11375.0, 10245.0, 9354.0, 8441.0, 7563.0, 6761.0, 5626.0, 4812.0, 3792.0, 3127.0, 2469.0, 1788.0, 1184.0, 847.0, 565.0, 652.0};
    double NoOscAD6_array[NB] = {5501.0, 6007.0, 8006.0, 9941.0, 11197.0, 12076.0, 12755.0, 13064.0, 12728.0, 12111.0, 11074.0, 9949.0, 9207.0, 8275.0, 7581.0, 6487.0, 5675.0, 4610.0, 3865.0, 3096.0, 2455.0, 1680.0, 1228.0, 795.0, 565.0, 622.0};

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
    double NoscTot[nAD] = {1494327.473145, 1520329.035156, 1354334.644531, 168308.216797, 169591.787720, 167289.893250};
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
    for (int j = 0 ; j < nEH ; j++) {
        ratio_histo[j] = (TH1F*)nosc_spect_histoPerMeV[j]->Clone(Form("ratio_histo_%d",j));
        ratio_histo[j]->Divide(our_noOsc_histo[j]);
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
    //BFit_spect_histoPerMeV[0]->Draw("same hist");
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
    //BFit_spect_histoPerMeV[1]->Draw("same hist");
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
    //BFit_spect_histoPerMeV[2]->Draw("same hist");
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
    TH2F *frame_ratio = new TH2F("frame_ratio","",NB,lo,hi,10,0.970,1.030);
    frame_ratio->GetXaxis()->SetTitle("Prompt Reconstructed Energy (MeV)");
    frame_ratio->GetXaxis()->SetTitleFont(ft);
    frame_ratio->GetXaxis()->SetTitleOffset(0.9);
    frame_ratio->GetXaxis()->SetTitleSize(1.4*sz);
    frame_ratio->GetXaxis()->SetLabelSize(1.4*sz);
    frame_ratio->GetXaxis()->SetLabelFont(ft);
    frame_ratio->GetYaxis()->SetTitle("(NoOsc DB)/(NoOsc Us)");
    frame_ratio->GetYaxis()->SetTitleFont(ft);
    frame_ratio->GetYaxis()->SetTitleOffset(0.7);
    frame_ratio->GetYaxis()->SetTitleSize(1.4*sz);
    frame_ratio->GetYaxis()->SetLabelSize(1.4*sz);
    frame_ratio->GetYaxis()->SetLabelFont(ft);

    TCanvas *ra = new TCanvas("ra", "ratios", 700, 600);
    TGaxis::SetMaxDigits(3);
    TLine* line1 = new TLine(0.7,1,12,1);
    line1->SetLineColor(kGray);

    TPad *padr1 = new TPad("padr1", "padr1", 0, 2./3., 1, 1.0);
    padr1->SetBottomMargin(0); // Upper and lower plot are joined
    padr1->Draw();             // Draw the upper pad: pad1
    padr1->cd();               // pad1 becomes the current pad
    frame_ratio->Draw();
    ratio_histo[0]->Draw("same h");
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
    frame_ratio->Draw();
    ratio_histo[1]->Draw("same h");
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
    frame_ratio->Draw();
    ratio_histo[2]->Draw("same h");
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

