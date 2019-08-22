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

    //---------------------------------------------------------

    //-- Our Best fit espectra - BEGIN------------------------------------------------------------//
    double bfOscAD1_array[NB] = {44906.769531, 48609.500000, 65805.437500, 80176.710938, 91900.109375, 100836.804688, 104935.085938, 106860.453125, 105016.460938, 100471.375000, 91899.437500, 84036.726562, 76400.031250, 70017.296875, 62114.125000, 54348.546875, 47022.136719, 39262.890625, 32155.544922, 25414.648438, 19850.259766, 14527.912109, 10496.083008, 7418.156250, 4629.794922, 6937.556641};
    double bfOscAD2_array[NB] = {45736.804688, 49759.835938, 66925.609375, 81274.960938, 92409.273438, 102589.648438, 106960.750000, 108730.562500, 106880.601562, 101899.742188, 94120.929688, 85332.101562, 78009.250000, 71388.085938, 63426.171875, 55402.843750, 48000.160156, 39843.972656, 32606.679688, 25826.085938, 20386.335938, 14810.592773, 10360.733398, 7593.277344, 4732.845703, 7040.208496};
    double bfOscAD3_array[NB] = {39504.968750, 43704.343750, 58559.621094, 71628.046875, 81242.210938, 90091.968750, 94247.351562, 96129.828125, 94338.851562, 90046.421875, 82609.617188, 75743.226562, 69294.718750, 62876.277344, 56105.203125, 50009.050781, 42603.187500, 35465.410156, 28597.164062, 23225.923828, 18378.447266, 13849.699219, 9753.091797, 6858.283203, 4728.520020, 16703.605469};
    double bfOscAD4_array[NB] = {5530.125000, 5579.536133, 7516.986328, 8939.312500, 10372.551758, 11131.632812, 11880.826172, 11853.464844, 11733.295898, 11054.114258, 10237.999023, 9256.074219, 8687.268555, 7926.119629, 7020.613281, 6353.561035, 5412.803711, 4447.867676, 3576.991943, 2935.256592, 2242.291016, 1685.780029, 1201.541748, 765.356995, 521.795105, 600.749268};
    double bfOscAD5_array[NB] = {5385.541504, 5662.872559, 7472.637207, 8968.308594, 10213.898438, 11181.324219, 11839.885742, 12069.756836, 11902.598633, 11362.651367, 10490.458984, 9491.621094, 8689.859375, 7872.372559, 7085.985352, 6370.272461, 5317.650879, 4556.246582, 3609.330566, 2985.270752, 2358.871582, 1718.446655, 1135.386475, 820.684875, 543.249939, 644.424927};
    double bfOscAD6_array[NB] = {5499.355957, 5702.975098, 7467.163086, 9126.190430, 10163.433594, 10965.411133, 11581.669922, 11876.322266, 11629.682617, 11096.867188, 10183.879883, 9227.711914, 8549.845703, 7726.782715, 7113.473145, 6094.327148, 5367.018066, 4369.830078, 3675.293701, 2953.593506, 2348.409912, 1608.205322, 1182.344971, 762.550171, 549.489807, 615.287292};
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
    double NoscTot[nAD] = {1496049.854492, 1522048.063965, 1356295.039551, 168463.915527, 169749.608154, 167437.114624};
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
    double IBD_candidates[nAD]  = { 101290,  102519,   92912,   13964,   13894,   13731};
    double emuem[nAD]           = { 0.7957,  0.7927,  0.8282,  0.9577,  0.9568,  0.9566};
    double daqTime[nAD]         = {191.001, 191.001, 189.645, 189.779, 189.779, 189.779};
    
    std::cout << std::endl;
    std::cout << "EH1 data/MeV Integral = " << data_spect_histoPerMeV[0]->Integral() << std::endl;
    std::cout << "EH2 data/MeV Integral = " << data_spect_histoPerMeV[1]->Integral() << std::endl;
    std::cout << "EH3 data/MeV Integral = " << data_spect_histoPerMeV[2]->Integral() << std::endl;
    std::cout << std::endl;
    std::cout << "EH1 data Integral = " << data_spect_histo[0]->Integral() << std::endl;
    std::cout << "EH2 data Integral = " << data_spect_histo[1]->Integral() << std::endl;
    std::cout << "EH3 data Integral = " << data_spect_histo[2]->Integral() << std::endl;
    std::cout << std::endl;
    data_spect_histo[0]-> Scale((IBD_candidates[0]+IBD_candidates[1])/data_spect_histo[0]->Integral());
    data_spect_histo[1]-> Scale((IBD_candidates[2])/data_spect_histo[1]->Integral());
    data_spect_histo[2]-> Scale((IBD_candidates[3]+IBD_candidates[4]+IBD_candidates[5])/data_spect_histo[2]->Integral());
    std::cout << "Mod EH1 data Integral = " << data_spect_histo[0]->Integral() << std::endl;
    std::cout << "Mod EH2 data Integral = " << data_spect_histo[1]->Integral() << std::endl;
    std::cout << "Mod EH3 data Integral = " << data_spect_histo[2]->Integral() << std::endl;
    std::cout << std::endl;

    
    double Dt = 0.0;
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
        NoO = NoOscAD1_array[iBIN]*(1.0*noOsc_IBDrate_perday[0]/noNoscTot[0])*emuem[0]*daqTime[0];
        NoOscAD1_histo->SetBinContent(iBIN+1,NoO);
        
        NoO = NoOscAD2_array[iBIN]*(1.0*noOsc_IBDrate_perday[1]/noNoscTot[1])*emuem[1]*daqTime[1];
        NoOscAD2_histo->SetBinContent(iBIN+1,NoO);
        
        NoO = NoOscAD3_array[iBIN]*(1.0*noOsc_IBDrate_perday[2]/noNoscTot[2])*emuem[2]*daqTime[2];
        NoOscAD3_histo->SetBinContent(iBIN+1,NoO);
        
        NoO = NoOscAD4_array[iBIN]*(1.0*noOsc_IBDrate_perday[3]/noNoscTot[3])*emuem[3]*daqTime[3];
        NoOscAD4_histo->SetBinContent(iBIN+1,NoO);
        
        NoO = NoOscAD5_array[iBIN]*(1.0*noOsc_IBDrate_perday[4]/noNoscTot[4])*emuem[4]*daqTime[4];
        NoOscAD5_histo->SetBinContent(iBIN+1,NoO);
        
        NoO = NoOscAD6_array[iBIN]*(1.0*noOsc_IBDrate_perday[5]/noNoscTot[5])*emuem[5]*daqTime[5];
        NoOscAD6_histo->SetBinContent(iBIN+1,NoO);
    }
    
    TH1F *our_BFit_histo[nEH];
    TH1F *our_noOsc_histo[nEH];
    for (int i = 0 ; i < nEH ; i++){
        our_BFit_histo[i] = new TH1F(Form("our_BFit_histo_%d",i),"",NB,xbins);
        our_BFit_histo[i]->SetLineWidth(2);
        our_BFit_histo[i]->SetLineColor(2);
        our_BFit_histo[i]->SetLineStyle(1);
        our_noOsc_histo[i] = new TH1F(Form("our_noOsc_histo_%d",i),"",NB,xbins);
        our_noOsc_histo[i]->SetLineWidth(1);
        our_noOsc_histo[i]->SetLineColor(4);
        our_noOsc_histo[i]->SetLineStyle(1);
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
            
            cont = data_spect_histo[j]->GetBinContent(i+1);
            data_spect_histo[j]->SetBinContent(i+1,(cont/wid));
        }
    }
    std::cout << std::endl;
    std::cout << "After dividing by the bin width:" << std::endl;
    std::cout << "Mod EH1 data Integral = " << data_spect_histo[0]->Integral() << std::endl;
    std::cout << "Mod EH2 data Integral = " << data_spect_histo[1]->Integral() << std::endl;
    std::cout << "Mod EH3 data Integral = " << data_spect_histo[2]->Integral() << std::endl;
    std::cout << std::endl;

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
    
    leg1->AddEntry(data_spect_histo[0],"Daya Bay Data","p");
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
    //data_spect_histoPerMeV[0]->Draw("P same");
    data_spect_histo[0]->Draw("P same");
    our_BFit_histo[0]->Draw("same hist");
    //nosc_spect_histoPerMeV[0]->Draw("same");
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
    //data_spect_histoPerMeV[1]->Draw("P same");
    data_spect_histo[1]->Draw("P same");
    our_BFit_histo[1]->Draw("same hist");
    //nosc_spect_histoPerMeV[1]->Draw("same");
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
    //data_spect_histoPerMeV[2]->Draw("P same");
    data_spect_histo[2]->Draw("P same");
    our_BFit_histo[2]->Draw("same hist");
    //nosc_spect_histoPerMeV[2]->Draw("same");
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

