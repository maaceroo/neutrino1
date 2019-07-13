//--------------------------------------------------------------------------------//
//--  db_CovaMatrix_6AD_6x26bins.C - By M.A. Acero O. & A.A. Aguilar-A. - 2018-01-27  --//
// This macro creates a plot of the 6x6 DB-AD Covariance matrix. This matrix was  //
// constructed as follows:                                                        //
//  1. Rebin the full correlation matrix from 37 bins as seen in                  //
//     "Daya Bay Oscillaton Analysis [Pure covariance approach]" presentation by  //
//     Henoch Wong (2016)                                                         //
//     to 26 bins as in the paper                                                 //
//     An et al. PRL112 061801 (2014) "Spectral Measurements of electron          //
//     antineutrino oscillation, amplitud and frecuency at Daya-Bay".             //
//  2. Rebin the full systematic errors from 25 (0 - 8 MeV) bins from             //
//     "Measurement of the Reactor Antineutrino Flux and Spectrum at Daya Bay",   //
//     arXiv:1508.04233,                                                          //
//     to 26 bin (0 - 12 MeV). The last bin has been put by hand.                 //
//  3. Calculate Covariance matrix as:                                            //
//     cov_ij = sigma_i*sigma_j*rho_ij                                            //
//  4. Replicate the Covariance Matrix for all the detector and make an array of  //
//  6x6 blocks.                                                                   //
//--------------------------------------------------------------------------------//

//Constants definition
//const int NB = 26;      //Correlation matrix Number of bins
//------------------------------------------------------
/*
void set_plot_style()
{//This function sets appropiate colors to the bar used for "COLZ" option
    const Int_t NRGBs = 5;
    const Int_t NCont = 255;
    
    Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
    Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
    Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    gStyle->SetNumberContours(NCont);
}
*/

void db_CovaMatrix_6AD_6x26bins()
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

    //---------------------------------------------------------------------------------------------------
    TFile *inFile = new TFile("./db_CovaMatrix_6AD_26bins.root");
    TH2F *covaMat_histo_1x1Det = (TH2F*)inFile->Get("covaMat_histo");
    covaMat_histo_1x1Det->SetName("covaMat_histo_1x1Det");
    
    TMatrixD covaMat_matrix_1x1Det(26,26);
    for (int i = 0 ; i < 26 ; i++) {
        for (int j = 0 ; j < 26 ; j++) {
            covaMat_matrix_1x1Det(i,j) = covaMat_histo_1x1Det->GetBinContent(i+1,j+1);
        }
    }
    cout << endl;
    TVectorD valuesTest;
    TMatrixD vectorsTest = covaMat_matrix_1x1Det.EigenVectors(valuesTest);
    for (Int_t i = 0; i < valuesTest.GetNrows(); ++i) {
        if (valuesTest(i) < 0)
            cout << "1x1Det-Matrix eigen-value " << i << " is " << valuesTest(i) << endl;
    }
    cout << endl;

    double rho_array[6][6]={{1.0,0.0,0.0,0.0,0.0,0.0},
                            {0.0,1.0,0.0,0.0,0.0,0.0},
                            {0.0,0.0,1.0,0.0,0.0,0.0},
                            {0.0,0.0,0.0,1.0,0.0,0.0},
                            {0.0,0.0,0.0,0.0,1.0,0.0},
                            {0.0,0.0,0.0,0.0,0.0,1.0}};
    
    //Construct correlation matrix upper triangle
    double corr_start = 0.99;
    for (int j = 0 ; j < 6 ; j++){
        for (int k = 0 ; k < 6 ; k++){
            for (int i = 0 ; i < 6 ; i++) {
                    int ij = i+j;
                    int ik = i+k;
                    if ( (ij)<6 && (ik<6) ) {
                        if (ij==ik)
                            rho_array[ij][ik]=1.0;
                        else if ( ik==(ij+1) )
                            rho_array[ij][ik]=corr_start;
                        else if ( ij<(6-1) )
                            rho_array[ij][ik] = rho_array[ij][ik-1]*rho_array[ij+1][ik];
                }//if
            }//for i
        }//for k
    }//for j
    
    //Assign correlation matrixlower triangle
    double corrMat_6x6Block[6][6];
    for(int ii = 0 ; ii < 6 ; ii++){
        for(int jj = ii+1 ; jj < 6 ; jj++){
            rho_array[jj][ii] = rho_array[ii][jj];
        }
    }
    
    cout << "The AD-correlation Matrix is:" << endl;
    for(int ii = 0 ; ii < 6 ; ii++){
        for(int jj = 0 ; jj < 6 ; jj++){
            corrMat_6x6Block[ii][jj] =  rho_array[ii][jj];
            cout << corrMat_6x6Block[ii][jj] << "\t";
        }
        cout << endl;
    }
    
    //-- Histograms definition
    //X axis
    int NBx = covaMat_histo_1x1Det->GetXaxis()->GetNbins();
    double lox = covaMat_histo_1x1Det->GetXaxis()->GetBinLowEdge(1);
    double hix = covaMat_histo_1x1Det->GetXaxis()->GetBinUpEdge(NBx);
    double xbins[NB+1];
    covaMat_histo_1x1Det->GetXaxis()->GetLowEdge(xbins);
    xbins[NB] = hix;
    //Y axis
    int NBy = covaMat_histo_1x1Det->GetYaxis()->GetNbins();
    double loy = covaMat_histo_1x1Det->GetYaxis()->GetBinLowEdge(1);
    double hiy = covaMat_histo_1x1Det->GetYaxis()->GetBinUpEdge(NBy);
    double ybins[NB+1];
    covaMat_histo_1x1Det->GetYaxis()->GetLowEdge(ybins);
    ybins[NB] = hiy;

    cout << "NBx = " << NBx << " NBy = " << NBy << endl;
    
    //-- Number of detectors
    int nDet = 6;
    int    NBx_6x6 = NBx*nDet;
    double lox_6x6 = 0;
    double hix_6x6 = NBx*nDet;
    int    NBy_6x6 = NBy*nDet;
    double loy_6x6 = 0;
    double hiy_6x6 = NBy*nDet;
    TH2F *covaMat_histo_6x6Det = new TH2F("covaMat_histo_6x6Det","",NBx_6x6,lox_6x6,hix_6x6,NBy_6x6,loy_6x6,hiy_6x6);
    cout << NBx_6x6 << "  " << NBy_6x6 << endl;
    for (int i = 0 ; i < NBx_6x6 ; i++) {
        for (int j = i ; j < NBy_6x6 ; j++) {
            int ii = i%NBx;
            int jj = j%NBy;
            int iB = int(i/NBx);
            int jB = int(j/NBy);
            double block_corrFact_ij = corrMat_6x6Block[iB][jB];
            //cout << i << " " << j << "   -   " << ii << " " << jj << "   -   " << iB << " " << jB << endl;
            double value = covaMat_histo_1x1Det->GetBinContent(ii+1,jj+1)*block_corrFact_ij;
            covaMat_histo_6x6Det->SetBinContent(i+1,j+1,value);
            covaMat_histo_6x6Det->SetBinContent(j+1,i+1,value);
        }
    }
    double Td_vector[156] = {3094.27,
        3328.37,
        4400.68,
        5464.73,
        6239.45,
        6825.7,
        7107.09,
        7254.89,
        7102.13,
        6719.56,
        6176.3,
        5648.9,
        5154.588,
        4625.2,
        4218.162,
        3616.7,
        3148.497,
        2574.81,
        2175.713,
        1713.45,
        1348.74,
        924.511,
        693.7145,
        469.199,
        298.912,
        460.933,
        3147.79,
        3435.51,
        4509.68,
        5495.3,
        6244.22,
        6933.23,
        7178.81,
        7262.43,
        7217.66,
        6817.68,
        6288.996,
        5647.067,
        5199.029,
        4701.6,
        4231.84,
        3681.09,
        3218.577,
        2602.97,
        2224.395,
        1696.686,
        1342.14,
        990.8852,
        720.6146,
        498.076,
        307.141,
        447.05,
        2824.71,
        3081.1,
        4183.73,
        4990.59,
        5700.55,
        6328.75,
        6595.58,
        6641.67,
        6542.75,
        6152.1,
        5745.987,
        5171.283,
        4730.165,
        4343.223,
        3923.978,
        3381.918,
        2891.33,
        2387.62,
        1933.05,
        1543.751,
        1237.223,
        884.256,
        662.1829,
        434.3892,
        296.7714,
        351.28,
        445.208,
        456.644,
        644.744,
        766.2,
        897.515,
        940.892,
        1000.83,
        1023.31,
        1014.24,
        952.328,
        865.1798,
        827.717,
        738.596,
        684.5727,
        586.382,
        526.8371,
        424.308,
        388.0292,
        301.2751,
        248.8286,
        191.254,
        138.8071,
        98.979,
        64.277203,
        38.6452,
        47.3206,
        430.331,
        470.371,
        635.584,
        793.8,
        861.44,
        919.362,
        1014.21,
        1024.71,
        1014.6,
        948.129,
        882.4329,
        803.9076,
        695.839,
        644.1373,
        625.0897,
        528.6827,
        443.16,
        380.5733,
        311.7677,
        232.4658,
        188.926,
        141.1113,
        98.73912,
        63.36419,
        42.3723,
        51.702,
        422.049,
        462.337,
        636.007,
        749.831,
        848.009,
        945.405,
        1017.38,
        1024.42,
        1005.25,
        948.534,
        856.223,
        787.3812,
        722.4516,
        644.6125,
        585.158,
        506.1464,
        440.042,
        378.6312,
        294.9266,
        240.9473,
        195.183,
        142.7697,
        87.2261,
        59.4545,
        37.5502,
        46.1555
    };
    double Md_vector[156] = {3039.94,
        3272.28,
        4351.97,
        5322.1,
        6034.94,
        6713.46,
        6941.25,
        7101.41,
        6954.63,
        6635.88,
        6099.08,
        5546.61,
        5064.89,
        4668.08,
        4200.52,
        3649.12,
        3110.81,
        2614.94,
        2133.24,
        1608.7,
        1340.27,
        929.638,
        660.483,
        504.835,
        307.041,
        477.479,
        3099.9,
        3318.7,
        4409.64,
        5389.58,
        6110.03,
        6796.79,
        7028.15,
        7189.67,
        7039.97,
        6717.12,
        6173.75,
        5614.58,
        5127.01,
        4725.38,
        4252.14,
        3694.02,
        3149.16,
        2647.27,
        2159.72,
        1628.82,
        1357.11,
        941.477,
        669.044,
        511.487,
        311.263,
        488.216,
        2807.21,
        2923.86,
        4008.58,
        4881.98,
        5553.97,
        6129.57,
        6447.77,
        6534,
        6553.21,
        6190.11,
        5637.32,
        5110.49,
        4596.41,
        4187.65,
        3831.99,
        3304.82,
        2948.64,
        2447.71,
        1934.18,
        1512.22,
        1234.88,
        970.804,
        641.312,
        417.044,
        285.022,
        318.935,
        409.609,
        446.861,
        619.205,
        703.954,
        795.304,
        927.254,
        916.734,
        935.412,
        910.737,
        891.562,
        784.429,
        749.902,
        695.07,
        621.849,
        560.751,
        491.587,
        432.521,
        361.262,
        290.007,
        228.839,
        173.707,
        130.898,
        90.1882,
        63.9215,
        43.5773,
        49.0334,
        409.743,
        445.277,
        616.549,
        700.684,
        791.456,
        922.742,
        912.336,
        930.872,
        906.204,
        887.101,
        780.502,
        746.148,
        691.59,
        618.737,
        557.946,
        489.129,
        430.36,
        359.459,
        288.563,
        227.703,
        172.85,
        130.255,
        89.7499,
        63.6132,
        43.3705,
        48.8809,
        405.773,
        440.301,
        609.483,
        692.557,
        782.217,
        911.961,
        901.7,
        919.999,
        895.577,
        876.689,
        771.34,
        737.389,
        683.472,
        611.474,
        551.397,
        483.389,
        425.31,
        355.242,
        285.179,
        225.034,
        170.825,
        128.731,
        88.7014,
        62.8709,
        42.8657,
        48.3429
    }; //tiene los valores Md
    
    TMatrixD delta_vector(156,1);
    TMatrixD transp_delta_vector(1,156);
    for (int i = 0 ; i < 156 ; i++) {
        delta_vector(i,0) = Md_vector[i] - Td_vector[i];
        transp_delta_vector(0,i) = delta_vector(i,0);
    }
    //delta_vector.Print();

    TH2F *fullnorm_covaMat_histo_6x6Det = new TH2F("fullnorm_covaMat_histo_6x6Det","",NBx_6x6,lox_6x6,hix_6x6,NBy_6x6,loy_6x6,hiy_6x6);
    for (int i = 0 ; i < NBx_6x6 ; i++) {
        for (int j = i ; j < NBy_6x6 ; j++) {
            int ii = i%NBx;
            int jj = j%NBy;
            int iB = int(i/NBx);
            int jB = int(j/NBy);
            double block_corrFact_ij = corrMat_6x6Block[iB][jB];
            //cout << iB << ", " << jB << " " << block_corrFact_ij << endl;
            double value = Td_vector[i]*Td_vector[j]*covaMat_histo_1x1Det->GetBinContent(ii+1,jj+1)*block_corrFact_ij;
            if (ii != jj) {
                value *= 1.0;
            }
            if (i == j) {
                value += Md_vector[i]*(1.0);
            }
            fullnorm_covaMat_histo_6x6Det->SetBinContent(i+1,j+1,value);
            fullnorm_covaMat_histo_6x6Det->SetBinContent(j+1,i+1,value);
            //verificar que norm_covaMat_histo_6x6Det se puede invertir (hacer igual que lo que está abajo)
        }
    }
    //-- Test matrix inversion
    //-- Creating the Rebinned Matrix with root tools.
    TMatrixD covaMat_matrix(NBx_6x6,NBy_6x6);
    TMatrixD fullnorm_covaMat_matrix(NBx_6x6,NBy_6x6);
    TMatrixD fullnorm_covaMatBlock_matrix(NBx,NBy);
    covaMat_matrix.Zero();
    fullnorm_covaMat_matrix.Zero();
    fullnorm_covaMatBlock_matrix.Zero();
    for (int i = 0 ; i < NBx_6x6 ; i++) {
        for (int j = 0 ; j < NBy_6x6; j++) {
            covaMat_matrix(i,j)      = covaMat_histo_6x6Det->GetBinContent(i+1,j+1);
            fullnorm_covaMat_matrix(i,j) = fullnorm_covaMat_histo_6x6Det->GetBinContent(i+1,j+1);
            int c1 = 4;
            int c2 = 5;
            if ((c1*NBx <= i && i < c2*NBx) && (c1*NBy <= j && j < c2*NBy)) {
                fullnorm_covaMatBlock_matrix(i-c1*NBx,j-c1*NBy) = fullnorm_covaMat_matrix(i,j);
            }
        }
    }
    cout << endl;
    //fullnorm_covaMatBlock_matrix.Print();
    TVectorD valuesBL;
    TMatrixD vectors = fullnorm_covaMatBlock_matrix.EigenVectors(valuesBL);
    for (Int_t i = 0; i < valuesBL.GetNrows(); ++i) {
        TVectorD vector(26);
        if (valuesBL(i) < 0)
            cout << "Block-Matrix eigen-value " << i << " is " << valuesBL(i) << endl;
    }
    cout << endl;
    
    TMatrixD inv_mat(NBx_6x6,NBy_6x6);
    TMatrixD uno_mat(NBx_6x6,NBy_6x6);
    inv_mat = covaMat_matrix;
    inv_mat.Invert();
    uno_mat = covaMat_matrix*inv_mat;
    TMatrixD inv_fullnormmat(NBx_6x6,NBy_6x6);
    TMatrixD uno_fullnormmat(NBx_6x6,NBy_6x6);
    inv_fullnormmat = fullnorm_covaMat_matrix;
    inv_fullnormmat.Invert();
    uno_fullnormmat = fullnorm_covaMat_matrix*inv_fullnormmat;

    double sqr_chi = 0.0;
    for (int i = 0 ; i < 156 ; i++) {
        for (int j = 0 ; j < 156 ; j++) {
            double delta_i = transp_delta_vector(0,i);
            double delta_j = delta_vector(j,0);
            double invMat_ij = inv_fullnormmat(i,j);
            sqr_chi += delta_i*invMat_ij*delta_j;
        }
    }
    
    cout << "chi2 = " << sqr_chi << endl;
/*
    TVectorD values;
    //TMatrixD vectors = fullnorm_covaMat_matrix.EigenVectors(values);
    TMatrixD vectors = fullnorm_covaMat_matrix.EigenVectors(values);
    for (Int_t i = 0; i < values.GetNrows(); ++i) {
        //TVectorD vector(TMatrixTColumn_const<double>(vectors, i));
        TVectorD vector(156);
        if (values(i) < 0)
            cout << "eigen-value " << i << " is " << values(i) << " with eigen-vector" << endl;
        //vector.Print();
    }
*/
    if(fullnorm_covaMat_matrix.Determinant() == 0){
    //    inv_mat = fullnorm_covaMat_matrix;
    //    inv_mat.Invert();
    //}
        //else {
        cout << "Determinant = " << fullnorm_covaMat_matrix.Determinant() << endl;
        cout << "Singular fullCovaMatrix!" << endl;
        cout << "Execution aborted!" << endl;
        //break;
    }
    
    
    TH2F *unoMat_histo = new TH2F("unoMat_histo","",NBx_6x6,lox_6x6,hix_6x6,NBy_6x6,loy_6x6,hiy_6x6);
    TH2F *norm_unoMat_histo = new TH2F("norm_unoMat_histo","",NBx_6x6,lox_6x6,hix_6x6,NBy_6x6,loy_6x6,hiy_6x6);
    for (int k = 0 ; k < NBx_6x6 ; k++) {
        for (int l = 0 ; l < NBy_6x6 ; l++) {
            unoMat_histo->SetBinContent(k+1,l+1,uno_mat(k,l));
            norm_unoMat_histo->SetBinContent(k+1,l+1,uno_fullnormmat(k,l));
        }
    }
    
    TH2F *invMat_histo = new TH2F("invMat_histo","",NBx_6x6,lox_6x6,hix_6x6,NBy_6x6,loy_6x6,hiy_6x6);
    TH2F *norm_invMat_histo = new TH2F("norm_invMat_histo","",NBx_6x6,lox_6x6,hix_6x6,NBy_6x6,loy_6x6,hiy_6x6);
    for (int k = 0 ; k < NBx_6x6 ; k++) {
        for (int l = 0 ; l < NBy_6x6 ; l++) {
            invMat_histo->SetBinContent(k+1,l+1,inv_mat(k,l));
            norm_invMat_histo->SetBinContent(k+1,l+1,inv_fullnormmat(k,l));
        }
    }

    // Drawng section
    set_plot_style();
    //---------------------------------------------------------
    //TH2F *frame = new TH2F("frame","",NBx_6x6,lox_6x6,hix_6x6,NBy_6x6,loy_6x6,hiy_6x6);
    TCanvas *canv0 = new TCanvas("canv0","",740,740);
    //frame->Draw();
    covaMat_histo_6x6Det->Draw("COLZ");
    gPad->SetTicks(1,1);
    //---------------------------------------------------------
    TCanvas *canv2 = new TCanvas("canv2","",3*400,400);
    canv2->Divide(3,1);
    canv2->cd(1);
    covaMat_histo_6x6Det->Draw("COLZ");
    canv2->cd(2);
    invMat_histo->Draw("COLZ");
    canv2->cd(3);
    unoMat_histo->Draw("COLZ");

    //canv0->Print("files_plots/canv_DB.pdf");
    //---------------------------------------------------------
    // write to output file
    TFile *fout = new TFile("../files_data/db_CovaMatrix_6AD_6x26bins.root","recreate");
    fout->cd();
    
    covaMat_histo_6x6Det->Write();

    fout->Close();

    //Testing Norm matrix
    TCanvas *canv01 = new TCanvas("canv01","",740,740);
    fullnorm_covaMat_histo_6x6Det->Draw("COLZ");
    gPad->SetTicks(1,1);
    //---------------------------------------------------------
    TCanvas *canv02 = new TCanvas("canv02","",3*400,400);
    canv02->Divide(3,1);
    canv02->cd(1);
    fullnorm_covaMat_histo_6x6Det->Draw("COLZ");
    canv02->cd(2);
    norm_invMat_histo->Draw("COLZ");
    canv02->cd(3);
    norm_unoMat_histo->Draw("COLZ");

/*    //Testing 1x1 Det matrix
    TCanvas *canv03 = new TCanvas("canv03","",740,740);
    covaMat_histo_1x1Det->Draw("COLZ");
    gPad->SetTicks(1,1);
*/
}// end
