//Constants definition
const int NB = 35;      //Correlation matrix Number of bins
//------------------------------------------------------
void set_plot_style()
{
    const Int_t NRGBs = 5;
    const Int_t NCont = 255;
    
    Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
    Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
    Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    gStyle->SetNumberContours(NCont);
}
//------------------------------------------------------

void db_CorrMatrix_6AD8ADx35bins()
{// begin

    //------------- Style --------------
    gROOT->SetStyle("Plain");
    gStyle->SetOptStat(0);
    //------------- Style ----------------------------------------------
    //----------  Text Style  ---------
    //ft = 10 * fontID + precision
    Int_t ft = 10 * 4 + 2;
    Double_t sz = 0.04;
    //------------------------------------------------------------------

    //--6AD 8x8 correlation matrix
    TFile *fmat6 = new TFile("db_CorrMatrix_6AD_35bins.root");
    TH2F *corrMat_histo_6Det_1x1 = (TH2F*)fmat6->Get("corrMatRebHisto");
    cout << "6AD Bins: " << corrMat_histo_6Det_1x1->GetXaxis()->GetNbins() << endl;
    //--8AD 1x1 correlation matrix
    TFile *fmat8 = new TFile("db_CorrMatrix_8AD_35bins.root");
    TH2F *corrMat_histo_8Det_1x1 = (TH2F*)fmat8->Get("corrMatRebHisto");
    cout << "8AD Bins: " << corrMat_histo_8Det_1x1->GetXaxis()->GetNbins() << endl;
    //matrix definition
    
    //---------------------------------------------------------------------------
    //---------------------------------------------------------------------------
    // 8x8 correlation matrix (constructed by hand ... and by "eye")
    
    double rho_array[8][8]={{1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},
                            {0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0},
        {0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0},
        {0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0},
        {0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0},
        {0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0},
        {0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0},
        {0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0}};
    
    //Construct correlation matrix upper triangle
    double corr_start = 0.999;
    for (int j = 0 ; j < 8 ; j++){
        for (int k = 0 ; k < 8 ; k++){
            for (int i = 0 ; i < 8 ; i++) {
                int ij = i+j;
                int ik = i+k;
                if ( (ij)<8 && (ik<8) ) {
                    if (ij==ik)
                        rho_array[ij][ik]=1.0;
                    else if ( ik==(ij+1) )
                        rho_array[ij][ik]=corr_start;
                    else if ( ij<(8-1) )
                        rho_array[ij][ik] = rho_array[ij][ik-1]*rho_array[ij+1][ik];
                }//if
            }//for i
        }//for k
    }//for j
    
    //Assign correlation matrixlower triangle
    double corrMat_8x8Block[8][8];
    for(int ii = 0 ; ii < 8 ; ii++){
        for(int jj = ii+1 ; jj < 8 ; jj++){
            rho_array[jj][ii] = rho_array[ii][jj];
        }
    }
    
    cout << "The AD-correlation Matrix is:" << endl;
    for(int ii = 0 ; ii < 8 ; ii++){
        for(int jj = 0 ; jj < 8 ; jj++){
            corrMat_8x8Block[ii][jj] =  rho_array[ii][jj];
            cout << setprecision(4) << corrMat_8x8Block[ii][jj] << "\t";
        }
        cout << endl;
    }
    //-- building the large (8*35)x(8*35) 6AD Correlation matrix
    //-- Number of detectors
    int NBx = NB;
    int NBy = NB;
    int nDet = 6+2;
    int    NBx_8x8 = NBx*nDet;
    double lox_8x8 = 0;
    double hix_8x8 = NBx*nDet;
    int    NBy_8x8 = NBy*nDet;
    double loy_8x8 = 0;
    double hiy_8x8 = NBy*nDet;
    TMatrixD corrMat_mat_8x8_6Det(NBx_8x8,NBy_8x8);
    TH2F *corrMat_histo_8x8_6Det = new TH2F("corrMat_histo_8x8_6Det","",NBx_8x8,lox_8x8,hix_8x8,NBy_8x8,loy_8x8,hiy_8x8);
    cout << NBx_8x8 << "  " << NBy_8x8 << endl;
    for (int i = 0 ; i < NBx_8x8 ; i++) {
        for (int j = i ; j < NBy_8x8 ; j++) {
            int ii = i%NBx;
            int jj = j%NBy;
            int iB = int(i/NBx);
            int jB = int(j/NBy);
            double block_corrFact_ij = corrMat_8x8Block[iB][jB];
            //cout << i << " " << j << "   -   " << ii << " " << jj << "   -   " << iB << " " << jB << endl;
            double value = corrMat_histo_6Det_1x1->GetBinContent(ii+1,jj+1)*block_corrFact_ij;
            if (iB == 3 || iB == 7)
                value = 0.0;
            if (jB == 3 || jB == 7)
                value = 0.0;
            corrMat_histo_8x8_6Det->SetBinContent(i+1,j+1,value);
            corrMat_histo_8x8_6Det->SetBinContent(j+1,i+1,value);
            corrMat_mat_8x8_6Det(i,j) = value;
            corrMat_mat_8x8_6Det(j,i) = value;
        }
    }
    //---------------------------------------------------------------------------
    //-- building the large (8*35)x(8*35) 8AD Correlation matrix
    TMatrixD corrMat_mat_8x8_8Det(NBx_8x8,NBy_8x8);
    TH2F *corrMat_histo_8x8_8Det = new TH2F("corrMat_histo_8x8_8Det","",NBx_8x8,lox_8x8,hix_8x8,NBy_8x8,loy_8x8,hiy_8x8);
    cout << NBx_8x8 << "  " << NBy_8x8 << endl;
    for (int i = 0 ; i < NBx_8x8 ; i++) {
        for (int j = i ; j < NBy_8x8 ; j++) {
            int ii = i%NBx;
            int jj = j%NBy;
            int iB = int(i/NBx);
            int jB = int(j/NBy);
            double block_corrFact_ij = corrMat_8x8Block[iB][jB];
            //cout << i << " " << j << "   -   " << ii << " " << jj << "   -   " << iB << " " << jB << endl;
            double value = corrMat_histo_8Det_1x1->GetBinContent(ii+1,jj+1)*block_corrFact_ij;
            corrMat_histo_8x8_8Det->SetBinContent(i+1,j+1,value);
            corrMat_histo_8x8_8Det->SetBinContent(j+1,i+1,value);
            corrMat_mat_8x8_8Det(i,j) = value;
            corrMat_mat_8x8_8Det(j,i) = value;
        }
    }
    TMatrixD offDiagBlock_6DetX8Det(NBx_8x8,NBy_8x8);
    offDiagBlock_6DetX8Det.Mult(corrMat_mat_8x8_6Det,corrMat_mat_8x8_8Det);
    //---------------------------------------------------------------------------
    double corrMat_2x2BigBlock[2][2] = {{1.0,0.99},{0.99,1.0}};
    double value;
    TH2F *corrMat_histo_6Det8Det = new TH2F("corrMat_histo_6Det8Det","",2*NBx_8x8,0,2*NBx_8x8,2*NBy_8x8,0,2*NBy_8x8);
    //-- there is something worng here, with the matrix filling - 2019.05.03
    for (int i = 0 ; i < 2*NBx_8x8 ; i++) {
        for (int j = 0 ; j < 2*NBy_8x8 ; j++) {
            if ( (i < NBx_8x8) && (j < NBy_8x8) ) {
                value = corrMat_mat_8x8_6Det(i,j);
                corrMat_histo_6Det8Det->SetBinContent(i+1,j+1,value);
                corrMat_histo_6Det8Det->SetBinContent(j+1,i+1,value);
            }
            if ( (i >= NBx_8x8) && (j >= NBy_8x8) ) {
                value = corrMat_mat_8x8_8Det(i-NBx_8x8,j-NBy_8x8);
                corrMat_histo_6Det8Det->SetBinContent(i+1,j+1,value);
                corrMat_histo_6Det8Det->SetBinContent(j+1,i+1,value);
            }/*
            if ( (i < NBx_8x8) && (j >= NBy_8x8) ) {
                value = offDiagBlock_6DetX8Det(i,j-NBy_8x8);
                corrMat_histo_6Det8Det->SetBinContent(i+1,j+1,value);
                corrMat_histo_6Det8Det->SetBinContent(j+1,i+1,value);
            }*/
        }
    }
    
    //---------------------------------------------------------
    set_plot_style();
    /*
    TCanvas *canv0 = new TCanvas("canv0","",2*500,500);
    canv0->Divide(2,1);
    canv0->cd(1);
    corrMat_histo_8x8_6Det->Draw("COLZ");
    canv0->cd(2);
    corrMat_histo_8x8_8Det->Draw("COLZ");*/
    //---------------------------------------------------------
    TCanvas *canv1 = new TCanvas("canv1","",900,100,500,500);
    canv1->cd(1);
    corrMat_histo_6Det8Det->Draw("COLZ");
    

    /*
    //---------------------------------------------------------------------------

    //--- Graph with errors from https://arxiv.org/pdf/1508.04233.pdf  Fig. 2
    const int Ngr = 26;
    double erxx[Ngr] = {0.975,1.375,1.625,1.875,2.125,2.375,2.625,2.875,3.125,3.375,
        3.625,3.875,4.125,4.375,4.625,4.875,5.125,5.375,5.625,5.875,
        6.125,6.375,6.625,6.875,7.500,10.0};
    //double eryy[Ngr] =
    //{
    //1.12448,1.04615,1.03497,1.03030,1.02751,1.02657,1.02471,1.02378,1.02284,1.02378,
    //                1.02471,1.02657,1.02751,1.02937,1.03217,1.03497,1.03776,1.04149,1.04336,1.04615,
    //                1.05175,1.06014,1.06853,1.07786,1.24009,1.5
    //};
    double eryy[Ngr] =
    {
        1.0115,1.0064,1.0053,1.0046,1.0043,1.0039,1.0031,1.0031,1.0027,1.0027,
        1.0035,1.0027,1.0035,1.0039,1.0039,1.0042,1.0050,1.0058,1.0065,1.0080,
        1.0099,1.0122,1.0152,1.0186,1.0227,1.0642
    };
    
    for(int ii = 0 ; ii < Ngr ; ii++)
        eryy[ii] = eryy[ii]*1.0; //2019.03.19 - Testing a larger systematics error
    
    TGraph *errgr = new TGraph(Ngr,erxx,eryy);
    printf("f(3.5)=%f\n",errgr->Eval(3.5));
    errgr->Draw("AC");
    //---------------
    

    //define histogram
    double    lo = 0.7;
    double    hi = 12.0;

    //Bins for the rebinned matrix
    double xbins[NB+1];
    xbins[0] = 0.7;
    double delta_bins2 = (7.9 - 1.3)/33; // 0.2 MeV/bin
    for (int i = 0 ; i < (NB-1) ; i++)
        xbins[i+1] = 1.3 + delta_bins2*i;
    xbins[NB] = hi;

    
    
    TH1F *e_syst_histo1 = new TH1F("e_syst_histo1","",NB,xbins);
    e_syst_histo1->SetLineColor(kRed);
    //e_syst_histo1->SetFillColor(kRed);
    TH1F *e_syst_histo2 = new TH1F("e_syst_histo2","",NB,xbins);
    e_syst_histo2->SetLineColor(kRed);
    //e_syst_histo2->SetFillColor(10);

    //- array with total systematic error interpolated from errgr TGraph
    double e_syst[NB];
    for (int i=0;i<NB;i++){
        double xpoint = e_syst_histo1->GetBinCenter(i+1);
        e_syst[i] = errgr->Eval(xpoint);
        //cout << e_syst[i] << endl; //checking
    }
    //cout << endl; //checking
    //cout << "end Check" << endl; //checking
    
    //---Fill e_syst_histo1 and e_syst_histo2
    for (int i = 0 ; i < NB ; i++) {
        e_syst_histo1->SetBinContent(i+1,e_syst[i]);
        e_syst_histo2->SetBinContent(i+1,2.0-e_syst[i]);
    }

    //- assign values to covariance matrix elements
    TH2F *covaMat_histo = new TH2F("covaMat_histo","",NB,xbins,NB,xbins);
    for (int i = 0 ; i < NB ; i++) {
        double sigma_i = e_syst_histo1->GetBinContent(i+1) - 1.0;
        for (int j = 0 ; j < NB ; j++) {
            double sigma_j = (e_syst_histo1->GetBinContent(j+1)) - 1.0;
            double rho_ij = corrMat_histo->GetBinContent(i+1,j+1);
            double cont = sigma_i*sigma_j*rho_ij;
            covaMat_histo->SetBinContent(i+1,j+1,cont);
        }
    }
    
    //covaMat_histo->SetMinimum(-1.0);
    //covaMat_histo->SetMaximum(+1.0);

    //Creating the Rebinned Matrix with root tools.
    TMatrixD covaMat_matrix(NB,NB);
    covaMat_matrix.Zero();
    for (int i = 0 ; i < NB ; i++) {
        for (int j = 0 ; j < NB; j++) {
            covaMat_matrix(i,j) = covaMat_histo->GetBinContent(i+1,j+1);
        }
    }
    //------------------------------------------------------
    // Positive definiteness correction
    
    TVectorD values;
    TMatrixD vectors = covaMat_matrix.EigenVectors(values);
    for (Int_t i = 0; i < values.GetNrows(); ++i) {
        cout << "eigen-value " << i << " is " << values(i) << " with eigen-vector" << endl;
    }
    //vectors.Print();
    
    TMatrixD vectors_I = vectors;
    vectors_I.T();
    
    TMatrixD temp(NB,NB);
    TMatrixD diagMat(NB,NB);
    diagMat.Zero();
    for (int i = 0 ; i < NB ; i++) {
        //if(values(i) < 0.0 && abs(values(i)) < 1.0e-4)
        if(values(i) < 0.0)
            diagMat(i,i) = 0.0;
        else
            diagMat(i,i) = values(i);
    }
    temp = vectors*diagMat;
    
    TMatrixD fixedMat(NB,NB);
    TMatrixD fixedMatScaled(NB,NB);
    fixedMat = temp*vectors_I;
    
    for (int i = 0 ; i < NB ; i++)
        for (int j = 0 ; j < NB ; j++)
            fixedMatScaled(i,j) = fixedMat(i,j)/sqrt(fixedMat(i,i)*fixedMat(j,j));
    
    TVectorD values2;
    vectors = fixedMatScaled.EigenVectors(values2);
    cout << "First iteration" << endl;
    for (Int_t i = 0 ; i < values2.GetNrows(); ++i) {
        cout << "eigen-value " << i << " is " << values(i) << " " << values2(i) << endl;
    }
    cout << "Printing Fixed Matrix Scaled..." << endl;
    //fixedMatScaled.Print();
    fixedMat.Print();
    
    //covaMat_matrix = fixedMatScaled; // replace covaMat_Matrix with fixed matrix
    covaMat_matrix = fixedMat; // replace covaMat_Matrix with fixed matrix
    //Replace histo
    for (int k = 0 ; k < NB ; k++) {
        for (int l = 0 ; l < NB ; l++) {
            //covaMat_histo->SetBinContent(k+1,l+1,fixedMatScaled(k,l));
            covaMat_histo->SetBinContent(k+1,l+1,fixedMat(k,l)); // 08.03.2019 - Fixing fixedMat.
        }
    }
    
    //------------------------------------------------------
    cout << "Let's invert the matrix..." << endl;
    TMatrixD inv_mat(NB,NB);
    TMatrixD uno_mat(NB,NB);
    inv_mat = covaMat_matrix;
    inv_mat.Invert();
    uno_mat = covaMat_matrix*inv_mat;
    //printf("\n Testing Unitarity \n");
    //invReb_mat.Print();
    TH2F *unoMat_histo = new TH2F("unoMat_histo","",NB,xbins,NB,xbins);
    for (int k = 0 ; k < NB ; k++) {
        for (int l = 0 ; l < NB ; l++) {
            unoMat_histo->SetBinContent(k+1,l+1,uno_mat(k,l));
        }
    }
    
    TH2F *invMat_histo = new TH2F("invMat_histo","",NB,xbins,NB,xbins);
    for (int k = 0 ; k < NB ; k++) {
        for (int l = 0 ; l < NB ; l++) {
            invMat_histo->SetBinContent(k+1,l+1,inv_mat(k,l));
        }
    }

    TFile *outf = new TFile("./db_CovaMatrix_6AD_35bins.root","recreate");
    outf->cd();
    covaMat_histo->Write();
    covaMat_matrix.Write();
    outf->Close();
    //---------------------------------------------------------
    set_plot_style();

    TCanvas *canv0 = new TCanvas("canv0","",2*500,500);
    canv0->Divide(2,1);
    canv0->cd(1);
    corrMat_histo->Draw("COLZ");
    canv0->cd(2);
    covaMat_histo->Draw("COLZ");
    //---------------------------------------------------------
    TCanvas *canv2 = new TCanvas("canv2","",3*400,400);
    canv2->Divide(3,1);
    canv2->cd(1);
    covaMat_histo->Draw("COLZ");
    canv2->cd(2);
    invMat_histo->Draw("COLZ");
    canv2->cd(3);
    unoMat_histo->Draw("COLZ");
    //---------------------------------------------------------
    TH2F *e_frame = new TH2F("e_frame","",10,0.7,12,10,0.5,1.5);
    TCanvas *canv1 = new TCanvas("canv1","",500,300);
    canv1->cd();
    e_frame->Draw();
    e_syst_histo1->Draw("same hist");
    e_syst_histo2->Draw("same hist");
    gPad->Update();
    TLine* line1 = new TLine(gPad->GetUxmin(),1,gPad->GetUxmax(),1);
    line1->Draw();
    gPad->SetTicks(1,1);
*/
}// end

