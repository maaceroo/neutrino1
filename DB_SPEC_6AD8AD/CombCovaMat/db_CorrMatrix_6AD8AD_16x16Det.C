//-- db_CorrMatrix_6AD8AD_16x16Det.C - By M.A. Acero and A.A. Aguilar-Arevao - 24.05.2019 --//
// This macro creates the correlation matrix necessary to perform the analysis
// of the Daya-Bay data from the 6AD and 8AD periods combined.
// The obtained correlation matrix is an approximate reproduction of the one
// presented by Henoch Wong's Slides at Joint DC-RENO-DYB workshop (2016)
// https://indico.snu.ac.kr/indico/event/4/session/25/contribution/32/material/slides/0.pdf

//Constants definition
const int NB = 35;      //Single detector Correlation matrix Number of bins
//------------------------------------------------------
/*
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
*/
//------------------------------------------------------

void db_CorrMatrix_6AD8AD_16x16Det()
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
    corrMat_histo_6Det_1x1->SetTitle("6AD Block");
    std::cout << "6AD Bins: " << corrMat_histo_6Det_1x1->GetXaxis()->GetNbins() << std::endl;
    //--8AD 1x1 correlation matrix
    TFile *fmat8 = new TFile("db_CorrMatrix_8AD_35bins.root");
    TH2F *corrMat_histo_8Det_1x1 = (TH2F*)fmat8->Get("corrMatRebHisto");
    corrMat_histo_8Det_1x1->SetTitle("8AD Block");
    std::cout << "8AD Bins: " << corrMat_histo_8Det_1x1->GetXaxis()->GetNbins() << std::endl;
    //matrix definition
    TMatrixD DiagBlock6AD_mat(NB,NB);
    TMatrixD DiagBlock8AD_mat(NB,NB);
    for (int i = 0 ; i < NB ; i++) {
        for (int j = 0 ; j < NB ; j++) {
            DiagBlock6AD_mat(i,j) = corrMat_histo_6Det_1x1->GetBinContent(i+1,j+1);
            DiagBlock8AD_mat(i,j) = corrMat_histo_8Det_1x1->GetBinContent(i+1,j+1);
        }
    }

    TMatrixD offDiagBlock_6DetX8Det(NB,NB);
    offDiagBlock_6DetX8Det.Mult(DiagBlock6AD_mat,DiagBlock8AD_mat);
    //offDiagBlock_6DetX8Det.Mult(DiagBlock8AD_mat,DiagBlock6AD_mat);
    TH2F *offDiagBlock_6DetX8Det_histo = new TH2F("offDiagBlock_6DetX8Det_histo","6ADBlk x 8ADBlk",NB,0,NB, NB,0,NB);
    for (int i = 0 ; i < NB ; i++) {
        for (int j = 0 ; j < NB ; j++) {
            double val = offDiagBlock_6DetX8Det(i,j);
            offDiagBlock_6DetX8Det_histo->SetBinContent(i+1,j+1,val);
        }
    }

    TVectorD valEigen;
    TMatrixD vecEigen = offDiagBlock_6DetX8Det.EigenVectors(valEigen);
    valEigen.Print();
    //vecEigen.Print();
    
    TMatrixD vectors_I = vecEigen;
    vectors_I.T();
    
    TMatrixD temp(NB,NB);
    TMatrixD diagMat(NB,NB);
    diagMat.Zero();
    for (int i = 0 ; i < NB ; i++) {
        //if(values(i) < 0.0 && abs(values(i)) < 1.0e-4)
        if(valEigen(i) < 0.0)
            diagMat(i,i) = 0.0;
        else
            diagMat(i,i) = valEigen(i);
    }
    temp = vecEigen*diagMat;
    
    TMatrixD fixedMat(NB,NB);
    fixedMat = temp*vectors_I;

    TVectorD valEigenFix;
    TMatrixD vecEigenFix = fixedMat.EigenVectors(valEigenFix);
    valEigenFix.Print();

    TMatrixD matP(NB,NB);
    matP.Zero();
    TMatrixD sqrtOffDiag_Mat(NB,NB);
    TMatrixD Q = vecEigen;
    TMatrixD QInv = Q;
    QInv.Invert();
    TMatrixD one(NB,NB);
    one.Mult(Q,QInv);
    //one.Print();
    TH2F *one_histo = new TH2F("one_histo","",NB,0,NB,NB,0,NB);
    for (int k = 0 ; k < NB ; k++) {
        for (int l = 0 ; l < NB ; l++) {
            one_histo->SetBinContent(k+1,l+1,one(k,l));
        }
        matP(k,k) = sqrt(abs(valEigenFix(k)));
    }
    
    TMatrixD diag(NB,NB);
    TMatrixD test(NB,NB);
    test = QInv*offDiagBlock_6DetX8Det;
    diag = test*Q;
    //diag.Print();
    TH2F *diagQ_histo = new TH2F("diagQ_histo","",NB,0,NB,NB,0,NB);

    TMatrixD matTemp(NB,NB);
    matTemp         = Q*matP;
    sqrtOffDiag_Mat = matTemp*QInv;
    
    TMatrixD sqrtOffDiag_MatScaled(NB,NB);
    //-- Forced correlation matrix to have ones in the diagonal and scales off-diagonal elements.
    for (int i = 0 ; i < NB ; i++)
        for (int j = 0 ; j < NB ; j++)
            sqrtOffDiag_MatScaled(i,j) = sqrtOffDiag_Mat(i,j)/sqrt(sqrtOffDiag_Mat(i,i)*sqrtOffDiag_Mat(j,j));
    
    TVectorD valEigenFixScaled;
    TMatrixD vecEigenFixScaled = sqrtOffDiag_MatScaled.EigenVectors(valEigenFixScaled);
    valEigenFixScaled.Print();
    
    TH2F *sqrtOffDiag_histo = new TH2F("sqrtOffDiag_histo","sqrt(6ADBlk x 8ADBlk)",NB,0,NB,NB,0,NB);

    for (int k = 0 ; k < NB ; k++) {
        for (int l = 0 ; l < NB ; l++) {
            diagQ_histo      ->SetBinContent(k+1,l+1,diag(k,l));
            sqrtOffDiag_histo->SetBinContent(k+1,l+1,sqrtOffDiag_MatScaled(k,l));
        }
    }
    sqrtOffDiag_histo->SetMinimum(-1.0);
    
    
    //---------------------------------------------------------------------------
    //---------------------------------------------------------------------------
    // 8x8 correlation matrix (constructed by hand ... and by "eye")
    
    double rho_array[16][16]={{1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},
                              {0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},
                              {0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},
                              {0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},
                              {0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},
                              {0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},
                              {0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},
                              {0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},
                              {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},
                              {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0},
                              {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0},
                              {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0},
                              {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0},
                              {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0},
                              {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0},
                              {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0}};
    
    //Construct correlation matrix upper triangle
    double corr_start = 0.999999;
    for (int j = 0 ; j < 16 ; j++){
        for (int k = 0 ; k < 16 ; k++){
            for (int i = 0 ; i < 16 ; i++) {
                int ij = i+j;
                int ik = i+k;
                if ( (ij)<16 && (ik<16) ) {
                    if (ij==ik)
                        rho_array[ij][ik]=1.0;
                    else if ( ik==(ij+1) )
                        rho_array[ij][ik]=corr_start;
                    else if ( ij<(16-1) )
                        rho_array[ij][ik] = rho_array[ij][ik-1]*rho_array[ij+1][ik];
                }//if
            }//for i
        }//for k
    }//for j
    
    //Assign correlation matrixlower triangle
    double corrMat_16x16Block[16][16];
    for(int ii = 0 ; ii < 16 ; ii++){
        for(int jj = ii+1 ; jj < 16 ; jj++){
            rho_array[jj][ii] = rho_array[ii][jj];
        }
    }
    
    std::cout << "The AD-correlation Matrix is:" << std::endl;
    for(int ii = 0 ; ii < 16 ; ii++){
        for(int jj = 0 ; jj < 16 ; jj++){
            corrMat_16x16Block[ii][jj] =  rho_array[ii][jj];
            std::cout << setprecision(7) << corrMat_16x16Block[ii][jj] << "   ";
        }
        std::cout << std::endl;
    }
    
    //-- building the large (8*35)x(8*35) 6AD Correlation matrix
    //-- Number of detectors
    int NBx = NB;
    int NBy = NB;
    int nDet = 16;
    int    NBx_16x16 = NBx*nDet;
    double lox_16x16 = 0;
    double hix_16x16 = NBx*nDet;
    int    NBy_16x16 = NBy*nDet;
    double loy_16x16 = 0;
    double hiy_16x16 = NBy*nDet;
    double value;
    TMatrixD corrMat_mat_16x16(NBx_16x16,NBy_16x16);
    TH2F *corrMat_histo_16x16 = new TH2F("corrMat_histo_16x16","",NBx_16x16,lox_16x16,hix_16x16,NBy_16x16,loy_16x16,hiy_16x16);
    std::cout << NBx_16x16 << "  " << NBy_16x16 << std::endl;
    for (int i = 0 ; i < NBx_16x16 ; i++) {
        for (int j = i ; j < NBy_16x16 ; j++) {
            int ii = i%NBx;
            int jj = j%NBy;
            int iB = int(i/NBx);
            int jB = int(j/NBy);
            double block_corrFact_ij = corrMat_16x16Block[iB][jB];
            //cout << i << " " << j << "   -   " << ii << " " << jj << "   -   " << iB << " " << jB << endl;
            if (iB < 8 && jB < 8) {
                value = corrMat_histo_6Det_1x1->GetBinContent(ii+1,jj+1)*block_corrFact_ij;
            }
            if (iB >= 8 && jB >= 8) {
                value = corrMat_histo_8Det_1x1->GetBinContent(ii+1,jj+1)*block_corrFact_ij;
            }
            if ((iB < 8 && jB >= 8) || (iB >= 8 && jB < 8)) {
                value = sqrtOffDiag_histo->GetBinContent(ii+1,jj+1)*block_corrFact_ij;
            }
            if (iB == 3 || iB == 7)
                value = 0.0;
            if (jB == 3 || jB == 7)
                value = 0.0;
            corrMat_histo_16x16->SetBinContent(i+1,j+1,value);
            corrMat_histo_16x16->SetBinContent(j+1,i+1,value);
            corrMat_mat_16x16(i,j) = value;
            corrMat_mat_16x16(j,i) = value;
        }
    }
    //-- Forced large-correlation matrix to have ones in the diagonal and scales off-diagonal elements.
    TMatrixD corrMat_mat_16x16_Scaled(NBx_16x16,NBy_16x16);
    TH2F *corrMat_histo_16x16_Scaled = new TH2F("corrMat_histo_16x16_Scaled","",NBx_16x16,lox_16x16,hix_16x16,NBy_16x16,loy_16x16,hiy_16x16);
    for (int i = 0 ; i < NBx_16x16 ; i++)
        for (int j = 0 ; j < NBy_16x16 ; j++){
            int ii = i%NBx;
            int jj = j%NBy;
            int iB = int(i/NBx);
            int jB = int(j/NBy);
            //corrMat_mat_16x16_Scaled(i,j) = corrMat_mat_16x16(i,j)/sqrt(corrMat_mat_16x16(i,i)*corrMat_mat_16x16(j,j));
            corrMat_mat_16x16_Scaled(i,j) = min(1.0,corrMat_mat_16x16(i,j)/sqrt(corrMat_mat_16x16(i,i)*corrMat_mat_16x16(j,j)));

            value = corrMat_mat_16x16_Scaled(i,j);
            if ((iB == 3 || iB == 7) && i != j)
                value = 0.0;
            if ((jB == 3 || jB == 7) && i != j)
                value = 0.0;
            corrMat_histo_16x16_Scaled->SetBinContent(i+1,j+1,value);
        }
    corrMat_histo_16x16->SetMinimum(-1.0);
    corrMat_histo_16x16_Scaled->SetMinimum(-1.0);
    
    std::cout << "MaxVal = " << corrMat_histo_16x16_Scaled->GetMaximum() << std::endl;

    //---------------------------------------------------------
    set_plot_style();
    
    //---------------------------------------------------------
    
    TCanvas *canv1 = new TCanvas("canv1","Blocks",100,100,2*500,2*500);
    canv1->Divide(2,2);
    canv1->cd(1);
    corrMat_histo_6Det_1x1->Draw("COLZ");
    canv1->cd(2);
    corrMat_histo_8Det_1x1->Draw("COLZ");
    canv1->cd(3);
    offDiagBlock_6DetX8Det_histo->Draw("COLZ");
    canv1->cd(4);
    sqrtOffDiag_histo->Draw("COLZ");

    
    TCanvas *canv2 = new TCanvas("canv2","one",200,100,500,500);
    canv2->cd(1);
    one_histo->Draw("COLZ");
    
    TCanvas *canv3 = new TCanvas("canv3","diagQ",200,100,700,700);
    canv3->cd(1);
    diagQ_histo->Draw("COLZ");
    
    //---------------------------------------------------------
    
    TCanvas *canv4 = new TCanvas("canv4","",900,100,700,700);
    canv4->cd(1);
    corrMat_histo_16x16_Scaled->Draw("COLZ");
    
    //---------------------------------------------------------
    // Write the rebinned correlation matrix to output file
    TFile *fout = new TFile("./db_CorrMatrix_6AD8AD_16x16Det.root","recreate");
    fout->cd();
    corrMat_histo_16x16_Scaled->Write();
    corrMat_mat_16x16_Scaled.Write();
    
    fout->Close();

     
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
*/
}// end

