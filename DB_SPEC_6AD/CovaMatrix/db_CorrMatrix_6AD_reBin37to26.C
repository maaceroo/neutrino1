// This macro implements the Cholesky descomposition to obtain the
// covariance starting from the Correlation matrix a test matrix

//Constants definition
const int NBC = 37;         //Correlation matrix Number of bins
const int NBCreb = 26;      //Rebinned Correlation matrix Number of bins
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

void db_CorrMatrix_6AD_reBin37to26()
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

    TFile *fmat = new TFile("db_CorrMatrix_6AD_37bins.root");
    TH2F *CorMat_histo = (TH2F*)fmat->Get("CorMat_histo");
    TH2F *RedCorMat_histo = new TH2F("RedCorMat_histo","",NBC,0,NBC,NBC,0,NBC);
    cout << "Bins: " << CorMat_histo->GetXaxis()->GetNbins() << endl;
    //matrix definition

    cout << "P1. Initial matrix" << endl;
    double matrix[NBC][NBC];
    double scaleFactor = 1.0 - 0.10;

    for (int i = 0 ; i < NBC ; i++) {
        for (int j = 0 ; j < NBC ; j++) {
            matrix[i][j] = CorMat_histo->GetBinContent(i+1,j+1);
            if(i != j)
                matrix[i][j] = scaleFactor*matrix[i][j];
            RedCorMat_histo->SetBinContent(i+1,j+1,matrix[i][j]);
            //printf("%7.4f ",matrix[i][j]);
        }
        //printf("\n");
    }

    //-------------------------------------------
    //-- Cholesky: BEGIN-------------------------
    cout << "P2. Cholesky function BEGINS" << endl;
    double epsilon = 1.0e-8;
    double sumh2;
    double h[NBC][NBC];
    for (int i = 0 ; i < NBC ; i++) {
        for (int j = 0 ; j < NBC ; j++) {
            h[i][j] = 0.0;
        }
    }

    // Calculate the lower triangular matrix such that:
    // h x ( h^T) = matrix
    for ( int i = 0 ; i < NBC ; i++){
        //if ( fabs(matrix[0][0]) > epsilon ) {
        if ( matrix[0][0] != 0.0 ) {//like Alexis' code
            h[i][0] = matrix[i][0]/sqrt(matrix[0][0]);
        }
        else
            h[i][0] = 0.0;
    }
    
    for (int j = 1 ; j < NBC ; j++){
        sumh2 = 0.0;
        for (int jj = 0 ; jj < j ; jj++){
            sumh2 += h[j][jj]*h[j][jj];
        }
        if(matrix[j][j] - sumh2 < 0)
            h[j][j] = 0.0;
        else
            h[j][j] = sqrt(matrix[j][j] - sumh2);
        
        for (int i = j+1 ; i < NBC ; i++){
            sumh2 = 0.0;
            for ( int ii = 0 ; ii < j ; ii++){
                sumh2 += h[i][ii]*h[j][ii];
            }
            //if (fabs(h[j][j]) > epsilon) {
            if ( h[j][j] != 0.0) {//like Alexis' code
                h[i][j] = (matrix[i][j] - sumh2)/h[j][j];
            }
            else
                h[i][j] = 0.0;
        }
    }
    cout << "P3. Cholesky function ENDS" << endl;
    //-- Cholesky: END-------------------------
    //-----------------------------------------

/*    cout << "\n This is the result of the Cholesky Function: " << endl;
    for (int i = 0 ; i < NBC ; i++){
        for (int j = 0 ; j < NBC; j++)
            printf("%7.4f ",h[i][j]);
        printf("\n");
    }
    printf("\n");
*/
    //Triangular matrix (result of the Cholesky Function)
    TH2F *triangHisto = new TH2F("triangHisto","",NBC,0,NBC,NBC,0,NBC);
    for (int i = 0 ; i < NBC ; i++) {
        for (int j = 0 ; j < NBC; j++) {
            double val = 0.0;
            val = h[i][j];
            triangHisto->SetBinContent(i+1,j+1,val);
        }
    }

    // h x ( h^T) = matrix (where h is the Triangular matrix)
    TH2F *newHisto = new TH2F("newHisto","",NBC,0,NBC,NBC,0,NBC);
    for (int i = 0 ; i < NBC ; i++) {
        for (int j = 0 ; j < NBC; j++) {
            double val = 0.0;
            for (int k = 0 ; k < NBC ; k++) {
                val += h[i][k] * h[j][k];
            }
            if (i != j) {
                //val = val/scaleFactor;
            }
            //cout << "val(" << i << "," << j << ") = " << val << endl;
            newHisto->SetBinContent(i+1,j+1,val);
        }
    }
    
    //-----------------------------------------------------------------------------------
    //-- Building the new (rebinned) covariance matrix
    //-----------------------------------------------------------------------------------
    cout << "P4. Defining histograms" << endl;
    //define histogram
    double    lo = 0.7;
    double    hi = 12.0;
    //Bins for the rebinned matrix
    double xbins[NBCreb+1];
    xbins[0] = 0.7;
    double delta_bins2 = (7.3 - 1.3)/24; // 0.25 MeV/bin
    for (int i = 0 ; i < (NBCreb-1) ; i++)
        xbins[i+1] = 1.3 + delta_bins2*i;
    xbins[NBCreb] = hi;
    //Bins for the original matrix (as a test)
    double xbins1[NBC+1];
    xbins1[0] = 0.7;
    double delta_bins1 = (8.0 - 1.0)/35; // 0.2 MeV/bin
    for (int i = 0 ; i < (NBC-1) ; i++)
        xbins1[i+1] = 1.0 + delta_bins1*i;
    xbins1[NBC] = hi;

    //Normal distribution centered at cero and with mean 1
    cout << "P5. Normal distribution definition" << endl;
    TF1 *fGauss = new TF1("fGauss","exp(-0.5*x*x)",-8.0,8.0);
    double xvec[NBC]; //-Random-filled vector
    double vecFluc[NBC]; //-Fluctuated vector (from the modified Correlation Matrix)
    double vecFlucReb[NBCreb]; //-Rebinned fluctuated vector

    const int NumSamples = 5000;
    TH1F *histoArray[NumSamples]; //-Array of histograms to stores the rebined vecors
    TH1F *histoArray1[NumSamples]; //-Array of histograms to stores the original vecors

    cout << "P6. Random-filling of vectors/matrices" << endl;
    for (int m = 0 ; m < NumSamples ; m++) {
        histoArray[m] = new TH1F(Form("histoArray_%d",m),"",NBCreb,xbins);
        histoArray1[m] = new TH1F(Form("histoArray1_%d",m),"",NBC,xbins1);
        // Find n vectors drawn from Multivariate Gaussian distribution With Mean 0 and variance 1 in all dimensions.
        for (int i = 0 ; i < NBC ; i++) {
            xvec[i] = fGauss->GetRandom();
        }
        //Transform the Vector xvec into a fluctuated one
        for (int i = 0 ; i < NBC ; i++) {
            vecFluc[i] = 0.0;
            for (int j = 0 ; j < NBC ; j++) {
                vecFluc[i] += h[i][j]*xvec[j];
            }
        }
        //Filling the fluctuated-rebinned vectors
        vecFlucReb[0] = (vecFluc[0] + vecFluc[1] + 0.5*vecFluc[2])/2.5;
        for (int n = 0 ; n < 6 ; n++) {
            vecFlucReb[1+4*n] = (0.50*vecFluc[2+5*n] + 0.75*vecFluc[3+5*n])/1.25;
            vecFlucReb[2+4*n] = (0.25*vecFluc[3+5*n] + 1.00*vecFluc[4+5*n])/1.25;
            vecFlucReb[3+4*n] = (1.00*vecFluc[5+5*n] + 0.25*vecFluc[6+5*n])/1.25;
            vecFlucReb[4+4*n] = (0.75*vecFluc[6+5*n] + 0.50*vecFluc[7+5*n])/1.25;
        }
        vecFlucReb[25] = (0.5*vecFluc[32] + vecFluc[33] + vecFluc[34] + vecFluc[35] + vecFluc[36])/4.5;
        for (int ii = 0 ; ii < NBCreb ; ii++) {
            histoArray[m]->SetBinContent(ii+1,vecFlucReb[ii]);
        }

        for (int ii = 0 ; ii < NBC ; ii++) {
            histoArray1[m]->SetBinContent(ii+1,vecFluc[ii]);
        }
    }
    
    cout << "P7. Test-histograms to be drawn" << endl;
    double corrMatReb[NBCreb][NBCreb];
    TH2F *corrMatRebHisto = new TH2F("corrMatRebHisto","",NBCreb,0,NBCreb,NBCreb,0,NBCreb);
    for (int k = 0 ; k < NBCreb ; k++) {
        for (int l = 0 ; l < NBCreb; l++) {
            corrMatReb[k][l] = 0.0;
            for (int m = 0 ; m < NumSamples ; m++) {
                corrMatReb[k][l] += (histoArray[m]->GetBinContent(k+1) - 0.0)*(histoArray[m]->GetBinContent(l+1) - 0.0)/(NumSamples-1);
            }
            corrMatRebHisto->SetBinContent(k+1,l+1,corrMatReb[k][l]/scaleFactor);
            if (k == l) corrMatRebHisto->SetBinContent(k+1,l+1,1.0);
        }
    }
    
    //Creating the Covariance Matrix with root tools.
    TMatrixD corrMat_matrix(NBC,NBC);
    corrMat_matrix.Zero();
    for (int i = 0 ; i < NBC ; i++) {
        for (int j = 0 ; j < NBC; j++) {
            corrMat_matrix(i,j) = CorMat_histo->GetBinContent(i+1,j+1);
        }
    }
    
    //Creating the Rebinned Matrix with root tools.
    TMatrixD corrMatReb_matrix(NBCreb,NBCreb);
    corrMatReb_matrix.Zero();
    for (int i = 0 ; i < NBCreb ; i++) {
        for (int j = 0 ; j < NBCreb; j++) {
            corrMatReb_matrix(i,j) = corrMatRebHisto->GetBinContent(i+1,j+1);
        }
    }
    
    //Eigenvalues of rebined Matrix
    cout << "P8. Checking Eigenvalues of Rebined Matrix" << endl;
    TVectorD values;
    TMatrixD vectors = corrMatReb_matrix.EigenVectors(values);
    for (Int_t i = 0; i < values.GetNrows(); ++i) {
        cout << "eigen-value " << i << " is " << values(i) << " with eigen-vector" << endl;
    }

    TMatrixD invReb_mat(NBCreb,NBCreb);
    TMatrixD unoReb_mat(NBCreb,NBCreb);
    invReb_mat = corrMatReb_matrix;
    invReb_mat.Invert();
    unoReb_mat = corrMatReb_matrix*invReb_mat;
    //printf("\n Testing Unitarity \n");
    //invReb_mat.Print();
    TH2F *unoRebMat_histo = new TH2F("unoRebMat_histo","",NBCreb,0,NBCreb,NBCreb,0,NBCreb);
    for (int k = 0 ; k < NBCreb ; k++) {
        for (int l = 0 ; l < NBCreb ; l++) {
            unoRebMat_histo->SetBinContent(k+1,l+1,unoReb_mat(k,l));
        }
    }

    TH2F *invRebMat_histo = new TH2F("invRebMat_histo","",NBCreb,0,NBCreb,NBCreb,0,NBCreb);
    for (int k = 0 ; k < NBCreb ; k++) {
        for (int l = 0 ; l < NBCreb ; l++) {
            invRebMat_histo->SetBinContent(k+1,l+1,invReb_mat(k,l));
        }
    }

    double corrMatChk[NBC][NBC];
    TH2F *corrMatChkHisto = new TH2F("corrMatChkHisto","",NBC,0,NBC,NBC,0,NBC);
    for (int k = 0 ; k < NBC ; k++) {
        for (int l = 0 ; l < NBC; l++) {
            corrMatChk[k][l] = 0.0;
            for (int m = 0 ; m < NumSamples ; m++) {
                corrMatChk[k][l] += (histoArray1[m]->GetBinContent(k+1) - 0.0)*(histoArray1[m]->GetBinContent(l+1) - 0.0)/(NumSamples-1);
            }
            corrMatChkHisto->SetBinContent(k+1,l+1,corrMatChk[k][l]);
        }
    }

    cout << "P9. Drawing Section" << endl;
    // Drawing section
    //---------------------------------------------------------
    //- Set color scale of 2D histograms
    triangHisto->SetMaximum(1.0);
    triangHisto->SetMinimum(-1.0);

    newHisto->SetMaximum(1.0);
    newHisto->SetMinimum(-1.0);

    RedCorMat_histo->SetMaximum(1.0);
    RedCorMat_histo->SetMinimum(-1.0);

    corrMatRebHisto->SetMaximum(1.0);
    corrMatRebHisto->SetMinimum(-1.0);

    corrMatChkHisto->SetMaximum(1.0);
    corrMatChkHisto->SetMinimum(-1.0);
    
    invRebMat_histo->SetMaximum(1000.0);
    invRebMat_histo->SetMinimum(-1000.0);
    
    set_plot_style();


    TCanvas *canv0 = new TCanvas("canv0","",3*500,500);
    canv0->Divide(3,1);
    canv0->cd(1);
    RedCorMat_histo->Draw("COLZ");
    //CorMat_histo->Draw("COLZ");
    //newHisto->Draw("COLZ TEXT");
    canv0->cd(2);
    newHisto->Draw("COLZ");
    canv0->cd(3);
    corrMatChkHisto->Draw("COLZ");
    //---------------------------------------------------------
    //---------------------------------------------------------
    TCanvas *canv1 = new TCanvas("canv1","",500,500);
    triangHisto->Draw("COLZ");

    //Draw the new matrix!!
    //---------------------------------------------------------
    TCanvas *canv2 = new TCanvas("canv2","",500,500);
    corrMatRebHisto->Draw("COLZ");
    //corrMatChkHisto->Draw("COLZ");
    //------------------------------------------------------------------

    //Draw the original matrix!!
    //---------------------------------------------------------
    TCanvas *canv3 = new TCanvas("canv3","",500,500);
    CorMat_histo->Draw("COLZ");
    //------------------------------------------------------------------

    //Draw the identity matrix!!
    //---------------------------------------------------------
    TCanvas *canv4 = new TCanvas("canv4","Identity",3*500,500);
    canv4->Divide(3,1);
    canv4->cd(1);
    corrMatRebHisto->Draw("COLZ");
    canv4->cd(2);
    invRebMat_histo->Draw("COLZ");
    canv4->cd(3);
    unoRebMat_histo->Draw("COLZ");
    //------------------------------------------------------------------
    //canv0->Print("canv0.pdf");
    //canv1->Print("canv1.pdf");
    //canv2->Print("canv2.pdf");
    
    //---------------------------------------------------------
    // Write the rebinned correlation matrix to output file
    TFile *fout = new TFile("./db_CorrMatrix_6AD_26bins.root","recreate");
    fout->cd();
    corrMatRebHisto->Write();
    
    fout->Close();
    
    cout << "P10. Successfull end!" << endl;

}// end

