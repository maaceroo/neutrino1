//Constants definition
const int NB = 26;      //Correlation matrix Number of bins
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

void db_CovaMatrix_6AD_26bins()
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

    TFile *fmat = new TFile("db_CorrMatrix_6AD_26bins.root");
    TH2F *corrMat_histo = (TH2F*)fmat->Get("corrMatRebHisto");
    cout << "Bins: " << corrMat_histo->GetXaxis()->GetNbins() << endl;
    //6x6 Det. Correlation matrix - BEGIN
    double rho_array[6][6]={{1.0,0.0,0.0,0.0,0.0,0.0},
        {0.0,1.0,0.0,0.0,0.0,0.0},
        {0.0,0.0,1.0,0.0,0.0,0.0},
        {0.0,0.0,0.0,1.0,0.0,0.0},
        {0.0,0.0,0.0,0.0,1.0,0.0},
        {0.0,0.0,0.0,0.0,0.0,1.0}};
    
    //Construct correlation matrix upper triangle
    double corr_start = 0.9999;
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
        }
    }

    //-- Number of detectors
    int nDet = 6;
    int    NBx_6x6 = NB*nDet;
    double lox_6x6 = 0;
    double hix_6x6 = NB*nDet;
    int    NBy_6x6 = NB*nDet;
    double loy_6x6 = 0;
    double hiy_6x6 = NB*nDet;
    TH2F *corrMat_histo_6x6Det = new TH2F("corrMat_histo_6x6Det","",NBx_6x6,lox_6x6,hix_6x6,NBy_6x6,loy_6x6,hiy_6x6);
    cout << NBx_6x6 << "  " << NBy_6x6 << endl;
    for (int i = 0 ; i < NBx_6x6 ; i++) {
        for (int j = i ; j < NBy_6x6 ; j++) {
            int ii = i%NB;
            int jj = j%NB;
            int iB = int(i/NB);
            int jB = int(j/NB);
            double block_corrFact_ij = corrMat_6x6Block[iB][jB];
            //cout << i << " " << j << "   -   " << ii << " " << jj << "   -   " << iB << " " << jB << endl;
            double value = corrMat_histo->GetBinContent(ii+1,jj+1)*block_corrFact_ij;
            //double value = corrMat_histo->GetBinContent(ii+1,jj+1);
            //std::cout << value << " ";
            //if (jj==25) {
                //std::cout << std::endl;
            //}
            corrMat_histo_6x6Det->SetBinContent(i+1,j+1,value);
            corrMat_histo_6x6Det->SetBinContent(j+1,i+1,value);
        }
    }
    //6x6 Det. Correlation matrix - END

    //-- Information obtained from Fig. 2 of PRL 116, 061801 (2016)
    double e_syst[NB] = {1.124, 1.046, 1.035, 1.030, 1.028,
                        1.026, 1.025, 1.024, 1.023, 1.024,
                        1.025, 1.026, 1.028, 1.029, 1.032,
                        1.033, 1.038, 1.041, 1.043, 1.046,
                        1.052, 1.060, 1.068, 1.078, 1.240,
                        1.440};

    //define histogram
    double    lo = 0.7;
    double    hi = 12.0;
    //Bins for the rebinned matrix
    double xbins[NB+1];
    xbins[0] = 0.7;
    double delta_bins2 = (7.3 - 1.3)/24; // 0.25 MeV/bin
    for (int i = 0 ; i < (NB-1) ; i++)
        xbins[i+1] = 1.3 + delta_bins2*i;
    xbins[NB] = hi;

    TH1F *e_syst_histo1 = new TH1F("e_syst_histo1","",NB,xbins);
    e_syst_histo1->SetLineColor(kRed);
    e_syst_histo1->SetFillColorAlpha(kRed, 0.55);
    //e_syst_histo1->SetFillColor(kRed);
    TH1F *e_syst_histo2 = new TH1F("e_syst_histo2","",NB,xbins);
    e_syst_histo2->SetLineColor(kRed);
    e_syst_histo2->SetFillColorAlpha(kWhite, 1.0);
    //e_syst_histo2->SetFillColor(10);
    for (int i = 0 ; i < NB ; i++) {
        e_syst_histo1->SetBinContent(i+1,e_syst[i]);
        e_syst_histo2->SetBinContent(i+1,2.0-e_syst[i]);
    }

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

    TFile *outf = new TFile("./db_CovaMatrix_6AD_26bins.root","recreate");
    outf->cd();
    covaMat_histo->Write();
    covaMat_matrix.Write();
    outf->Close();
    //---------------------------------------------------------
    set_plot_style();

    TCanvas *canv0 = new TCanvas("canv0","canv0",2*500,500);
    canv0->Divide(2,1);
    canv0->cd(1);
    corrMat_histo->Draw("COLZ");
    canv0->cd(2);
    covaMat_histo->Draw("COLZ");
    //---------------------------------------------------------
    TCanvas *canv2 = new TCanvas("canv2","canv2",3*400,400);
    canv2->Divide(3,1);
    canv2->cd(1);
    covaMat_histo->Draw("COLZ");
    canv2->cd(2);
    invMat_histo->Draw("COLZ");
    canv2->cd(3);
    unoMat_histo->Draw("COLZ");
    //---------------------------------------------------------
    TH2F *e_frame = new TH2F("e_frame","",10,0.7,12,10,0.5,1.5);
    e_frame->GetXaxis()->SetTitle("Prompt Energy (MeV)");
    //e_frame->GetXaxis()->SetTitleFont(ft);
    //e_frame->GetXaxis()->SetTitleOffset(0.9);
    //e_frame->GetXaxis()->SetTitleSize(1.2*sz);
    //e_frame->GetXaxis()->SetLabelSize(1.2*sz);
    //e_frame->GetXaxis()->SetLabelFont(ft);
    //e_frame->GetYaxis()->SetLabelSize(1.2*sz);
    //e_frame->GetYaxis()->SetLabelFont(ft);
    TCanvas *canv1 = new TCanvas("canv1","canv1",700,400);
    canv1->cd();
    e_frame->Draw();
    e_syst_histo1->Draw("same hist");
    e_syst_histo2->Draw("same hist");
    gPad->Update();
    gPad->RedrawAxis();
    TLine* line1 = new TLine(gPad->GetUxmin(),1,gPad->GetUxmax(),1);
    line1->Draw();
    gPad->SetTicks(1,1);

    //---------------------------------------------------------
    corrMat_histo_6x6Det->SetMaximum(1.0);
    corrMat_histo_6x6Det->SetMinimum(-1.0);
    corrMat_histo_6x6Det->GetXaxis()->SetLabelSize(0.8*sz);
    corrMat_histo_6x6Det->GetXaxis()->SetLabelFont(ft);
    corrMat_histo_6x6Det->GetYaxis()->SetLabelSize(0.8*sz);
    corrMat_histo_6x6Det->GetYaxis()->SetLabelFont(ft);
    corrMat_histo_6x6Det->GetZaxis()->SetLabelSize(0.6*sz);
    corrMat_histo_6x6Det->GetZaxis()->SetLabelFont(ft);
    TCanvas *canv6x6 = new TCanvas("canv6x6","canv6x6",700,700);
    corrMat_histo_6x6Det->Draw("COLZ");
    gPad->SetTicks(1,1);
    
    canv6x6->Print("canv6x6.eps");

}// end

