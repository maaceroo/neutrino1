//-- db_CovaMatrix_6AD8AD_16x16Det.C - By M.A. Acero and A.A. Aguilar-Arevao - 24.05.2019 --//
// This macro creates the correlation matrix necessary to perform the analysis
// of the Daya-Bay data from the 6AD and 8AD periods combined.
// The obtained correlation matrix is an approximate reproduction of the one
// presented by Henoch Wong's Slides at Joint DC-RENO-DYB workshop (2016)
// https://indico.snu.ac.kr/indico/event/4/session/25/contribution/32/material/slides/0.pdf
// We then use this correlation matrix together with the fraction error to generate
// the corresponding covariance matrix

//Constants definition
//const int NB = 35;      //Single detector Correlation matrix Number of bins
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

void db_CovaMatrix_6AD8AD_16x16Det()
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

    //--Correlation matrix
    TFile *fmat = new TFile("./db_CorrMatrix_6AD8AD_16x16Det.root");
    //-- Extract Scaled histo with diagonal elements scaled to one
    TH2F *corrMat_histo_16x16 = (TH2F*)fmat->Get("corrMat_histo_16x16_Scaled");
    corrMat_histo_16x16->SetTitle("(6AD+8AD) Blocks");
    std::cout << "Bins: " << corrMat_histo_16x16->GetXaxis()->GetNbins() << std::endl;

    const int Ngr = 26;
    double erxx[Ngr] = {0.975,1.375,1.625,1.875,2.125,2.375,2.625,2.875,3.125,3.375,
                        3.625,3.875,4.125,4.375,4.625,4.875,5.125,5.375,5.625,5.875,
                        6.125,6.375,6.625,6.875,7.500,10.0};
    
    //-- Information obtained from Fig. 2 of PRL 116, 061801 (2016) (last value is made up by us!)
    double eryy6AD[Ngr] =
    {
        1.12448,1.04615,1.03497,1.03030,1.02751,1.02657,1.02471,1.02378,1.02284,
        1.02378,1.02471,1.02657,1.02751,1.02937,1.03217,1.03497,1.03776,1.04149,
        1.04336,1.04615,1.05175,1.06014,1.06853,1.07786,1.24009,1.5
    };
    
    //--- Graph with errors from https://arxiv.org/pdf/1508.04233.pdf Fig. 2
    double eryy8AD[Ngr] =
    {
        1.0115,1.0064,1.0053,1.0046,1.0043,1.0039,1.0031,1.0031,1.0027,1.0027,
        1.0035,1.0027,1.0035,1.0039,1.0039,1.0042,1.0050,1.0058,1.0065,1.0080,
        1.0099,1.0122,1.0152,1.0186,1.0227,1.0642
    };

    //for(int ii = 0 ; ii < Ngr ; ii++)
    //eryy[ii] = eryy[ii]*1.04; //2019.03.19 - Testing a larger systematics error
    
    TGraph *errgr6AD = new TGraph(Ngr,erxx,eryy6AD);
    TGraph *errgr8AD = new TGraph(Ngr,erxx,eryy8AD);
    printf("f(3.5)=%f\n",errgr6AD->Eval(3.5));
    printf("f(3.5)=%f\n",errgr8AD->Eval(3.5));
    errgr6AD->SetLineStyle(1);
    errgr8AD->SetLineStyle(2);
    errgr6AD->Draw("AC");
    errgr8AD->Draw("SAME");

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
    
    TH1F *e_syst_histo6AD = new TH1F("e_syst_histo6AD","",NB,xbins);
    e_syst_histo6AD->SetLineColor(kRed);
    TH1F *e_syst_histo8AD = new TH1F("e_syst_histo8AD","",NB,xbins);
    e_syst_histo8AD->SetLineColor(kBlue);
    //---Fill e_syst_histo6AD and e_syst_histo8AD
    for (int i = 0 ; i < NB ; i++) {
        double xpoint6AD = e_syst_histo6AD->GetBinCenter(i+1);
        e_syst_histo6AD->SetBinContent(i+1,errgr6AD->Eval(xpoint6AD)-1.0);
        double xpoint8AD = e_syst_histo8AD->GetBinCenter(i+1);
        e_syst_histo8AD->SetBinContent(i+1,errgr8AD->Eval(xpoint8AD)-1.0);
    }
    TCanvas *canv1 = new TCanvas("canv1","one",200,100,700,500);
    canv1->cd();
    e_syst_histo6AD->Draw("h");
    e_syst_histo8AD->Draw("h same");
    

    //-- building the large (8*35)x(8*35) 6AD Covariance matrix
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
    double sigmai, sigmaj;
    TMatrixD covaMat_mat_16x16(NBx_16x16,NBy_16x16);
    TH2F *covaMat_histo_16x16 = new TH2F("covaMat_histo_16x16","",NBx_16x16,lox_16x16,hix_16x16,NBy_16x16,loy_16x16,hiy_16x16);
    std::cout << NBx_16x16 << "  " << NBy_16x16 << std::endl;
    for (int i = 0 ; i < NBx_16x16 ; i++) {
        for (int j = i ; j < NBy_16x16 ; j++) {
            int ii = i%NBx; // from 0 to 34
            int jj = j%NBy; // from 0 to 34
            int iB = int(i/NBx); // from 0 to 7
            int jB = int(j/NBy); // from 0 to 7
            //cout << i << " " << j << "   -   " << ii << " " << jj << "   -   " << iB << " " << jB << endl;
            if (iB < 8 && jB < 8) {
                sigmai = e_syst_histo6AD->GetBinContent(ii+1);
                sigmaj = e_syst_histo6AD->GetBinContent(jj+1);
            }
            if (iB >= 8 && jB >= 8) {
                sigmai = e_syst_histo8AD->GetBinContent(ii+1);
                sigmaj = e_syst_histo8AD->GetBinContent(jj+1);
            }
            if (iB < 8 && jB >= 8) {
                sigmai = e_syst_histo6AD->GetBinContent(ii+1);
                sigmaj = e_syst_histo8AD->GetBinContent(jj+1);
            }
            if (iB >= 8 && jB < 8) {
                sigmai = e_syst_histo8AD->GetBinContent(ii+1);
                sigmaj = e_syst_histo6AD->GetBinContent(jj+1);
            }
            double rho_ij = corrMat_histo_16x16->GetBinContent(i+1,j+1);
            value = rho_ij*sigmai*sigmaj;
            //if(value < 0)
                //std::cout << "val = " << value << std::endl;
            covaMat_histo_16x16->SetBinContent(i+1,j+1,value);
            covaMat_histo_16x16->SetBinContent(j+1,i+1,value);
            //corrMat_mat_16x16(i,j) = value;
            //corrMat_mat_16x16(j,i) = value;
        }
    }
    //set_plot_style();

    //covaMat_histo_16x16->SetMaximum( 1e-4);
    //covaMat_histo_16x16->SetMinimum(1e-8);

    TCanvas *canv2 = new TCanvas("canv2","one",200,100,2*700,700);
    canv2->Divide(2,1);
    canv2->cd(1);
    corrMat_histo_16x16->Draw("COLZ");
    canv2->cd(2);
    gPad->SetLogz(1);
    covaMat_histo_16x16->Draw("COLZ");
    
    // Write the correlation and covariant matrices to output file
    TFile *fout = new TFile("./db_CovaMatrix_6AD8AD_16x16Det.root","recreate");
    fout->cd();
    covaMat_histo_16x16->Write();
    corrMat_histo_16x16->Write();
    
    fout->Close();

}// end

