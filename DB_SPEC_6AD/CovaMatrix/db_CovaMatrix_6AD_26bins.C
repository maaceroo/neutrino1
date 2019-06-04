//Constants definition
const int NB = 26;      //Correlation matrix Number of bins
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
    //matrix definition
    
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
    //e_syst_histo1->SetFillColor(kRed);
    TH1F *e_syst_histo2 = new TH1F("e_syst_histo2","",NB,xbins);
    e_syst_histo2->SetLineColor(kRed);
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

}// end

