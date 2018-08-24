// This macro creates a plot of the DB-AD Correlation matrix as seen in
// "Daya Bay Oscilla+on Analysis [Pure covariance approach]"
// presentation by Henoch Wong (2016)
//
// The input color digitized matrix suffers from being not
// sem ipositive-definite, hence a process is imposed on it
// to force it to become PosDef 
// Refs.  
// [1] "Fixing a broken correlation matrix" Phil Joubert and Stephen Langdell
// https://www.nag.co.uk/IndustryArticles/fixing-a-broken-correlation-matrix.pdf
// [2] "Correcting non-positive definite correlation matrices" by Stef Marée, Netherlands)
    

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

void db_CorrMatrix_8AD_FixPosDef_37bins()
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
    //get data from digitized file
    
    TFile *fmat = new TFile("./db_CorrMatrix_8AD_DiagBlk1_37bins.root");
    TH2F *CorMat_histo = (TH2F*)fmat->Get("CorMat_histo");
    cout << "Bins: " << CorMat_histo->GetXaxis()->GetNbins() << endl;
     

    //histogram binnung

    const int NBx = 37;
    const int NBy = 37;

/*
    double lox = 0.0;
    double hix = 37.0;
    double xbins[38];
    for (int i = 0 ; i < NBx ; i++)
        xbins[i] = lox + delta_bins*i;
    xbins[37] = hix;
    //Y axis
    double    loy = 0.0;
    double    hiy = 37.0;
    double ybins[38];
    for (int i = 0 ; i < NBy ; i++)
        ybins[i] = loy + delta_bins*i;
    ybins[37] = hiy;

*/

/*
    int cols = 37;
    int rows = 37;

    double ** mat_array;
    mat_array = new double*[rows];
    for(int i = 0 ; i < rows ; i++)
        mat_array[i] = new double[cols];
    for(int j = 0 ; j < rows ; j++)
        for(int k = 0 ; k < cols ; k++)
            mat_array[j][k] = CorMat_histo->GetBinContent(j+1,k+1);

*/

//---------------------------------------------------------------------------------------------------

/*  may remove

    //define histogram
    double delta_bins = 1.0;
    //X axis
    double NBx = 37.0;
    double lox = 0.0;
    double hix = 37.0;
    double xbins[38];
    for (int i = 0 ; i < NBx ; i++)
        xbins[i] = lox + delta_bins*i;
    xbins[37] = hix;
    //Y axis
    double    NBy = 37.0;
    double    loy = 0.0;
    double    hiy = 37.0;
    double ybins[38];
    for (int i = 0 ; i < NBy ; i++)
        ybins[i] = loy + delta_bins*i;
    ybins[37] = hiy;



    //cout << "Histogram definition " << endl;
    TH2F *CorMat_histo = new TH2F("CorMat_histo","",NBx,xbins,NBy,ybins);
    CorMat_histo->SetMaximum(+1.0);
    CorMat_histo->SetMinimum(-1.0);
    set_plot_style();
    
    // set histogram contents
    for (int i = 0 ; i < NBx ; i++)
        for (int j = 0 ; j < NBy ; j++)
        {
            CorMat_histo->SetBinContent(i+1,j+1,mat_array[36-i][j]);
        }
    //cout << "Histogram has been filled " << endl;

*/

    CorMat_histo->SetMaximum(+1.0);
    CorMat_histo->SetMinimum(-1.0);
    set_plot_style();

    TMatrixD CorMat_mat(NBx,NBy);
    //TMatrixDSym CorMat_mat(NBx);
    CorMat_mat.Zero();
    for (int i = 0 ; i < NBx ; i++)
        for (int j = 0 ; j < NBy ; j++)
            CorMat_mat(i,j) = CorMat_histo->GetBinContent(i+1,j+1);
            //CorMat_mat(i) = CorMat_histo->GetBinContent(i+1,j+1);

    TVectorD values;
    TMatrixD vectors = CorMat_mat.EigenVectors(values);
    for (Int_t i = 0; i < values.GetNrows(); ++i) {
        cout << "eigen-value " << i << " is " << values(i) << " with eigen-vector" << endl;
    }
    //vectors.Print();
    
    TMatrixD vectors_I = vectors;
    vectors_I.T();
    
    TMatrixD temp(NBx,NBy);
    TMatrixD diagMat(NBx,NBy);
    diagMat.Zero();
    for (int i = 0 ; i < NBx ; i++) {
        //if(values(i) < 0.0 && abs(values(i)) < 1.0e-4)
        if(values(i) < 0.0)
            diagMat(i,i) = 0.0;
        else
            diagMat(i,i) = values(i);
    }
    temp = vectors*diagMat;

    TMatrixD fixedMat(NBx,NBy);
    TMatrixD fixedMatScaled(NBx,NBy);
    fixedMat = temp*vectors_I;

    for (int i = 0 ; i < NBx ; i++)
        for (int j = 0 ; j < NBy ; j++)
            fixedMatScaled(i,j) = fixedMat(i,j)/sqrt(fixedMat(i,i)*fixedMat(j,j));
    
    TVectorD values2;
    vectors = fixedMatScaled.EigenVectors(values2);
    cout << "First iteration" << endl;
    for (Int_t i = 0 ; i < values2.GetNrows(); ++i) {
        cout << "eigen-value " << i << " is " << values(i) << " " << values2(i) << endl;
    }
    fixedMatScaled.Print();

    //TH2F *fixedMat_histo = new TH2F("fixedMat_histo","",NBx,xbins,NBy,ybins);
    //fixedMat_histo->SetMaximum(+1.0);
    //fixedMat_histo->SetMinimum(-1.0);

   TH2F *fixedMat_histo = (TH2F*)CorMat_histo->Clone("fixedMat_histo");
   //fixedMat_histo->SetName("fixedMat_histo");
    for (int i = 0 ; i < NBx ; i++)
        for (int j = 0 ; j < NBy ; j++){
            fixedMat_histo->SetBinContent(i+1,j+1,fixedMatScaled(i,j));
            //CorMat_histo->SetBinContent(i+1,j+1,fixedMatScaled(i,j));
        }
    // Drawng section
    //---------------------------------------------------------

    //TH2F *frame = new TH2F("frame","",NBx,xbins,NBy,ybins);
    
    TCanvas *canv0 = new TCanvas("canv0","",2*740,740);
    //frame->Draw();
    canv0->Divide(2,1);
    canv0->cd(1);
    CorMat_histo->Draw("COLZ");

    gPad->SetTicks(1,1);
    canv0->cd(2);
    fixedMat_histo->Draw("COLZ");
    gPad->SetTicks(1,1);

    //canv0->Print("files_plots/canv_DB.pdf");
    //---------------------------------------------------------
    // write to output file
    TFile *fout = new TFile("./db_CorrMatrix_8AD_FixPosDef_37bins.root","recreate");
    fout->cd();
    
    fixedMat_histo->Write();

    fout->Close();

}// end

