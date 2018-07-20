// This macro creates a plot of the DB-AD Correlation matrix as seen in
// "Daya Bay Oscilla+on Analysis [Pure covariance approach]"
// presentation by Henoch Wong (2016)
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

void db_CorrMatrix_8AD_DiagBlk1_37bins()
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
    //get information from Color-Scale file
    string fileScale = "./db_CorrMatrix_8AD_ColorScale_RGB_37bins.txt";
    ifstream sclfile((fileScale).c_str());
    int colSc = 4;
    int rowSc = 599;

    double ** scl_array;
    scl_array = new double*[rowSc];
    for(int i = 0 ; i < rowSc ; i++)
        scl_array[i] = new double[colSc];
    for(int j = 0 ; j < rowSc ; j++)
        for(int k = 0 ; k < colSc ; k++)
            sclfile >> scl_array[j][k];

    //get information from Color-Matrix file
    string fileMatrix = "./db_CorrMatrix_8AD_lowTriangle_RGB_37bins.txt";
    ifstream mtxfile((fileMatrix).c_str());
    int colMt = 37*3;
    int rowMt = 37;

    double ** mtx_array;
    mtx_array = new double*[rowMt];
    for(int i = 0 ; i < rowMt ; i++)
        mtx_array[i] = new double[colMt];
    for(int j = 0 ; j < rowMt ; j++)
        for(int k = 0 ; k < colMt ; k++)
            mtxfile >> mtx_array[j][k];
    int count = 0;
    /*
    for (int i = 0 ; i < colMt ; i++) {
        if (i%3 == 0){
            cout << "-----" << count << "-----" << endl;
            count++;
        }
        cout << mtx_array[28][i] <<  endl;
    }
     */
    //---------------------------------------------------------------------------------------------------

    //define histogram to store the resulting correlation matrix
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
    double NBy = 37.0;
    double loy = 0.0;
    double hiy = 37.0;
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
    double Red = 0.0;
    double Gre = 0.0;
    double Blu = 0.0;
    double diff = 0.0;
    int indx = 0;
    int theRow;
    for (int row = 0 ; row < rowMt ; row++) {
        for (int col = 0 ; col < 37 ; col++) {
            double mtxVal = 0.0;
            double diff_min = 1e6; // initialize to a large number for every pixel
            int col2 = col + (2*col);
            Red = mtx_array[row][col2];
            Gre = mtx_array[row][col2+1];
            Blu = mtx_array[row][col2+2];
            //cout << "R G B -> " << Red << " " << Blu << " " << Gre << endl;
            //if (Red == R_s && Gre == G_s && Blu == B_s) {
            //  indx += 1;
            //cout << indx << " The Correlation value is " << scl_array[sclR][3] << endl;
            //}
            for (int sclR = 0 ; sclR < rowSc ; sclR++) {
                double R_s = scl_array[sclR][0];
                double G_s = scl_array[sclR][1];
                double B_s = scl_array[sclR][2];
                //cout << sclR << "\t" << R_i << "\t" << G_i << "\t" << B_i << endl;
                diff = sqrt(pow(R_s-Red,2) + pow(G_s-Gre,2) + pow(B_s-Blu,2));
                if (diff < diff_min) {
                    diff_min = diff;
                    theRow = sclR;
                }// end if diff
            }//end for sclR
            indx += 1;
            mtxVal     = scl_array[theRow][3];
            if (Red==0 && Gre==0 && Blu==0) {
                mtxVal = 0.0;
                //cout << setprecision(4);
                //cout << "(row,col) = " << row << ", " << col << "\t color " << mtxVal << endl;
            }
            CorMat_histo->SetBinContent(37-row,col+1,mtxVal);
            //CorMat_histo->SetBinContent(row+1,37-col,mtxVal);
        }//end for col
        //cout << endl;
    }//end for row
    //break;

    for (int i = 1 ; i <= 37 ; i++) {
        CorMat_histo->SetBinContent(i,i,1.0);
    }
    for (int i = 0 ; i < rowMt ; i++) {
        for (int j = i ; j < rowMt ; j++) {
            //cout << i << " " << j << "   -   " << ii << " " << jj << endl;
            double value = CorMat_histo->GetBinContent(i+1,j+1);
            CorMat_histo->SetBinContent(i+1,j+1,value);
            CorMat_histo->SetBinContent(j+1,i+1,value);
        }
    }

    cout << "Congratulations! Histogram has been successfully filled!! " << endl;
    
    
    // Drawng section
    //---------------------------------------------------------
    TH2F *frame = new TH2F("frame","",NBx,xbins,NBy,ybins);
    
    int nCanv = 1;
    TCanvas *canv0 = new TCanvas("canv0","",nCanv*700,700);
    //frame->Draw();
    canv0->Divide(nCanv,1);
    canv0->cd(1);
    CorMat_histo->Draw("COLZ");
    gPad->SetTicks(1,1);

    //canv0->Print("files_plots/canv_DB.pdf");
    //---------------------------------------------------------
    // write to output file
    TFile *fout = new TFile("./db_CorrMatrix_8AD_DiagBlk1_37bins.root","recreate");
    fout->cd();
    
    CorMat_histo->Write();

    fout->Close();

}// end

