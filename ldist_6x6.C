//---------------------------------------------------------------------//
//--  ldist_6x6.C - By M.A. Acero O. & A.A. Aguilar-A. - 2016-01-02  --//
//---------------------------------------------------------------------//
// This macro can be executed under ROOT typing                        //
// "root[0] .x ldist_6x6.C"                                            //
// For the rate-only analysis, using information from:                 //
// - F.P. An et al., PRL 112 061801 (2014)                             //
// - F.P. An et al., Chin.Phys.C37 011001 (2013)                       //
//---------------------------------------------------------------------//

void ldist_6x6()
{ // begin

    //------------- Style --------------
    gROOT->SetStyle("Plain");
    gStyle->SetOptStat(0);
    //------------- Style --------------


    const int nDet = 6; //Number of detectorr = 6
    const int nRea = 6; //Number of reactors  = 6
    
    //Baseline Distances (cm)
    char  *detNames[nDet] = {"EH1-AD1", "EH1-AD2", "EH2-AD1", /*"EH2-AD2",*/
                             "EH3-AD1", "EH3-AD2", "EH3-AD3"/*, "EH3-AD4"*/};

    //From Nucl.Inst.Meth.Phys.Research A 811 (2016) 133–161 (Table 2)
    double baselines[nDet][nRea] = {
		{362.380,371.763,903.466,817.158,1353.618,1265.315},
        {357.940,368.414,903.347,816.896,1354.229,1265.886},
		{1332.479,1358.148,467.574,489.577,557.579,499.207},
//		{1337.429,1362.876,472.971,495.346,558.707,501.071},
		{1919.632,1894.337,1533.180,1533.628,1551.384,1524.940},
		{1917.519,1891.977,1534.919,1535.032,1554.767,1528.079},
		{1925.255,1899.861,1538.930,1539.468,1556.344,1530.078}//,
//		{1923.149,1897.507,1540.667,1540.872,1559.721,1533.179}
		};
/*
    for (int i = 0; i < nDet; i++) {
        cout << detNames[i];
        for (int j = 0; j < nRea; j++) {
            cout << " " << baselines[i][j] << "  ";
        }
        cout << endl;
    }
*/
    //Number of protons in the GdLS per kg. From JINST 8 P09015 (2013) - arXiv:1307.1089
    double protperKg = 7.169e25;

    //GdLS Mass in each detector (kg). From Nucl.Inst.Meth.Phys.Research A 811 (2016) 133–161 (Table 10)
    double massesDet[nDet] = {19941.0,19967.0,19891.0,/*19944.0,*/19917.0,19989.0,19892.0/*,19931.0*/};

    double Nprot[nDet];
    for (int i = 0 ; i < nDet ; i++) {
        Nprot[i] = massesDet[i]*protperKg;
    }
    
    //Individual Reactor (maximum) Thermal Power (GW)
    double th_pow[nRea] = {2.9,2.9,2.9,2.9,2.9,2.9};
    //Mean Fission fractions (U235,U238,Pu239,Pu241). From PRL 112, 061801 (2014)
    double fissFrac[4] = {0.573,0.076,0.301,0.050};
    //Thermal Fission Energies in MeV (U235,U238,Pu239,Pu241). From Phys.Rev.C88, 014605 (2013)
    double ThFisEn[4] = {202.36,205.99,211.12,214.26};
    double AvMeVperFiss = 0.0;
    for (int i = 0 ; i < 4 ; i++) {
        AvMeVperFiss += fissFrac[i]*ThFisEn[i];
    }
    double EperFis = AvMeVperFiss*1.0e6*1.6e-19*1.0e-9; //(Energy [GJ] per fission)
    //Individual Reactor fission rate
    double fisrate[nRea];
    for (int i = 0 ; i < nRea ; i++) {
        fisrate[i] = th_pow[i]/EperFis; //Number of fission per second
        //cout << "Fission rate reactor " << i << ": " << fisrate[i] << " s^-1" << endl;
    }
    //break;

//    int nb = 48;
    int nb = 36;
    double lo = 0;
    double hi = nb;
    TH1F *histo_ldist_6x6 = new TH1F("histo_ldist_6x6","",nb,lo,hi);
    histo_ldist_6x6->SetXTitle("EH-AD index");
    histo_ldist_6x6->SetYTitle("a. u.");

    for (int id=0; id<nDet; id++){
        for (int ir=0; ir<nRea; ir++){

            int ii = id*nRea+ir;

            double wgt = massesDet[id]*th_pow[ir]/baselines[id][ir]**2;

            histo_ldist_6x6->SetBinContent(ii+1,wgt);
            
        } //for ir
    } //for id
    
    double integ_6x6 = histo_ldist_6x6->Integral();
    histo_ldist_6x6->Scale(1.0/integ_6x6);
    
    TFile *fout = new TFile("daya-bay-ldist_6x6.root","recreate");
    fout->cd();
    histo_ldist_6x6->Write();

    printf("\n");

   // Drawing section

    TCanvas *canv2 = new TCanvas("canv2","",600,470);
    canv2->cd();
    
    histo_ldist_6x6->Draw("hist");
    
    canv2->Print("ldist_6x6.pdf");
    canv2->Print("ldist_6x6.eps");
    
} //end
