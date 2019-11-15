void ldist()
{ // begin

    //------------- Style --------------
    gROOT->SetStyle("Plain");
    gStyle->SetOptStat(0);
    //------------- Style --------------


    const int nDet = 8; //Number of detectorr = 8
    const int nRea = 6; //Number of reactors  = 6
    
    //Baseline Distances (cm)
    const char *detNames[nDet] = {"EH1-AD1", "EH1-AD2", "EH2-AD1", "EH2-AD2",
          	                  "EH3-AD1", "EH3-AD2", "EH3-AD3", "EH3-AD4"};

    //From Nucl.Inst.Meth.Phys.Research A 811 (2016) 133–161 (Table 2)
    double baselines[nDet][nRea] = {
		{ 362.380, 371.763, 903.466, 817.158,1353.618,1265.315},
        { 357.940, 368.414, 903.347, 816.896,1354.229,1265.886},
		{1332.479,1358.148, 467.574, 489.577, 557.579, 499.207},
		{1337.429,1362.876, 472.971, 495.346, 558.707, 501.071},
		{1919.632,1894.337,1533.180,1533.628,1551.384,1524.940},
		{1917.519,1891.977,1534.919,1535.032,1554.767,1528.079},
		{1925.255,1899.861,1538.930,1539.468,1556.344,1530.078},
		{1923.149,1897.507,1540.667,1540.872,1559.721,1533.179}
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
    double massesDet[nDet] = {19941.0,19967.0,19891.0,19944.0,19917.0,19989.0,19892.0,19931.0};

    double Nprot[nDet];
    for (int i = 0 ; i < nDet ; i++) {
        Nprot[i] = massesDet[i]*protperKg;
    }
    
    //Individual Reactor (maximum) Thermal Power (GW)
    //double th_pow[nRea] = {2.9,2.9,2.9,2.9,2.9,2.9}; //This was used for the 6AD analysis
    //- Average thermal power as in Table I, PRD 95, 072006 (2017) for the 8AD period
    double th_pow[nRea] = {2.514,2.447,2.566,2.519,2.519,2.550};
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

    int nb = 48;
    double lo = 0;
    double hi = 48;
    TH1F *histo_ldist = new TH1F("histo_ldist","",nb,lo,hi);
    histo_ldist->Sumw2();
    TH1F *histo_ldist_eh3 = new TH1F("histo_ldist_eh3","",nb,lo,hi);
    TH1F *histo_ldist_6Det = new TH1F("histo_ldist_6Det","",nb,lo,hi);

    for (int id=0; id<nDet; id++){
        for (int ir=0; ir<nRea; ir++){

            int ii = id*nRea+ir;
            double wgt = massesDet[id]*th_pow[ir]/(pow(baselines[id][ir],2));

            //printf("%2d \t %2d %2d \t %s: %7.3f m %f \n",ii,id,ir,detNames[id],baselines[id][ir], wgt);
      
            histo_ldist->SetBinContent(ii+1,wgt);
            //printf("%2d \t %2d %2d \t %s: %7.3f m %f \n",ii,id,ir,detNames[id],baselines[id][ir], wgt);

            if (id>=4){
                histo_ldist_eh3->SetBinContent(ii+1,wgt);
            }
            if (id != 3 && id != 7){
                histo_ldist_6Det->SetBinContent(ii+1,wgt);
            }
        } //for ir
    } //for id
    
    double integ = histo_ldist->Integral();
    histo_ldist->Scale(1.0/integ);

    double integ_eh3 = histo_ldist_eh3->Integral();
    histo_ldist_eh3->Scale(1.0/integ_eh3);
    
    double integ_6Det = histo_ldist_6Det->Integral();
    histo_ldist_6Det->Scale(1.0/integ_6Det);
    
    TFile *fout = new TFile("files_data/daya-bay-ldist.root","recreate");
    fout->cd();
    histo_ldist->Write();
    histo_ldist_eh3->Write();
    histo_ldist_6Det->Write();

    //Test generation of baselines
    TH1F *histo_ldist_gen = new TH1F("histo_ldist_gen","",nb,lo,hi);
    histo_ldist_gen->SetMarkerStyle(8);
    histo_ldist_gen->SetMarkerSize(1.0);

    printf("\n");

    int Nevt = 10000000;
    for (int i=0 ; i<Nevt ; i++){
        int bl_idx = histo_ldist->GetRandom();
        histo_ldist_gen->Fill(bl_idx);

        int idet =  (bl_idx/nRea);
        int irea =  (bl_idx- idet*nRea);

        //printf("%2d \t %2d %2d\n",bl_idx, idet, irea);
    }
    double integ_gen = histo_ldist_gen->Integral();
    histo_ldist->Scale(integ_gen);

    TH1F *histo_ldist_eh3_gen = new TH1F("histo_ldist_eh3_gen","",nb,lo,hi);
    histo_ldist_eh3_gen->SetMarkerStyle(8);
    histo_ldist_eh3_gen->SetMarkerSize(1.0);

    Nevt=10000000;
    for (int i=0;i<Nevt;i++){
        int bl_idx = histo_ldist_eh3->GetRandom();
        histo_ldist_eh3_gen->Fill(bl_idx);

        int idet =  (bl_idx/nRea);
        int irea =  (bl_idx- idet*nRea);

        //printf("%2d \t %2d %2d\n",bl_idx, idet, irea);
    }
    double integ_eh3_gen = histo_ldist_eh3_gen->Integral();
    histo_ldist_eh3->Scale(integ_eh3_gen);

   // Drawing section

    TCanvas *canv0 = new TCanvas("canv0","",600,470);
    canv0->cd();

    histo_ldist_gen->Draw("PE");
    histo_ldist->Draw("hist same");
    histo_ldist_gen->Draw("PE same");

    canv0->Print("files_plots/ldist.pdf");

/*
    TCanvas *canv1 = new TCanvas("canv1","",600,470);
    canv1->cd();
    
    histo_ldist_eh3_gen->Draw("PE");
    histo_ldist_eh3->Draw("hist same");
    histo_ldist_eh3_gen->Draw("PE same");
    
    canv1->Print("files_plots/ldist_eh3.pdf");
*/  
} //end
