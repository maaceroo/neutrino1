//-------------------------------------------------------------------------------------//
//--  ldist_2x6_RENO.C - By M.A. Acero O. & A.A. Aguilar-A. -D.J. Polo - 2017-10-02  --//
//-------------------------------------------------------------------------------------//
// This macro can be executed under ROOT typing                  	               //
// "root[0] .x ldist_2x6_RENO.C"                                                       //
// For the rate-only analysis, using information from:              		       //
// - F.P. An et al., arXiv:1003.1391 (2010)                           		       //
// - F.P. An et al., arXiv:1610.04326 (2017)                       		       //
//-------------------------------------------------------------------------------------//

void ldist_2x6_RENO()
{ // begin
  
  //------------- Style --------------
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  //------------- Style --------------
  
  
  const int nDet = 2; //Number of detectors = 2
  const int nRea = 6; //Number of reactors  = 6
    
  //Baseline Distances (m)
  char  *detNames[nDet] = {"Near-AD1", "Far-AD2"}; 
  
  //Information from article "arXiv:1003.1391" 6 Mar 2010. (Table 1.2)
  double baselines[nDet][nRea] = {
    {667.9,451.8,304.8,336.1,513.9,739.1},
    {1556.5,1456.2,1395.9,1381.3,1413.8,1490.1}
  };
  
  // Number of free protons in the Target. From  (16 Aug 2017) - arXiv:1610.04326
  double protperKg = 1.189e30;
  
  //Target Mass in each detector (kg). From "arXiv:1003.1391" 6 Mar 2010.
  double massesDet[nDet] = {16080.0,16080.0};
  
  double Nprot[nDet];
  for (int i = 0 ; i < nDet ; i++) {
    Nprot[i] = massesDet[i]*protperKg;
  }
  
  // Information From (16 Aug 2017) - arXiv:1610.04326
  //Individual Reactor (maximum) Thermal Power (GW)
  double th_pow[nRea] = {2.775,2.775,2.815,2.815,2.815,2.815};
  //Mean Fission fractions (U235,U238,Pu239,Pu241).
  double fissFrac[4] = {0.571,0.073,0.300,0.056};
  
  //Thermal Fission Energies in MeV (U235,U238,Pu239,Pu241). Russian Research Centre ”Kurchatov Institute”, Moscow, 0410100v1 (7 oct 2004)
  double ThFisEn[4] = {201.92,205.52,209.99,213.60};
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
    int nb = 12;
    double lo = 0;
    double hi = nb;
    TH1F *histo_ldist_RENO_2x6 = new TH1F("histo_ldist_RENO_2x6","",nb,lo,hi);
    histo_ldist_RENO_2x6->SetXTitle("EH-AD index");
    histo_ldist_RENO_2x6->SetYTitle("a. u.");
    
    for (int id=0; id<nDet; id++){
      for (int ir=0; ir<nRea; ir++){
	
	int ii = id*nRea+ir;
	
	double wgt = massesDet[id]*th_pow[ir]/baselines[id][ir]**2;
	
	histo_ldist_RENO_2x6->SetBinContent(ii+1,wgt);
        
      } //for ir
    } //for id
    
    double integ_2x6 = histo_ldist_RENO_2x6->Integral();
    histo_ldist_RENO_2x6->Scale(1.0/integ_2x6);
    
    TFile *fout = new TFile("files_root/ldist_RENO_2x6.root","recreate");
    fout->cd();
    histo_ldist_RENO_2x6->Write();
    
    TH1F *histo_ldist_gen = new TH1F("histo_ldist_gen","",nb,lo,hi);
    histo_ldist_gen->SetMarkerStyle(8);
    histo_ldist_gen->SetMarkerSize(1.0);
    
    printf("\n");
    
    int Nevt=5000000;
    for (int i=0;i<Nevt;i++){
      int bl_idx = (int*)histo_ldist_RENO_2x6->GetRandom();
      histo_ldist_gen->Fill(bl_idx);
      
      int idet = (int*) (bl_idx/nRea);
      int irea = (int*) (bl_idx- idet*nRea);
      
      //printf("%2d \t %2d %2d\n",bl_idx, idet, irea);
    }
    
    double integ_gen = histo_ldist_gen->Integral();
    histo_ldist_RENO_2x6->Scale(integ_gen);
    
    TFile *fout = new TFile("files_root/ldist_RENO_gen.root","recreate");
    fout->cd();
    histo_ldist_gen->Write();
    
    // Drawing section
    
    TCanvas *canv0 = new TCanvas("canv0","",600,470);
    canv0->cd();
    
    histo_ldist_gen->Draw("PE");
    histo_ldist_RENO_2x6->Draw("hist same");
    histo_ldist_gen->Draw("PE same");

    canv0->Print("Plots/ldist.pdf");
    TCanvas *canv1 = new TCanvas("canv1","",600,470);
    histo_ldist_RENO_2x6->Draw("hist");

    
} //end
