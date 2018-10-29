//------------------------------------------------------------------------------------//
//- RENO_ntuple.C - By M.A. Acero O. & A.A. Aguilar-A.  && D. J. Polo - 2018-07-06  --//
//------------------------------------------------------------------------------------//
// This macro can be executed under ROOT typing                        		      //
// "root[0] .x RENO_ntuple.C"                                            	      //
// For the spectral analysis, using information from:                  		      //
// - F.P. An et al.,  RENO Coll. (2017) - arXiv:1610.04326                            //
// - F.P. An et al., arXiv:1003.1391" 6 Mar 2010. (Table 1.2)              	      //
//------------------------------------------------------------------------------------//
// Simulated events, following the RENO energy spectra from the  		      //
// article 1610.04326.                                                                //
//------------------------------------------------------------------------------------//

#include<constants.h>

void RENO_ntuple_spect()
{ // begin
  
  //-------------------
  // Energy Histograms
  //-------------------
  TFile *fenergy = new TFile("files_root/RENOplots.root","read");
  
  //The histogram of near and far data spectra
  const int nd = 2; // Number of detectors
  TH1F *data_spect_histo[nd];
  for(int n=0;n<nd;n++){
    data_spect_histo[n] = (TH1F*) fenergy->Get(Form("data_spect_histo_%d",n));  // Get data
  }
  
  //-------------------
  // Distance Histogram
  //-------------------
  TFile *fpathl = new TFile("files_root/ldist_RENO_2x6.root","read");
  // Get histogram - histo_ldist_RENO_2x6
  TH1F *histo_ldist_RENO_2x6 = (TH1F*) fpathl->Get("histo_ldist_RENO_2x6");
  
  const int nDet = 2;   // number of ad at Reno 
  const int nRea = 6;   // number of reactors
  //Baseline Distances (m)
  char  *detNames[nDet] = {"near-AD1", "far-AD2"};
  
  //The 12 baselines in RENO ad*nr = 12 ()
  double baselines[nDet][nRea] =
    {
      {667.9,451.8,304.8,336.1,513.9,739.1},
      {1556.5,1456.2,1395.9,1381.3,1413.8,1490.1}
    };
  
  //make ntuple
  TFile *fout = new TFile("files_root/RENO-ntuple.root","RECREATE");
  TTree *T = new TTree("T","Monte Carlo neutrino events");
  
  float Ep, En, Ln;
  int   blid,ir,id,ad;
  T->Branch("Ep"  ,&Ep  ,"Ep/F");	  //prompt energy
  T->Branch("En"  ,&En  ,"En/F");	  //neutrino energy
  T->Branch("Ln"  ,&Ln  ,"Ln/F");	  //neutrino baseline
  T->Branch("blid",&blid,"blid/s"); //Baseline id (0-11)
  
  T->Branch("ir", &ir, "ir/s"); //reactor
  T->Branch("id", &id, "id/s"); //detector
  
  //int Nevents = 5000000; // CAUTION!! This must be commented out when using the script
  int Nevents = atoi(getenv("NTUPLE_EVENTS")); // This must be uncommented when using the script
  for (int i = 0 ; i < Nevents ; i++)
    {
      // generate a baseline (blid uniquely identifies the baseline)
      blid = (int*) histo_ldist_RENO_2x6->GetRandom();
      id =   (int*) (blid/nRea);
      ir =   (int*) (blid - id*nRea);
      Ln = baselines[id][ir];
      
      if(id==0)  ad=0;
      else if (id==1) ad=1;
      Ep = data_spect_histo[ad]->GetRandom();
      En = Ep + avg_nRecoilE + avg_constE;
      //En = Mn + Ep - Mp ; // Neutrino energy. Where Mp and Mn are the proton and neutron masses
      
      
      T->Fill();
    }
    
  fout->Write();
  
} // end
