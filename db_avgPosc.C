//For the Daya Bay rate-only analysis, using information from
//F.P. An et al., PRL 112 061801 (2014)
#include <constants.h>

void db_avgPosc()
{ //begin

    //------------- Style --------------
    gROOT->SetStyle("Plain");
    gStyle->SetOptStat(0);
    //------------- Style --------------
    //----------  Text Style  ---------
    //ft = 10 * fontID + precision
    Int_t ft = 10 * 4 + 2;
    Double_t sz = 0.04;
    //----------  Text Style  ---------

    //---------------------------------------------------
    // Open ntuple file to read simulated data
    TFile *fntuple = new TFile("db-ntuple.root","READ");
    TTree *T = (TTree*)fntuple->Get("T");
    //---------------------------------------------------
    // histogram binning for the simulated data
    double    NB = 26;
    double    lo = 0.7;
    double    hi = 12.0;
    double    xbins[27];
    xbins[0] = lo;
    double delta_bins2 = (7.3 - 1.3)/24.;// = 0.25 MeV/bin
    for (int i=0;i<(NB-1);i++)
    {
        xbins[i+1] = 1.3 + delta_bins2*i;
    }
    xbins[26] = hi;
    //---------------------------------------------------
    const int nAD = 6; //Number os Antineutrino detectors
    TH1F *Posc_AD_BF[nAD];
    for (int i = 0 ; i < nAD ; i++)
    {
        Posc_AD_BF[i] = new TH1F(Form("Posc_AD_BF_%d",i),"",1000,0,1);//to compute <POsc(BF)>
        Posc_AD_BF[i]->SetLineColor(i+1);
    }
    //---------------------------------------------------
    
    //Computing <POsc(s2t_BF,dm2_31)> for each AD
    // AD1 -> id = 0; AD2 -> id = 1; AD3 -> id = 2; AD4 -> id = 4; AD5 -> id = 5; AD6 -> id = 6
    /*
     const double s22th13_BF = 0.089;                   //PRL 112 061801 (2014)
     const double th13 = 0.5*asin(sqrt(s22th13_BF));
     const double c4th13 = (cos(th13))**4;
     const double dm2_31 = 2.32e-3; //eV^2,             //PRL 108 171803 (2012)
     const double dm2_21 = 7.59e-5; //eV^2,             //PRL 108 171803 (2012)
     const double s22th12 = 0.861;  //                  //PRL 108 171803 (2012)
     */
    
    //-- Best Fit values from PRL2014 (initial rate-only analysis)
    double s2t = 0.090;
    double dm2 = 2.32e-3;
    
    int sel;
    double avgPosc_AD[nAD]; //<POsc(s2t_BF,dm2_31)>
    double noOsc_IBDrate_perday[nAD];
    double integ;
    for (int iAD = 0 ; iAD < nAD ; iAD++)
    {
        if (iAD < 3) sel = iAD;
        else if (iAD >= 3) sel = iAD +1;
        
        T->Draw(Form("(1.0 - %e*((sin( 1.267 * %e * Ln/En ))**2) - (0.25*pow(1.0 + sqrt(1.0-%e),2)) * 0.861 * pow(sin( 1.267 * 7.59e-5 * Ln/En ),2)) >> Posc_AD_BF_%d",s2t,dm2,s2t,iAD),Form("id==%d",sel));
    
        avgPosc_AD[iAD] = Posc_AD_BF[iAD]->GetMean();

        cout << "(avgPosc_AD)_" << iAD+1 << " = (" << avgPosc_AD[iAD] << ") " << endl;
        
        integ = Posc_AD_BF[iAD]->Integral();
        Posc_AD_BF[iAD]->Scale(1.0/integ);
    }

    


} //end
