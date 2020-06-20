//Constants to be used in the macro RENO_test_EnergyS_ntuple.C

const int nDet = 2;   // number of ad at Reno
const int nRea = 6;   // number of reactors

// histogram binning for the simulated data
const int    NB = 26; 
const double lo = 1.2;  
const double hi = 8.4;

//BF oscillation parameters from PRL121,201801 (2018) by RENO Coll.
double     ssq2th13RENO = 0.0896; // +/- 0.0048(stat) +/- 0.0047(syst)
double     dmsqeeRENO   = 2.68e-3;// +/- 0.12(stat)   +/- 0.07(syst) x 10^{-3} eV^2
double     ssq2th12RENO = 0.851;
double     dmsq21RENO   = 7.53e-5;// +/- 0.18 x 10^{-5} eV^2

//IBD rate (per day)from Table I, PRL121 2018
//double IBDrate_data[nDet][2] = { {470.53,0.51},{47.06,0.15} };
double IBDrate_data[nDet][2] = { {461.00,0.58},{44.82,0.18} };
