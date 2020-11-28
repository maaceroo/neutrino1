//Physical constants

double pi = 3.141592653589793238;
double Mn = 939.565;        //Neutron mass in MeV
double Mp = 938.272;        //Proton mass in MeV
double DeltaMnp = 1.2933;   //MeV
double mel = 0.510998928;   //Electron mass in(MeV)

double avg_nRecoilE = 10.0e-3; //MeV
double avg_constE = 0.78; //MeV

const int nDet = 2;   // number of ad at Reno
const int nRea = 6;   // number of reactors

// histogram binning for the simulated data
const int    NB = 26; 
const double lo = 1.2;  
const double hi = 8.4;

//Grid of oscillation parameters.
const int     N_s2t = atoi(getenv("NS2T"));
const int     N_dm2 = atoi(getenv("NDM2"));

double       lo_s2t = atof(getenv("LO_S2T"));
double       hi_s2t = atof(getenv("HI_S2T"));

double       lo_dm2 = atof(getenv("LO_DM2"));
double       hi_dm2 = atof(getenv("HI_DM2"));

//BF oscillation parameters from PRL121,201801 (2018) by RENO Coll.
double     ssq2th13RENO = 0.0896; // +/- 0.0048(stat) +/- 0.0047(syst)
double     dmsqeeRENO   = 2.68e-3;// +/- 0.12(stat)   +/- 0.07(syst) x 10^{-3} eV^2
double     ssq2th12RENO = 0.851;
double     dmsq21RENO   = 7.53e-5;// +/- 0.18 x 10^{-5} eV^2

double      fudge   = atof(getenv("fudge"));  //Adjust to Far/Near relative normalization
double      fFac1   = atof(getenv("fFac1"));  //Energy scale factor - Vertical displacement
double      fFac2   = atof(getenv("fFac2"));  //Energy scale factor - Vertical widening
double    resFac    = 1.00 ;  //Energy Resolution factor
double      eEscl   = 0.00 ;  //Energy scale factor for near det (alternative)

double      BFtoObs[2] = {1.0,1.0};  //Correction to the BF IBD Rate per day

//IBD rate (per day)from Table I, PRL121 2018
//double IBDrate_data[nDet][2] = { {470.53,0.51},{47.06,0.15} };
double IBDrate_data[nDet][2] = { {461.00,0.58},{44.82,0.18} };

//Path for output
//Directory to save files
const char *dirVar = getenv("JOBID");
string dirName = dirVar;
