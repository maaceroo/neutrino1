//Physical constants

double pi = 3.141592653589793238;
double Mn = 939.565;        //Neutron mass in MeV
double Mp = 938.272;        //Proton mass in MeV
double DeltaMnp = 1.2933;   //MeV
double mel = 0.510998928;   //Electron mass in(MeV)

double avg_nRecoilE = 10.0e-3; //MeV
double avg_constE = 0.78; //MeV

//Number of Antineutrino Detectors
const int nAD = 8;
//Number of Nuclear Reactors
const int nNR = 6;

// histogram binning for the simulated data
const int  NB = 35;
double  lo = 0.7;
double  hi = 12.0;

//Grid of oscillation parameters.
const int     N_s2t = atoi(getenv("NS2T"));
const int     N_dm2 = atoi(getenv("NDM2"));

double       lo_s2t = atof(getenv("LO_S2T"));
double       hi_s2t = atof(getenv("HI_S2T"));

double       lo_dm2 = atof(getenv("LO_DM2"));
double       hi_dm2 = atof(getenv("HI_DM2"));

//Fixed neutrino oscillations parameters
//double   dm2_21 = 7.59e-5; //eV^2,       //PRL 108 171803 (2012)
//double s22th_12 = 0.861;
double   dm2_21 = 7.53e-5; //eV^2,         //PRD 95 072006 (2017) - DB1230Days
double s22th_12 = 0.846;

//double      fudge   = 1.00 ; //Adjust to Far/Near relative normalization
//double      fFac6AD   = 1.00 ;  //Energy scale factor -
//double      fFac8AD   = 1.00 ;  //Energy scale factor - 

double covMatAct = 1.0; // Set the CovMatirx On (1.0) or Off (0.0)
