//Physical constants

double pi = 3.141592653589793238;
double Mn = 939.565;        //Neutron mass in MeV
double Mp = 938.272;        //Proton mass in MeV
double DeltaMnp = 1.2933;   //MeV
double mel = 0.510998928;   //Electron mass in(MeV)

double avg_nRecoilE = 10.0e-3; //MeV
double avg_constE = 0.78; //MeV

// histogram binning for the simulated data
const int  NB = 26;
double  lo = 0.7;
double  hi = 12.0;

//Number of Antineutrino Detectors
const int nAD = 6;

//Number of Nuclear Reactors
const int nNR = 6;

//Fixed neutrino oscillations parameters
double   dm2_21 = 7.59e-5; //eV^2,                //PRL 108 171803 (2012)
double s22th_12 = 0.861;

//Grid of oscillation parameters
const int     N_s2t = atoi(getenv("NS2T"));
const int     N_dm2 = atoi(getenv("NDM2"));

double       lo_s2t = atof(getenv("LO_S2T"));
double       hi_s2t = atof(getenv("HI_S2T"));

double       lo_dm2 = atof(getenv("LO_DM2"));
double       hi_dm2 = atof(getenv("HI_DM2"));


