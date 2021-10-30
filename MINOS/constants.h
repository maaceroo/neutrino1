//Physical constants

double pi = 3.141592653589793238;
double Mn = 939.565;        //Neutron mass in MeV
double Mp = 938.272;        //Proton mass in MeV
double DeltaMnp = 1.2933;   //MeV
double mel = 0.510998928;   //Electron mass in(MeV)

// histogram binning for the data (PRL 110, 251801 (2013))
const int NB_numu110    = 23;
double    lo110         =  0.5;
double    hi110         = 14.0;
const int NB_numubar110 = 12;
double    lob110        = 1.0;
double    hib110        = 14.0;

//Fixed neutrino oscillations parameters
double  dm2_21 = 7.54e-5; //eV^2, //PRL 112 191801 (2014)
double s2th_12 = 0.307;
double s2th_13 = 0.0242; //+-0.0025 Gaussian Penalty

//Grid of oscillation parameters
const int     N_s2t = atoi(getenv("NS2T"));
const int     N_dm2 = atoi(getenv("NDM2"));

double       lo_s2t = atof(getenv("LO_S2T"));
double       hi_s2t = atof(getenv("HI_S2T"));

double       lo_dm2 = atof(getenv("LO_DM2"));
double       hi_dm2 = atof(getenv("HI_DM2"));


//Path for output
//Directory to save files
const char *dirVar = getenv("JOBID");
string dirName = dirVar;
