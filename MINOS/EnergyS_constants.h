//Constants to be used in the macro MINOS_test_EnergyS_ntuple.C

// histogram binning for the data (PRL 110, 251801 (2013))  
const int NB_numu110      = 23; 
const double lo110        =  0.5;  
const double hi110        = 14.0 ;
const int NB_numubar110   = 12;
const double lob110       =  1.0;
const double hib110       = 14.0;
const int NB_numubarWS110 = 11;
const double lobWS110     =  0.0;
const double hibWS110     = 25.0;

//Path for output
//Directory to save files
const char *dirVar = getenv("JOBID");
string dirName = dirVar;
