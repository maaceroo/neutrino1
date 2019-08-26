//---------------------------------------------------------------------------------//
//--         combined_deltachi2_min.cpp - By: M.A. Acero-O. - 20.07-2019         --//
//---------------------------------------------------------------------------------//
//-- This code computes the minimum value of the chi^2 function for the combined --//
//-- DB+RENO data analysis, defining the combined chi^2-function as              --//
//--    chi^2_comb = (chi2_BD - chi2_DB_MIN) + (chi2_RENO - chi2_RENO_MIN)       --//
//-- The code also finds de minimum of this function and the corresponding values--//
//-- of the oscillation parameters.  At the same time, the marginalization to    --//
//-- those parameters is perfomed, and four txt files are created, as follows:   --//
//-- +files_data/chi2_minumum_COMB.txt -> the chi2_min and parameters            --//
//-- +files_data/chi2_s2t-dm2_surface_COMBINED.txt -> the chi2_combined surface  --//
//-- +files_data/db_s2t_chi2_COMBINED.txt -> marginalization wrt dm2             --//
//-- +files_data/db_dm2_chi2_COMBINED.txt -> marginalization wrt s2t             --//
//---------------------------------------------------------------------------------//
//-- To execute, compiled it and then run the .exe file:                         --//
//-- >g++ -o comb_chi2.exe combined_deltachi2_min.cpp                            --//
//-- >./comb_chi2.exe Ns2t Sdm2 ./                                               --//
//-- replacing Ns2t and Ndm2 for the numbers of the grid used for the anlyses.   --//
//---------------------------------------------------------------------------------//

#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <string>

using namespace std;

#define COLS 3
#define C_MIN 50000.0

//***************************************************************************//

//int main(void)
int main(int argc, char* argv[])
{
    if(argc < 4)
    {
        cout << "use: " << endl;
        cout << argv[0] << "  s2T(int) dm2(int) path_to_data_dir" << endl;
        exit(1);
    }
    // path to read/write data
    string path1 = argv[3];
    string min_chi   = "files_data/chi2_minumum_COMB_rate.txt";
    string delta_chi = "files_data/chi2_s2t-dm2_surface_COMBINED_rate.txt";

    ofstream min_file, comb_file;
    min_file.open((path1 + min_chi).c_str());
    min_file << fixed << setprecision(10);
    comb_file.open((path1 + delta_chi).c_str());
    comb_file << setprecision(6);

    int sn = atoi(argv[1]);
    int dmn= atoi(argv[2]);
    int cols, val;
    string db_dat_f = "files_data/chi2_s2t-eps_surface_RATE.txt";
    string rn_dat_f = "../RENO/files/chi2_s2t-a_surface_RATE.txt";
    ifstream db_data((path1 + db_dat_f).c_str());
    ifstream rn_data((path1 + rn_dat_f).c_str());
    string db_chi2M = "files_data/chi2_minumum_RATE.txt";
    string rn_chi2M = "../RENO/files/chi2_minimum_RATE.txt";
    ifstream db_min((path1 + db_chi2M).c_str());
    ifstream rn_min((path1 + rn_chi2M).c_str());
    
    cols = COLS;
    int fils = (sn)*(dmn);

    double db_chi2min[COLS] = {0.0};
    double rn_chi2min[COLS] = {0.0};
    for (int i = 0 ; i < COLS ; i++) {
        db_min >> db_chi2min[i];
        rn_min >> rn_chi2min[i];
    }
    
    std::cout << db_chi2min[0] << "\t" << db_chi2min[1] << "\t" << db_chi2min[2] << std::endl;
    std::cout << rn_chi2min[2] << "\t" << rn_chi2min[0] << "\t" << rn_chi2min[1] << std::endl;

    //double chi[fils][cols];
    double ** db_chi;
    double ** rn_chi;
    double ** comb_chi;
    db_chi = new double*[fils];
    rn_chi = new double*[fils];
    comb_chi = new double*[fils];
    for(int i = 0 ; i < fils ; i++){
        db_chi[i] = new double[cols];
        rn_chi[i] = new double[cols];
        comb_chi[i] = new double[cols];
    }

    double x,min=C_MIN;
  
    //-- BEGIN - Getteing the DB and RENO results and building the combined chi2-function
    for(int j = 0 ; j < fils ; j++){
        for(int k = 0 ; k < cols ; k ++){
            db_data >> db_chi[j][k];
            rn_data >> rn_chi[j][k];
            if (k < 2) {
                comb_chi[j][k] = db_chi[j][k];
            }
            else
                comb_chi[j][k] = (db_chi[j][k] - db_chi2min[0]) + (rn_chi[j][k] - rn_chi2min[2]);
        }
        //std::cout << comb_chi[j][0] << "\t" << comb_chi[j][1] << "\t" << comb_chi[j][2] << endl;
    }
    //-- END - Getteing the DB and RENO results and building the combined chi2-function

    //-- BEGIN - Finding the minimum of the combined chi2-function. Printing files
    for(int i=0;i<fils;i++)
        {
            x = comb_chi[i][2];
            if(x<min)
            {
                min = x;
                val = i;
            }
        }
    cout << "The minimum Chi^2 value is   " << setprecision(6) << min << endl;
    cout << "For s2t =  " << comb_chi[val][0] << ",  eps =  " << comb_chi[val][1] << endl;

    min_file << min << "  " << comb_chi[val][0] << "  " << comb_chi[val][1] << endl;
    
    for (int k = 0 ; k < fils ; k++) {
        for (int j = 0 ; j < cols ; j++) {
            comb_file << comb_chi[k][j] << "\t";
        }
        if (comb_chi[k][1] == 3.5e-3) {
            comb_file << endl;
        }
        comb_file << endl;
    }

    db_data.close();
    rn_data.close();
    min_file.close();
    //-- END - Finding the minimum of the combined chi2-function. Printing files

    //-- BEGIN - Marginalization of the combined chi2
    ofstream fileS,fileM;
    string s2t_file = "files_data/db_s2t_chi2_COMBINED_rate.txt";
    string dmn_file = "files_data/db_eps_chi2_COMBINED_rate.txt";
    fileS.open((path1 + s2t_file).c_str());
    fileM.open((path1 + dmn_file).c_str());
    
    double minimum = C_MIN;
    double compare_s, compare_m;
    int fix_s,fix_m;
    
    const int lim = atoi(argv[1]);
    const int limi = atoi(argv[2]);
    for(int k = 0 ; k < lim ; k++)
    {
        for(int i = (limi*k) ; i < limi*(k+1) ; i++)
        {
            compare_s = comb_chi[i][2];
            if(compare_s < minimum)
            {
                minimum = compare_s;
                fix_s = i;
            }
        }
        
        fileS << comb_chi[fix_s][0] << "  " << minimum - min << endl;
        minimum = C_MIN;
    }
    
    for(int k = 0 ; k < limi ; k++)
    {
        for(int i = k ; i < limi*lim ; i = i+limi)
        {
            compare_m = comb_chi[i][2];
            if(compare_m < minimum)
            {
                minimum = compare_m;
                fix_m = i;
            }
        }
        fileM << comb_chi[fix_m][1] << "  " << minimum - min << endl;
        minimum = C_MIN;
    }
    
    fileS.close();
    fileM.close();
    //-- END - Marginalization of the combined chi2

    return 0;
}
