/********************************************************************************/
/****** RENO_margin_spect.cpp -D. Polo T. M.A. Acero O. & - 09-27-2018     ******/
/*     Program to get the marginalization analysis for the oscillation          */
/*     parameters. It works for the anylisis of DayaBay data.                   */
/*     The progrma produces a pair of files (for each experiment) which         */
/*     have the data to make the maginalized Chi^2 plots for each of the        */
/*     oscillation parameters.                                                  */
/* ***              Executable file:  marg.exe ***                              */
/*                                                                              */
/********************************************************************************/

#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <vector>
#include <string>

using namespace std;

#define COL 3
#define MIN 2.0e6

//int main(void)
int main(int argc, char* argv[])
{
    if(argc < 4)
    {
        cout << "use:  " << endl;
        cout <<  argv[0] << " s2t(int) dm2(int) path_to_data_dir" << endl;
        exit(1);
    }
    // path to read/write data
    string path1 = argv[3];
    //string filePath = dirName;
    string chi_file = "/files/chi2_s2t-dm2_surface_spect-noFL.txt";
    string min_chi = "/files/chi2_minimun_spect.txt";
  
    ofstream file,filem;
    ifstream min1((path1 + chi_file).c_str());
    ifstream min2((path1 + min_chi).c_str());
  
    double * c_min;                    //Best Fit point
    c_min = new double[3];
    for(int j = 0 ; j < 3 ; j++)
      min2 >> c_min[j];

    int mal = atoi(argv[1]);
    int sal = atoi(argv[2]);
    
    //const int row = (sal+1)*(mal+1);
    const int row = (sal)*(mal);
    //double matrix[row][COL];
    double ** matrix;
    matrix = new double*[row];
    for(int i = 0 ; i < row ; i++)
      matrix[i] = new double[COL];
    
    double minchi = c_min[2];

    cout << minchi <<endl;
    cout << c_min[0] << " " << c_min[1] << " " << c_min[2] << endl;
    
    for(int l=0;l<row;l++)             //Chi2 data
      for(int m=0;m<COL;m++)
	min1 >> matrix[l][m];
    
    string s2t_file = "/files/RENO_s2t_chi2_SPEC.txt";
    string dmn_file = "/files/RENO_dm2_chi2_SPEC.txt";
    file.open((path1 + s2t_file).c_str());
    filem.open((path1 + dmn_file).c_str());
    
    double minimum = MIN;
    double compare, compare_m;
    int fix,fix_m;
  
    //const int lim = sal + 1;
    //const int limi = mal + 1;
    const int lim = sal;
    const int limi = mal;
    for(int k=0;k<lim;k++)
    {
        for(int i=(limi*k);i<limi*(k+1);i++)
        {
            compare = matrix[i][2];
            if(compare < minimum)
            {
                minimum = compare;
                fix = i;
            }
        }

    file << matrix[fix][0] << "  " << minimum-minchi << endl;
    minimum = MIN;
    }
  
    for(int k=0;k<limi;k++)
    {
        for(int i=k;i<limi*lim;i=i+limi)
        {
            compare_m = matrix[i][2];
            if(compare_m < minimum)
            {
                minimum = compare_m;
                fix_m = i;
            }
        }
        filem << matrix[fix_m][1] << "  " << minimum-minchi << endl;
        minimum = MIN;
    }
  
    min1.close();
    min2.close();
    file.close();
    filem.close();
    return 0;
    
}
