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
  string min_chi = "/data/numu_chi2_minumum.txt";

  ofstream min_file;
  min_file.open((path1 + min_chi).c_str());
  min_file << fixed << setprecision(10);
  
  int sn = atoi(argv[1]);
  int dmn= atoi(argv[2]);
  int cols, val;
  string dat_f = "/data/numu_chi2_s2t-dm2_surface-noFL.txt";
  ifstream data((path1 + dat_f).c_str());
  
  cols = COLS;
  //int fils = (sn+1)*(dmn+1);
  int fils = (sn)*(dmn);
  
  //double chi[fils][cols];
  double ** chi;
  chi = new double*[fils];
  for(int i = 0 ; i < fils ; i++)
    chi[i] = new double[cols];

  double x,min=C_MIN;
  
  for(int j=0;j<fils;j++)
    for(int k=0;k<cols;k++)
      data >> chi[j][k];
  
  for(int i=0;i<fils;i++)
    {
      x = chi[i][2];
      if(x<min)
	{
	  min = x;
	  val = i;
	}
    }
  cout << "The minimum Chi^2 value is   " << setprecision(6) << min << endl;
  cout << "For s2t =  " << chi[val][0] << ",  dmn =  " << chi[val][1] << endl;

  min_file << min << "  " << chi[val][0] << "  " << chi[val][1] << endl;

  data.close();
  min_file.close();

  return 0;
}
