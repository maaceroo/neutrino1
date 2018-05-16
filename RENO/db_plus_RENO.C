void db_plus_RENO(){
  
  ofstream chi2_s2t_RDB;
  string chi2_s2t_RENO_plus_DB = "files/chi2_s2t_RENO_plus_DB.txt";
  chi2_s2t_RDB.open((chi2_s2t_RENO_plus_DB).c_str());
  
  TGraph *db_s2t_chi2_RATE  = new TGraph("../DB_RATE/files_data/db_s2t_chi2_RATE.txt","%lg  %lg","");
  TGraph *RENO_s2t_chi2_RATE  = new TGraph("files/RENO_s2t_chi2_RATE.txt","%lg  %lg","");

  //db_s2t_chi2_RATE->Draw("AP");
  //RENO_s2t_chi2_RATE->Draw("PE");
  double xmin_db = 0.00413567;
  double xmin_RENO = 0.00948247;
  double a = 0;
  double b = 0;
  double s2t;
  double chi2_RENO = 0,chi2_db = 0;
  double chi2_sum = 0;
  for (int j = 0 ; j < 200 ; j++)
    {
      
      chi2_RENO = (RENO_s2t_chi2_RATE->GetY()[j])+ xmin_RENO;
      chi2_db = (db_s2t_chi2_RATE->GetY()[j])+ xmin_db;
      s2t = db_s2t_chi2_RATE->GetX()[j];
      
      chi2_sum = chi2_RENO + chi2_db;
      chi2_s2t_RDB << s2t << "\t" << chi2_sum << endl;
      //cout << " chi2_RENO = " <<  chi2_RENO  << " chi2_db = " <<  chi2_db  << " chi2_sum = " << chi2_sum <<endl; 
      
    }

      int n;
    const int rows = 2;
    const int columns = 200;
    ifstream matriz("files/chi2_s2t_RENO_plus_DB.txt"); 
    double ** matr;  
    double minimo[rows];
    matr = new double*[rows];
    for(int k = 0 ; k < rows ; k++){
      matr[k] = new double[columns];
    }		
    
    for(int l = 0 ; l < columns ; l++){		
      for(int j = 0 ; j < rows ; j++){
	matriz >> matr[j][l];
      }
    }
    
    for(int i = 0 ; i < rows ; i++){
      minimo [i] = matr[i][0];
      for(int ll = 0 ; ll < columns ; ll++){
	if(matr[i][ll] < minimo[i]){
	  minimo[i] = matr[i][ll];
	  n =ll;
	}
      }
    }

    cout << "sin2t_13 = " << matr[0][n]   << " x2_min = " << matr[1][n]  << " n = " << n  << endl;
    /////////////////////////////////////////////
    

  
}
