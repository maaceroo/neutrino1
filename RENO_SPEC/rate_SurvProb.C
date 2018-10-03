double rate_SurvProb(double *x, double *par)
{//begin

    double s2t = x[0];
    
    double avgD21 = par[0];
    double avgDee = par[1];
    
    double s22th_12 = 0.846;
    
    double SProb = 1.0 - 0.25*((1.0 + sqrt(1.0 - s2t))**2)*s22th_12*avgD21 - s2t*avgDee;
    
    return SProb;
    
}//end
