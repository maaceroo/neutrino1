double rate_SurvProb(double *x, double *par)
{//begin

    double s2t = x[0];
    
    double avgD21 = par[0];
    double avgD31 = par[1];
    
    double s22th_12 = 0.861;
    
    double SProb = 1.0 - 0.25*(pow(1.0 + sqrt(1.0 - s2t),2))*s22th_12*avgD21 - s2t*avgD31;
    
    return SProb;
    
}//end
