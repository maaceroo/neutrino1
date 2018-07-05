void db_CovaMatrix_6AD_6x6Det_build()
{
    gROOT->ProcessLine(".x db_CorrMatrix_6AD_37bins.C");
    
    gROOT->ProcessLine(".x db_CorrMatrix_6AD_reBin37to26.C");
    
    gROOT->ProcessLine(".x db_CovaMatrix_6AD_26bins.C");
    
    gROOT->ProcessLine(".x db_CovaMatrix_6AD_6x26bins.C");

}
