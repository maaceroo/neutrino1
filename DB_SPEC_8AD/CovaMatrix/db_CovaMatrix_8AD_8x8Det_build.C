void db_CovaMatrix_8AD_8x8Det_build()
{

    gROOT->ProcessLine(".x db_CorrMatrix_8AD_DiagBlk1_37bins.C");

    gROOT->ProcessLine(".x db_CorrMatrix_8AD_FixPosDef_37bins.C");
    
    gROOT->ProcessLine(".x db_CorrMatrix_8AD_reBin37to35.C");
    
    gROOT->ProcessLine(".x db_CovaMatrix_8AD_35bins.C");
    
    gROOT->ProcessLine(".x db_CovaMatrix_8AD_8x35bins.C");

}
