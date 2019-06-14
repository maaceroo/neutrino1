void db_CovaMatrix_6AD8AD_16x16Det_build()
{

    gROOT->ProcessLine(".x db_CorrMatrix_6AD_reBin37to35.C");

    gROOT->ProcessLine(".x db_CorrMatrix_6AD8AD_16x16Det.C");
    
    gROOT->ProcessLine(".x db_CovaMatrix_6AD8AD_16x16Det.C");
    
}
