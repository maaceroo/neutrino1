# neutrino1
Oscillation analysis of Daya-Bay data.
Files in this package are intended to be used to analyze data of the Daya Bay reactor antineutrino experiment. There are codes to do a rate-only analysis and others to make a spectral analysis.

# Quick Reference (How-To):

## Building the Covariance Matrix

### Runnig with a single command line
Runnig a single root macro will execute all the necessary scripts to create the Covariance Matrix needed for the 6AD spectral Analysis:

    > root -l -n -b -q db_CovaMatrix_6AD_6x6Det_build.C

The `db_CorrMatrix_6AD_DiagBlk3_37bins.txt` file contains the digitalised correlation matris extracted from slides "Daya Bay Oscilla+on Analysis" by Henoch Wong (2016).

The final product ir the `root` file `db_CovaMatrix_6AD_6x26bins.root`, automatically saved to `../files_data/`, from it is called by the analysis macro.

### Runnig step-by-step
**1.** Execute (_requires `db_CorrMatrix_6AD_DiagBlk3_37bins.txt`_):

    > root -b -l -n -q db_CorrMatrix_6AD_37bins.C

_Output_:  
- db_CorrMatrix_6AD_37bins.root

**2.** Execute

    > root -b -l -n -q db_CorrMatrix_6AD_reBin37to26.C

_Output_:  
- db_CorrMatrix_6AD_26bins.root

**3.** Execute

    > root -b -l -n -q db_CovaMatrix_6AD_26bins.C

_Output_:
- db_CovaMatrix_6AD_26bins.root

**4.** Execute

    > root -b -l -n -q db_CovaMatrix_6AD_6x26bins.C

_Output_:  
- ../files_data/db_CovaMatrix_6AD_6x26bins.root

