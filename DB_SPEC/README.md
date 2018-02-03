# neutrino1
Oscillation analysis of Daya-Bay data.
Files in this package are intended to be used to analyze data of the Daya Bay reactor antineutrino experiment. There are codes to do a rate-only analysis and others to make a spectral analysis.

# Quick Reference (How-To):
## Rate-Only analysis (1D results)

**1.** Execute root macro ldist.C:

    > root -b -l -n -q ldist.C

_Output_:  
- daya-bay-ldist.root
- ldist.pdf
- ldist_6Det.pdf

**2.** Execute root macro db_ntuple.C (_requires `PRL112_data.root`_)

    > root -b -l -n -q db_ntuple.C

_Output_:  
- files_data/db-ntuple.root

**3.** Execute macro db_osc_spec.C (_requires `db_ntuple.root`_)

    > root -b -l -n -q db_osc_spec.C

_Output_:
- Declaration of arrays (needed in Step 4): 
    - noOsc_IBDrate_perday[nAD] 
    - avgSinDelta21[nAD] 
    - avgSinDelta31[nAD] 

**4.** Execute root macro db_minuit_spec.C [_NOTE: the number of grid points are hard-coded in the statements at Line 43 (`#define N_s2t  200`) and Line 48 (`#define N_dm2  200`). Change those numbers at you prefer._]

    > root -b -l -n -q db_minuit_spec.C

_Output_:  
- files_data/chi2_s2t-eps_surface_SPEC.txt (_contains three columns: `sin^2(2th)`, `Delta M^2`, `chi2`_)

**5.** Compile and execute db_chi2_min.cpp and db_margin.cpp

5.1. Compile

    > g++ -o db_chi2_min.exe db_chi2_min.cpp
    > g++ -o db_margin.exe db_margin.cpp

5.2. Execute: arguments are the number of grid points and the path to the file `chi2_s2t-dm2_surface_SPEC.txt`

    > ./db_chi2_min.exe 200 200 ./
    > ./db_margin.exe 200 200 ./

_Output_:
- files_data/chi2_minumum_SPEC.txt (_File containig three tab-separated values: `chi2_min  s2t_BF  epsilon_BF`_)
- files_data/db_s2t_chi2_SPEC.txt (_chi2 vs. s2th, where chi2 is marginalized over all the pull terms and dm2_)
- files_data/db_dm2_chi2_SPEC.txt (_chi2 vs. dm2, where chi2 is marginalized over all the pull terms and s2t_)

**6.** Execute the gnuplot macro multi_plot_margin_SPEC.gnu (_requires the three `.txt` output files from step **4**_)

    > gnuplot multi_plot_margin_SPEC.gnu

_Output_:  
- files_data/db_plots_SPEC.eps (_Contour plot of s2t vs. dm2 and marginalized chi2 plots_)
