# neutrino1
Oscillation analysis of Daya-Bay data.
Files in this package are intended to be used to analyze data of the Daya Bay reactor antineutrino experiment. There are codes to do a rate-only analysis and others to make a spectral analysis.

# Quick Reference (How-To):
## Rate-Only analysis (1D results)

**1.** Execute root macro ldist.C:

    > root -b -l -n -q ldist.C

_Output_:  
- daya-bay-ldist.root
- ldist.pdf (this is not really needed)
- ldist_6Det.pdf (this is not really needed)

**2.** Execute root macro db_ntuple.C (requires "PRL112_data.root")

    > root -b -l -n -q db_ntuple.C

_Output_:  
- files_data/db-ntuple.root

**3.** Execute root macro db_minuit.C [_NOTE_: the number of grid points are hard-coded in the statements at Line 43 (`#define N_s2t  200`) and Line 48 (`#define N_eps  200`)]

    > root -b -l -n -q db_minuit.C

_Output_:  
- files_data/chi2_s2t-eps_surface_RATE.txt (contains three columns: sin^2(2th), epsilon, chi^2).

**4.** Compile and execute db_chi2_min.cpp and db_margin.cpp

4.1. Compile

    > g++ -o db_chi2_min.exe db_chi2_min.cpp
    > g++ -o db_margin.exe db_margin.cpp

4.2. Execute: arguments are the number of grid points and the path to the file chi2_s2t-eps_surface_RATE.txt

    > ./db_chi2_min.exe 200 200 ./
    > ./db_margin.exe 200 200 ./

_Output_:
- files_data/chi2_minumum_RATE.txt (File containig three tab-separated values: chi2_min  s2t_BF  epsilon_BF) 
- files_data/db_s2t_chi2_RATE.txt (chi2 vs. s2th, where chi2 is marginalized over all the pull terms and epsilon)
- files_data/db_eps_chi2_RATE.txt (chi2 vs. epsilon, where chi2 is marginalized over all the pull terms and s2t)

**5.** Execute the gnuplot macro multi_plot_margin_RATE.gnu (requiresthe three .txt output files from step 4)

    > gnuplot multi_plot_margin_RATE.gnu

_Output_:  
- files_data/db_plots_RATE.eps (Contour plot of s2t vs. epsilon and marginalized chi2 plots)
