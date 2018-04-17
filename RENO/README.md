# neutrino1
Oscillation analysis of RENO data.
Files in this package are intended to be used to analyze data of the Reactor Experiment of Neutrinos Oscillation (RENO). There are codes to do a rate-only analysis.

# Quick Reference (How-To):
## Rate-Only analysis (1D results)

**1.** Execute root macro renograph.C:

    > root -l -s -q renograph.C

_Output_:  
- files_root/RENOplots.root
- Plots/fardetector.pdf
- Plots/fardetectorbg.pdf 
- Plots/Neardetector.pdf 
- Plots/Neardetectorbg.pdf 
- Plots/Predictions.pdf 
- Plots/ratio_obs.pdf 
- Plots/ratio_spect.pdf 

**2.** Execute root macro ldist_2x6_RENO.C:

    > root -l -s -q ldist_2x6_RENO.C

_Output_:  
- file_root/ldist_RENO_2x6.root
- file_root/ldist_RENO_gen.root
- Plots/ldist_near_det.pdf 
- Plots/ldist_far_det.pdf 
- Plots/ldist.pdf

**3.** Execute root macro RENO_ntuple.C (_requires `renograph.root`and constant.h_)

    > root -l -s -q RENO_ntuple.C

_Output_:  
- files_root/RENO-ntuple.root

**4.** Execute macro RENO_osc_rate.C (_requires `RENO_ntuple.root`_)

    > root -b -l -n -q RENO_osc_rate.C

_Output_:
- Declaration of arrays (needed in Step 5): 
    - noOsc_IBDrate_perday[nAD][nNR] 
    - avgSinDelta21[nAD][nNR] 
    - avgSinDeltaee[nAD][nNR] 

**5.** Execute root macro RENO_minuit2.C [_NOTE: the number of grid points are hard-coded in the statements at Line 42 (`#define N_s2t  200`) and Line 47 (`#define N_eps  200`). Change those numbers at you prefer._]

    > root -b -l -n -q RENO_minuit2.C

_Output_:  
- files/chi2_s2t-a_surface_RATE.txt (_contains three columns: `sin^2(2th_13)`, `a`, `chi2`_)
- files/chi2_minimun_RATE.txt (_File containig three tab-separated values: `chi2_min  s2t_BF  a_BF`_)
- files/chi2_pullTerms_RATE.txt 

**6.** Compile and execute RENO_margin.cpp

6.1. Compile

     > g++ -o RENO_margin.exe db_margin.cpp

6.2. Execute: arguments are the number of grid points and the path to the file `chi2_s2t-a_surface_RATE.txt`

    > ./RENO_margin.exe 200 200 ./

_Output_:
- files/RENO_s2t_chi2_RATE.txt (_chi2 vs. s2th, where chi2 is marginalized over all the pull terms and a_)
- files/RENO_a_chi2_RATE.txt (_chi2 vs. a, where chi2 is marginalized over all the pull terms and s2t_)

**7.** Execute the gnuplot macro multi_plot_margin_RATE.gnu (_requires the three `.txt` second output files from step **5** and output files from step **6** _)

    > gnuplot multi_plot_margin_RATE_RENO.gnu

_Output_:  
- Plots/RENO_plots_RATE.eps (_Contour plot of s2t vs. a and marginalized chi2 plots_)


**8.** Execute root macro db_plus_RENO.C (_requires the two files. The first is took from Daya Bay analisys`db_s2t_chi2_RATE.txt` and the second one is took from RENO_s2t_chi2_RATE.txt _)

    > root -b -l -n -q db_plus_RENO.C

_Output_:  
- files/chi2_s2t_RENO_plus_DB.txt (_contains three columns: `sin^2(2th_13)`, `chi2_DB + chi2_RENO`__)

**9.** Execute the gnuplot macro plots_RATE_db_and_RENO.gnu (_requires three files `chi2_s2t_RENO_plus_DB.txt` `db_s2t_chi2_RATE.txt` and `db_s2t_chi2_RATE.txt`_)

    > gnuplot plots_RATE_db_and_RENO.gnu

_Output_:  
- Plots/plots_RATE_db_and_RENO.eps (_Plot of s2t vs. (`chi2_DB` , `chi2_RENO`and `chi2_DB + chi2_RENO`) in the same plot_)

