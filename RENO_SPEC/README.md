# neutrino1
Oscillation analysis of RENO data.
Files in this are package intended to be used to analyze data from the Reactor Experiment of Neutrinos Oscillation (RENO). These are codes to do a spectral analysis.

# Quick Reference (How-To):
## Spectral analysis (1D results)

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

**3.** Execute root macro RENO_ntuple_spect.C (_requires `renograph.root`and constant.h_)

    > root -l -s -q RENO_ntuple_spect.C

_Output_:  
- files_root/RENO-ntuple.root

**4.** Execute macro RENO_osc_spect.C (_requires `RENO_ntuple.root`_)  [_NOTE: the number of grid points are hard-coded in the statements at Line 185 (`#define N_s2t  100`) and Line 186 (`#define N_dm2  100`). Change those numbers at you prefer._]

    > root -b -l -n -q RENO_osc_spect.C

_Output_:
- RENO_gridOscSpectra_test.txt
- RENO_gridOscSpectra_Ratio.txt

**5.** Execute root macro RENO_minuit_spect.C  (_requires `RENOplots.root`, `RENO_ntuple.root` and `RENO_gridOscSpectra_test.txt`_) 

    > root -b -l -n -q RENO_minuit_spect.C

_Output_:  
- files/chi2_s2t-dm2_surface_spect.txt (_contains three columns: `sin^2(2th_13)`, `dm^2_ee`, `chi2`_)
- files/chi2_minimun_spect.txt (_File containig three tab-separated values: `s2t_BF   chi2_min   dm2_BF`_)
- files/chi2_pullTerms_spect.txt 

**6.** Compile and execute RENO_margin_spect.cpp

6.1. Compile

     > g++ -o RENO_margin_spect.exe RENO_margin_spect.cpp

6.2. Execute: arguments are the number of grid points and the path to the file `chi2_s2t-dm2_surface_SPEC.txt`

    > ./RENO_margin_spect.exe 100 100 ./

_Output_:
- files/RENO_s2t_chi2_SPEC.txt (_chi2 vs. s2th, where chi2 is marginalized over all the pull terms and a_)
- files/RENO_dm2_chi2_SPEC.txt (_chi2 vs. dm2, where chi2 is marginalized over all the pull terms and s2t_)

**7.** Execute the gnuplot macro multi_plot_margin_spect_RENO.gnu (_requires the three `.txt` first output files from step **5** and output files from step **6** _)

    > gnuplot multi_plot_margin_spect_RENO.gnu

_Output_:  
- Plots/RENO_plots_SPEC.eps (_Contour plot of s2t vs. dm2 and marginalized chi2 plots_)


