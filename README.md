# neutrino1
Oscillation analysis of Daya-Bay data.
Files in this package are intended to be used to analyze data of the Daya Bay reactor antineutrino experiment. There are codes to do a rate-only analysis and others to make a spectral analysis.

############################
# quick reference (How-To):
############################
## Rate-Only analysis (1D results)
############################

1. Execute root macro ldist.C:

  > root -b -l -n -q ldist.C

  Output: daya-bay-ldist.root
          ldist.pdf (this is not really needed)
          ldist_eh3.pdf (this is not really needed)
          ldist_6Det.pdf (this is not really needed)

2. Execute root macro db_ntuple.C (requires "PRL112_data.root")

  > root -b -l -n -q db_ntuple.C

   Output: db-ntuple.root

3. Execute root macro db_osc_rate.C (requires "rate_SurvProb.C")

  > root -b -l -n -q db_osc_rate.C

  Output: db_chi2_rate.txt
          POsc_avg.pdf (this is not really necessary)

  'db_chi2_rate.txt' contains three columns: sin^2(2th) Surv_Prob chi^2
  whic allows to plot the marginalized plot sin^2(2th) vs. chi^2, showing the BF for this oscillation parameter. You can use your prefered plotting software.

############################
## Rate-Only analysis (2D results)
############################

1. Execute root macro ldist.C:

   > root -b -l -n -q ldist_6x6.C

   Output: daya-bay-ldist_6x6.root
        ldist_6x6.pdf (this is not really needed)
        ldist_6x6.eps (this is not really needed)

2. Execute root macro db_ntuple.C (requires "PRL112_data.root")

   > root -b -l -n -q db_ntuple.C

   Output: db-ntuple.root

