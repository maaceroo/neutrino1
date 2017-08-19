# neutrino1
Oscillation analysis of Daya-Bay data.
Files in this package are intended to be used to analyze data of the Daya Bay reactor antineutrino experiment. There are codes to do a rate-only analysis and others to make a spectral analysis.

############################
## Rate-Only analysis
# quick reference (How-To):
############################

1. Execute root macro ldist.C:

  > root -b -l -n -q ldist.C

  Output: daya-bay-ldist.root
          ldist.pdf (this is not really necessary)
          ldist_eh3.pdf (this is not really necessary)
          ldist_6Det.pdf (this is not really necessary)

2. Execute root macro db_ntuple.C (requires "PRL112_data.root")

  > root -b -l -n -q db_ntuple.C

   Output: db-ntuple.root

3. Execute root macro db_osc_rate.C (requires "rate_SurvProb.C")

  > root -b -l -n -q db_osc_rate.C

  Output: db_chi2_rate.txt
          POsc_avg.pdf (this is not really necessary)

