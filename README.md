# neutrino1
Oscillation analysis of Daya-Bay data.
Files in this package are intended to be used to analyze data of the Daya Bay reactor antineutrino experiment. There are condes to do a rate-only analysis and others to make a spectral analysis.


# quick reference (How-To):

1. Execute root macro ldist.C:

  > root -b -n -q ldist.C

  Output: daya-bay-ldist.root
          ldist.pdf
          ldist_eh3.pdf
          ldist_6Det.pdf

2. Execute root macro db_ntuple.C (requires "PRL112_data.root")

  > root -b -n -q db_ntuple.C

   Output: db-ntuple.root

3. Execute root macro db_osc_rate.C (requires "rate_SurvProb.C")

  > root -b -n -q db_osc_rate.C

  Output: 
