# Daya Bay Spectral Analysis - 1230 days
Oscillation analysis of Daya-Bay experiment.
Files in this package are intended to be used to analyze data of the Daya Bay reactor antineutrino experiment.
Data are taken from the 6AD (PRL 112, 061801 (2014)) and 8AD (PRD 95, 072006 (2017)) periods in order to perform a full data analysis

Th first step to be done is taking the data an put them in an apporpiate format (histograms), separated for the different DAQ time periods

To do so, execute
> root -b -n -l -q db_PRL112_26to35_rebin.C

This macro should create the file _histos_6AD217Days_8AD1913Days_data.root_ which contains the necessary information to be used by the other macros.

This analysis is based on a Covariance Matrix which is to be build inside the directory _CombCovaMat_ by executing 
> root -b -n -l -q db_CovaMatrix_6AD8AD_16x16Det_build.C

The process will create the necessary _root_ files containing the Covariance Matrix.

### Runnig with a single command line
Then, a simplified way to run the complete analysis can be performed by executing the shell script _db_osc_spec.sh_:

> source db_osc_spec.sh

This will go step-by-step running each macro (for 6AD and 8AD) to generate the necessary information for the analysis, as described bellow. You may change some of the parameters in this file:

- NTUPLE_EVENTS: this is the number of MC simulated events. Note that 1e6 is good enough, but it is set to 5e6.
- Grid Parameters:
- (NS2T,NDM2) defines the number of points to be used in the parameter space (mixing angle, mass-squared difference). The larger the number, the finer the grid (more points in the parameter space), the longer the time to execute.
- (LO_S2T, HI_S2T) and (LO_DM2, HI_DM2) difine the lower and higher values to be used for the two mixing parameters to be fitted.
