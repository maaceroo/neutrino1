#!/bin/bash
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
#Define grid 
echo '=========================================='
echo '0) Define Grid'
echo '=========================================='
echo 

export NS2T=20
export NDM2=20

export LO_S2T=0.01
export HI_S2T=0.20

export LO_DM2=1.2e-3
export HI_DM2=3.5e-3

echo 'Grid in to be used:'
echo 'NS2T='$NS2T
echo 'LO_S2T='$LO_S2T
echo 'HI_S2T='$HI_S2T
echo
echo 'NDM2='$NDM2
echo 'LO_DM2='$LO_DM2
echo 'HI_DM2='$HI_DM2
echo


#-----------------------------------------------------------------------------
echo

# Construct L distribution
echo '=========================================='
echo '1) Running renograph.C'
echo '=========================================='
echo 
time root -b -l -n -q renograph.C

echo

#-----------------------------------------------------------------------------
echo

# Construct L distribution
echo '=========================================='
echo '2) Running ldist_2x6_RENO.C'
echo '=========================================='
echo 
time root -b -l -n -q ldist_2x6_RENO.C

echo

#-----------------------------------------------------------------------------
# Construct ntuple
echo '=========================================='
echo '3) Running RENO_ntuple_spect.C'
echo '=========================================='
echo 
export NTUPLE_EVENTS=1000000
echo $NTUPLE_EVENTS ntuple events
time root -b -l -n -q RENO_ntuple_spect.C

echo

#-----------------------------------------------------------------------------
# construct oscillated spectra for all points in the grid
echo '=========================================='
echo '4) Running RENO_osc_spect.C'
echo '=========================================='
echo
time root -b -l -n -q RENO_osc_spect.C

#-----------------------------------------------------------------------------
# run minimization
echo '=========================================='
echo '5) Running RENO_minuit_spect.C'
echo '=========================================='
echo
time root -b -l -n -q RENO_minuit_spect.C

echo

#-----------------------------------------------------------------------------
#Remove first line from file 
tail -n +2 files/chi2_s2t-dm2_surface_spect.txt > files/chi2_s2t-dm2_surface_spect-noFL.txt

#-----------------------------------------------------------------------------
#compile routines for minimization and marginalization
echo '=========================================='
echo 'compiling  RENO_margin_spect.cpp'
echo '=========================================='
echo
g++ -o RENO_margin_spect.exe RENO_margin_spect.cpp
#clang++ -o RENO_margin_spect.exe RENO_margin_spect.cpp

echo
#-----------------------------------------------------------------------------
echo '=========================================='
echo 'executing RENO_margin_spect.exe'
echo '=========================================='
echo
time ./RENO_margin_spect.exe $NS2T $NDM2 ./

echo

#-----------------------------------------------------------------------------
#Extract BF_CHI2, BF_S2T, BF_DM2 from chi2_minumum_SPEC.txt

read BF_S2T BF_DM2 BF_CHI2 <<< `cat files/chi2_minimun_spect.txt`

#----------------------------------------------------------------------------
#----------------------------------------------------------------------------
# Form gnuplot script
echo '=========================================='
echo 'Editting gnu plot script ...'
echo
sed -i "136s/.*/set label 35 '+' at $BF_S2T,$BF_DM2 center font 'CharterBT-Roman,15'/" multi_plot_margin_spect_RENO.gnu

sed -i "138s/.*/min = $BF_CHI2/" multi_plot_margin_spect_RENO.gnu

echo

#----------------------------------------------------------------------------
#Execute gnuplot script
echo '=========================================='
echo 'Runnign gnuplot macro'
echo '=========================================='
echo
gnuplot multi_plot_margin_spect_RENO.gnu

echo

#----------------------------------------------------------------------------
#Open in ghostview
gv Plots/RENO_plots_SPEC.eps &

#----------------------------------------------------------------------------
echo Done!
