#!/bin/bash
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
echo

# Construct L distribution
echo '=========================================='
echo '1) Running ldist.C'
echo '=========================================='
echo 
root -b -l -n -q ldist.C

echo

#-----------------------------------------------------------------------------
# Construct ntuple
echo '=========================================='
echo '2) Running ntuple.C'
echo '=========================================='
echo 
export NTUPLE_EVENTS=5000000
echo $NTUPLE_EVENTS ntuple events
root -b -l -n -q db_ntuple.C

echo

#-----------------------------------------------------------------------------
#Define grid 

export NS2T=60
export NEPS=60

export LO_S2T=0.01
export HI_S2T=0.25

export LO_EPS=-1.0e-2
export HI_EPS=+1.0e-2

echo 'Grid in to be used:'
echo 'NS2T='$NS2T
echo 'LO_S2T='$LO_S2T
echo 'HI_S2T='$HI_S2T
echo
echo 'NEPS='$NEPS
echo 'LO_EPS='$LO_EPS
echo 'HI_EPS='$HI_EPS
echo

#-----------------------------------------------------------------------------
# construct oscillated spectra for all points in the grid
echo '=========================================='
echo '3) Running db_osc_rate.C'
echo '=========================================='
echo
root -b -l -n -q db_osc_rate.C

#-----------------------------------------------------------------------------
# run minimization
echo '=========================================='
echo '3) Running db_minuit.C'
echo '=========================================='
echo
root -b -l -n -q db_minuit.C

echo

#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
#compile routines for minimization and marginalization
echo '=========================================='
echo 'compiling  db_chi2_min.cpp and db_margin.cpp'
echo '=========================================='
echo
g++ -o db_chi2_min.exe db_chi2_min.cpp
g++ -o db_margin.exe db_margin.cpp
##clang++ -o db_chi2_min.exe db_chi2_min.cpp
##clang++ -o db_margin.exe db_margin.cpp

echo
#-----------------------------------------------------------------------------
echo '=========================================='
echo 'executing  db_chi2_min.exe and db_margin.exe'
echo '=========================================='
echo
./db_chi2_min.exe $NS2T $NEPS ./
./db_margin.exe $NS2T $NEPS ./

echo

#-----------------------------------------------------------------------------
#Extract BF_CHI2, BF_S2T, BF_DM2 from chi2_minumum_SPEC.txt

read BF_CHI2 BF_S2T BF_DM2 <<< `cat files_data/chi2_minumum_RATE.txt`

#----------------------------------------------------------------------------
#----------------------------------------------------------------------------
# Edit gnuplot script
echo '=========================================='
echo 'Editting gnu plot script ...'
echo
#Create temporary file with new line in place
sed -i '' '125s/.*/set label 35 "+" at '$BF_S2T','$BF_DM2' center font "CharterBT-Roman,15"/' multi_plot_margin_RATE.gnu

sed -i '' '127s/.*/min = '$BF_CHI2'/' multi_plot_margin_RATE.gnu

echo

#----------------------------------------------------------------------------
#Execute gnuplot script
echo '=========================================='
echo 'Runnign gnuplot macro'
echo '=========================================='
echo
gnuplot multi_plot_margin_RATE.gnu

echo

#----------------------------------------------------------------------------
#Open in ghostview
gv files_plots/db_plots_RATE.eps &

#----------------------------------------------------------------------------
echo Done!
