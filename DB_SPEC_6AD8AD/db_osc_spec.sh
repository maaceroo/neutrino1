#!/bin/bash
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
#Define grid 
echo '=========================================='
echo '0) Define Grid'
echo '=========================================='
echo 

export NS2T=10
export NDM2=10

export LO_S2T=0.05
export HI_S2T=0.12

export LO_DM2=2.0e-3
export HI_DM2=3.0e-3

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
echo '1) Running ldist.C'
echo '=========================================='
echo
#time root -b -l -n -q ldist.C

echo

#-----------------------------------------------------------------------------
# Construct ntuple
echo '=========================================='
echo '2) Running ntuple.C'
echo '=========================================='
echo
export NTUPLE_EVENTS=1000000
echo $NTUPLE_EVENTS ntuple events
#time root -b -l -n -q db_ntuple.C

echo

#-----------------------------------------------------------------------------
# construct oscillated spectra for all points in the grid
echo '=========================================='
echo '3) Running db_osc_spec.C'
echo '=========================================='
echo
time root -b -l -n -q db_osc_spec.C

#-----------------------------------------------------------------------------
# run minimization
echo '=========================================='
echo '3) Running db_minuit.C'
echo '=========================================='
echo
#time root -b -l -n -q db_minuit_spec_CovMat_9pars.C
time root -b -l -n -q db_minuit_spec_CovMat_1par.C
#time root -b -l -n -q db_minuit_spec_CovMat.C

echo

#-----------------------------------------------------------------------------
#Remove first line from file
tail -n +2 files_data/chi2_s2t-dm2_surface_SPEC.txt > files_data/chi2_s2t-dm2_surface_SPEC-noFL.txt

#-----------------------------------------------------------------------------
#compile routines for minimization and marginalization
echo '=========================================='
echo 'compiling  db_chi2_min.cpp and db_margin.cpp'
echo '=========================================='
echo
g++ -o db_chi2_min.exe db_chi2_min.cpp
g++ -o db_margin.exe db_margin.cpp
#clang++ -o db_chi2_min.exe db_chi2_min.cpp
#clang++ -o db_margin.exe db_margin.cpp

echo
#-----------------------------------------------------------------------------
echo '=========================================='
echo 'executing  db_chi2_min.exe and db_margin.exe'
echo '=========================================='
echo
./db_chi2_min.exe $NS2T $NDM2 ./
./db_margin.exe $NS2T $NDM2 ./

echo

#-----------------------------------------------------------------------------
#Extract BF_CHI2, BF_S2T, BF_DM2 from chi2_minumum_SPEC.txt

read BF_CHI2 BF_S2T BF_DM2 <<< `cat files_data/chi2_minumum_SPEC.txt`

##Extract fudge, fFac1 and fFac2 from constants.h
#fFac1=$(awk 'NR == 39 {print $4}' constants.h)
#fFac2=$(awk 'NR == 40 {print $4}' constants.h)
#echo 'fFac1 = ' $fFac1
#echo 'fFac2 = ' $fFac2

#----------------------------------------------------------------------------
#----------------------------------------------------------------------------
# Form gnuplot script
echo '=========================================='
echo 'Editting gnu plot script ...'
echo '=========================================='
echo 'Multiplot Script... Done!'
echo '-------------------------'
echo
#Editing the multiplot file - Contours and marginalizatins
sed -i'' -e "132s/.*/set label 35 '+' at $BF_S2T,$BF_DM2*1e3 center font 'CharterBT-Roman,15'/" multi_plot_margin_SPEC.gnu

sed -i'' -e "134s/.*/min = $BF_CHI2/" multi_plot_margin_SPEC.gnu
echo 'Comparisson plot Script... Done!'
echo '--------------------------------'

#Editing the plot file - Contours comparison
sed -i'' -e "42s/.*/set label 5 '+' at $BF_S2T,$BF_DM2*1e3 center font 'CharterBT-Roman,15'/" plot.gnu

sed -i'' -e "44s/.*/min = $BF_CHI2/" plot.gnu

#sed -i'' -e "9s/.*/set output \"files_plots\/db_plot_COMPARE_1stbin_in_0.7_noER.pdf\"/" plot.gnu
sed -i'' -e "9s/.*/set output \"files_plots\/db_plot_COMPARE_1stbin_in_mel_noER_1M.pdf\"/" plot.gnu
#sed -i'' -e "9s/.*/set output \"files_plots\/db_plot_COMPARE_1stbin_unPhys_noER.pdf\"/" plot.gnu

echo

#----------------------------------------------------------------------------
#Execute gnuplot script
echo '=========================================='
echo 'Runnign gnuplot macro'
echo '=========================================='
echo
#gnuplot multi_plot_margin_SPEC.gnu
gnuplot plot.gnu
#rm *.gnu-e

echo
echo '=========================================='
echo 'A new plot has been generated. Look at it '
echo 'in the files_plots/ directory.'
echo '=========================================='
echo
#----------------------------------------------------------------------------
#Open in ghostview
#gv files_plots/db_plots_SPEC.eps &
#gv files_plots/db_plot_COMPARE.eps &

#----------------------------------------------------------------------------
echo Done!
