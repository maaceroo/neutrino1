#!/bin/bash
#-----------------------------------------------------------------------------
export JOBID=`(echo $PBS_JOBID | cut -d. -f1)`
echo 'JOBID='${JOBID}
#-----------------------------------------------------------------------------
#Define grid
echo '=========================================='
echo '0) Define Grid'
echo '=========================================='
echo

export NS2T=40
export NDM2=40

export LO_S2T=0.01
export HI_S2T=0.20

#export LO_DM2=1.5e-3 #to use for the combined (DB+RENO) ana.
export LO_DM2=1.7e-3
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

#Correction factors
export fudge=$FUDGE
export fFac1=$FFAC1
export fFac2=$FFAC2

echo 'Correction factors used in this run:'
echo '------------------------------------'
echo 'fudge='$fudge
echo 'fFac1='$fFac1
echo 'fFac2='$fFac2
echo

#-----------------------------------------------------------------------------
echo
#
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
echo '3) Running RENO_ntuple_noosc_spect.C'
echo '=========================================='
echo
export NTUPLE_EVENTS=1000000
echo $NTUPLE_EVENTS ntuple events
time root -b -l -n -q RENO_ntuple_noosc_spect.C

echo

#-----------------------------------------------------------------------------
# construct oscillated spectra for all points in the grid
echo '=========================================='
echo '4) Running RENO_osc_spect.C'
echo '=========================================='
echo
time root -b -l -n -q RENO_osc_spect.C

#-----------------------------------------------------------------------------
# Construct ntuple
echo '=========================================='
echo '5) Running RENO_EScaleDeriv_ntuple.C'
echo '=========================================='
echo
time root -b -l -n -q RENO_EScaleDeriv_ntuple.C

echo

#-----------------------------------------------------------------------------
# run minimization
echo '=========================================='
echo '6) Running RENO_minuit_spect.C'
echo '=========================================='
echo

#sel=1 #RENO_minuit_spect.C
sel=2 #RENO_minuit_spect_EScale.C
if [ $sel = 1 ]
then
    echo "sel = " $sel
    time root -b -l -n -q RENO_minuit_spect.C
else    
     echo "sel = " $sel
    time root -b -l -n -q RENO_minuit_spect_EScale.C
fi

echo

#-----------------------------------------------------------------------------
#Remove first line from file
tail -n +2 ${JOBID}/files/chi2_s2t-dm2_surface_spect.txt > ${JOBID}/files/chi2_s2t-dm2_surface_spect-noFL.txt

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
time ./RENO_margin_spect.exe $NS2T $NDM2 ./${JOBID}

echo

#-----------------------------------------------------------------------------
#Extract BF_CHI2, BF_S2T, BF_DM2 from chi2_minumum_SPEC.txt
read BF_S2T BF_DM2 BF_CHI2 <<< `cat ${JOBID}/files/chi2_minimun_spect.txt`

#Extract fudge, fFac1 and fFac2 from constants.h
#fudge=$(awk 'NR == 36 {print $4}' constants.h)
#fFac1=$(awk 'NR == 37 {print $4}' constants.h)
#fFac2=$(awk 'NR == 38 {print $4}' constants.h)
#echo 'fudge = ' $fudge
#echo 'fFac1 = ' $fFac1

#----------------------------------------------------------------------------
#----------------------------------------------------------------------------
# Form gnuplot script
echo '=========================================='
echo 'Editting gnu plot script ...'
echo '=========================================='
#echo '-------------------------'
#sed -i'' -e "12s/.*/set output \"${JOBID}\/Plots\/RENO_plots_SPEC_fudge_$fudge\_fFac1_$fFac1\_fFac2_$fFac2.pdf\"/" multi_plot_margin_spect_RENO.gnu

#sed -i'' -e "136s/.*/set label 35 '+' at $BF_S2T,$BF_DM2*1e3 center font 'CharterBT-Roman,15'/" multi_plot_margin_spect_RENO.gnu

#sed -i'' -e "138s/.*/min = $BF_CHI2/" multi_plot_margin_spect_RENO.gnu
#echo 'Multiplot Script... Done!'

echo '--------------------------------'
sed -i'' -e "9s/.*/set output \"${JOBID}\/Plots\/plot_SPEC_fudge_$fudge\_fFac1_$fFac1\_fFac2_$fFac2\_EScale.pdf\"/" plot.gnu

sed -i'' -e "51s/.*/set label 35 '+' at $BF_S2T,$BF_DM2*1e3 center font 'CharterBT-Roman,15'/" plot.gnu

sed -i'' -e "53s/.*/min = $BF_CHI2/" plot.gnu

sed -i'' -e "61s/.*/splot '${JOBID}\/files\/chi2_s2t-dm2_surface_spect-noFL.txt' u 1:((\$2)*1e3):((\$3)-min) w l lw 3/" plot.gnu

if [ $sel -eq 1 ]
then
    echo "sel = " $sel
    sed -i'' -e "75s/.*/set label 55 'Minuit Normal' at 0.145,1.75 font 'CharterBT-Roman,15'/" plot.gnu
else
    echo "sel = " $sel
    sed -i'' -e "75s/.*/set label 55 'Minuit EScale' at 0.145,1.75 font 'CharterBT-Roman,15'/" plot.gnu
fi

echo 'Comparisson plot Script... Done!'
echo

#----------------------------------------------------------------------------
#Execute gnuplot script
echo '=========================================='
echo 'Runnign gnuplot macro'
echo '=========================================='
echo
#gnuplot multi_plot_margin_spect_RENO.gnu
gnuplot plot.gnu
#rm *.gnu-e

echo

##----------------------------------------------------------------------------
##Open in ghostview
##gv Plots/RENO_plots_SPEC.eps &
##gv Plots/plot_SPEC.eps &
##----------------------------------------------------------------------------

echo Done!
