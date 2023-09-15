#!/bin/bash
#---------------------------------------------------------------
export JOBID=`(echo $SLURM_JOB_ID | cut -d. -f1)`
echo 'JOBID='${JOBID}
#---------------------------------------------------------------
#Define grid 
echo '=========================================='
echo '0) Define Grid'
echo '=========================================='
echo 

export NS2T=80
export NDM2=80

export LO_S2T=0.75
export HI_S2T=1.00

#export LO_DM2=1.0e-3
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


##---------------------------------------------------------------
#echo
## Construct ntuple
echo '=========================================='
echo '1) Running ntuple.C'
echo '=========================================='
echo
export NTUPLE_EVENTS=10000000
echo $NTUPLE_EVENTS ntuple events
time root -b -l -n -q MINOS_ntupla.C

echo

FILE=${JOBID}/data/minos_EScaleDerivative.root
if [[ -f "$FILE" ]]; then
echo "$FILE exists."
echo "The analysis will continue."
echo
else
echo "$FILE does not exist."
echo "Executing MINOS_EnergyS_ntuple.C"
time root -b -l -n -q MINOS_EnergyS_ntuple.C
echo "$FILE should exist now!"
echo "The analysis will continue."
echo
fi 
##---------------------------------------------------------------
## construct oscillated spectra for all points in the grid
echo '=========================================='
echo '2) Running MINOS_osc_2nu.C'
echo '=========================================='
echo
#time root -b -l -n -q MINOS_osc_2nu.C
#
#---------------------------------------------------------------
# run minimization
echo '=========================================='
echo '3) Running MINOS_minuit.C'
echo '=========================================='
echo
#root -b -l -n -q MINOS_minuit.C

echo

#---------------------------------------------------------------
#Remove first line from file
#tail -n +2 ${JOBID}/data/numu_chi2_s2t-dm2_surface.txt > ${JOBID}/data/numu_chi2_s2t-dm2_surface-noFL.txt
#tail -n +2 ${JOBID}/data/numub_chi2_s2t-dm2_surface.txt > ${JOBID}/data/numub_chi2_s2t-dm2_surface-noFL.txt
#tail -n +2 ${JOBID}/data/numubWS_chi2_s2t-dm2_surface.txt > ${JOBID}/data/numubWS_chi2_s2t-dm2_surface-noFL.txt

#---------------------------------------------------------------
#compile routines for minimization and marginalization
#echo '=========================================='
#echo 'compiling  minos_numu_chi2_min.cpp and minos_numu_margin.cpp'
#echo '=========================================='
#g++ -o minos_numu_chi2_min.exe minos_numu_chi2_min.cpp
#g++ -o minos_numu_margin.exe minos_numu_margin.cpp
#echo 'compiling  minos_numub_chi2_min.cpp and minos_numub_margin.cpp'
#echo '=========================================='
#g++ -o minos_numub_chi2_min.exe minos_numub_chi2_min.cpp
#g++ -o minos_numub_margin.exe minos_numub_margin.cpp
#echo 'compiling  minos_numubWS_chi2_min.cpp and minosWS_numub_margin.cpp'
#echo '=========================================='
#g++ -o minos_numubWS_chi2_min.exe minos_numubWS_chi2_min.cpp
#g++ -o minos_numubWS_margin.exe minos_numubWS_margin.cpp

echo
#---------------------------------------------------------------
#echo '=========================================='
#echo 'executing  minos_numu_chi2_numu_min.exe and minos_numu_margin.exe'
#./minos_numu_chi2_min.exe $NS2T $NDM2 ./${JOBID}
#./minos_numu_margin.exe $NS2T $NDM2 ./${JOBID}
#echo 'executing  minos_numub_chi2_min.exe and minos_numub_margin.exe'
#echo '=========================================='
#./minos_numub_chi2_min.exe $NS2T $NDM2 ./${JOBID}
#./minos_numub_margin.exe $NS2T $NDM2 ./${JOBID}
#echo 'executing  minos_numubWS_chi2_min.exe and minos_numubWS_margin.exe'
#echo '=========================================='
#./minos_numubWS_chi2_min.exe $NS2T $NDM2 ./${JOBID}
#./minos_numubWS_margin.exe $NS2T $NDM2 ./${JOBID}
#echo

#---------------------------------------------------------------
#Extract BF_CHI2, BF_S2T, BF_DM2 from chi2_minumum_SPEC.txt

#read BF_NUMU_CHI2   BF_NUMU_S2T   BF_NUMU_DM2 <<< `cat ${JOBID}/data/numu_chi2_minumum.txt`
#read BF_NUMUB_CHI2  BF_NUMUB_S2T  BF_NUMUB_DM2 <<< `cat ${JOBID}/data/numub_chi2_minumum.txt`
#read BF_NUMUBWS_CHI2  BF_NUMUBWS_S2T  BF_NUMUBWS_DM2 <<< `cat ${JOBID}/data/numubWS_chi2_minumum.txt`

#---------------------------------------------------------------
#---------------------------------------------------------------
# Edit gnuplot script
#echo '=========================================='
#echo 'Editting gnu plot scripts ...'
#echo

#echo '--------------- numu --------------'

#sed -i'' -e "12s/.*/set output \"${JOBID}\/files_plots\/MINOS_plots_numu.pdf\"/" multi_plot_margin_numu.gnu

#sed -i'' -e "117s/.*/set label 35 '+' at $BF_NUMU_S2T,$BF_NUMU_DM2*1e3 center font 'CharterBT-Roman,15'/" multi_plot_margin_numu.gnu

#sed -i'' -e "119s/.*/min = $BF_NUMU_CHI2/" multi_plot_margin_numu.gnu

#sed -i'' -e "48s/.*/plot '${JOBID}\/data\/numu_dm2_chi2.txt' u 2:(10**3*(\$1)) w l lw 2/" multi_plot_margin_numu.gnu

#sed -i'' -e "76s/.*/plot '${JOBID}\/data\/numu_s2t_chi2.txt' u 1:2 w l lw 2 t \"\", 16.0 lt 2 lw 2 t \"99.99% C.L. (4{\/Symbol s})\", 9.0 lt 3 lw 2 t \"99.73% C.L. (3{\/Symbol s})\", 4.0 lt 4 lw 2 t \"95.45% C.L. (2{\/Symbol s})\", 2.71 lt 5 lw 2 t \"90.00% C.L.\", 1.0 lt 6 lw 2 t \"68.27% C.L. (1{\/Symbol s})\"/" multi_plot_margin_numu.gnu

#sed -i'' -e "127s/.*/splot '${JOBID}\/data\/numu_chi2_s2t-dm2_surface-noFL.txt' u 1:(10**3*(\$2)):((\$3)-min) w l lw 2/" multi_plot_margin_numu.gnu

#echo
#echo '--------------- numubar ---------------'

#sed -i'' -e "12s/.*/set output \"${JOBID}\/files_plots\/MINOS_plots_numub.pdf\"/" multi_plot_margin_numub.gnu

#sed -i'' -e "117s/.*/set label 35 '+' at $BF_NUMUB_S2T,$BF_NUMUB_DM2*1e3 center font 'CharterBT-Roman,15'/" multi_plot_margin_numub.gnu

#sed -i'' -e "119s/.*/min = $BF_NUMUB_CHI2/" multi_plot_margin_numub.gnu

#sed -i'' -e "48s/.*/plot '${JOBID}\/data\/numub_dm2_chi2.txt' u 2:(10**3*(\$1)) w l lw 2/" multi_plot_margin_numub.gnu

#sed -i'' -e "76s/.*/plot '${JOBID}\/data\/numub_s2t_chi2.txt' u 1:2 w l lw 2 t \"\", 16.0 lt 6 lw 2 t \"99.99% C.L. (4{\/Symbol s})\", 9.0 lt 2 lw 2 t \"99.73% C.L. (3{\/Symbol s})\", 4.0 lt 3 lw 2 t \"95.45% C.L. (2{\/Symbol s})\", 2.71 lt 5 lw 2 t \"90.00% C.L.\", 1.0 lt 4 lw 2 t \"68.27% C.L. (1{\/Symbol s})\"/" multi_plot_margin_numub.gnu

#sed -i'' -e "127s/.*/splot '${JOBID}\/data\/numub_chi2_s2t-dm2_surface-noFL.txt' u 1:(10**3*(\$2)):((\$3)-min) w l lw 2/" multi_plot_margin_numub.gnu

#echo '--------------- numubarWS ---------------'

#sed -i'' -e "12s/.*/set output \"${JOBID}\/files_plots\/MINOS_plots_numubWS.pdf\"/" multi_plot_margin_numubWS.gnu

#sed -i'' -e "117s/.*/set label 35 '+' at $BF_NUMUBWS_S2T,$BF_NUMUBWS_DM2*1e3 center font 'CharterBT-Roman,15'/" multi_plot_margin_numubWS.gnu

#sed -i'' -e "119s/.*/min = $BF_NUMUBWS_CHI2/" multi_plot_margin_numubWS.gnu

#sed -i'' -e "48s/.*/plot '${JOBID}\/data\/numubWS_dm2_chi2.txt' u 2:(10**3*(\$1)) w l lw 2/" multi_plot_margin_numubWS.gnu

#sed -i'' -e "76s/.*/plot '${JOBID}\/data\/numubWS_s2t_chi2.txt' u 1:2 w l lw 2 t \"\", 16.0 lt 6 lw 2 t \"99.99% C.L. (4{\/Symbol s})\", 9.0 lt 2 lw 2 t \"99.73% C.L. (3{\/Symbol s})\", 4.0 lt 3 lw 2 t \"95.45% C.L. (2{\/Symbol s})\", 2.71 lt 5 lw 2 t \"90.00% C.L.\", 1.0 lt 4 lw 2 t \"68.27% C.L. (1{\/Symbol s})\"/" multi_plot_margin_numubWS.gnu

#sed -i'' -e "127s/.*/splot '${JOBID}\/data\/numubWS_chi2_s2t-dm2_surface-noFL.txt' u 1:(10**3*(\$2)):((\$3)-min) w l lw 2/" multi_plot_margin_numubWS.gnu

#echo

#---------------------------------------------------------------
#Execute gnuplot script
echo '=========================================='
echo 'Runnign gnuplot macros '
echo '=========================================='
echo
#gnuplot multi_plot_margin_numu.gnu
#gnuplot multi_plot_margin_numub.gnu
#gnuplot multi_plot_margin_numubWS.gnu

#gnuplot multi_plot_margin_compare.gnu
#gnuplot multi_plot_margin_compare.gnu

echo

echo Done!
