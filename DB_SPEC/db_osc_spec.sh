#!/bin/bash
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------

# Construct L distribution
root -b -l -n -q ldist.C

#-----------------------------------------------------------------------------
# Construct ntuple
export NTUPLE_EVENTS=5000000
echo $NTUPLE_EVENTS ntuple events
root -b -l -n -q db_ntuple.C


#-----------------------------------------------------------------------------
#Define grid size

export NS2T=100
export NDM2=100

#-----------------------------------------------------------------------------
# construct oscillated spectra for all points in the grid
root -b -l -n -q db_osc_spec.C

#-----------------------------------------------------------------------------
# run minimization
root -b -l -n -q db_minuit_spec.C

#-----------------------------------------------------------------------------
#Remove first line from file 
tail -n +2 files_data/chi2_s2t-dm2_surface_SPEC.txt > files_data/chi2_s2t-dm2_surface_SPEC-noFL.txt


#-----------------------------------------------------------------------------
#compile routines

g++ -o db_chi2_min.exe db_chi2_min.cpp
g++ -o db_margin.exe db_margin.cpp

#-----------------------------------------------------------------------------

./db_chi2_min.exe $NS2T $NDM2 ./
./db_margin.exe $NS2T $NDM2 ./

#-----------------------------------------------------------------------------
#Extract BF_CHI2, BF_S2T, BF_DM2 from chi2_minumum_SPEC.txt

read BF_CHI2 BF_S2T BF_DM2 <<< `cat files_data/chi2_minumum_SPEC.txt`

#----------------------------------------------------------------------------
#----------------------------------------------------------------------------
# Form gnuplot script

cat > multi_plot_margin_SPEC.gnu <<EOF
##  multi_plot_margin.gnu - By: M.A.AceroO. - 23.Sep.2017  ##
## Script to put multple plots in one eps file. This code  ##
## generates the contour and the marginalization of the    ##
## oscillation params.                                     ##

########################################################################################

set terminal postscript color "CharterBT-Roman" 12 enhanced size 7,7

########################################################################################
## Figure File
set output "files_plots/db_plots_SPEC.eps"
########################################################################################

set multiplot
set origin 0,0

#######################################
##  Marginalization of delta m^2_31  ##
#######################################
## Location
set lmargin at screen 0.65
set rmargin at screen 0.95
set bmargin at screen 0.10
set tmargin at screen 0.65

## x-axis settings
xmin =  0
xmax =  17
set xrange[xmin:xmax]
set xlabel "{/Symbol D}{/Symbol c}^{2}"

## y-axis settings
#set logscale y
#ymin = +1e-4
#ymax = +1e-2
ymin = +1.5
ymax = +3.5
set yrange[ymin:ymax]
unset ytics
#set format y "10^{%T}"

unset key

## Vertical lines ad 1,2,3 sigma (1D)
set arrow 1 from  1.00,ymin to 1.00,ymax nohead lt 4 lw 2
set arrow 3 from  4.00,ymin to 4.00,ymax nohead lt 3 lw 2
set arrow 5 from  9.00,ymin to 9.00,ymax nohead lt 2 lw 2
set arrow 7 from 16.00,ymin to 16.00,ymax nohead lt 6 lw 2

#plot 'files_data/db_dm2_chi2_SPEC.txt' u 2:1 w l lw 1
#plot 'files_data/db_dm2_chi2_SPEC.txt' u 2:(10**3*(\$1)) w l lw 1
#plot 'files_data/db_dm2_chi2_SPEC.txt' u 2:(10**3*(\$1)) smooth csplines w l lw 1
plot 'files_data/db_dm2_chi2_SPEC.txt' u 2:(10**3*(\$1)) w l lw 1

reset
######################################
## Marginalization of sin^2(2Theta) ##
######################################
## Location
set lmargin at screen 0.10
set rmargin at screen 0.65
set bmargin at screen 0.65
set tmargin at screen 0.95

## x-axis settings
#set logscale x
xmin = 0.01
#xmax = 0.3
xmax = 0.2
set xrange[xmin:xmax]
unset xtics
set format x "10^{%T}"

## y-axis settings
ymin =  0
ymax =  17
set yrange[ymin:ymax]
set ylabel "{/Symbol D}{/Symbol c}^{2}"

set key at 0.31,3.8

plot 'files_data/db_s2t_chi2_SPEC.txt' u 1:2 w l lw 1 t "", \
    16.0 lt 6 lw 2 t "99.99% C.L. (4{/Symbol s})", \
     9.0 lt 2 lw 2 t "99.73% C.L. (3{/Symbol s})", \
     4.0 lt 3 lw 2 t "95.45% C.L. (2{/Symbol s})", \
     1.0 lt 4 lw 2 t "68.27% C.L. (1{/Symbol s})"

reset


#####################################
## Contour plot for allowed region ##
#####################################
## Location
set lmargin at screen 0.10
set rmargin at screen 0.65
set bmargin at screen 0.10
set tmargin at screen 0.65

## Contour settings
unset surface
set view 0,0
set size 1.0,1.0
set contour base
set cntrparam bspline
set cntrparam order 10
set cntrparam levels discret 2.30,6.18,11.83

## x-axis settings
#set logscale x
xmin = 0.01
#xmax = 0.3
xmax = 0.2
set xrange[xmin:xmax]
set xtics 0.0,0.05,0.19
#set mxtics
#set format x "10^{%T}"
set label 2 "sin^{2}2{/Symbol q}_{13}" at 0.1,1.25 center

## y-axis settings
#set logscale y
set ytics offset -43.5
#ymin = +1e-4
#ymax = +1e-2
ymin = 1.5
ymax = 3.5
set yrange[ymin:ymax]
set ytics 1.5,0.5,3.4
#set format y "10^{%T}"
set label 4 "{/Symbol D}m^{2}_{31} (eV^2)" at -0.0175,2.5 center rotate by 90
#set label 4 "{/Symbol D}m^{2}_{31} (10^{-3} eV^2)" at 0.006,2.5 center rotate by 90

## Mark at the BF
set label 35 "+" at $BF_S2T,$BF_DM2*1e3 center font "CharterBT-Roman,15"
#set label 35 "+" at 0.0965,2.15 center font "CharterBT-Roman,15"
## Minimum chi2 value
min = $BF_CHI2

unset ztics
set clabel
unset key
set grid ytics lc rgb "#bbbbbb" lw 1 lt 0
set grid xtics lc rgb "#bbbbbb" lw 1 lt 0

#splot 'files_data/chi2_s2t-dm2_surface_SPEC.txt' u 1:2:(($3)-min) w l lw 2
#splot 'files_data/chi2_s2t-dm2_surface_SPEC.txt' u 1:(10**3*(\$2)):((\$3)-min) w l lw 2
splot 'files_data/chi2_s2t-dm2_surface_SPEC-noFL.txt' u 1:(10**3*(\$2)):((\$3)-min) w l lw 2

########################################################################################

unset multiplot

set output
set terminal x11

EOF
#----------------------------------------------------------------------------
#----------------------------------------------------------------------------

#----------------------------------------------------------------------------
#Execute gnuplot script
gnuplot multi_plot_margin_SPEC.gnu


#----------------------------------------------------------------------------
#Open in ghostview
gv files_plots/db_plots_SPEC.eps &
