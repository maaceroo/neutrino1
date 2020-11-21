##  multi_plot_margin.gnu - By: D. Polo T.M.A.AceroO. - 29.Sep.2018  ##
## Script to put multple plots in one eps file. This code  ##
## generates the contour and the marginalization of the    ##
## oscillation params.                                     ##

########################################################################################

set terminal pdfcairo enhanced color font "CharterBT-Roman,14" size 7,7

########################################################################################
## Figure File
set output "27406148/Plots/RENO_plots_SPEC_fudge_1.000_fFac1_0.990_fFac2_1.00.pdf"
########################################################################################

set multiplot
set origin 0,0

#######################################
##  Marginalization of delta m^2_ee  ##
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
#ymin = +1.2e-3
#ymax = +1e-2
ymin = +1.7
ymax = +3.5
set yrange[ymin:ymax]
unset ytics
#set format y "10^{%T}"

unset key

## Vertical lines ad 1,2,3,4 sigma (1D)
set arrow 1 from  1.00,ymin to 1.00,ymax nohead lt 4 lw 2
set arrow 3 from  4.00,ymin to 4.00,ymax nohead lt 3 lw 2
set arrow 5 from  9.00,ymin to 9.00,ymax nohead lt 2 lw 2
set arrow 7 from 16.00,ymin to 16.00,ymax nohead lt 6 lw 2

plot 'files/RENO_dm2_chi2_SPEC.txt' u 2:(10**3*($1)) w l lw 2

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
xmin = 0.0
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

set key at 0.29,7.0

plot 'files/RENO_s2t_chi2_SPEC.txt' u 1:2 w l lw 2 t "Spectral", \
16.0 lt 6 lw 2 t "99.99% C.L. (4{/Symbol s})", \
9.0 lt 2 lw 2 t "99.73% C.L. (3{/Symbol s})", \
4.0 lt 3 lw 2 t "95.45% C.L. (2{/Symbol s})", \
1.0 lt 4 lw 2 t "68.27% C.L. (1{/Symbol s})"
#'../RENO/files/RENO_s2t_chi2_RATE.txt' u 1:2 w l lw 1 lc black dt 4 t "Rate-only", \

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
xmin = 0.0
#xmax = 0.3
xmax = 0.2
set xrange[xmin:xmax]
set xtics 0.0,0.05,0.19
set mxtics
#set format x "10^{%T}"
set label 2 "sin^{2}2{/Symbol q}_{13}" at 0.1,1.50 center

## y-axis settings
#set logscale y
set ytics offset -58.0
#ymin = +1.7e-3
#ymax = +3.5e-3
ymin = 1.7
ymax = 3.5
set yrange[ymin:ymax]
set ytics 1.5,0.5,3.4
#set format y "10^{%T}"
set label 4 "{/Symbol D}m^{2}_{ee} (10^{-3} eV^2)" at -0.021,2.6 center rotate by 90
#set label 4 "{/Symbol D}m^{2}_{31} (10^{-3} eV^2)" at 0.006,2.5 center rotate by 90

#set label 5 "{+ Best-fit}" at 0.12,0.0038 center
#set label 6 "{* Best-fit}" at 0.12,0.0037 center 

## Mark at the BF
set label 35 '+' at 0.105,0.00305*1e3 center font 'CharterBT-Roman,15'
## Minimum chi2 value
min = 425.419

unset ztics
set clabel
unset key
set grid ytics lc rgb "#bbbbbb" lw 1 lt 0
set grid xtics lc rgb "#bbbbbb" lw 1 lt 0

splot 'files/chi2_s2t-dm2_surface_spect-noFL.txt' u 1:(($2)*1e3):(($3)-min) w l lw 3

unset xtics
unset ytics
unset label 4
unset label 2

set arrow 11 from 0.11,3.35 to 0.135,3.35 nohead lw 3
set label 11 'This work ({/Symbol \053} BF)' at 0.14,3.35 font 'CharterBT-Roman,13'
set arrow 22 from 0.11,3.2 to 0.135,3.2 nohead lw 2 dt (2,8,2,8)
set label 22 'RENO ({/Symbol \264} BF)' at 0.14,3.2 font 'CharterBT-Roman,13'
set label 25 '{/Symbol \264}' at 0.0896,2.68 center font 'CharterBT-Roman,15'

plot  'files/2018_RENO_68CL.txt' u 1:2 w l lt 4 lw 2 dt (2,8,2,8), \
'files/2018_RENO_95CL.txt' u 1:2 w l lt 3 lw 2 dt (2,8,2,8), \
'files/2018_RENO_99CL.txt' u 1:2 w l lt 2 lw 2 dt (2,8,2,8)

########################################################################################

unset multiplot

set output
set terminal x11


