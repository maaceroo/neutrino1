##  plot.gnu - By: M.A.AceroO. - 21.Nov.2018  ##

########################################################################################

set terminal pdfcairo enhanced color font "CharterBT-Roman,12" size 6,6

########################################################################################
## Figure File
set output "27430293/Plots/plot_SPEC_fudge_1.002_fFac1_0.992_fFac2_1.000_EScale.pdf"
########################################################################################

set multiplot
set origin 0,0

#####################################
## Contour plot for allowed region ##
#####################################
## Location
set lmargin at screen 0.12
set rmargin at screen 0.95
set bmargin at screen 0.12
set tmargin at screen 0.95

## Contour settings
unset surface
set view 0,0
set size 1.0,1.0
set contour base
set cntrparam bspline
set cntrparam order 10
set cntrparam levels discret 2.30,6.18,11.83

## x-axis settings
xmin = 0.0
xmax = 0.2
set xrange[xmin:xmax]
set xtics 0.0,0.05,0.2
set mxtics
set label 2 "sin^{2}2{/Symbol q}_{13}" at 0.1,1.55 center

## y-axis settings
set ytics offset -75.0
ymin = +1.7
ymax = +3.5
set yrange[ymin:ymax]
set ytics 1.5,0.5,3.5
set mytics
set label 4 "{/Symbol D}m^{2}_{ee} (10^{-3} eV^2)" at -0.022,2.6 center rotate by 90

## Mark at the BF
set label 35 '+' at 0.09,0.002553*1e3 center font 'CharterBT-Roman,15'
## Minimum chi2 value
min = 29.3818

unset ztics
set clabel
unset key
set grid ytics lc rgb "#bbbbbb" lw 1 lt 0
set grid xtics lc rgb "#bbbbbb" lw 1 lt 0

splot '27430293/files/chi2_s2t-dm2_surface_spect-noFL.txt' u 1:(($2)*1e3):(($3)-min) w l lw 3
#splot 'files/test.txt' u 1:2:(($3)-min) w l lw 2

unset xtics
unset ytics
unset label 4
unset label 2

set arrow 11 from 0.11,3.35 to 0.14,3.35 nohead lw 3
set label 11 'This work ({/Symbol \053} BF)' at 0.145,3.35 font 'CharterBT-Roman,11'
set arrow 22 from 0.11,3.2 to 0.14,3.2 nohead lw 1 dt 4
set label 22 'RENO ({/Symbol \264} BF)' at 0.145,3.2 font 'CharterBT-Roman,11'
set label 25 '{/Symbol \264}' at 0.0896,2.68 center font 'CharterBT-Roman,15'

set label 55 'Minuit EScale' at 0.145,1.75 font 'CharterBT-Roman,15'

plot  'files/2018_RENO_68CL.txt' u 1:2 w l lt 4 lw 1 dt 4, \
      'files/2018_RENO_95CL.txt' u 1:2 w l lt 3 lw 1 dt 4, \
      'files/2018_RENO_99CL.txt' u 1:2 w l lt 2 lw 1 dt 4

########################################################################################

unset multiplot

set output
set terminal x11
