##  plot.gnu - By: M.A.AceroO. - 21.Nov.2018  ##

########################################################################################

set terminal postscript color "CharterBT-Roman" 12 enhanced size 6,6

########################################################################################
## Figure File
set output "Plots/plot_SPEC.eps"
########################################################################################

set multiplot
set origin 0,0

#####################################
## Contour plot for allowed region ##
#####################################
## Location
set lmargin at screen 0.17
set rmargin at screen 0.97
set bmargin at screen 0.15
set tmargin at screen 0.95

## Contour settings
unset surface
set view 0,0
set size 2.0,2.0
set contour base
set cntrparam bspline
set cntrparam order 10
set cntrparam levels discret 2.30,6.18,11.83

## x-axis settings
xmin = 0.00
xmax = 0.2
set xrange[xmin:xmax]
set xtics 0.0,0.02,0.19
set mxtics
set label 2 "sin^{2}2{/Symbol q}_{13}" at 0.1,0.0009 center

## y-axis settings
set ytics offset -55
ymin = +1.2e-3
ymax = +3.52e-3
set yrange[ymin:ymax]
set ytics 0.0,0.0005,0.0039
set mytics
set label 4 "{/Symbol |D}m^{2}_{ee}| (eV^2)" at -0.032,0.0024 center rotate by 90

#set label 5 "{+ Best-fit}" at 0.12,0.0038 center

## Mark at the BF
set label 35 '+' at 0.073333,0.002611 center font 'CharterBT-Roman,15'
## Minimum chi2 value
min = 21.9017

unset ztics
set clabel
unset key
set grid ytics lc rgb "#bbbbbb" lw 1 lt 0
set grid xtics lc rgb "#bbbbbb" lw 1 lt 0

#splot 'files/chi2_s2t-dm2_surface_spect-noFL.txt' u 1:2:(($3)-min) w l lw 2
splot 'files/test.txt' u 1:2:(($3)-min) w l lw 2

unset xtics
unset ytics
unset label 4
unset label 2

set arrow 11 from 0.11,3.4e-3 to 0.14,3.4e-3 nohead lw 2
set label 11 'Our Ana. (+ BF)' at 0.145,3.4e-3 font 'CharterBT-Roman,11'
set arrow 22 from 0.11,3.2e-3 to 0.14,3.2e-3 nohead lw 2 dt 4
set label 22 'RENO Ana. (* BF)' at 0.145,3.2e-3 font 'CharterBT-Roman,11'
set label 25 '*' at 0.082,0.00262 center font 'CharterBT-Roman,15'

plot  'files/RENO_68CL.txt' u 1:(($2)*1e-3) w l lt 4 lw 2 dt 4, \
      'files/RENO_95CL.txt' u 1:(($2)*1e-3) w l lt 3 lw 2 dt 4, \
      'files/RENO_99CL.txt' u 1:(($2)*1e-3) w l lt 2 lw 2 dt 4

########################################################################################

unset multiplot

set output
set terminal x11


