##  multi_plot_margin.gnu - By: M.A.AceroO. - 23.Sep.2017  ##
## Script to put multple plots in one eps file. This code  ##
## generates the contour and the marginalization of the    ##
## oscillation params.                                     ##

########################################################################################

set terminal postscript color "CharterBT-Roman" 12 enhanced size 7,7

########################################################################################
## Figure File
set output "files_plots/db_plots_COMPARE.eps"
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
set arrow 1 from 1.00,ymin to 1.00,ymax nohead lt 4 lw 2
set arrow 3 from 4.00,ymin to 4.00,ymax nohead lt 3 lw 2
set arrow 5 from 9.00,ymin to 9.00,ymax nohead lt 2 lw 2
set arrow 7 from 16.0,ymin to 16.0,ymax nohead lt 6 lw 2

#plot 'files_data/db_dm2_chi2_SPEC.txt' u 2:1 w l lw 1
#plot 'files_data/db_dm2_chi2_SPEC.txt' u 2:(10**3*($1)) w l lw 1
plot 'files_data/db_dm2_chi2_SPEC.txt' u 2:(10**3*($1)) w l lw 1#,\
     'files_data/db_dm2_DB_PRL112.txt' u 2:1 w l lt 0 lw 2

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

set key at 0.305,9.0

plot 'files_data/db_s2t_chi2_SPEC.txt' u 1:2 w l lw 1 t "",      16.0 lt 6 lw 2 t "99.99% C.L. (4{/Symbol s})",      9.0 lt 2 lw 2 t "99.73% C.L. (3{/Symbol s})",      4.0 lt 3 lw 2 t "95.45% C.L. (2{/Symbol s})",      1.0 lt 4 lw 2 t "68.27% C.L. (1{/Symbol s})"#,\
    'files_data/db_s2t_DB_PRL112.txt' u 1:2 w l lt 0 lw 2 t "Daya Bay PRL112"

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
set ytics 1.5,0.2,3.4
#set format y "10^{%T}"
set label 4 "{/Symbol D}m^{2}_{ee} (10^{-3} eV^2)" at -0.0175,2.5 center rotate by 90
#set label 4 "{/Symbol D}m^{2}_{31} (10^{-3} eV^2)" at 0.006,2.5 center rotate by 90

## Mark at the BF
set label 35 "+" at 0.092821,0.002556*1e3 center font "CharterBT-Roman,15"
## Minimum chi2 value
min = 67.6181

unset ztics
set clabel
unset key

#splot 'files_data/chi2_s2t-dm2_surface_SPEC.txt' u 1:2:(()-min) w l lw 2
splot 'files_data/chi2_s2t-dm2_surface_SPEC-noFL.txt' u 1:(10**3*($2)):(($3)-min) w l lw 3

set arrow 11 from 0.08,1.7 to 0.11,1.7 nohead lw 3
set label 11 'Our Ana. (+ BF)' at 0.115,1.7 font 'CharterBT-Roman,11'
set arrow 22 from 0.08,1.6 to 0.11,1.6 nohead lw 1
set label 22 'DB PRL112 Ana. (x BF)' at 0.115,1.6 font 'CharterBT-Roman,11'
set label 25 'x' at 0.09,2.59 center font 'CharterBT-Roman,15'

splot 'files_data/chi2_scan.dat' u 1:2:3 w l lw 1 dt 4

########################################################################################

unset multiplot

set output
set terminal x11

