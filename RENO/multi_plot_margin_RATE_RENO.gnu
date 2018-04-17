##  multi_plot_margin.gnu - By: M.A.AceroO. and D.J. Polo- 18.03.2018  ##
##     Script to put multple plots in one eps file. This code          ##
##      generates the contour and the marginalization of the           ##
##                         oscillation params.                         ##
#########################################################################

set terminal postscript color "CharterBT-Roman" 12 enhanced size 7,7

########################################################################################
## Figure File
set output "Plots/RENO_plots_RATE.eps"
########################################################################################

set multiplot
set origin 0,0

#######################################
##  Marginalization of epsilon       ##
#######################################
## Location
set lmargin at screen 0.65
set rmargin at screen 0.95
set bmargin at screen 0.10
set tmargin at screen 0.65

## x-axis settings
xmin =  0
xmax =  4.0
set xrange[xmin:xmax]
set xlabel "{/Symbol D}{/Symbol c}^{2}"

## y-axis settings
#set logscale y
ymin = -1e-2
ymax = 1e-2
set yrange[ymin:ymax]
unset ytics
#set format y "5^{%T}"

unset key

## Vertical lines ad 1,2,3 sigma (1D)
set arrow 1 from 1.00,ymin to 1.00,ymax nohead lt 4 lw 2
set arrow 3 from 4.00,ymin to 4.00,ymax nohead lt 3 lw 2
set arrow 5 from 9.00,ymin to 9.00,ymax nohead lt 2 lw 2

plot 'files/RENO_a_chi2_RATE.txt' u 2:1 w l lw 1 
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
xmax = 0.2
set xrange[xmin:xmax]
unset xtics
#set format x "10^{%T}"

## y-axis settings
ymin =  0.0
ymax =  10.0
set yrange[ymin:ymax]
set ylabel "{/Symbol D}{/Symbol c}^{2}"

set key at 0.31,5

plot 'files/RENO_s2t_chi2_RATE.txt' u 1:2 w l lw 1 t "",  \
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
xmin = 0.0
xmax = 0.2
set xrange[xmin:xmax]
set xtics 0,0.1,0.29
set mxtics
#set format x "10^{%T}"
set label 2 "sin^{2}2{/Symbol q}_{13}" at 0.090,-0.025 center

## y-axis settings
#set logscale y
set ytics offset -44.5
#set ytics -0.001,0.0005,0.0009
set mytics
ymin = -1e-2
ymax = +1e-2
set yrange[ymin:ymax]
set label 3 "a" at -0.03,-0.000125 center rotate by 90

set label 4 "{+ Rate-only}" at 0.16,0.008 center 

set label 35 "+" at 0.0879397, 0.000251256 center font "CharterBT-Roman,15"
## Mark at the BF
## Minimum chi2 value
min = 0.00948247

unset ztics
set clabel
unset key

splot 'files/chi2_s2t-a_surface_RATE.txt' u 1:2:(($3)-min) w l lw 2

########################################################################################

unset multiplot

set output
set terminal x11
