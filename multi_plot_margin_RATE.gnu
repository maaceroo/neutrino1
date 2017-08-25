##  multi_plot_margin.gnu - By: M.A.AceroO. - 29.May.2017  ##
## Script to put multple plots in one eps file. This code  ##
## generates the contour and the marginalization of the    ##
## oscillation params.                                     ##

########################################################################################

set terminal postscript color "CharterBT-Roman" 12 enhanced size 7,7

########################################################################################
## Figure File
set output "db_plots_RATE.eps"
########################################################################################

set multiplot
set origin 0,0

#######################################
##  Marginalization of epsilon       ##
#######################################

set lmargin at screen 0.65
set rmargin at screen 0.95
set bmargin at screen 0.10
set tmargin at screen 0.65

xmin =  0
xmax =  10
set xrange[xmin:xmax]
set xlabel "{/Sym D}{/Sym c}^{2}"

#set logscale y
ymin = -1e-2
ymax = +1e-2
set yrange[ymin:ymax]
unset ytics
#set format y "10^{%T}"

unset key
set arrow 1 from 1.00,ymin to 1.00,ymax nohead lt 4 lw 2
set arrow 3 from 4.00,ymin to 4.00,ymax nohead lt 3 lw 2
set arrow 5 from 9.00,ymin to 9.00,ymax nohead lt 2 lw 2

plot 'db_eps_chi2_RATE.txt' u 2:1 w l lw 1

reset
######################################
## Marginalization of sin^2(2Theta) ##
######################################

set lmargin at screen 0.10
set rmargin at screen 0.65
set bmargin at screen 0.65
set tmargin at screen 0.95

#set logscale x
xmin = 0.0
xmax = 0.5
set xrange[xmin:xmax]
unset xtics
#set format x "10^{%T}"

ymin =  0
ymax =  10
set yrange[ymin:ymax]
set ylabel "{/Sym D}{/Sym c}^{2}"

set key at 0.77,5

plot 'db_s2t_chi2_RATE.txt' u 1:2 w l lw 1 t "", \
     9.0 lt 2 lw 2 t "99.73% C.L. (3{/Sym s})", \
     4.0 lt 3 lw 2 t "95.45% C.L. (2{/Sym s})", \
     1.0 lt 4 lw 2 t "68.27% C.L. (1{/Sym s})"

reset


#####################################
## Contour plot for allowed region ##
#####################################

set lmargin at screen 0.10
set rmargin at screen 0.65
set bmargin at screen 0.10
set tmargin at screen 0.65

unset surface
set view 0,0
set size 1.0,1.0

set contour base
set cntrparam bspline
set cntrparam order 10
set cntrparam levels discret 2.30,6.18,11.83

#set logscale x
xmin = 0.0
xmax = 0.5
set xrange[xmin:xmax]
#set format x "10^{%T}"
set label 2 "sin^{2}2{/Sym q}_{13}" at 0.25,-0.0125 center

#set logscale y
#set format y "10^{%T}"
ymin = -1e-2
ymax = +1e-2
set yrange[ymin:ymax]
set label 4 "{/Sym D}m^{2}_{ee} (eV^2)" at -0.07,0.0 center rotate by 90

set label 35 "+" at 0.0904523,0.000150754 center font "CharterBT-Roman,15"

unset ztics
set clabel
unset key

min = 0.00484571

splot 'chi2_s2t-eps_surface_RATE.txt' u 1:2:(($3)-min) w l lw 2

########################################################################################

unset multiplot

set output
set terminal x11
