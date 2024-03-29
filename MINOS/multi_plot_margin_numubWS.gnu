##  multi_plot_margin.gnu - By: M.A.AceroO. - 27.Oct.2021  ##
## Script to put multple plots in one pdf file. This code  ##
## generates the contour and the marginalization of the    ##
## oscillation parameters.                                 ##

########################################################################################

set terminal pdfcairo enhanced color font "CharterBT-Roman,14" size 7,7

########################################################################################
## Figure File
set output "2827952/files_plots/MINOS_plots_numubWS.pdf"
########################################################################################

set multiplot
set origin 0,0

#######################################
##  Marginalization of delta m^2_32  ##
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
ymin = +1.7
ymax = +3.5
set yrange[ymin:ymax]
unset ytics

unset key

## Vertical lines ad 1,2,3 sigma (1D)
set arrow 1 from  1.00,ymin to 1.00,ymax nohead lt 6 lw 2
set arrow 2 from  2.71,ymin to 2.71,ymax nohead lt 5 lw 2
set arrow 3 from  4.00,ymin to 4.00,ymax nohead lt 4 lw 2
set arrow 5 from  9.00,ymin to 9.00,ymax nohead lt 3 lw 2
set arrow 7 from 16.00,ymin to 16.00,ymax nohead lt 2 lw 2

plot '2827952/data/numubWS_dm2_chi2.txt' u 2:(10**3*($1)) w l lw 2

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
xmin = 0.75
xmax = 1.00
set xrange[xmin:xmax]
unset xtics
set format x "10^{%T}"

## y-axis settings
ymin =  0
ymax =  17
set yrange[ymin:ymax]
set ylabel "{/Symbol D}{/Symbol c}^{2}"

set key at 1.14,9.0

plot '2827952/data/numubWS_s2t_chi2.txt' u 1:2 w l lw 2 t "", 16.0 lt 6 lw 2 t "99.99% C.L. (4{/Symbol s})", 9.0 lt 2 lw 2 t "99.73% C.L. (3{/Symbol s})", 4.0 lt 3 lw 2 t "95.45% C.L. (2{/Symbol s})", 2.71 lt 5 lw 2 t "90.00% C.L.", 1.0 lt 4 lw 2 t "68.27% C.L. (1{/Symbol s})"

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
set cntrparam levels discret 2.30,4.61,6.18,11.83,19.33

## x-axis settings
xmin = 0.75
xmax = 1.00
set xrange[xmin:xmax]
#set xtics 0.75,0.05,0.96
#set mxtics
#set format x "10^{%T}"
set label 2 "sin^{2}2{/Symbol q}_{23}" at 0.875,1.55 center

## y-axis settings
set ytics offset -51.1
ymin = 1.7
ymax = 3.5
set yrange[ymin:ymax]
#set ytics 2.1,0.1,2.8
#set format y "10^{%T}"
set label 4 "{/Symbol D}m^{2}_{32} (10^{-3} eV^2)" at 0.725,2.5 center rotate by 90

## Mark at the BF
set label 35 '+' at 0.0000000000,0.0000000000*1e3 center font 'CharterBT-Roman,15'
## Minimum chi2 value
min = 0.0000000000

unset ztics
set clabel
unset key
set grid ytics lc rgb "#bbbbbb" lw 1 lt 0
set grid xtics lc rgb "#bbbbbb" lw 1 lt 0

splot '2827952/data/numubWS_chi2_s2t-dm2_surface-noFL.txt' u 1:(10**3*($2)):(($3)-min) w l lw 2

unset xtics
unset ytics
unset label 4
unset label 2
unset label 35

##set label 35 'X' at 0.955,0.00238*1e3 center font 'CharterBT-Roman,15'

##plot './data/contour90CL_NuB_A_tab.txt' w l lt rgb "black" lw 3, './data/contour90CL_NuB_B_tab.txt' w l lt rgb "black" lw 3
########################################################################################

unset multiplot

set output
set terminal x11
