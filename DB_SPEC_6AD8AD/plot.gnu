##  plot.gnu - By: M.A.AceroO. - 27.Dec.2018  ##

########################################################################################

set terminal pdfcairo enhanced color font "CharterBT-Roman,14" size 6,6

########################################################################################
## Figure File
set output "files_plots/db_plot_COMPARE_1stbin_in_mel_noER_1M.pdf"
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

## x-axis settings
xmin = 0.05
xmax = 0.12
set xrange[xmin:xmax]
set xtics 0.05,0.01,0.12
set mxtics
set label 2 "sin^{2}2{/Symbol q}_{13}" at 0.085,2.03 center

## y-axis settings
set ytics offset -72
ymin = +2.1
ymax = +2.92
set yrange[ymin:ymax]
set ytics 2.1,0.1,2.9
set mytics
set label 4 "{/Symbol |D}m^{2}_{ee}| (10^{-3} eV^2)" at 0.042,2.5 center rotate by 90

## Mark at the BF
set label 5 '+' at 0.0811110000,0.0023330000*1e3 center font 'CharterBT-Roman,15'
## Minimum chi2 value
min = 1462.5800000000

unset ztics
set clabel
unset key
set grid ytics lc rgb "#bbbbbb" lw 1 lt 0
set grid xtics lc rgb "#bbbbbb" lw 1 lt 0

## Set Legend information
set arrow 11 from 0.09,2.85 to 0.1,2.85 nohead lw 3
set label 11 'Our Ana. (+ BF)' at 0.102,2.85 font 'CharterBT-Roman,13'
set arrow 22 from 0.09,2.8 to 0.1,2.8 nohead lw 3 dt 4
set label 22 'DB Ana. (* BF)' at 0.102,2.8 font 'CharterBT-Roman,13'
set label 35 '*' at 0.0841,2.50 center font 'CharterBT-Roman,15'

## Contour settings - Our Results
unset surface
set view 0,0
set size 2.0,2.0
set contour base
set cntrparam bspline
set cntrparam order 10
set cntrparam levels discret 2.30,6.18,11.83

## Our Result - Surface
splot 'files_data/chi2_s2t-dm2_surface_SPEC-noFL.txt' u 1:(1e3*($2)):(($3)-min) w l lw 3

## Contour settings - DB 1230-Days Results
#--Contour color and dashtype--
set linetype 101 lc 2 dt 4
set linetype 102 lc 3 dt 4
set linetype 103 lc 4 dt 4
#------------------------------
unset surface
set view 0,0
set size 2.0,2.0
set contour base
set cntrparam bspline
set cntrparam order 10
set cntrparam levels discret 2.30,6.18,11.83
set cntrparam firstlinetype 101

## DB1230 Result - Surface
splot 'files_data/DB_DeltaChiSq_1230days.txt' u 1:(1e3*($2)):3 w l lw 3

########################################################################################

unset multiplot

set output
set terminal x11
