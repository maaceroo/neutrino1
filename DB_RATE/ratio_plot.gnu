### Script to plot the ratio Ndetected/Nexpected ###
### for the Daya Bay data analysis.              ###
### Based on data from PRL112 021801 (2014)      ###

set terminal postscript color "CharterBT-Roman" 12 enhanced size 9,5 \
fontfile add "~/Dropbox/UAtlantico/Docs_Research_UdelA/gnuplot/Fonts/bchr8a.pfa" \
fontfile add "~/Dropbox/UAtlantico/Docs_Research_UdelA/gnuplot/Fonts/mathcal.pfa" \
fontfile add "~/Dropbox/UAtlantico/Docs_Research_UdelA/gnuplot/Fonts/sym.pfa" 

set output "files_plots/db_ratio.eps"

set xlabel "Weighted Baseline (m)"
set ylabel "N_{detected}/N_{expected}"

set xrange [0:2000]
set yrange[0.90:1.01]

s2th = 0.09
dm2_31 = 2.32e-3
dm2_21 = 7.59e-5
NuE = 4.23

plot 1-(s2th*(sin(1.267*dm2_31*x/NuE))**2)-(0.25*(1+sqrt(1-s2th))**2)*0.861*(sin(1.267*dm2_21*x/NuE))**2 t "SurvP", \
     "survProb.txt" u 1:3:4 w errorbar t ""

set output
set terminal x11
