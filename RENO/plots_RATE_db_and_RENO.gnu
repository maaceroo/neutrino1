set terminal postscript color "CharterBT-Roman" 12 enhanced size 7,7

########################################################################################
## Figure File
set output "Plots/COMB_plots_RATE.eps"
########################################################################################

## x-axis settings
#set logscale x
xmin = 0.0
xmax = 0.2
set xrange[xmin:xmax]
set xtics 0.0,0.05,0.19
#set grid xtics lc rgb "#bbbbbb" lw 1 lt 0
set xlabel "sin^{2}2{/Symbol q}_{13}" 
#set format x "10^{%T}"

## y-axis settings
ymin =  0.0
ymax =  17.0
set yrange[ymin:ymax]
set ylabel "{/Symbol D}{/Symbol c}^{2}"

min = 0.0
set key at 0.198,14

plot 'files/RENO_s2t_chi2_RATE.txt' u 1:2 w l lw 1 lc black dt 5 t "RENO", \
     '../DB_RATE_6AD/files_data/db_s2t_chi2_RATE.txt' u 1:2 w l lw 1 dt 3 black t "Daya Bay", \
     '../DB_RATE_6AD/files_data/db_s2t_chi2_COMBINED_rate.txt' u 1:(($2)-min)  w l lw 2 lc 1 t "Daya Bay + RENO", \
     16.0 lt 6 lw 2 t "99.97% C.L. (4{/Symbol s})", \
      9.0 lt 2 lw 2 t "99.73% C.L. (3{/Symbol s})", \
      4.0 lt 3 lw 2 t "95.45% C.L. (2{/Symbol s})", \
      1.0 lt 4 lw 2 t "68.27% C.L. (1{/Symbol s})"

