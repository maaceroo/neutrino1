set terminal postscript color "CharterBT-Roman" 12 enhanced size 7,7

########################################################################################
## Figure File
set output "Plots/plots_RATE_db_and_RENO.eps"
########################################################################################

##set origin 0,0

#set lmargin at screen 0.10
#set rmargin at screen 0.65
#set bmargin at screen 0.65
#set tmargin at screen 0.95

## x-axis settings
#set logscale x
xmin = 0.0
xmax = 0.25
set xrange[xmin:xmax]
set xlabel "sin^{2}2{/Symbol q}_{13}" 
#unset xtics
#set format x "10^{%T}"

## y-axis settings
ymin =  0.0
ymax =  18.0
set yrange[ymin:ymax]
set ylabel "{/Symbol D}{/Symbol c}^{2}"

set label " 68.27% C.L. (1{/Symbol s}) " at 0.18,1.3
set label " 95.45% C.L. (2{/Symbol s}) " at 0.18,4.3
set label " 99.73% C.L. (3{/Symbol s}) " at 0.18,9.3
set label " 99.98% C.L. (4{/Symbol s}) " at 0.18,15.7

#set key at 0.31,5

plot 'files/RENO_s2t_chi2_RATE.txt' u 1:2 w l lw 1 lc 1 t "RENO", '../DB_RATE/files_data/db_s2t_chi2_RATE.txt' u 1:2 w l lw 1 lc 10 t "Daya Bay",'files/chi2_s2t_RENO_plus_DB.txt' u 1:2 w l lw 1 lc 7 t "Daya Bay + RENO", \
     9.0 lt 1 lw 2 t "", \
     4.0 lt 2 lw 2 t "", \
     1.0 lt 3 lw 2 t "", \
     16.0 lt 4 lw 2 t "" 
  

    

