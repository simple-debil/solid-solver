# Размер картинки относительно 512х512
S = 1

# Названия осей
set xlabel "time"
set ylabel "{/Symbol e}"
set xlabel font "Helvetica,14"
set ylabel font "Helvetica,18"
set xrange [-1:1000]

# Форматы чисел
set format y "%.1tE%+003T"

# Легенда
#set key off
#set label "релаксация" at graph 0.5, 0.85 offset 0, 1 center font "Helvetica,14"
#set key on box
set key outside box top center #offset 0, 5

set key reverse Left 
#set key at graph 0.5, 0.5 center center
set key spacing 1
set key samplen 1
set key height 0.2
set key width 4
set key font "Helvetica Bold,14"


# Терминал
set terminal pngcairo enhanced fontscale S size 512*S, 512*S
set output fn_out
set size ratio 1

plot \
fn_eps11 using 1:2 with points pointtype 6 ps S*1 lw S*1 lc rgb 'black' title '{/Symbol e}_{xx}',\
fn_eps11 using 1:3 with points pointtype 6 ps S*1 lw S*1 lc rgb 'black' notitle,\
fn_eps11_analit using 1:2 with lines lw S lt 1 dt 1 lc rgb 'black' notitle,\
fn_eps22 using 1:2 with points pointtype 4 ps S*1 lw S*1 lc rgb 'black' title '{/Symbol e}_{yy}',\
fn_eps22 using 1:3 with points pointtype 4 ps S*1 lw S*1 lc rgb 'black' notitle,\
fn_eps22_analit using 1:2 with lines lw S lt 1 dt 1 lc rgb 'black' notitle,\
fn_eps33 using 1:2 with points pointtype 8 ps S*1 lw S*1 lc rgb 'black' title '{/Symbol e}_{zz}',\
fn_eps33 using 1:3 with points pointtype 8 ps S*1 lw S*1 lc rgb 'black' notitle,\
fn_eps33_analit using 1:2 with lines lw S lt 1 dt 1 lc rgb 'black' notitle
