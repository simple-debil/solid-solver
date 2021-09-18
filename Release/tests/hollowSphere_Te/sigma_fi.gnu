# Размер картинки относительно 512х512
S = 1

# Названия осей
set xlabel "r"
set ylabel "{/Symbol s}_{/Symbol f}"
set xlabel font "Helvetica,14"
set ylabel font "Helvetica,18"
set xrange [:4]

# Форматы чисел
set format y "%.1tE%+003T"

# Легенда
set key off
#set label "нагружение" at graph 0.3, 0.8 offset 0, 1 center font "Helvetica,14"
#set label "разгрузка" at graph 0.3, 0.5 offset 0, 1 center font "Helvetica,14"


# Терминал
set terminal pngcairo enhanced fontscale S size 512*S, 512*S
set output fn_out
set size ratio 1

plot \
fn_sigma_fi with points pointtype 7 ps S*.2 lw S*1 lc rgb 'black' notitle, \
fn_sigma_fi_analit with lines lw S lt 1 dt 1 lc rgb 'black' notitle


