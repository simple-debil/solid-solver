# Размер картинки относительно 512х512
S = 1

# Названия осей
set xlabel "time"
set ylabel "{/Symbol s}_{экв}"
set xlabel font "Helvetica,14"
set ylabel font "Helvetica,18"
set xrange [-1:]

# Форматы чисел
set format y "%.1tE%+003T"

# Легенда
set key off
#set label "релаксация" at graph 0.5, 0.85 offset 0, 1 center font "Helvetica,14"


# Терминал
set terminal pngcairo enhanced fontscale S size 512*S, 512*S
set output fn_out
set size ratio 1

plot \
fn_sigmaEqv using 1:2 with points pointtype 7 ps S*0.5 lw S*1 lc rgb 'black' notitle,\
fn_sigmaEqv using 1:3 with points pointtype 7 ps S*0.5 lw S*1 lc rgb 'black' notitle
