# Размер картинки относительно 512х512
S = 1

# Названия осей
set xlabel "r"
set ylabel "T"
set xlabel font "Helvetica,14"
set ylabel font "Helvetica,18"
set xrange [:4]

# Форматы чисел
set format y "%.1tE%+003T"

# Легенда
set key off


# Терминал
set terminal pngcairo enhanced fontscale S size 512*S, 512*S
set output "./tests/hollowSphere/T.png"
set size ratio 1

plot\
"./tests/hollowSphere/T_ch.txt" with points pointtype 7 ps S*.2 lw S*1 lc rgb 'black' notitle, \
"./tests/hollowSphere/T.txt" with lines lw S lt 1 dt 1 lc rgb 'black' notitle
