# Размер картинки относительно 512х512
S = 4

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
set label "нагружение" at graph 0.3, 0.8 offset 0, 1 center font "Helvetica,14"
set label "разгрузка" at graph 0.3, 0.5 offset 0, 1 center font "Helvetica,14"


# Терминал
set terminal pngcairo enhanced fontscale S size 512*S, 512*S
set output "./tests/hollowSphere/sigma_fi.png"
set size ratio 1

plot\
"./tests/hollowSphere/load_sigma_fi_ch.txt" with points pointtype 7 ps S*.2 lw S*1 lc rgb 'black' notitle, \
"./tests/hollowSphere/unload_sigma_fi_ch.txt" with points pointtype 7 ps S*.2 lw S*1 lc rgb 'black' notitle, \
"./tests/hollowSphere/load_sigma_fi.txt" with lines lw S lt 1 dt 1 lc rgb 'black' notitle, \
"./tests/hollowSphere/unload_sigma_fi.txt" with lines lw S lt 1 dt 1 lc rgb 'black' notitle
