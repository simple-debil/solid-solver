# Размер картинки относительно 512х512
S = 1

# Названия осей
set xlabel "r"
set ylabel "T"
set xlabel font "Helvetica,14"
set ylabel font "Helvetica,18"
set xrange [:4]
#set logscale y

# Форматы чисел
set format y "%.1tE%+003T"

# Легенда
set key off


# Терминал
set terminal pngcairo enhanced fontscale S size 512*S, 512*S
set output fn_out
set size ratio 1

plot \
fn_T_pogr with points pointtype 7 ps S*.2 lw S*1 lc rgb 'black' notitle

