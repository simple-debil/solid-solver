# Размер картинки относительно 512х512
S = 1

# Названия осей
set xlabel "r"
set ylabel "F_{contact}"
set xlabel font "Helvetica,14"
set ylabel font "Helvetica,14"
set yrange [0:]

# Форматы чисел
set format y "%.1tE%+003T"

# Легенда
set key box at graph 0.5, 1.01 bot center horizontal
set key reverse Left 
set key spacing 1
set key samplen 1
set key height -0.1
set key width 1
set key font "Helvetica,14"

# Терминал
set terminal pngcairo enhanced fontscale S size 512*S, 512*S
set output fn_out
set size ratio 1

plot \
fn_F with points pointtype 6 ps S lw S lc rgb 'black' title 'F_{contact}'
#fn_F with linespoints pointtype 6 ps S lw S lt 1 dt 1 lc rgb 'black' title 'F '
