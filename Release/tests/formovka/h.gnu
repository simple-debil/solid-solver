# Размер картинки относительно 512х512
S = 1

# Названия осей
set xlabel "X"
set ylabel "h"
set xlabel font "Helvetica,14"
set ylabel font "Helvetica,16"

# Форматы чисел
set format y "%.1tE%+003T"

# Легенда
set key box at graph 0.5, 1.01 bot center horizontal
set key reverse Left 
set key spacing 1
set key samplen 1
set key height -0.3
set key width 0
set key font "Helvetica,14"

#set key on box
#set key reverse Left 
#set key at graph 0.5, 0.05 bot center
#set key spacing 1
#set key samplen 1

# Терминал
set terminal pngcairo enhanced fontscale S size 512*S, 512*S
set output fn_out
set size ratio 1

plot \
fn_h with linespoints pointtype 6 ps S lw S lt 1 dt 1 lc rgb 'black' title 'h '
