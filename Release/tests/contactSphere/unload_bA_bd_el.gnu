# Размер картинки относительно 512х512
S = 1

# Названия осей
set xlabel "unload\\_bd"
set ylabel "unload\\_ba"
set xlabel font "Helvetica,14"
set ylabel font "Helvetica,14"

# Форматы чисел
#set format y "%.1tE%+003T"

# Легенда
set key box at graph 0.5, 1.01 bot center horizontal
set key reverse Left 
set key spacing 1
set key samplen 1
set key height -0.1
set key width 0
set key font "Helvetica,14"

# Терминал
set terminal pngcairo enhanced fontscale S size 512*S, 512*S
set output fn_out
set size ratio 1

plot \
fn_unload_ba_bd        with points pointtype 6 ps S*1 lw S*1.0 lc rgb 'black' title 'unload\_ba_{numb}', \
fn_unload_ba_bd_analit with lines lw S lt 1 dt 1 lc rgb 'black' title 'unload\_ba_{analit}'

