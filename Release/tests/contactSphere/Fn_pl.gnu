# Размер картинки относительно 512х512
S = 1

# Названия осей
set xlabel "d"
set ylabel "P"
set xlabel font "Helvetica,14"
set ylabel font "Helvetica,14"

# Форматы чисел
set format y "%.1tE%+003T"

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
fn_Fn        with points pointtype 6 ps S*1 lw S*1.0 lc rgb 'black' title 'P_{numb}', \
fn_Fn_analit with lines lw S lt 1 dt 2 lc rgb 'black' title 'P@^{elastic}_{analit}'

