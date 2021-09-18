# Размер картинки относительно 512х512
S = 1

# Названия осей
set xlabel "r"
set ylabel "эквивалентная пластическая деформация"
#set ylabel "{~{/Symbol e}{0.6\\~}}^{  p}_& "
#set ylabel "{/Symbol e}@_{eqv}^c"
set xlabel font "Helvetica,14"
set ylabel font "Helvetica,14"
set xrange [1:4]

# Форматы чисел
set format y "%.1tE%+003T"

# Легенда
#set key off
set key box at graph 0.5, 1.01 bot center horizontal
set key reverse Left 
set key spacing 1
set key samplen 1
set key height -0.1
set key width 2
set key font "Helvetica,14"


# Терминал
set terminal pngcairo enhanced fontscale S size 512*S, 512*S
set output fn_out
set size ratio 1

plot \
fn_eps_pl_eqv with points pointtype 7 ps S*.2 lw S*1 lc rgb 'black' title '{~{/Symbol e}{0.6\~}}^{  p}'
