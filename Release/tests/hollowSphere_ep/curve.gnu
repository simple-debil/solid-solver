# Размер картинки относительно 512х512
S = 1

# Названия осей
set xlabel "q"
set ylabel "{/Symbol s}"
set xlabel font "Helvetica,14"
set ylabel font "Helvetica,14"

# Форматы чисел
set format y "%.1tE%+003T"

# Количество делений по x
ntics = 4
stats fn_curve using 1 name 'x' nooutput
set xtics x_max/ntics

# Легенда
set key box at graph 0.5, 1.01 bot center horizontal
set key reverse Left 
set key spacing 1
set key samplen 1
set key height 0
set key width 1
set key font "Helvetica,14"


# Терминал
set terminal pngcairo enhanced fontscale S size 512*S, 512*S
set output fn_out
set size ratio 1

plot \
fn_curve with lines lw S lt 1 dt 1 lc rgb 'black' title 'curve'

