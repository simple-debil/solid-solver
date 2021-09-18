# Размер картинки относительно 512х512
S = 1

# Названия осей
set xlabel "r"
set ylabel "p"
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
set key width 2
set key font "Helvetica,14"

# Терминал
set terminal pngcairo enhanced fontscale S size 512*S, 512*S
set output fn_out
set size ratio 1


plot \
fn_P_analit_d with lines lw S lt 1 dt 1 lc rgb 'black' title 'p_{analit}', \
fn_P          with points pointtype 6 ps S*1 lw S*1.0 lc rgb 'black' title 'p_{contact}', \
fn_sigma1     with points pointtype 5 ps S*0.5 lw S*1.0 lc rgb 'black' title '{/Symbol s}_{z}', \

#plot \
#fn_P          with points pointtype 6 ps S*1 lw S*1.0 lc rgb 'black' title 'P', \
#fn_P_analit_d with lines lw S lt 1 dt 1 lc rgb 'black' title 'P_{a}', \
#fn_sigma1     with points pointtype 4 ps S*1 lw S*1.0 lc rgb 'black' title '{/Symbol s}_{1}', \
#fn_sigma2     with points pointtype 8 ps S*1 lw S*1.0 lc rgb 'black' title '{/Symbol s}_{2}', \
#fn_sigma3     with points pointtype 10 ps S*1 lw S*1.0 lc rgb 'black' title '{/Symbol s}_{3}', \

