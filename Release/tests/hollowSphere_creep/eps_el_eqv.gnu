# Размер картинки относительно 512х512
S = 1

# Названия осей
set xlabel "r"
set ylabel "эквивалентная упругая деформация"
#set ylabel "{~{/Symbol e}{0.6\~}}^{  e}_& "
#set ylabel "~{/Symbol e}{0.6\\~}_& "
set xlabel font "Helvetica,14"
set ylabel font "Helvetica,14"
set xrange [1:4]

# Форматы чисел
set format y "%.1tE%+003T"

# Легенда
#set key off
set key box at graph 0.5, 1.01 bot center horizontal
#set key box outside top center horizontal
set key reverse Left 
set key spacing 1
set key samplen 1
set key height -0.1
set key width 2
set key font "Helvetica,14"
#set label "нагружение" at graph 0.3, 0.8 offset 0, 1 center font "Helvetica,14"
#set label "разгрузка" at graph 0.3, 0.5 offset 0, 1 center font "Helvetica,14"


# Терминал
set terminal pngcairo enhanced fontscale S size 512*S, 512*S
set output fn_out
set size ratio 1

plot \
fn_eps_el_eqv with points pointtype 7 ps S*.2 lw S*1 lc rgb 'black' title '{~{/Symbol e}{0.6\~}}^{  e}'
