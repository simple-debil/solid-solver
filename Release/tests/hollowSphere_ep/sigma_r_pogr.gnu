# Размер картинки относительно 512х512
S = 1

# Названия осей
set xlabel "r"
set ylabel "%"
set xlabel font "Helvetica,14"
set ylabel font "Helvetica,18"
set logscale y

# Форматы чисел
set format y "%.1tE%+003T"

# Легенда
set label "нагружение" at graph 0.5, 0.95 offset 0, 0 center font "Helvetica,14"
set key on box
set key reverse Left 
set key at graph 0.5, 0.8 bot center
set key spacing 1
set key samplen 1
set key height 0.2
set key width 4
set key font "Helvetica Bold,14"


# Терминал
set terminal pngcairo enhanced fontscale S size 512*S, 512*S
set output fn_out
set size ratio 1

plot \
fn_sigma_r_pogr with points pointtype 7 ps S*.2 lw S*1 lc rgb 'black' title '{/Symbol D}{/Symbol s}_r / {/Symbol s}_r{/Symbol \264}100%'
#"./tests/hollowSphere_ep/unload_sigma_r_pogr.txt" with points pointtype 7 ps S*.2 lw S*1 lc rgb 'black' notitle, \
