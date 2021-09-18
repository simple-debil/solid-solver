# Размер картинки относительно 512х512
S = 1

# Названия осей
set xlabel "номер шага"
set ylabel "невязка"
set xlabel font "Helvetica,14"
set ylabel font "Helvetica,14"
#set grid back xtics lw S*1 lt 1 lc rgb '#009900'
#set xtics 1

# Форматы чисел
set format y "%.1tE%+003T"
set logscale y

# Легенда
set key box at graph 0.5, 1.01 bot center horizontal
#set key box outside top center horizontal
set key reverse Left 
set key spacing 1
set key samplen 1
set key height -0.3
set key width 2
set key font "Helvetica,14"

#set label "(Initial {/Symbol s})" at 0, 1 offset -10, 2 left


# Терминал
# set terminal pngcairo  transparent enhanced - прозрачный фон
set terminal pngcairo enhanced fontscale S size 512*S, 512*S
set output f_out
set size ratio 1

a = 1
b = 1
plot \
f_in using 1:12  with linespoints pointtype 4 ps S*a lw S*b lc rgb 'black' title '{/Symbol D}{/Symbol s}/{/Symbol s}', \
f_in using 1:8  with linespoints pointtype 2 ps S*a lw S*b lc rgb 'black' title '{/Symbol D}F/F', \
f_in using 1:9  with linespoints pointtype 3 ps S*a lw S*b lc rgb 'black' title '{/Symbol D}r'


#f_in using 1:11  with linespoints pointtype 6 ps S*a lw S*b lc rgb 'black' title '{/Symbol D}{/Symbol e}/{/Symbol e}', \
#f_in using 1:7  with linespoints pointtype 8 ps S*a lw S*b lc rgb 'black' title '{/Symbol D}R/R'
#|{/Symbol D}{/:Bold r}|




#f_in using 1:11 with linespoints pointtype 6 ps S*a lw S*b lc rgb 'black' title 'max\({/Symbol D}~{/Symbol e}{0.6\~}  / ~{/Symbol e}{0.6\~} \)'
#'deps_{creep}/eps_{creep}' 

