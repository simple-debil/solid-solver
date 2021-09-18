# Размер картинки относительно 512х512
S = 1

# Названия осей
set xlabel "номер шага"
set ylabel "количество итераций"
set xlabel font "Helvetica,14"
set ylabel font "Helvetica,14"
set yrange [0:]

# Форматы чисел
#set format y "%.1tE%+003T"
#set logscale y

# Легенда
#set key off
#set key outside box top center
set key box at graph 0.5, 1.005 bot center horizontal
set key reverse Left 
set key spacing 1
set key samplen 1
set key height -0.15
set key width -12
set key font "Helvetica,12"

# Терминал
# set terminal pngcairo  transparent enhanced - прозрачный фон
set terminal pngcairo enhanced fontscale S size 512*S, 512*S
set output f_out
set size ratio 1

a = 1
b = 1
plot \
f_in using 1:2 with linespoints pointtype 4 ps S*a lw S*b lc rgb 'black' title 'Итерации', \
f_in using 1:3 with linespoints pointtype 2 ps S*a lw S*b lc rgb 'black' title 'Изм. контакта', \
f_in using 1:4 with linespoints pointtype 6 ps S*a lw S*b lc rgb 'black' title 'Изм. ~{C}{0.7\~}'



#notitle  #, \
#f_in using 1:5 with linespoints pointtype 6 ps S*a lw S*b lc rgb 'black' title 'max\({/Symbol D}~{/Symbol s}{0.6\~}  / ~{/Symbol s}{0.6\~}  \)'


