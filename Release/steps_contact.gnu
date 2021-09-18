# Размер картинки относительно 512х512
S = 1

# Названия осей
set xlabel "номер шага"
set ylabel "невязка"
set xlabel font "Helvetica,14"
set ylabel font "Helvetica,14"
set grid back xtics lw S*1 lt 1 lc rgb '#009900'
set xtics 1

# Форматы чисел
set format y "%.1tE%+003T"
set logscale y

# Легенда
set key box at graph 0.5, 1.01 bot center horizontal
#set key outside box top center
set key reverse Left 
set key spacing 1
set key samplen 1
set key height -0.3
set key width 1
set key font "Helvetica,14"


# Терминал
# set terminal pngcairo  transparent enhanced - прозрачный фон
set terminal pngcairo enhanced fontscale S size 512*S, 512*S
set output f_out
set size ratio 1

a = 0.7
b = 0.7
plot \
f_in using 1:5 with linespoints pointtype 4 ps S*a lw S*b lc rgb 'black' title '{/Symbol D}F/F' , \
f_in using 1:6 with linespoints pointtype 6 ps S*a lw S*b lc rgb 'black' title '{/Symbol D}r' 
