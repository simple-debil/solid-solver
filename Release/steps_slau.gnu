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
#set key box at graph 0.5, 1.01 bot center vertical
#set key outside box top center #offset 0, 5
set key reverse Left
#set key at graph 0.5, -0.3 bot center
set key spacing 1
set key samplen 1
#set key height -5
set key height -0
set key width -18
set key font "Helvetica Bold,14"


# Терминал
# set terminal pngcairo  transparent enhanced - прозрачный фон
set terminal pngcairo enhanced fontscale S size 512*S, 512*S
set output f_out
set size ratio 1

a = 0.7
b = 0.7
plot \
f_in using 1:4 with linespoints pointtype 4 ps S*a lw S*b lc rgb 'black' title '|{/:Bold A q}_{ prev} ‒ {/:Bold b}| / |{/:Bold b}|', \
f_in using 1:2 with linespoints pointtype 6 ps S*a lw S*b lc rgb 'black' title '|{/:Bold A q} ‒ {/:Bold b}| / |{/:Bold b}|'


#f_in using 1:3 with linespoints pointtype 8 ps S*a lw S*b lc rgb 'black' title 'slau rel' , \
