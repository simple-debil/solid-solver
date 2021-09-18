size_r = 4
size_z = 4

# Размер картинки относительно 512х512
S = 1

# Названия осей
set xrange [ * : * ] noreverse writeback
set x2range [ * : * ] noreverse writeback
set yrange [ * : * ] noreverse writeback
set y2range [ * : * ] noreverse writeback
set zrange [ * : * ] noreverse writeback
set cbrange [ * : * ] noreverse writeback
set rrange [ * : * ] noreverse writeback
set xlabel "r/a"
set ylabel "z/a" 
#set zlabel  offset character 1, 0, 0 font "" textcolor lt -1 norotate
set ytics 1
set xtics 1

# Форматы чисел
#set format y "%.1tE%+003T"

# Легенда
set key nobox at graph 0.5, 1.1 bot center horizontal
#set key reverse Left 
set key spacing 1
set key samplen 1
set key height 1
set key width 1
set key font "Helvetica,14"

# Терминал
set terminal pngcairo enhanced fontscale S size 512*S, 512*S
set size ratio 1

# переменные
set dummy r, z


# сплайны(численное решение)
set dgrid3d  32, 32,  4
set table $Gridded_sigma1
splot fn_2d_sigma1
unset table
unset dgrid3d


# контур
set style increment default
set style textbox opaque margins  0.5,  0.5 fc  bgnd noborder linewidth  1.0
set view map
set samples 50, 50
set isosamples 50, 50

set contour base
set cntrparam order 8
set cntrparam bspline #cubicspline
set cntrparam levels auto 10 
#set cntrparam levels discrete -0.4, -0.1, -0.02, 0.0, 0.008
set cntrlabel format '%8.3g' font ',7' start 1 interval -1 onecolor
set style data lines


unset surface
# построение графиков
set output fn_out
set multiplot
splot $Gridded_sigma1 with lines dt 1 lc rgb 'black' notitle
splot $Gridded_sigma1 with labels boxed notitle
unset multiplot
