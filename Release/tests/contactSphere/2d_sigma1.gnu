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

# аналитическое решение
SQR(x) = x*x
pow(x, n) = x**n

w(r,z) = sqrt(0.5*(SQR(r) + SQR(z) - SQR(a) + sqrt(SQR(SQR(r) + SQR(z) - SQR(a)) + 4*SQR(a*z))))
sigma_fi(r,z) = -((1.-2.*NU)/2.*SQR(a/r)*(1 - pow(z/w(r,z), 3))\
                      + 3./2.*z/w(r,z)*(2*NU + (1 - NU)*SQR(w(r,z))/(SQR(a) + SQR(w(r,z))) - (1 + NU)*w(r,z)/a*atan(a/w(r,z))))
sigma_z(r,z) = -3./2.*pow(z/w(r,z), 3)*SQR(a*w(r,z))/(pow(w(r,z), 4) + SQR(a*z))
sigma_r(r,z) = -(sigma_fi(r,z) + sigma_z(r,z)) - 3.*z/w(r,z)*(1 + NU)*(1 - w(r,z)/a*atan(a/w(r,z)));
sigma_rz(r,z) = -3./2.*r*w(r,z)*SQR(z)/(pow(w(r,z), 4) + SQR(a*z))*SQR(a)/(SQR(a) + SQR(w(r,z)));

sigma1(r,z) = (sigma_r(r,z) + sigma_z(r,z))/2. + sqrt(SQR((sigma_r(r,z) - sigma_z(r,z))/2.) + SQR(sigma_rz(r,z)));
sigma2(r,z) = sigma_fi(r,z);
sigma3(r,z) = (sigma_r(r,z) + sigma_z(r,z))/2. - sqrt(SQR((sigma_r(r,z) - sigma_z(r,z))/2.) + SQR(sigma_rz(r,z)));

set xrange [0.0001 : size_r] noreverse writeback      #r
set yrange [-size_z : -0.0001] noreverse writeback    #z



# контур
set style increment default
set style textbox opaque margins  0.5,  0.5 fc  bgnd noborder linewidth  1.0
set view map
set samples 50, 50
set isosamples 50, 50

set contour base
set cntrparam order 8
set cntrparam bspline #cubicspline
#set cntrparam levels auto 10 
set cntrparam levels discrete -0.4, -0.1, -0.02, 0.0, 0.008
set cntrlabel format '%8.3g' font ',7' start 1 interval -1 onecolor
set style data lines


unset surface
# построение графиков
set output fn_out
set multiplot
splot $Gridded_sigma1 with lines dt 2 lc rgb 'black' notitle
#splot $Gridded_sigma1 with labels boxed notitle
splot sigma1(r*a,-z*a) with lines lc rgb 'black' notitle #title '{/Symbol s}_{}'
splot sigma1(r*a,-z*a) with labels boxed notitle
unset multiplot
