set term pngcairo lw 2 font "AR PL UKai CN,20" size 640*4, 600*2
set output "result.png"
set grid
set sample 10000

set yrange [0:2]

set multiplot layout 2, 3

set title "Interpolation curve"

set xlabel "x"
set ylabel "y"
plot 'lagrange.dat' using 1:3:4 with line lw 1.5 lc rgb "#2B60DE" t 'Lagrange', 'samples.dat' with points pt 5 ps 1.5 lc black t 'Samples', 1/(1+x**2) with line ls 1 lw 1 lc black t 'Origin'

set xlabel "x"
set ylabel "y"
plot 'newton.dat' using 1:3:4 with line lw 1.5 lc rgb "#2B60DE" t 'Newton', 'samples.dat' with points pt 5 ps 1.5 lc black t 'Samples', 1/(1+x**2) with line ls 1 lw 1 lc black t 'Origin'

set xlabel "x"
set ylabel "y"
plot 'cubic_spline.dat' using 1:3:4 with line lw 1.5 lc rgb "#2B60DE" t 'Cubic spline', 'samples.dat' with points pt 5 ps 1.5 lc black t 'Samples', 1/(1+x**2) with line ls 1 lw 1 lc black t 'Origin'

set title "Error distribution"

set xlabel "x"
set ylabel "y"
plot 'lagrange.dat' using 1:3:4 with yerrorlines lw 0.1 ps 0 lc rgb "#2B60DE" t 'Lagrange', 'samples.dat' with points pt 5 ps 1.5 lc black t 'Samples', 1/(1+x**2) with line ls 1 lw 1 lc black t 'Origin'

set xlabel "x"
set ylabel "y"
plot 'newton.dat' using 1:3:4 with yerrorlines lw 0.1 ps 0 lc rgb "#2B60DE" t 'Newton', 'samples.dat' with points pt 5 ps 1.5 lc black t 'Samples', 1/(1+x**2) with line ls 1 lw 1 lc black t 'Origin'

set xlabel "x"
set ylabel "y"
plot 'cubic_spline.dat' using 1:3:4 with yerrorlines lw 0.1 ps 0 lc rgb "#2B60DE" t 'Cubic spline', 'samples.dat' with points pt 5 ps 1.5 lc black t 'Samples', 1/(1+x**2) with line ls 1 lw 1 lc black t 'Origin'

unset multiplot