set terminal postscript eps color solid linewidth 2 "Helvetica" 20
set title "Convergence comparation"
set xlabel "iter"
set ylabel "x"
set output "compare.eps"
set xrange [0:10]
plot "history_bisection.dat" with linespoints lw 2 t "bisection", "history_downhill.dat" with linespoints lw 2 t "downhill", "history_post.dat" lw 2 t "post", "history_jacobi.dat" with linespoints lw 2 t "jacobi", "history_post.dat" with linespoints lw 2 t "post", "history_atiken.dat" with linespoints lw 2 t "atiken"
set output