set terminal postscript eps
set output "history.eps"
plot "history.dat" with linespoints
set output