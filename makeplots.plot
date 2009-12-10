set terminal png truecolor enhanced
set output "function.png"
plot "functionvalues.out"
set output "delta.png"
plot "delta.out" with lines
set output "backtracks.png"
plot "line.out" with lines
