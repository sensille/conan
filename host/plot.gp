set term png size 800,800
set output "path.png"
plot "path.csv" using 2:3 title "Path" with lines
set term png size 600,260
set output "v.png"
plot "path.csv" using 1:4 title "Velocity" with lines
set output "a.png"
plot "path.csv" using 1:5 title "Acceleration" with lines
set output "j.png"
plot "path.csv" using 1:6 title "Jerk" with lines

set output "xy-err.png"
plot "path.csv" using 1:8 title "xy error" with lines
set output "v-err.png"
plot "path.csv" using 1:9 title "v error" with lines
set output "a-err.png"
plot "path.csv" using 1:10 title "a error" with lines

set term png size 600,260
set output "x-fn.png"
plot "path.csv" using 1:12 title "x target" with lines
set output "y-fn.png"
plot "path.csv" using 1:13 title "y target" with lines
set output "vx-fn.png"
plot "path.csv" using 1:14 title "vx target" with lines
set output "vy-fn.png"
plot "path.csv" using 1:15 title "vy target" with lines
set output "ax-fn.png"
plot "path.csv" using 1:16 title "ax target" with lines
set output "ay-fn.png"
plot "path.csv" using 1:17 title "ay target" with lines

set output "x-sim.png"
plot "path.csv" using 1:2 title "x actual" with lines
set output "y-sim.png"
plot "path.csv" using 1:3 title "y actual" with lines
set output "vx-sim.png"
plot "path.csv" using 1:19 title "vx actual" with lines
set output "vy-sim.png"
plot "path.csv" using 1:20 title "vy actual" with lines
set output "ax-sim.png"
plot "path.csv" using 1:21 title "ax actual" with lines
set output "ay-sim.png"
plot "path.csv" using 1:22 title "ay actual" with lines

set output "x-5th.png"
plot "path.csv" using 1:24 title "x 5th" with lines
set output "y-5th.png"
plot "path.csv" using 1:25 title "y 5th" with lines
set output "vx-5th.png"
plot "path.csv" using 1:26 title "vx 5th" with lines
set output "vy-5th.png"
plot "path.csv" using 1:27 title "vy 5th" with lines
set output "ax-5th.png"
plot "path.csv" using 1:28 title "ax 5th" with lines
set output "ay-5th.png"
plot "path.csv" using 1:29 title "ay 5th" with lines
