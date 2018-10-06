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

set term png size 1200,800
set output "vx.png"
plot "path.csv" using 1:12 title "vx" with points
