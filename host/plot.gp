set term png size 800,800
set output "path.png"
plot "path.csv" using 2:3 title "Path" with points
set term png size 600,260
set output "v.png"
plot "path.csv" using 1:4 title "Velocity"
set output "a.png"
plot "path.csv" using 1:5 title "Acceleration"
set output "j.png"
plot "path.csv" using 1:6 title "Jerk"
