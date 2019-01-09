set term png size 6000,6000
set output "ba/ba.png"
plot "ba.data" using 2:3 title "B-Spline" with points pointtype 8 pointsize 0.05, \
     "orig.data" using 2:3 title "Original" with lines lt rgb "blue", \
     "orig.data" using 2:3 title "Original" with points lt rgb "red", \
     "extra.data" using 2:3 title "Extra" with points lt rgb "green", \
     "ctrlp.data" using 1:2 title "hull" with points lt rgb "grey", \
     "ctrlp.data" using 1:2 title "hull" with lines lt rgb "grey", \
     "knots.data" using 1:2 title "Knots" with points lt rgb "black" pt 8

set term png size 6000,400
set output "ba/v.png"
plot "vajcs.data" using 1:2 title "velocity" with lines lt rgb "blue", \
     "vajcs.data" using 1:9 title "velocity(r)" with lines lt rgb "red"
set term png size 6000,400
set output "ba/a.png"
plot "vajcs.data" using 1:3 title "acceleration" with lines lt rgb "blue", \
     "vajcs.data" using 1:10 title "acceleration(r)" with lines lt rgb "red"
set term png size 6000,400
set output "ba/j.png"
plot "vajcs.data" using 1:4 title "jerk" with lines lt rgb "blue", \
     "vajcs.data" using 1:11 title "jerk" with lines lt rgb "red"
set term png size 6000,400
set output "ba/c.png"
plot "vajcs.data" using 1:5 title "curvature" with lines lt rgb "blue", \
     "vajcs.data" using 1:12 title "curvature" with lines lt rgb "red"
set term png size 6000,400
set output "ba/s.png"
plot "vajcs.data" using 1:6 title "sharpness" with lines lt rgb "blue", \
     "vajcs.data" using 1:13 title "sharpness" with lines lt rgb "red"
