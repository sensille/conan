all: conan serial motion b-approx

VER=28
SRC = serial.c crc16.c main.c
HDR = checksum.h serial.h
C_SRC = motion.c clothoid.c circle.c planner.c line.c conan.c gcode.c b-approx.c
C_HDR = conan.h

serial: $(SRC) $(HDR) Makefile
	gcc -g -Wall $(SRC) -lm -o serial

conan: $(C_SRC) $(C_HDR) Makefile
	gcc -g -Wall $(C_SRC) -lm  -lgmp -lmpfr -o conan

motion: motion.c Makefile circle.c
	gcc -g -Wall motion.c circle.c -DTEST -lm  -lgmp -lmpfr -o motion

b-approx: b-approx.c Makefile
	gcc -g -Wall b-approx.c -DBATEST -lm  -lgmp -lmpfr -o b-approx

plot:
	gnuplot plot.gp
	mkdir -p /var/www/html/3dpfs/$(VER)
	cp index.html *.png /var/www/html/3dpfs/$(VER)/

gcode:
	gnuplot gc/data.gp
	mkdir -p /var/www/html/3dpfs/$(VER)
	cp -R gc /var/www/html/3dpfs/$(VER)/

clean:
	rm serial *.png gmp
