all: conan serial

VER=4
SRC = serial.c crc16.c main.c
HDR = checksum.h serial.h
C_SRC = motion.c clothoid.c circle.c planner.c line.c conan.c
C_HDR = conan.h

serial: $(SRC) $(HDR) Makefile
	gcc -g -Wall $(SRC) -lm -o serial

conan: $(C_SRC) $(C_HDR) Makefile
	gcc -g -Wall $(C_SRC) -lm  -lgmp -lmpfr -o conan

plot:
	gnuplot plot.gp
	-mkdir /var/www/html/3dpfs/$(VER)
	cp index.html *.png /var/www/html/3dpfs/$(VER)/

clean:
	rm serial *.png gmp
