all: motion serial gmp

SRC = serial.c crc16.c main.c
HDR = checksum.h serial.h
serial: $(SRC) $(HDR)
	gcc -g -Wall $(SRC) -lm -o serial

gmp: gmp.c Makefile
	gcc -Wall -g gmp.c -o gmp -lm -lgmp -lmpfr

motion: motion.c Makefile
	gcc -Wall -g motion.c -o motion -lm -lgmp -lmpfr

plot:
	gnuplot plot.gp
	cp index.html *.png /var/www/html/3dpfs/

clean:
	rm serial *.png gmp
