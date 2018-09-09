all: serial

SRC = serial.c crc16.c main.c
HDR = checksum.h serial.h
serial: $(SRC) $(HDR)
	gcc -g -Wall $(SRC) -lm -o serial

clean:
	rm serial
