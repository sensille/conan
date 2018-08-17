all: serial

SRC = serial.c crc16.c
HDR = checksum.h
serial: $(SRC) $(HDR)
	gcc -g -Wall $(SRC) -o serial

clean:
	rm serial
