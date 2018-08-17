#ifndef __CRC16__
#define __CRC16__

#include <stdint.h>

#define		CRC_START_16		0x0000

#define		CRC_POLY_16		0xA001
#define		CRC_POLY_CCITT		0x1021

uint16_t crc_16(const unsigned char *input_str, size_t num_bytes);
uint16_t update_crc_16( uint16_t crc, unsigned char c );

#endif
