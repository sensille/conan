#ifndef __SERIAL_H__
#define __SERIAL_H__

#include <stdint.h>
#include <termios.h>

typedef struct _comm {
	int	cm_fd;
	int	cm_timeout;
	uint8_t	cm_seq;
	uint8_t cm_recv_seq;
} comm_t;

const char *bytes(void *data, int len);

int open_terminal(const char *name, speed_t baud, int timeout, comm_t *cp);
int read_frame(comm_t *cm, unsigned char *outbuf, int *outlen, int maxlen);

int
sendrecv(comm_t *cm, void *in, int inlen, void *out, int *outlen, int maxout);

void send0(comm_t *cm, void *in, int inlen);
void resync(comm_t *cm);

void
spi(comm_t *cm, int chan,
	uint8_t cmd1, uint32_t data1,
	uint8_t cmd2, uint32_t data2,
	uint8_t cmd3, uint32_t data3,
	uint8_t *status1, uint32_t *read1,
	uint8_t *status2, uint32_t *read2,
	uint8_t *status3, uint32_t *read3);

void
spi_read(comm_t *cm, int chan, uint8_t cmd,
	uint8_t *status1, uint32_t *data1,
	uint8_t *status2, uint32_t *data2,
	uint8_t *status3, uint32_t *data3);

void
spi_write(comm_t *cm, int chan,
	uint8_t cmd1, uint32_t data1,
	uint8_t cmd2, uint32_t data2,
	uint8_t cmd3, uint32_t data3);

void tmcw(comm_t *cm, int chanmask, uint8_t reg, uint32_t data);
int print_status(comm_t *cm);

#endif
