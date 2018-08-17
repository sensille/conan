#include <stdio.h>
#include <unistd.h>
#include <fcntl.h>
#include <termios.h>
#include <string.h>
#include <errno.h>
#include <stdlib.h>
#include <poll.h>
#include <assert.h>

#include "checksum.h"
#include "tmc2130.h"

typedef struct _comm {
	int	cm_fd;
	int	cm_timeout;
	uint8_t	cm_seq;
	uint8_t cm_recv_seq;
} comm_t;

int
open_terminal(const char *name, speed_t baud, int timeout, comm_t *cp)
{
	int fd;

	memset(cp, 0, sizeof(*cp));

	fd = open(name, O_RDWR); /* connect to port */

	/* set the other settings (in this case, 9600 8N1) */
	struct termios settings;
	tcgetattr(fd, &settings);
	cfmakeraw(&settings);
	cfsetspeed(&settings, baud); /* baud rate */
	settings.c_cflag |= CRTSCTS;
	tcsetattr(fd, TCSANOW, &settings); /* apply the settings */
	tcflush(fd, TCOFLUSH);

	cp->cm_fd = fd;
	cp->cm_timeout = timeout;
	cp->cm_seq = 0;

	return fd;
}

#define MAXBYTES 300

const char *
bytes(void *data, int len)
{
	int i;
	unsigned char *p = data;
	static char outbuf[MAXBYTES * 3 + 1];

	outbuf[0] = 0;
	for (i = 0; i < len; ++i)
		sprintf(outbuf + i * 3, "%02x ", p[i]);

	return outbuf;
}

int
escape_byte(unsigned char b, unsigned char *buf)
{
	if (b == 0x7d || b == 0x7e) {
		buf[0] = 0x7d;
		buf[1] = b ^ 0x20;
		return 2;
	}
	buf[0] = b;
	return 1;
}
int
send_frame(comm_t *cm, unsigned char *data, int datalen)
{
	unsigned char *src = data;
	/* worst case: seq, data, crc (possibly doubled by escaping), EOF */
	unsigned char buf[(datalen + 3) * 2 + 1];
	uint16_t crc;
	int i;
	unsigned char *p = buf;

	cm->cm_recv_seq = cm->cm_seq;
	cm->cm_seq = (cm->cm_seq + 1) & 0x7f;

	crc = update_crc_16(0, cm->cm_recv_seq);
	for (i = 0; i < datalen; ++i)
		crc = update_crc_16(crc, data[i]);

	p += escape_byte(cm->cm_recv_seq, p);
	while (datalen--)
		p += escape_byte(*src++, p);
	p += escape_byte(crc >> 8, p);
	p += escape_byte(crc & 0xff, p);
	*p++ = 0x7e;

printf("send %s\n", bytes(buf, p - buf));
	return write(cm->cm_fd, buf, p - buf);
}

int
read_all(comm_t *cm, unsigned char *buf, int *len, int maxlen)
{
	int ret;
	unsigned char c;
	struct pollfd pfd;

	*len = 0;

	pfd.fd = cm->cm_fd;
	pfd.events = POLLIN;
printf("read");
	while (1) {
		ret = poll(&pfd, 1, cm->cm_timeout);
		if (ret == 0)
{
printf("poll timed out, len %d\n", *len);
			return -1;
}
		if (ret != 1) {
			printf("unexpected return value from poll: %d\n", ret);
			exit(1);
		}
#if 0
		printf("revents: %x\n", pfd.revents);
#endif
		ret = read(cm->cm_fd, &c, 1);
		if (ret != 1) {
			printf("read failed with %d=%s\n", errno,
				strerror(errno));
			exit(1);
		}
		printf(" %02x", c);
		if (c == 0x7e) {
			printf("\n");
			return 1;
		}
		*buf++ = c;
		++*len;
		if (*len == maxlen) {
			printf("\n");
			return 0;
		}
	}
}

/*
 * timeout is the timeout per character
 * maxlen has to be large enough to hold the escaped packet including crc
 */
int
read_frame(comm_t *cm, unsigned char *outbuf, int *outlen, int maxlen)
{
	int ret;
	int in_overlong_packet = 0;
	unsigned char recvbuf[maxlen * 2 + 2];
	int recvlen = 0;
	int i;

	while (1) {
		ret = read_all(cm, recvbuf, &recvlen, maxlen * 2 + 2);
		if (ret == -1)	/* timeout */
			return -1;
		if (ret == 0) {	/* max len reached */
			in_overlong_packet = 1;
			goto next_packet; /* just keep on reading packets */
		}
		if (in_overlong_packet == 1) {
			/* end of packet reached, discard and start over */
			in_overlong_packet = 0;
			goto next_packet;
		}
		/* de-escape packet */
		int in_escape = 0;
		unsigned char *src = recvbuf;
		unsigned char *dst = recvbuf;
		for (i = 0; i < recvlen; ++i) {
			assert(*src != 0x7e);
			if (*src == 0x7d) {
				if (in_escape) {
					printf("bad esc char sequence\n");
					goto next_packet;
				}
				in_escape = 1;
				++src;
			} else if (in_escape) {
				*dst++ = *src++ ^ 0x20;
				in_escape = 0;
			} else { 
				*dst++ = *src++;
			}
		}
		*outlen = dst - recvbuf - 2;
		if (*outlen > maxlen + 1) {
			printf("oversized packet received\n");
			goto next_packet;
		}
		if (*outlen < 1) {
			printf("short packet received\n");
			goto next_packet;
		}
		/* calc crc */
		uint16_t crc = crc_16(recvbuf, *outlen);
		if (recvbuf[*outlen    ] != (crc >> 8) ||
		    recvbuf[*outlen + 1] != (crc & 0xff)) {
			printf("bad crc received\n");
			goto next_packet;
		}
		if (recvbuf[0] & 0x80) {
			printf("target reports receive error\n");
			goto next_packet;
		}
		if (recvbuf[0] != (cm->cm_recv_seq & 0x7f)) {
			/* XXX TODO retransmit */
			printf("unexpected seq: 0x%02x != 0x%02x\n",
				recvbuf[0], cm->cm_recv_seq);
			goto next_packet;
		}
		--*outlen; /* don't copy out seq */
		memcpy(outbuf, recvbuf + 1, *outlen);
		return 0;
next_packet:;
	}
}

void
check(comm_t *cm, void *in, int inlen, void *exp, int explen)
{
	int ret;
	int len_out;
	unsigned char read_buf[explen];

	printf("send %s ", bytes(in, inlen));
	printf("expect %s ", bytes(exp, explen));
	ret = send_frame(cm, in, inlen);
	if (ret < 0) {
		printf("sending data failed: %d=%s\n", errno, strerror(errno));
		exit(1);
	}
	ret = read_all(cm, read_buf, &len_out, explen);
	printf("received %s\n", bytes(read_buf, len_out));
	if (ret < 0) {
		printf("read_all failed\n");
		exit(1);
	}
	if (memcmp(read_buf, exp, explen) != 0) {
		printf("failed to read expected result\n");
		exit(1);
	}
}

int
sendrecv(comm_t *cm, void *in, int inlen, void *out, int *outlen, int maxout)
{
	int ret;

	printf("send %s ", bytes(in, inlen));
	ret = send_frame(cm, in, inlen);
	if (ret < 0) {
		printf("sending data failed: %d=%s\n", errno, strerror(errno));
		exit(1);
	}
	ret = read_frame(cm, out, outlen, maxout);
	printf("received %s\n", bytes(out, *outlen));
	if (ret < 0) {
		printf("read_frame failed\n");
		exit(1);
	}
	return 0;
}

/*
 * send packet with expected response len 0
 */
void
send0(comm_t *cm, void *in, int inlen)
{
	unsigned char recv_buf[1];
	int recv_len;

	sendrecv(cm, in, inlen, recv_buf, &recv_len, sizeof(recv_buf));
	assert(recv_len == 0);
}

void
resync(comm_t *cm)
{
	int ret;
	uint8_t recvbuf[100];
	int recvlen;

	/* terminate eventually running frame */
	ret = write(cm->cm_fd, "\x7e", 1);
	if (ret != 1) {
		printf("resync: sending EOF failed\n");
		exit(1);
	}

	/* read all pending data */
	while (1) {
		ret = read_all(cm, recvbuf, &recvlen, sizeof(recvbuf));
		if (ret == -1)	/* timeout */
			break;
	}
	cm->cm_seq = 0x80;
	send0(cm, "", 0);
	cm->cm_seq = 0x01;
}

static void
w32(uint8_t *out, uint32_t d)
{
	out[0] = (d >> 24) & 0xff;
	out[1] = (d >> 16) & 0xff;
	out[2] = (d >>  8) & 0xff;
	out[3] = (d >>  0) & 0xff;
}

static uint32_t
r32(uint8_t *d)
{
	return (d[0] << 24) + (d[1] << 16) + (d[2] << 8) + d[3];
}

void
spi(comm_t *cm, int chan,
	uint8_t cmd1, uint32_t data1,
	uint8_t cmd2, uint32_t data2,
	uint8_t cmd3, uint32_t data3,
	uint8_t *status1, uint32_t *read1,
	uint8_t *status2, uint32_t *read2,
	uint8_t *status3, uint32_t *read3)
{
	uint8_t wbuf[1 + 3 * 5];
	uint8_t rbuf[3 * 5] = { 0 };
	int rlen;
	int ret;

	wbuf[0] = 0x80 + chan;
	wbuf[1] = cmd1;
	w32(wbuf + 2, data1);
	wbuf[6] = cmd2;
	w32(wbuf + 7, data2);
	wbuf[11] = cmd3;
	w32(wbuf + 12, data3);
	
	send0(cm, wbuf, 16);

	/* recv response data */
	ret = read_frame(cm, rbuf, &rlen, sizeof(rbuf));
	printf("received %s\n", bytes(rbuf, rlen));
	assert(ret == 0);
	assert(rlen == 15);

	if (status1)
		*status1 = rbuf[0];
	if (read1)
		*read1 = r32(rbuf + 1);
	if (status2)
		*status2 = rbuf[5];
	if (read2)
		*read2 = r32(rbuf + 6);
	if (status3)
		*status3 = rbuf[10];
	if (read3)
		*read3 = r32(rbuf + 11);
}

void
spi_read(comm_t *cm, int chan, uint8_t cmd,
	uint8_t *status1, uint32_t *data1,
	uint8_t *status2, uint32_t *data2,
	uint8_t *status3, uint32_t *data3)
{
	assert((cmd & 0x80) == 0);
	/* send command */
	spi(cm, chan, cmd, 0, cmd, 0, cmd, 0,
		NULL, NULL, NULL, NULL, NULL, NULL);
	/* read result */
	spi(cm, chan, 0, 0, 0, 0, 0, 0, status1, data1, status2, data2,
		status3, data3);
}

void
spi_write(comm_t *cm, int chan,
	uint8_t cmd1, uint32_t data1,
	uint8_t cmd2, uint32_t data2,
	uint8_t cmd3, uint32_t data3)
{
	spi(cm, chan, cmd1, data1, cmd2, data2, cmd3, data3,
		NULL, NULL, NULL, NULL, NULL, NULL);
}

void
tmcw(comm_t *cm, int chanmask, uint8_t reg, uint32_t data)
{
	if (chanmask & 0x07) {
		spi_write(cm, 0,
			(chanmask & 0x01) ? (reg | 0x80) : 0, data,
			(chanmask & 0x02) ? (reg | 0x80) : 0, data,
			(chanmask & 0x04) ? (reg | 0x80) : 0, data);
	}
	if (chanmask & 0x38) {
		spi_write(cm, 1,
			(chanmask & 0x08) ? (reg | 0x80) : 0, data,
			(chanmask & 0x10) ? (reg | 0x80) : 0, data,
			(chanmask & 0x20) ? (reg | 0x80) : 0, data);
	}
}

void
test(comm_t *cm)
{
	check(cm, "\x80\x55", 2, "\x00\x00", 2);
	check(cm, "\x01zsdf", 5, "\x00\x01", 2);
	check(cm, "\x01zsdf", 5, "\x03\x01", 2);
	check(cm, "\x02zfddffddffdd", 13, "\x00\x02", 2);
}

void
print_status(comm_t *cm)
{
	uint8_t s[6];	/* STATUS */
	uint32_t v[6];	/* DRV_STATUS */
	uint32_t g[6];	/* GSTAT */
	int i;

	spi_read(cm, 0, 0x6f, &s[0], &v[0], &s[1], &v[1], &s[2], &v[2]);
	spi_read(cm, 1, 0x6f, &s[3], &v[3], &s[4], &v[4], &s[5], &v[5]);
	spi_read(cm, 0, 0x01, NULL, &g[0], NULL, &g[1], NULL, &g[2]);
	spi_read(cm, 1, 0x01, NULL, &g[3], NULL, &g[4], NULL, &g[5]);

	printf("status:");
	for (i = 0; i < 6; ++i)
		printf(" %02x:%08x-%x", s[i], v[i], g[i]);
	printf("\n");
}

int
main(int argc, char **argv)
{
	int ret;
	comm_t cc;	/* control connection */
	comm_t cd;	/* data connection */

	ret = open_terminal("/dev/ttyS1", B115200, 100, &cd);
	if (ret < 0) {
		printf("failed to open serial port S1: %s\n", strerror(errno));
		exit(1);
	}
	ret = open_terminal("/dev/ttyS2", B115200, 100, &cc);
	if (ret < 0) {
		printf("failed to open serial port S2: %s\n", strerror(errno));
		exit(1);
	}

#if 0
	test(fd);
#endif

#if 0
	for (int i = 0; i < 10000; ++i) {
		printf("sending block %d\n", i);
		int ret = send_frame(fd, "\x80\x01\x00\x00\x00", 5);
		if (ret < 0) {
			printf("sending data failed: %d=%s\n", errno,
				strerror(errno));
			exit(1);
		}
	}
#endif

printf("control channel\n");
	resync(&cc);

	/* read tmc2130 versions */
	uint32_t v[6];
	int i;
	spi_read(&cc, 0, 0x04, NULL, &v[0], NULL, &v[1], NULL, &v[2]);
	spi_read(&cc, 1, 0x04, NULL, &v[3], NULL, &v[4], NULL, &v[5]);

	for (i = 0; i < 6; ++i)
		v[i] >>= 24;

	printf("versions:");
	for (i = 0; i < 6; ++i)
		printf(" %02x", v[i]);
	printf("\n");

	for (i = 0; i < 6; ++i) {
		if (v[i] != 0x11) {
			printf("at least one version is bad\n");
			exit(1);
		}
	}

#if 0
	/* disable all drivers upfront */
	tmcw(&cc, 0x3f, TMCR_GCONF, TMC_GCONF_STOP_ENABLE);
#endif

	int cmask = 0x3f;
	/* configuration without StealthChop */
	tmcw(&cc, cmask, TMCR_CHOPCONF,
		TMC_CHOPCONF_DEDGE |
		(2 << TMC_CHOPCONF_TBL_SHIFT) |
		(7 << TMC_CHOPCONF_HEND_SHIFT) |
		(2 << TMC_CHOPCONF_HSTRT_SHIFT) |
		(5 << TMC_CHOPCONF_TOFF_SHIFT));
	tmcw(&cc, cmask, TMCR_IHOLD_IRUN,
		( 6 << TMC_IHOLD_IRUN_IHOLDDELAY_SHIFT) |
#if 0
		(24 << TMC_IHOLD_IRUN_IRUN_SHIFT) |
#else
		(10 << TMC_IHOLD_IRUN_IRUN_SHIFT) |	/* XXX differ for E */
#endif
		( 3 << TMC_IHOLD_IRUN_IHOLD_SHIFT));
	tmcw(&cc, cmask, TMCR_TPOWER_DOWN, 0x0a);
	tmcw(&cc, cmask, TMCR_GCONF,
		TMC_GCONF_DIAG0_INT_PUSHPULL |
		TMC_GCONF_DIAG1_PUSHPULL |
		TMC_GCONF_DIAG1_STALL |
		TMC_GCONF_DIAG0_OTPW
#if 1
		| TMC_GCONF_EN_PWM_MODE
#endif
	);
	tmcw(&cc, cmask, TMCR_TPWMTHRS, 0);
	tmcw(&cc, cmask, TMCR_PWMCONF,
		TMC_PWMCONF_PWM_AUTOSCALE |
		(  1 << TMC_PWMCONF_PWM_GRAD_SHIFT) |
		(200 << TMC_PWMCONF_PWM_AMPL_SHIFT));

	/* enable, GPOUT 0 to Lo */
	send0(&cc, "\x71\x00", 2);

	/* stop motion controller */
	send0(&cc, "\x61", 1);

printf("data channel\n");
	resync(&cd);
	/* reset motion controller */
	send0(&cd, "\x40", 1);
	/* routing: controller 1 to outputs 0, 1, 2 */
#if 0
	send0(&cd, "\x60\x00\x01", 3);
	send0(&cd, "\x60\x01\x01", 3);
	send0(&cd, "\x60\x02\x01", 3);
	send0(&cd, "\x60\x03\x01", 3);
	send0(&cd, "\x60\x04\x01", 3);
	send0(&cd, "\x60\x05\x01", 3);
#else
	send0(&cd, "\x60\x01\x01", 3);
#endif

	send0(&cd, "\x6a\x00", 2);	/* endstop mask */
	send0(&cd, "\x6b\x00", 2);	/* endstop polarity */

	/* notify */
	send0(&cd, "\x61\x55\xaa", 3);

	/* positive jerk */
	send0(&cd, "\x80\x55\xaa\x00", 4);
	send0(&cd, "\x71\x02\x00\x00", 4);

	/* negative jerk */
	send0(&cd, "\x80\xaa\x55\x01", 4);
	send0(&cd, "\x71\x02\x00\x00", 4);

	/* constant velocity */
	send0(&cd, "\x80\x00", 1);
	send0(&cd, "\x71\x08\x00\x00\x00", 5);

#if 1
	/*
	 * set wait mask so the currently running command can be interrupted
	 * by an endstop
	 */
	send0(&cd, "\x69\x01", 2); 
#endif

	/* negative jerk */
	send0(&cd, "\x80\xaa\x55\x01", 4);
	send0(&cd, "\x71\x02\x00\x00", 4);

	/* positive jerk again */
	send0(&cd, "\x80\x55\xaa\x00", 4);
	send0(&cd, "\x71\x02\x00\x00", 4);

	/* notify */
	send0(&cd, "\x61\x7a\xbb", 3);
	/* stop */
	send0(&cd, "\x41", 1);

	/* start motion controller */
	send0(&cc, "\x60", 1);

	/* wait for notifications */
	while (1) {
		print_status(&cc);

		unsigned char recv_buf[20];
		int recv_len = 0;
		ret = read_frame(&cd, recv_buf, &recv_len, sizeof(recv_buf));
		if (ret == 0) {
			printf("received %s\n", bytes(recv_buf, recv_len));
			if (recv_len == 4 &&
			    memcmp(recv_buf, "\x00\x00\x7a\xbb", 4) == 0)
				break;
		}
	}
	sleep(1);
	send0(&cc, "\x70\x00", 2);

	return 0;
}

