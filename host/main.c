#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <errno.h>
#include <math.h>

#include "serial.h"
#include "tmc2130.h"

#define C_OFF	0
#define C_E	1
#define C_XY1	2
#define C_XY2	3
#define C_Z	4

#define R_E	0
#define R_XY1	1
#define R_XY2	2
#define R_Z1	3
#define R_Z2	4
#define R_Z3	5

#define E_Y	1
#define E_X	2
#define E_Z	4

#define CMD_LOADJERK	0x80
#define CMD_LOADSNAP	0x90
#define CMD_LOADCRACKLE	0xa0

#define MAX_PARAM_STR 100
#define PARAM_SIZE (208 / 8)

static uint64_t hz;
static int steps_per_mm;
static double step;
static double f_spm;

static void
print_tmc2130_versions(comm_t *cm)
{
	/* read tmc2130 versions */
	uint32_t v[6];
	int i;
	spi_read(cm, 0, 0x04, NULL, &v[0], NULL, &v[1], NULL, &v[2]);
	spi_read(cm, 1, 0x04, NULL, &v[3], NULL, &v[4], NULL, &v[5]);

	for (i = 0; i < 6; ++i)
		v[i] >>= 24;

	printf("versions:");
	for (i = 0; i < 6; ++i)
		printf(" %02x", v[i]);
	printf("\n");

	for (i = 0; i < 3; ++i) { /* XXX */
		if (v[i] != 0x11) {
			printf("at least one version is bad\n");
			exit(1);
		}
	}
}

static void
init_drivers(comm_t *cm)
{
	int cmask = 0x3f;
	int cmask_xy = 0x06;
	int ret;

	/* disable drivers, TOFF=0 */
	tmcw(cm, cmask, TMCR_CHOPCONF,
		TMC_CHOPCONF_DEDGE |
		(2 << TMC_CHOPCONF_TBL_SHIFT) |
		(7 << TMC_CHOPCONF_HEND_SHIFT) |
		(2 << TMC_CHOPCONF_HSTRT_SHIFT) |
		(0 << TMC_CHOPCONF_TOFF_SHIFT));

	/* clear GSTAT by reading */
	print_status(cm);

	/* configuration without StealthChop */
	tmcw(cm, cmask, TMCR_CHOPCONF,
		TMC_CHOPCONF_DEDGE |
#if 0 /* bad idea, let magic smoke out */
		TMC_CHOPCONF_DISS2G |
#endif
		(2 << TMC_CHOPCONF_TBL_SHIFT) |
		(7 << TMC_CHOPCONF_HEND_SHIFT) |
		(2 << TMC_CHOPCONF_HSTRT_SHIFT) |
		(5 << TMC_CHOPCONF_TOFF_SHIFT));
	/* all but X/Y */
	tmcw(cm, cmask & ~cmask_xy, TMCR_IHOLD_IRUN,
		( 6 << TMC_IHOLD_IRUN_IHOLDDELAY_SHIFT) |
		(10 << TMC_IHOLD_IRUN_IRUN_SHIFT) |
		( 3 << TMC_IHOLD_IRUN_IHOLD_SHIFT));
	/* X/Y */
	tmcw(cm, cmask_xy, TMCR_IHOLD_IRUN,
		( 6 << TMC_IHOLD_IRUN_IHOLDDELAY_SHIFT) |
		(24 << TMC_IHOLD_IRUN_IRUN_SHIFT) |
		( 3 << TMC_IHOLD_IRUN_IHOLD_SHIFT));
	tmcw(cm, cmask, TMCR_TPOWER_DOWN, 0x0a);
	tmcw(cm, cmask, TMCR_GCONF,
		TMC_GCONF_DIAG0_INT_PUSHPULL |
		TMC_GCONF_DIAG1_PUSHPULL |
		TMC_GCONF_DIAG1_STALL |
		TMC_GCONF_DIAG0_OTPW
#if 0
		| TMC_GCONF_EN_PWM_MODE
#endif
	);
	tmcw(cm, cmask, TMCR_TPWMTHRS, 0);
	tmcw(cm, cmask, TMCR_PWMCONF,
		TMC_PWMCONF_PWM_AUTOSCALE |
		(  1 << TMC_PWMCONF_PWM_GRAD_SHIFT) |
		(200 << TMC_PWMCONF_PWM_AMPL_SHIFT));

	/* check if any driver reports an error */
	ret = print_status(cm);
	if (ret != 0)
		exit(1);
}

static void
send_data_cmd(comm_t *cd, uint8_t cmd, int64_t val)
{
	uint8_t data[9];
	uint8_t *p = data;
	int i;

	for (i = 0; i < 7; ++i) {
		uint16_t msh = val >> 48;

		if ((msh & 0xff80) != 0xff80 &&
		    (msh & 0xff80) != 0)
			break;

		val <<= 8;
	}
	*p++ = cmd;
	for (; i < 8; ++i) {
		*p++ = val >> 56;
		val <<= 8;
	}

	send0(cd, data, p - data);
}

/*
 * route channel chan to controller cntrl
 */
static void
m_set_routing(comm_t *cd, int cntrl, int chan)
{
	send_data_cmd(cd, 0x60, (chan << 8) | cntrl);
}

static void
m_preload_jerk(comm_t *cd, int cntrl, int64_t jerk)
{
	send_data_cmd(cd, 0x80 + cntrl - 1, jerk);
}

static void
m_load_cnt(comm_t *cd, int64_t cnt)
{
	send_data_cmd(cd, 0x71, cnt);
}

static void
m_notify(comm_t *cd, int64_t code)
{
	send_data_cmd(cd, 0x61, code);
}

static void
m_wait_endstop(comm_t *cd, int endstop)
{
	send_data_cmd(cd, 0x69, endstop);
}

static void
send_longreg_cmd(comm_t *cd, uint8_t cmd, uint8_t val[PARAM_SIZE])
{
	uint8_t data[PARAM_SIZE + 1];
	int i;

	for (i = 0; i < PARAM_SIZE - 2; ++i) {
		if (((val[i] == 0xff) && ( val[i + 1] & 0x80)) ||
		    ((val[i] == 0x00) && (~val[i + 1] & 0x80)))
			continue;
		break;
	}
	data[0] = cmd;
	memcpy(data + 1, val + i, PARAM_SIZE - i);

	send0(cd, data, PARAM_SIZE - i + 1);
}

/*
 * we always use controller 0 for homing, routing to the respective steppers
 */
void
home_z(comm_t *cc, comm_t *cd)
{
	int cntrl = C_Z;

	/* disable all endstops XXX */
	send0(cd, "\x6a\x04", 2);	/* endstop mask */

	/* positive jerk */
	m_preload_jerk(cd, cntrl, 0x55aa00);
	m_load_cnt(cd, 0x20000);

	/* negative jerk */
	m_preload_jerk(cd, cntrl, -0x55aa00);
	m_load_cnt(cd, 0x20000);

	/* constant velocity */
	m_preload_jerk(cd, cntrl, 0);
	m_load_cnt(cd, 0x20000000);

	/*
	 * set wait mask so the currently running command can be interrupted
	 * by an endstop
	 */
	m_wait_endstop(cd, E_Z);

	/* negative jerk */
	m_preload_jerk(cd, cntrl, -0x55aa00);
	m_load_cnt(cd, 0x20000);

	/* positive jerk again */
	m_preload_jerk(cd, cntrl, 0x55aa00);
	m_load_cnt(cd, 0x20000);

	/* notify */
	m_notify(cd, 0x7abb);

	/* stop */
	send0(cd, "\x41", 1);

	/* start motion controller */
	send0(cc, "\x60", 1);

	/* wait for notifications */
	while (1) {
		int ret;
#if 1
		print_status(cc);
#endif

		unsigned char recv_buf[20];
		int recv_len = 0;
		ret = read_frame(cd, recv_buf, &recv_len, sizeof(recv_buf));
		if (ret == 0) {
			printf("received %s\n", bytes(recv_buf, recv_len));
			if (recv_len == 4 &&
			    memcmp(recv_buf, "\x00\x00\x7a\xbb", 4) == 0)
				break;
		}
	}

	m_preload_jerk(cd, cntrl, 0);
}

void
move_z(comm_t *cc, comm_t *cd, double stroke)
{
	int cntrl = C_Z;
#if 0
	int64_t v_jerk;
	int64_t t_jerk;
#endif

	/* enable all endstops XXX TODO */
	send0(cd, "\x6a\x00", 2);	/* endstop mask */

	/* positive jerk */
	m_preload_jerk(cd, cntrl, -0x55aa00);
	m_load_cnt(cd, 0x30000);

	/* negative jerk */
	m_preload_jerk(cd, cntrl, 0x55aa00);
	m_load_cnt(cd, 0x30000);

	/* constant velocity */
	m_preload_jerk(cd, cntrl, 0);
	m_load_cnt(cd, 0x2000000);

	/* negative jerk */
	m_preload_jerk(cd, cntrl, 0x55aa00);
	m_load_cnt(cd, 0x30000);

	/* positive jerk again */
	m_preload_jerk(cd, cntrl, -0x55aa00);
	m_load_cnt(cd, 0x30000);

	/* notify */
	m_notify(cd, 0x7abb);

	/* stop */
	send0(cd, "\x41", 1);

	/* start motion controller */
	send0(cc, "\x60", 1);

	/* wait for notifications */
	while (1) {
		int ret;
#if 1
		print_status(cc);
#endif

		unsigned char recv_buf[20];
		int recv_len = 0;
		ret = read_frame(cd, recv_buf, &recv_len, sizeof(recv_buf));
		if (ret == 0) {
			printf("received %s\n", bytes(recv_buf, recv_len));
			if (recv_len == 4 &&
			    memcmp(recv_buf, "\x00\x00\x7a\xbb", 4) == 0)
				break;
		}
	}

	m_preload_jerk(cd, cntrl, 0);
}

void
home_xy(comm_t *cc, comm_t *cd, int xy)
{
	int64_t jv = 0x35aa00;
	int64_t jerk1;
	int64_t jerk2;

	if (xy == 0) {
		jerk1 = -jv;
		jerk2 = jv;
	} else {
		jerk1 = -jv;
		jerk2 = -jv;
	}

	/* disable all endstops */
	send0(cd, "\x6a\x00", 2);	/* endstop mask */

	/* positive jerk */
	m_preload_jerk(cd, C_XY1, jerk1);
	m_preload_jerk(cd, C_XY2, jerk2);
	m_load_cnt(cd, 0x20000);

	/* negative jerk */
	m_preload_jerk(cd, C_XY1, -jerk1);
	m_preload_jerk(cd, C_XY2, -jerk2);
	m_load_cnt(cd, 0x20000);

	/* constant velocity */
	m_preload_jerk(cd, C_XY1, 0);
	m_preload_jerk(cd, C_XY2, 0);
	m_load_cnt(cd, 0x20000000);

	/*
	 * set wait mask so the currently running command can be interrupted
	 * by an endstop
	 */
	m_wait_endstop(cd, xy ? E_Y : E_X);

	/* negative jerk */
	m_preload_jerk(cd, C_XY1, -jerk1);
	m_preload_jerk(cd, C_XY2, -jerk2);
	m_load_cnt(cd, 0x20000);

	/* positive jerk again */
	m_preload_jerk(cd, C_XY1, jerk1);
	m_preload_jerk(cd, C_XY2, jerk2);
	m_load_cnt(cd, 0x20000);

	/* notify */
	m_notify(cd, 0x7abb);

	/* stop */
	send0(cd, "\x41", 1);

	/* start motion controller */
	send0(cc, "\x60", 1);

	/* wait for notifications */
	while (1) {
		int ret;
#if 1
		print_status(cc);
#endif

		unsigned char recv_buf[20];
		int recv_len = 0;
		ret = read_frame(cd, recv_buf, &recv_len, sizeof(recv_buf));
		if (ret == 0) {
			printf("received %s\n", bytes(recv_buf, recv_len));
			if (recv_len == 4 &&
			    memcmp(recv_buf, "\x00\x00\x7a\xbb", 4) == 0)
				break;
		}
	}

	m_preload_jerk(cd, C_XY1, 0);
	m_preload_jerk(cd, C_XY2, 0);
}

void
s_curve(double j, double a, double v, double s,
	uint64_t *ptj, uint64_t *pta, uint64_t *ptv)
{
	double tj = a / j;	/* time for jerk */
	double ta = v/a - tj;	/* time for acceleration */

	if (ta < 0) {
		printf("full acceleration can't be reached\n");
		tj = sqrt(v / j);
		ta = 0;
	}

	double tv = (s - j*tj * (2*tj*tj + 3*tj*ta + ta*ta)) / v;

	if (tv < 0) {
		printf("full velocity can't be reached\n");
		if (ta == 0) {
			tj = cbrt(s / (2 * j));
		} else {
			ta = -3.0 / 2 * tj + sqrt(tj*tj / 4 + s / (j * tj));
			if (ta < 0) {
				printf("need to reduce tj\n");
				tj = cbrt(s / (2 * j));
				ta = 0;
			}
		}
		tv = 0;
	}

	/* translate into controler values */
	uint64_t c_tj = tj * hz;
	uint64_t c_ta = ta * hz;
	uint64_t c_tv = tv * hz;

	printf("tj: %ld\n", c_tj);
	printf("ta: %ld\n", c_ta);
	printf("tv: %ld\n", c_tv);

	if (tj <= 0 || ta < 0 || tv < 0) {
		printf("invalid parameters, can't move\n");
		exit(1);
	}

	*ptj = c_tj;
	*pta = c_ta;
	*ptv = c_tv;
}

void
jerk(comm_t *cd, int64_t c_ja, int64_t c_jb, uint64_t c_tj)
{
	printf("jerk: ja %ld jb %ld tj %ld\n", c_ja, c_jb, c_tj);

	if (c_tj == 0)
		return;

	m_preload_jerk(cd, C_XY1, c_ja);
	m_preload_jerk(cd, C_XY2, c_jb);
	m_load_cnt(cd, c_tj);
}

void
check_status(comm_t *cc, comm_t *cd)
{
	int ret = print_status(cc);
	if (ret) {
		/* disable */
		send0(cd, "\x40", 1);
		send0(cc, "\x71\x01", 2);
		exit(1);
	}
}

void
move_xy(comm_t *cc, comm_t *cd, double x, double y,
	double v, double j, double a)
{
	uint64_t c_tj;
	uint64_t c_ta;
	uint64_t c_tv;

	/* total travel distance */
	double s = sqrt(x * x + y * y);

	/* project on corexy */
	double xy_a =  x + y;
	double xy_b = -x + y;
	
	s_curve(j, a, v, s, &c_tj, &c_ta, &c_tv);

	double ja = j * xy_a / s;
	double jb = j * xy_b / s;

	printf("xy_a %lf yz_b %lf\n", xy_a, xy_b);
	printf("s %lf v %lf\n", s, v);
	printf("ja %lf jb %lf\n", ja, jb);

	printf("steps_per_mm %d f_spm %lg\n", steps_per_mm, f_spm);
	int64_t c_ja = ja * f_spm / hz / hz;
	int64_t c_jb = jb * f_spm / hz / hz;

	printf("ja: %ld\n", c_ja);
	printf("jb: %ld\n", c_jb);
	printf("tj: %ld\n", c_tj);
	printf("ta: %ld\n", c_ta);
	printf("tv: %ld\n", c_tv);

	jerk(cd,  c_ja,  c_jb, c_tj);
	jerk(cd,     0,     0, c_ta);
	jerk(cd, -c_ja, -c_jb, c_tj);
	jerk(cd,     0,     0, c_tv);
	jerk(cd, -c_ja, -c_jb, c_tj);
	jerk(cd,     0,     0, c_ta);
	jerk(cd,  c_ja,  c_jb, c_tj);
}

void
move_square(comm_t *cc, comm_t *cd, double s, double v, double j, double a,
	int rounds)
{
	uint64_t c_tj;
	uint64_t c_ta;
	uint64_t c_tv;
	double x = s;
	double y = 0;
	int i;

	s_curve(j, a, v, s, &c_tj, &c_ta, &c_tv);

	/* project on corexy */
	int64_t c_ja = j * ( x + y) / s * f_spm / hz / hz;
	int64_t c_jb = j * (-x + y) / s * f_spm / hz / hz;

	/* lead in */
	jerk(cd,  c_ja,  c_jb, c_tj);
	jerk(cd,     0,     0, c_ta);
	jerk(cd, -c_ja, -c_jb, c_tj);
	jerk(cd,     0,     0, c_tv);

	/* start motion controller */
	send0(cc, "\x60", 1);

	for (i = 0; i < (rounds * 4 - 1); ++i) {
		switch (i % 4) {
			case 0: x =  0; y =  s; break;
			case 1: x = -s; y =  0; break;
			case 2: x =  0; y = -s; break;
			case 3: x =  s; y =  0; break;
		}

		int64_t n_ja = j * ( x + y) / s * f_spm / hz / hz;
		int64_t n_jb = j * (-x + y) / s * f_spm / hz / hz;

		jerk(cd, -c_ja + n_ja, -c_jb + n_jb, c_tj);
		jerk(cd,     0,     0, c_ta);
		jerk(cd,  c_ja - n_ja,  c_jb - n_jb, c_tj);

		jerk(cd,     0,     0, c_tv);

		c_ja = n_ja;
		c_jb = n_jb;

		check_status(cc, cd);
	}
	/* lead out */
	jerk(cd, -c_ja, -c_jb, c_tj);
	jerk(cd,     0,     0, c_ta);
	jerk(cd,  c_ja,  c_jb, c_tj);
}

void
s_curve_test_one(double j, double a, double v, double s)
{
	printf("j %lf a %lf v %lf s %lf\n", j, a, v, s);

	double tj = a / j;	/* time for jerk */
	double ta = v/a - tj;	/* time for acceleration */

	if (ta < 0) {
		printf("full acceleration can't be reached\n");
		tj = sqrt(v / j);
		ta = 0;
	}

	double tv = (s - j*tj * (2*tj*tj + 3*tj*ta + ta*ta)) / v;

	if (tv < 0) {
		printf("full velocity can't be reached: tv %lf ta %lf\n", tv, ta);
		if (ta == 0) {
			tj = cbrt(s / (2 * j));
		} else {
			ta = -3.0 / 2 * tj + sqrt(tj*tj / 4 + s / (j * tj));
			if (ta < 0) {
				printf("need to reduce tj\n");
				tj = cbrt(s / (2 * j));
				ta = 0;
			}
		}
		tv = 0;
	}
		
	printf("tj %lf ta %lf tv %lf\n", tj, ta, tv);

	if (tj <= 0 || ta < 0 || tv < 0) {
		printf("invalid parameters, can't move\n");
		exit(1);
	}

	double a_is = j * tj;
	double v_is = j * tj * (ta + tj);
	double s_is = j*tj*(2*tj*tj + 3*tj*ta + ta*ta) + tv*v_is;
	printf("is: a %lf v %lf s %lf\n", a_is, v_is, s_is);
}

void
s_curve_test(void)
{
	double j = 60000;	/* mm/s^3 */
	double a = 6000;	/* mm/s^2 */
	double s = 100;		/* mm */
	double v;

	for (v = 50; v < 1000; v += 50)
		s_curve_test_one(j, a, v, s);
}

void
square_test(comm_t *cc, comm_t *cd)
{
	double v = 150;
	double a = 20000;
	double j = 500000;
	int i;

	for (i = 0; i < 20; ++i) {
		move_xy(cc, cd,  120,    0, v, j, a);
		if (i == 0) {
			/* start motion controller */
			send0(cc, "\x60", 1);
		}
		move_xy(cc, cd,    0,  120, v, j, a);
		move_xy(cc, cd, -120,    0, v, j, a);
		move_xy(cc, cd,    0, -120, v, j, a);
		v += 50;
	}
}

void
line_test(comm_t *cc, comm_t *cd)
{
	double v = 400;
	double a = 10000;
	double j = 10000;
	int i;

	for (i = 0; i < 50; ++i) {
		if ((i % 6) < 3)
			j = 500000;
		else
			j = 200000000;

		move_xy(cc, cd,  20,    0, v, j, a);
		if (i == 0) {
			/* start motion controller */
			send0(cc, "\x60", 1);
		}
		move_xy(cc, cd,  -20,   0, v, j, a);
	}
}

void
star_test(comm_t *cc, comm_t *cd, double j)
{
	int i;
	int ret;

	double v = 400;
	double a = 10000;

	double r = 60;
	double step = 10;
	double angle = 0;
	double last_x = 0;
	double last_y = 0;

	for (i = 0; i < 300; ++i) {
		double next_x = sin(angle) * r;
		double next_y = cos(angle) * r;

		printf("move X %lf y %lf\n", next_x - last_x, next_y - last_y);

		move_xy(cc, cd,  next_x - last_x, next_y - last_y, v, j, a);
		if (i == 0) {
			/* start motion controller */
			send0(cc, "\x60", 1);
		}
		last_x = next_x;
		last_y = next_y;
		angle += (180 - step) * 2 * M_PI / 360;
		if ((i % 20) == 0)
			r *= 3.0 / 4;
		ret = print_status(cc);
		if (ret) {
			/* disable */
			send0(cd, "\x40", 1);
			send0(cc, "\x71\x01", 2);
			exit(1);
		}
	}
	move_xy(cc, cd,  -last_x, -last_y, v, j, a);
}

void
parse_param(uint8_t v[PARAM_SIZE], char *s_v)
{
	int l = strlen(s_v);
	uint8_t *dst = v + PARAM_SIZE - 1;
	int upper = 0;
	int neg = 0;
	int i;

	memset(v, 0, PARAM_SIZE);

	for (i = l - 1; i >= 0; --i) {
		char c = s_v[i];
		uint8_t d = c >= 'a' ? c - 'a' + 10 : c - '0';

		if (c == '-') {
			if (i != 0) {
				printf("invalid - in parameter %s\n", s_v);
				abort();
			}
			neg = 1;
		} else if (upper) {
			*dst |= d << 4;
			--dst;
		} else {
			if (dst < v) {
				printf("overlong parameter %s\n", s_v);
				abort();
			}
			*dst = d;
		}
		upper = !upper;
	}

	if (neg) {
		uint16_t result;

		for (i = 0; i < PARAM_SIZE; ++i)
			v[i] = v[i] ^ 255;
		for (i = PARAM_SIZE - 1; i >= 0; --i) {
			result = v[i] + 1;
			v[i] = result & 0xff;
			if (result <= 255)
				break;
		}
	}

	printf("parsed %53s\n to     ", s_v);
	for (i = 0; i < PARAM_SIZE; ++i)
		printf("%.02x", v[i]);
	printf("\n");
}

void
test_from_file(comm_t *cc, comm_t *cd, char *fname)
{
	FILE *fp = fopen(fname, "r");
	unsigned int steps;
	char s_jx[MAX_PARAM_STR];
	char s_sx[MAX_PARAM_STR];
	char s_cx[MAX_PARAM_STR];
	char s_jy[MAX_PARAM_STR];
	char s_sy[MAX_PARAM_STR];
	char s_cy[MAX_PARAM_STR];
	uint8_t jx[PARAM_SIZE];
	uint8_t sx[PARAM_SIZE];
	uint8_t cx[PARAM_SIZE];
	uint8_t jy[PARAM_SIZE];
	uint8_t sy[PARAM_SIZE];
	uint8_t cy[PARAM_SIZE];
	int ret;
	int i = 0;

	if (fp == NULL) {
		printf("failed to open input file\n");
		exit(1);
	}

	while(++i) {
		ret = fscanf(fp, "%d %s %s %s %s %s %s\n",
			&steps, s_jx, s_sx, s_cx, s_jy, s_sy, s_cy);
		if (ret == EOF)
			break;
		if (ret != 7) {
			printf("failed to match line from input file\n");
			exit(1);
		}
		parse_param(jx, s_jx);
		parse_param(sx, s_sx);
		parse_param(cx, s_cx);
		parse_param(jy, s_jy);
		parse_param(sy, s_sy);
		parse_param(cy, s_cy);

		send_longreg_cmd(cd, CMD_LOADJERK    + C_XY1 - 1, jx);
		send_longreg_cmd(cd, CMD_LOADSNAP    + C_XY1 - 1, sx);
		send_longreg_cmd(cd, CMD_LOADCRACKLE + C_XY1 - 1, cx);
		send_longreg_cmd(cd, CMD_LOADJERK    + C_XY2 - 1, jy);
		send_longreg_cmd(cd, CMD_LOADSNAP    + C_XY2 - 1, sy);
		send_longreg_cmd(cd, CMD_LOADCRACKLE + C_XY2 - 1, cy);
		m_load_cnt(cd, steps);

		ret = print_status(cc);
		if (ret) {
			/* disable */
			send0(cd, "\x40", 1);
			send0(cc, "\x71\x01", 2);
			exit(1);
		}
		if (i == 4) {
			/* start motion controller */
			send0(cc, "\x60", 1);
		}
	}
	if (i <= 4)
		send0(cc, "\x60", 1);
}

int
main(int argc, char **argv)
{
	int ret;
	comm_t cc;	/* control connection */
	comm_t cd;	/* data connection */
	char *fname = NULL;

	if (argc == 2)
		fname = argv[1];

/*
	256 microsteps
	400 steps / rev
	16 T / rev
	2 mm / T
	20000000 add /s
	2^^64 per step

	v: steps / s

	400 * 256  steps / 32 mm
*/
	hz = 20000000;
	steps_per_mm = 400 * 256 / 32;
	step = exp2(64) / hz;
	f_spm = (double)steps_per_mm * step;

#if 0
	s_curve_test();
	exit(1);
#endif

	ret = open_terminal("/dev/ttyS1", B115200, 10000, &cd);
	if (ret < 0) {
		printf("failed to open serial port S1: %s\n", strerror(errno));
		exit(1);
	}

	ret = open_terminal("/dev/ttyS2", B115200, 100, &cc);
	if (ret < 0) {
		printf("failed to open serial port S2: %s\n", strerror(errno));
		exit(1);
	}

	resync(&cc);

	print_tmc2130_versions(&cc);

	init_drivers(&cc);

	/* enable, GPOUT 0 to Lo */
	send0(&cc, "\x71\x00", 2);

	/*
	 * stop motion controller. This also clears a potential error condition
	 */
	send0(&cc, "\x61", 1);

	resync(&cd);

	/* reset motion controller */
	send0(&cd, "\x40", 1);

	send0(&cd, "\x6b\x00", 2);	/* endstop polarity, all NC */

	/* initialize routing */
	m_set_routing(&cd, C_E, R_E);
	m_set_routing(&cd, C_XY1, R_XY1);
	m_set_routing(&cd, C_XY2, R_XY2);
	m_set_routing(&cd, C_Z, R_Z1);
	m_set_routing(&cd, C_Z, R_Z2);
	m_set_routing(&cd, C_Z, R_Z3);

#if 0
	home_z(&cc, &cd);
	move_z(&cc, &cd, -100);
#endif
#if 0
	home_xy(&cc, &cd, 0);
	home_xy(&cc, &cd, 1);
#endif

	int n = 5;

	if (fname != NULL) {
		test_from_file(&cc, &cd, fname);
	} else if (n == 1) {
		square_test(&cc, &cd);
	} else if (n == 2) {
		line_test(&cc, &cd);
	} else if (n == 3) {
		star_test(&cc, &cd, 100000);
		star_test(&cc, &cd, 500000);
	} else if (n == 4) {
		int s = 100;
		int i;
		double rounds = 3;
		for (i = 0; i < 9; ++i) {
			/*                    s,   v,      j,     a */
			move_square(&cc, &cd, s, 400, 500000, 10000, rounds);
			s = s * 4.0 / 5;
			rounds = rounds * 5.0 / 4;
		}
	} else if (n == 5) {
		move_square(&cc, &cd, 150, 1300, 500000, 12000, 20);
	}

	/* notify */
	m_notify(&cd, 0x7abb);

	/* stop */
	send0(&cd, "\x41", 1);

	/* wait for notifications */
	while (1) {
		int ret;
#if 1
		print_status(&cc);
#endif

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

	m_preload_jerk(&cd, C_XY1, 0);
	m_preload_jerk(&cd, C_XY2, 0);

	/* disable drivers, GPOUT 0 to Lo */
	send0(&cc, "\x71\x01", 2);

	return 0;
}
