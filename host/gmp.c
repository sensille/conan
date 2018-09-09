#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <mpfr.h>

#if 0
XXX next up: model move-to
#endif

mpfr_rnd_t rnd = MPFR_RNDN;
int hz = 20000000;
mpfr_t scale;
mpz_t scale_z, maxval, minval, hz_z;

/*
 * calculate parameters for 5th oder polynomial
 */
void
fit5th(mpfr_t u, mpfr_t x0, mpfr_t v0, mpfr_t a0, mpfr_t x1, mpfr_t v1, mpfr_t a1,
	mpfr_t x, mpfr_t v, mpfr_t a, mpfr_t j, mpfr_t s, mpfr_t c)
{
	mpfr_t tmp, m, n, o;
	mpfr_inits(tmp, m, n, o, NULL);

	/* m = x1 - 1.0 / 2 * a0 * u * u - v0 * u - x0 */
	mpfr_set(m, x1, rnd);
	mpfr_set_ui(tmp, 1, rnd);
	mpfr_div_ui(tmp, tmp, 2, rnd);
	mpfr_mul(tmp, tmp, a0, rnd);
	mpfr_mul(tmp, tmp, u, rnd);
	mpfr_mul(tmp, tmp, u, rnd);
	mpfr_sub(m, m, tmp, rnd);
	mpfr_set(tmp, v0, rnd);
	mpfr_mul(tmp, tmp, u, rnd);
	mpfr_sub(m, m, tmp, rnd);
	mpfr_sub(m, m, x0, rnd);

	/* n = v1 - s0 * u - v0; */
	mpfr_set(n, v1, rnd);
	mpfr_set(tmp, a0, rnd);
	mpfr_mul(tmp, tmp, u, rnd);
	mpfr_sub(n, n, tmp, rnd);
	mpfr_sub(n, n, v0, rnd);

	/* o = a1 - a0 */
	mpfr_sub(o, a1, a0, rnd);

	mpfr_set(x, x0, rnd);
	mpfr_set(v, v0, rnd);
	mpfr_set(a, a0, rnd);

	/* j = (3 * o * u * u - 24 * n * u + 60 * m) / (u * u * u) */
	mpfr_mul_ui(j, o, 3, rnd);
	mpfr_mul(j, j, u, rnd);
	mpfr_mul(j, j, u, rnd);
	mpfr_mul_ui(tmp, n, 24, rnd);
	mpfr_mul(tmp, tmp, u, rnd);
	mpfr_sub(j, j, tmp, rnd);
	mpfr_mul_ui(tmp, m, 60, rnd);
	mpfr_add(j, j, tmp, rnd);
	mpfr_div(j, j, u, rnd);
	mpfr_div(j, j, u, rnd);
	mpfr_div(j, j, u, rnd);

	/* s = -(24 * o * u * u - 168 * n * u + 360 * m) /
		(u * u * u * u) */
	mpfr_mul_ui(s, o, 24, rnd);
	mpfr_mul(s, s, u, rnd);
	mpfr_mul(s, s, u, rnd);
	mpfr_mul_ui(tmp, n, 168, rnd);
	mpfr_mul(tmp, tmp, u, rnd);
	mpfr_sub(s, s, tmp, rnd);
	mpfr_mul_ui(tmp, m, 360, rnd);
	mpfr_add(s, s, tmp, rnd);
	mpfr_div(s, s, u, rnd);
	mpfr_div(s, s, u, rnd);
	mpfr_div(s, s, u, rnd);
	mpfr_div(s, s, u, rnd);
	mpfr_neg(s, s, rnd);

	/* c = (60 * o * u * u - 360 * n * u + 720 * m) /
		(u * u * u * u * u) */
	mpfr_mul_ui(c, o, 60, rnd);
	mpfr_mul(c, c, u, rnd);
	mpfr_mul(c, c, u, rnd);
	mpfr_mul_ui(tmp, n, 360, rnd);
	mpfr_mul(tmp, tmp, u, rnd);
	mpfr_sub(c, c, tmp, rnd);
	mpfr_mul_ui(tmp, m, 720, rnd);
	mpfr_add(c, c, tmp, rnd);
	mpfr_div(c, c, u, rnd);
	mpfr_div(c, c, u, rnd);
	mpfr_div(c, c, u, rnd);
	mpfr_div(c, c, u, rnd);
	mpfr_div(c, c, u, rnd);

	mpfr_printf("m %.10Re n %.10Re o %.10Re\n", m, n, o);
	mpfr_printf("%.10Re, %.10Re, %.10Re, %.10Re, %.10Re, %.10Re\n", x, v, a, j, s, c);
}

void
eval5th(mpfr_t t, mpfr_t x, mpfr_t v, mpfr_t a, mpfr_t j, mpfr_t s, mpfr_t c,
	mpfr_t rx, mpfr_t rv, mpfr_t ra, mpfr_t rj)
{
	mpfr_t tmp;

	mpfr_init(tmp);

	/*
	 * rx = x + v * t + 1.0/ 2 * a * t * t +
	 *	1.0 / 6 * j * t * t * t +
	 *	1.0 / 24 * s * t * t * t * t +
	 *	1.0 / 120 * c * t * t * t * t * t;
	 */
	mpfr_set(rx, x, rnd);
	mpfr_mul(tmp, v, t, rnd);
	mpfr_add(rx, rx, tmp, rnd);
	mpfr_set_ui(tmp, 1, rnd);
	mpfr_div_ui(tmp, tmp, 2, rnd);
	mpfr_mul(tmp, tmp, a, rnd);
	mpfr_mul(tmp, tmp, t, rnd);
	mpfr_mul(tmp, tmp, t, rnd);
	mpfr_add(rx, rx, tmp, rnd);
	mpfr_set_ui(tmp, 1, rnd);
	mpfr_div_ui(tmp, tmp, 6, rnd);
	mpfr_mul(tmp, tmp, j, rnd);
	mpfr_mul(tmp, tmp, t, rnd);
	mpfr_mul(tmp, tmp, t, rnd);
	mpfr_mul(tmp, tmp, t, rnd);
	mpfr_add(rx, rx, tmp, rnd);
	mpfr_set_ui(tmp, 1, rnd);
	mpfr_div_ui(tmp, tmp, 24, rnd);
	mpfr_mul(tmp, tmp, s, rnd);
	mpfr_mul(tmp, tmp, t, rnd);
	mpfr_mul(tmp, tmp, t, rnd);
	mpfr_mul(tmp, tmp, t, rnd);
	mpfr_mul(tmp, tmp, t, rnd);
	mpfr_add(rx, rx, tmp, rnd);
	mpfr_set_ui(tmp, 1, rnd);
	mpfr_div_ui(tmp, tmp, 120, rnd);
	mpfr_mul(tmp, tmp, c, rnd);
	mpfr_mul(tmp, tmp, t, rnd);
	mpfr_mul(tmp, tmp, t, rnd);
	mpfr_mul(tmp, tmp, t, rnd);
	mpfr_mul(tmp, tmp, t, rnd);
	mpfr_mul(tmp, tmp, t, rnd);
	mpfr_add(rx, rx, tmp, rnd);

	/*
	 *	rv = v + a * t + 1.0/ 2 * j * t * t +
	 * 	1.0 / 6 * s * t * t * t +
	 * 	1.0 / 24 * c * t * t * t * t;
	 */
	mpfr_set(rv, v, rnd);
	mpfr_mul(tmp, a, t, rnd);
	mpfr_add(rv, rv, tmp, rnd);
	mpfr_set_ui(tmp, 1, rnd);
	mpfr_div_ui(tmp, tmp, 2, rnd);
	mpfr_mul(tmp, tmp, j, rnd);
	mpfr_mul(tmp, tmp, t, rnd);
	mpfr_mul(tmp, tmp, t, rnd);
	mpfr_add(rv, rv, tmp, rnd);
	mpfr_set_ui(tmp, 1, rnd);
	mpfr_div_ui(tmp, tmp, 6, rnd);
	mpfr_mul(tmp, tmp, s, rnd);
	mpfr_mul(tmp, tmp, t, rnd);
	mpfr_mul(tmp, tmp, t, rnd);
	mpfr_mul(tmp, tmp, t, rnd);
	mpfr_add(rv, rv, tmp, rnd);
	mpfr_set_ui(tmp, 1, rnd);
	mpfr_div_ui(tmp, tmp, 24, rnd);
	mpfr_mul(tmp, tmp, c, rnd);
	mpfr_mul(tmp, tmp, t, rnd);
	mpfr_mul(tmp, tmp, t, rnd);
	mpfr_mul(tmp, tmp, t, rnd);
	mpfr_mul(tmp, tmp, t, rnd);
	mpfr_add(rv, rv, tmp, rnd);

	/*
	 * ra = a + j * t + 1.0/ 2 * s * t * t +
	 * 	1.0 / 6 * c * t * t * t;
	 */
	mpfr_set(ra, a, rnd);
	mpfr_mul(tmp, j, t, rnd);
	mpfr_add(ra, ra, tmp, rnd);
	mpfr_set_ui(tmp, 1, rnd);
	mpfr_div_ui(tmp, tmp, 2, rnd);
	mpfr_mul(tmp, tmp, s, rnd);
	mpfr_mul(tmp, tmp, t, rnd);
	mpfr_mul(tmp, tmp, t, rnd);
	mpfr_add(ra, ra, tmp, rnd);
	mpfr_set_ui(tmp, 1, rnd);
	mpfr_div_ui(tmp, tmp, 6, rnd);
	mpfr_mul(tmp, tmp, c, rnd);
	mpfr_mul(tmp, tmp, t, rnd);
	mpfr_mul(tmp, tmp, t, rnd);
	mpfr_mul(tmp, tmp, t, rnd);
	mpfr_add(ra, ra, tmp, rnd);

	/* rj = j + s * t + 1.0/ 2 * c * t * t; */
	mpfr_set(rj, j, rnd);
	mpfr_mul(tmp, s, t, rnd);
	mpfr_add(rj, rj, tmp, rnd);
	mpfr_set_ui(tmp, 1, rnd);
	mpfr_div_ui(tmp, tmp, 2, rnd);
	mpfr_mul(tmp, tmp, c, rnd);
	mpfr_mul(tmp, tmp, t, rnd);
	mpfr_mul(tmp, tmp, t, rnd);
	mpfr_add(rj, rj, tmp, rnd);

	mpfr_printf("[%Rf] rx %.10Re rv %.10Re ra %.10Re rj %.10Re\n", t, rx, rv, ra, rj);
}

void
calccircle(mpfr_t omega, mpfr_t t, mpfr_t r,
	mpfr_t x, mpfr_t vx, mpfr_t ax, mpfr_t jx,
	mpfr_t y, mpfr_t vy, mpfr_t ay, mpfr_t jy)
{
	mpfr_t sin_o, cos_o, omega2, omega3;

	mpfr_inits(sin_o, cos_o, omega2, omega3, NULL);

	mpfr_mul(omega2, omega, omega, rnd);
	mpfr_mul(omega3, omega2, omega, rnd);
	mpfr_mul(sin_o, omega, t, rnd);
	mpfr_sin(sin_o, sin_o, rnd);
	mpfr_mul(cos_o, omega, t, rnd);
	mpfr_cos(cos_o, cos_o, rnd);

	/*
	 *  x =  r *          sin(omega * t);
	 *  y =  r *          cos(omega * t);
	 * vx =  r * omega  * cos(omega * t);
	 * vy = -r * omega  * sin(omega * t);
	 * ax = -r * omega2 * sin(omega * t);
	 * ay = -r * omega2 * cos(omega * t);
	 * jx = -r * omega3 * cos(omega * t);
	 * jy =  r * omega3 * sin(omega * t);
	 */
	mpfr_mul(x, r, sin_o, rnd);
	mpfr_mul(y, r, cos_o, rnd);

	mpfr_mul(vx, omega, cos_o, rnd);
	mpfr_mul(vx, vx, r, rnd);
	mpfr_mul(vy, omega, sin_o, rnd);
	mpfr_mul(vy, vy, r, rnd);
	mpfr_neg(vy, vy, rnd);

	mpfr_mul(ax, omega2, sin_o, rnd);
	mpfr_mul(ax, ax, r, rnd);
	mpfr_neg(ax, ax, rnd);
	mpfr_mul(ay, omega2, cos_o, rnd);
	mpfr_mul(ay, ay, r, rnd);
	mpfr_neg(ay, ay, rnd);

	if (jx != NULL) {
		mpfr_mul(jx, omega3, cos_o, rnd);
		mpfr_mul(jx, jx, r, rnd);
		mpfr_neg(jx, jx, rnd);
	}
	if (jy != NULL) {
		mpfr_mul(jy, omega3, sin_o, rnd);
		mpfr_mul(jy, jy, r, rnd);
	}
}

/*
 * 256 * 400 microsteps make 32mm
 */
static int scale_inited = 0;
static mpfr_t scale_f;	/* steps per mm */
void
init_scale(void)
{
	if (scale_inited)
		return;

	mpfr_init(scale_f);
	mpfr_set_ui(scale_f, 256 * 400, rnd);
	mpfr_div_ui(scale_f, scale_f, 32, rnd);
	scale_inited = 1;
}

void
scale_one_to_dev(mpfr_t v, mpz_t dv, char *name)
{
	mpfr_t tmp;
	mpfr_init(tmp);

	init_scale();

	mpfr_mul(tmp, v, scale_f, rnd);
	mpfr_mul(tmp, tmp, scale, rnd);
	mpfr_get_z(dv, tmp, rnd);

#if 0
	mpfr_printf("%s %Zd\n", name, dv);
#endif
}

/*
 * scaling from mm to device steps, assuming a unit in t = 1 as 1/hz
 */
void
scale_to_dev(mpfr_t x, mpfr_t v, mpfr_t a, mpfr_t j, mpfr_t s, mpfr_t c,
	mpz_t dx, mpz_t dv, mpz_t da, mpz_t dj, mpz_t ds, mpz_t dc)
{
	scale_one_to_dev(c, dc, "c");
	scale_one_to_dev(s, ds, "s");
	scale_one_to_dev(j, dj, "j");
	scale_one_to_dev(a, da, "a");
	scale_one_to_dev(v, dv, "v");
	scale_one_to_dev(x, dx, "x");
}

void
scale_one_from_dev(mpz_t dx, mpfr_t x)
{
	init_scale();

	mpfr_set_z(x, dx, rnd);
	mpfr_div(x, x, scale, rnd);
	mpfr_div(x, x, scale_f, rnd);
}

void
check_limits(mpz_t z, char *name)
{
	if (mpz_cmp(z, maxval) > 0) {
		printf("%s exceeded maxval\n", name);
		exit(1);
	}
	if (mpz_cmp(z, minval) < 0) {
		printf("%s exceeded minval\n", name);
		exit(1);
	}
}

void
step_dev(int steps,
	mpz_t f_x, mpz_t f_vx, mpz_t f_ax, mpz_t f_jx, mpz_t f_sx, mpz_t f_cx,
	mpz_t f_y, mpz_t f_vy, mpz_t f_ay, mpz_t f_jy, mpz_t f_sy, mpz_t f_cy)
{
	int i;
	mpz_t oldx, oldy, newx, newy, d_x, d_y;
	mpz_inits(oldx, oldy, newx, newy, d_x, d_y, NULL);

	mpz_div(oldx, f_x, scale_z);
	mpz_div(oldy, f_y, scale_z);

	for (i = 0; i < steps; ++i) {

		mpz_add(f_x, f_x, f_vx);
		mpz_add(f_vx, f_vx, f_ax);
		mpz_add(f_ax, f_ax, f_jx);
		mpz_add(f_jx, f_jx, f_sx);
		mpz_add(f_sx, f_sx, f_cx);
		mpz_add(f_y, f_y, f_vy);
		mpz_add(f_vy, f_vy, f_ay);
		mpz_add(f_ay, f_ay, f_jy);
		mpz_add(f_jy, f_jy, f_sy);
		mpz_add(f_sy, f_sy, f_cy);

		check_limits(f_sx, "f_sx");
		check_limits(f_jx, "f_jx");
		check_limits(f_ax, "f_ax");
		check_limits(f_vx, "f_vx");
		check_limits(f_sy, "f_sy");
		check_limits(f_jy, "f_jy");
		check_limits(f_ay, "f_ay");
		check_limits(f_vy, "f_vy");

		mpz_div(newx, f_x, scale_z);
		mpz_div(newy, f_y, scale_z);

		if (mpz_cmp(oldx, newx)) {
			printf("step on x at %d\n", i);
		}
		if (mpz_cmp(oldy, newy)) {
			printf("step on y at %d\n", i);
		}

		mpz_set(oldx, newx);
		mpz_set(oldy, newy);
	}
}

static void
norm(mpfr_t n, mpfr_t x, mpfr_t y)
{
	mpfr_t tmp;
	mpfr_init(tmp);

	mpfr_mul(tmp, x, x, rnd);
	mpfr_mul(n, y, y, rnd);
	mpfr_add(n, n, tmp, rnd);
	mpfr_sqrt(n, n, rnd);
}

typedef struct _param {
	int	p_inited;
	unsigned long p_step;
	FILE	*p_fp;
	FILE	*p_fp2;
	mpz_t	p_x;
	mpz_t	p_vx;
	mpz_t	p_ax;
	mpz_t	p_y;
	mpz_t	p_vy;
	mpz_t	p_ay;
} param_t;
/*
 * store parameters for the next move. strictly speaking x/y/vx/vy/ax/ay
 * are not necessary as they are tracked inside the device, but they
 * can be used to verify that we're still on track
 */
void
push_parameters(param_t *pp, unsigned int cyc,
	mpz_t d_x, mpz_t d_vx, mpz_t d_ax, mpz_t d_jx, mpz_t d_sx, mpz_t d_cx,
	mpz_t d_y, mpz_t d_vy, mpz_t d_ay, mpz_t d_jy, mpz_t d_sy, mpz_t d_cy)
{
	if (pp->p_inited == 0) {
		mpz_inits(pp->p_x, pp->p_vx, pp->p_ax, NULL);
		mpz_inits(pp->p_y, pp->p_vy, pp->p_ay, NULL);
		mpz_set(pp->p_x, d_x);
		mpz_set(pp->p_vx, d_vx);
		mpz_set(pp->p_ax, d_ax);
		mpz_set(pp->p_y, d_y);
		mpz_set(pp->p_vy, d_vy);
		mpz_set(pp->p_ay, d_ay);
		pp->p_inited = 1;
		pp->p_fp = fopen("path.csv", "w");
		if (pp->p_fp == NULL) {
			printf("failed to open output file\n");
			exit(1);
		}
		pp->p_fp2 = fopen("path.params", "w");
		if (pp->p_fp2 == NULL) {
			printf("failed to open output file\n");
			exit(1);
		}
	}

	mpfr_fprintf(pp->p_fp2, "%d %Zx %Zx %Zx %Zx %Zx %Zx\n", cyc, d_jx, d_sx, d_cx, d_jy, d_sy, d_cy);

	/* simulate and plot */
	mpz_t p_jx, p_sx, p_cx, p_jy, p_sy, p_cy;
	mpz_inits(p_jx, p_sx, p_cx, p_jy, p_sy, p_cy, NULL);
	mpz_set(p_jx, d_jx);
	mpz_set(p_sx, d_sx);
	mpz_set(p_cx, d_cx);
	mpz_set(p_jy, d_jy);
	mpz_set(p_sy, d_sy);
	mpz_set(p_cy, d_cy);

	mpfr_t r_x, r_y, r_vx, r_vy, r_ax, r_ay, r_jx, r_jy, r_v, r_a, r_j;
	mpfr_inits(r_x, r_y, r_vx, r_vy, r_ax, r_ay, r_jx, r_jy, r_v, r_a, r_j, NULL);

	while (cyc > 0) {
		int steps = cyc > 1000 ? 1000 : cyc;

		step_dev(steps, pp->p_x, pp->p_vx, pp->p_ax, p_jx, p_sx, p_cx,
			pp->p_y, pp->p_vy, pp->p_ay, p_jy, p_sy, p_cy);

		pp->p_step += steps;

		scale_one_from_dev(pp->p_x, r_x);
		scale_one_from_dev(pp->p_vx, r_vx);
		scale_one_from_dev(pp->p_ax, r_ax);
		scale_one_from_dev(p_jx, r_jx);
		scale_one_from_dev(pp->p_y, r_y);
		scale_one_from_dev(pp->p_vy, r_vy);
		scale_one_from_dev(pp->p_ay, r_ay);
		scale_one_from_dev(p_jy, r_jy);
		norm(r_v, r_vx, r_vy);
		norm(r_a, r_ax, r_ay);
		norm(r_j, r_jx, r_jy);
		/* translate from dev time to seconds */
		mpfr_mul_z(r_v, r_v, hz_z, rnd);
		mpfr_mul_z(r_a, r_a, hz_z, rnd);
		mpfr_mul_z(r_a, r_a, hz_z, rnd);
		mpfr_mul_z(r_j, r_j, hz_z, rnd);
		mpfr_mul_z(r_j, r_j, hz_z, rnd);
		mpfr_mul_z(r_j, r_j, hz_z, rnd);

		mpfr_fprintf(pp->p_fp, "%.5f %.3Rf %.3Rf %.3Rf %.3Rf %.3Rf\n",
			(double)pp->p_step / hz, r_x, r_y, r_v, r_a, r_j);

		cyc -= steps;
	}
}

void
approx_circle(param_t *pp, double radius, double v, int divisions)
{
	mpfr_t s, t, cyc, omega, r, dt, t0, u;
	unsigned int cycles;
	int i;

	/*
	 * calculate clock cycles per division
	 *
	 * s = 2 * M_PI * r;
	 * t = s / v;
	 * cycles = t * hz;
	 */

	mpfr_inits(s, t, cyc, omega, r, dt, t0, u, NULL);
	mpfr_set_d(s, 2 * M_PI * radius, rnd);
	mpfr_div_ui(s, s, divisions, rnd);
	mpfr_div_d(t, s, v, rnd);
	mpfr_mul_z(cyc, t, hz_z, rnd);
	cycles = mpfr_get_ui(cyc, rnd);
	mpfr_set_d(r, radius, rnd);

	/* to avoid overflows, split into more segments */
	if (cycles > 16 * 1048576) {
		printf("XXX TODO need more divisions\n");
		exit(1);
	}
	mpfr_set(u, cyc, rnd);
	mpfr_set_d(omega, v, rnd);
	mpfr_div(omega, omega, r, rnd);
	mpfr_div_z(omega, omega, hz_z, rnd);

	mpfr_printf("approx_circle: cycles %d omega %.10Re t %.10Re\n", cycles, omega, t);

	mpfr_t e0_x, e0_y, e0_vx, e0_vy, e0_ax, e0_ay;
	mpfr_t e1_x, e1_y, e1_vx, e1_vy, e1_ax, e1_ay;
	mpfr_t x, vx, ax, jx, sx, cx;
	mpfr_t y, vy, ay, jy, sy, cy;
	mpfr_inits(e0_x, e0_y, e0_vx, e0_vy, e0_ax, e0_ay, NULL);
	mpfr_inits(e1_x, e1_y, e1_vx, e1_vy, e1_ax, e1_ay, NULL);
	mpfr_inits(x, vx, ax, jx, sx, cx, NULL);
	mpfr_inits(y, vy, ay, jy, sy, cy, NULL);

	mpz_t d_x, d_vx, d_ax, d_jx, d_sx, d_cx, d_y, d_vy, d_ay, d_jy, d_sy, d_cy;
	mpz_inits(d_x, d_vx, d_ax, d_jx, d_sx, d_cx, d_y, d_vy, d_ay, d_jy, d_sy, d_cy, NULL);

	mpfr_set_ui(t0, 0, rnd);
	calccircle(omega, t0, r, e0_x, e0_vx, e0_ax, NULL, e0_y, e0_vy, e0_ay, NULL);
	mpfr_set(u, cyc, rnd);
	for (i = 0; i <= divisions; ++i) {
		mpfr_mul_ui(dt, u, i, rnd);
		mpfr_add(dt, dt, u, rnd); /* dt = t * (i + 1) */

		calccircle(omega, dt, r, e1_x, e1_vx, e1_ax, NULL, e1_y, e1_vy, e1_ay, NULL);

		mpfr_printf("e0_xy %.10Re %.10Re e0_vxy %.10Re %.10Re e0_axy %.10Re %.10Re\n", e0_x, e0_y, e0_vx, e0_vy, e0_ax, e0_ay);
		mpfr_printf("e1_xy %.10Re %.10Re e1_vxy %.10Re %.10Re e1_axy %.10Re %.10Re\n", e1_x, e1_y, e1_vx, e1_vy, e1_ax, e1_ay);

		fit5th(u, e0_x, e0_vx, e0_ax, e1_x, e1_vx, e1_ax, x, vx, ax, jx, sx, cx);
		fit5th(u, e0_y, e0_vy, e0_ay, e1_y, e1_vy, e1_ay, y, vy, ay, jy, sy, cy);
		scale_to_dev(x, vx, ax, jx, sx, cx, d_x, d_vx, d_ax, d_jx, d_sx, d_cx);
		scale_to_dev(y, vy, ay, jy, sy, cy, d_y, d_vy, d_ay, d_jy, d_sy, d_cy);

		mpfr_printf("jx %Zd sx %Zd cx %Zd jy %Zd sy %Zd cy %Zd\n", d_jx, d_sx, d_cx, d_jy, d_sy, d_cy);

		push_parameters(pp, cycles, d_x, d_vx, d_ax, d_jx, d_sx, d_cx, d_y, d_vy, d_ay, d_jy, d_sy, d_cy);

		mpfr_set(e0_x, e1_x, rnd);
		mpfr_set(e0_vx, e1_vx, rnd);
		mpfr_set(e0_ax, e1_ax, rnd);
		mpfr_set(e0_y, e1_y, rnd);
		mpfr_set(e0_vy, e1_vy, rnd);
		mpfr_set(e0_ay, e1_ay, rnd);
	}
}

void
enter_circle(param_t *pp, double v)
{
	int cycles = hz;
	mpfr_t e0_x, e0_y, e0_vx, e0_vy, e0_ax, e0_ay;
	mpfr_t e1_x, e1_y, e1_vx, e1_vy, e1_ax, e1_ay;
	mpfr_t x, vx, ax, jx, sx, cx;
	mpfr_t y, vy, ay, jy, sy, cy;
	mpfr_inits(e0_x, e0_y, e0_vx, e0_vy, e0_ax, e0_ay, NULL);
	mpfr_inits(e1_x, e1_y, e1_vx, e1_vy, e1_ax, e1_ay, NULL);
	mpfr_inits(x, vx, ax, jx, sx, cx, NULL);
	mpfr_inits(y, vy, ay, jy, sy, cy, NULL);
	mpz_t d_x, d_vx, d_ax, d_jx, d_sx, d_cx, d_y, d_vy, d_ay, d_jy, d_sy, d_cy;
	mpz_inits(d_x, d_vx, d_ax, d_jx, d_sx, d_cx, d_y, d_vy, d_ay, d_jy, d_sy, d_cy, NULL);
	mpfr_t u;
	mpfr_init(u);
	mpfr_set_ui(u, cycles, rnd);

#if 0
approx_circle: cycles 8377580 omega 1.8750000000e-08 t 4.1887902048e-01
e0_xy 0.0000000000e+00 8.0000000000e+01 e0_vxy 1.5000000000e-06 -0.0000000000e+00 e0_axy -0.0000000000e+00 -2.8125000000e-14
e1_xy 1.2514757203e+01 7.9015067248e+01 e1_vxy 1.4815325109e-06 -2.3465169756e-07 e1_axy -4.3997193293e-15 -2.7778734579e-14
#endif

	mpfr_set_d(e0_x, -10, rnd);
	mpfr_set_d(e0_y, 70, rnd);
	mpfr_set_d(e0_vx, 0, rnd);
	mpfr_set_d(e0_vy, 0, rnd);
	mpfr_set_d(e0_ax, 0, rnd);
	mpfr_set_d(e0_ay, 0, rnd);

	mpfr_set_d(e1_x, 0, rnd);
	mpfr_set_d(e1_y, 80, rnd);
	mpfr_set_d(e1_vx, v / hz, rnd);
	mpfr_set_d(e1_vy, 0, rnd);
	mpfr_set_d(e1_ax, 0, rnd);
	mpfr_set_d(e1_ay, -2.8125e-14, rnd);

	fit5th(u, e0_x, e0_vx, e0_ax, e1_x, e1_vx, e1_ax, x, vx, ax, jx, sx, cx);
	fit5th(u, e0_y, e0_vy, e0_ay, e1_y, e1_vy, e1_ay, y, vy, ay, jy, sy, cy);
	scale_to_dev(x, vx, ax, jx, sx, cx, d_x, d_vx, d_ax, d_jx, d_sx, d_cx);
	scale_to_dev(y, vy, ay, jy, sy, cy, d_y, d_vy, d_ay, d_jy, d_sy, d_cy);

	mpfr_printf("jx %Zd sx %Zd cx %Zd jy %Zd sy %Zd cy %Zd\n", d_jx, d_sx, d_cx, d_jy, d_sy, d_cy);

	push_parameters(pp, cycles, d_x, d_vx, d_ax, d_jx, d_sx, d_cx, d_y, d_vy, d_ay, d_jy, d_sy, d_cy);
}

int
main(int argc, char **argv)
{
	param_t params = { 0 };

	mpfr_set_default_prec(240);

	mpfr_init(scale);
	mpz_inits(maxval, minval, hz_z, NULL);

	/* init global constants */
	mpz_set_ui(hz_z, hz);
	int bitshift = 208;
	mpfr_set_ui(scale, 1, rnd);
	mpfr_mul_2exp(scale, scale, bitshift, rnd);
	mpz_set_ui(maxval, 1);
	mpz_mul_2exp(maxval, maxval, bitshift);
	mpz_set(scale_z, maxval);
	mpz_neg(minval, maxval);
	mpz_sub_ui(maxval, maxval, 1);

	enter_circle(&params, 30);
	approx_circle(&params, 80, 30, 40);
	fclose(params.p_fp);

exit(0);
	mpfr_t r, div, omega, d, u;
	mpfr_inits(r, div, omega, d, u, NULL);

	int divisions = 20;
	mpfr_set_ui(r, 10, rnd);
	mpfr_set_ui(div, divisions, rnd);

	double stepsize = 10000;
	int steps = 20000000; /* one second, the maximum we want to support */
	int total_steps = steps * divisions;
	mpfr_set_d(omega, (2.0 * M_PI) / total_steps, rnd);	/* one round per second */
	mpfr_ui_div(u, total_steps, div, rnd);

	mpfr_printf("r %.10Re, div %.10Re, u %.10Re\n", r, div, u);
	for (mpfr_set_ui(d, 0, rnd); mpfr_cmp_ui(d, total_steps); mpfr_add(d, d, u, rnd)) {
		mpfr_t e0_x, e0_y, e0_vx, e0_vy, e0_ax, e0_ay;
		mpfr_t e1_x, e1_y, e1_vx, e1_vy, e1_ax, e1_ay;
		mpfr_inits(e0_x, e0_y, e0_vx, e0_vy, e0_ax, e0_ay, NULL);
		mpfr_inits(e1_x, e1_y, e1_vx, e1_vy, e1_ax, e1_ay, NULL);

		mpfr_printf("d is %.10Re\n", d);
		/*
		 * calculate start and end values of sin/cos
		 */
		mpfr_t du;
		mpfr_init(du);
		mpfr_add(du, d, u, rnd);

		calccircle(omega,  d, r, e0_x, e0_vx, e0_ax, NULL, e0_y, e0_vy, e0_ay, NULL);
		calccircle(omega, du, r, e1_x, e1_vx, e1_ax, NULL, e1_y, e1_vy, e1_ay, NULL);

mpfr_printf("e0_xy %.10Re %.10Re e0_vxy %.10Re %.10Re e0_axy %.10Re %.10Re\n", e0_x, e0_y, e0_vx, e0_vy, e0_ax, e0_ay);
mpfr_printf("e1_xy %.10Re %.10Re e1_vxy %.10Re %.10Re e1_axy %.10Re %.10Re\n", e1_x, e1_y, e1_vx, e1_vy, e1_ax, e1_ay);

		/*
		 * calculate parameters for 5th oder polynomial
		 */
		mpfr_t x, y, vx, vy, ax, ay, jx, jy, sx, sy, cx, cy;
		mpfr_inits(x, y, vx, vy, ax, ay, jx, jy, sx, sy, cx, cy, NULL);

		fit5th(u, e0_x, e0_vx, e0_ax, e1_x, e1_vx, e1_ax, x, vx, ax, jx, sx, cx);
		fit5th(u, e0_y, e0_vy, e0_ay, e1_y, e1_vy, e1_ay, y, vy, ay, jy, sy, cy);

		mpfr_t t, ug, tmp;
		mpfr_inits(t, ug, tmp, NULL);
		mpfr_add_d(ug, u, 0.00000001, rnd);

		mpz_t f_cx, f_sx, f_jx, f_ax, f_vx, f_x;
		mpz_t f_cy, f_sy, f_jy, f_ay, f_vy, f_y;
		mpz_inits(f_cx, f_sx, f_jx, f_ax, f_vx, f_x, NULL);
		mpz_inits(f_cy, f_sy, f_jy, f_ay, f_vy, f_y, NULL);

		scale_to_dev(x, vx, ax, jx, sx, cx, f_x, f_vx, f_ax, f_jx, f_sx, f_cx);
		scale_to_dev(y, vy, ay, jy, sy, cy, f_y, f_vy, f_ay, f_jy, f_sy, f_cy);

		for (mpfr_set_ui(t, 0, rnd); mpfr_cmp(t, ug) < 0; mpfr_add_d(t, t, stepsize, rnd)) {
			mpfr_t dt;
			mpfr_t e_x, e_vx, e_ax, e_jx;
			mpfr_t e_y, e_vy, e_ay, e_jy;

			mpfr_inits(dt, NULL);
			mpfr_inits(e_x, e_vx, e_ax, e_jx, NULL);
			mpfr_inits(e_y, e_vy, e_ay, e_jy, NULL);

			mpfr_add(dt, d, t, rnd);

			/*
			 * expected value from circle
			 */
			calccircle(omega, dt, r, e_x, e_vx, e_ax, e_jx, e_y, e_vy, e_ay, e_jy);

			mpfr_t r_x, r_vx, r_ax, r_jx;
			mpfr_t r_y, r_vy, r_ay, r_jy;

			mpfr_inits(r_x, r_vx, r_ax, r_jx, NULL);
			mpfr_inits(r_y, r_vy, r_ay, r_jy, NULL);

			/*
			 * approximated value by polynomial
			 */
			eval5th(t, x, vx, ax, jx, sx, cx, r_x, r_vx, r_ax, r_jx);
			eval5th(t, y, vy, ay, jy, sy, cy, r_y, r_vy, r_ay, r_jy);

			mpfr_t d_x, d_y, d_vx, d_vy, d_ax, d_ay, d_jx, d_jy;
			mpfr_inits(d_x, d_y, d_vx, d_vy, d_ax, d_ay, d_jx, d_jy, NULL);
			mpfr_sub(d_x, e_x, r_x, rnd);
			mpfr_sub(d_y, e_y, r_y, rnd);
			mpfr_sub(d_vx, e_vx, r_vx, rnd);
			mpfr_sub(d_vy, e_vy, r_vy, rnd);
			mpfr_sub(d_ax, e_ax, r_ax, rnd);
			mpfr_sub(d_ay, e_ay, r_ay, rnd);
			mpfr_sub(d_jx, e_jx, r_jx, rnd);
			mpfr_sub(d_jy, e_jy, r_jy, rnd);

			mpfr_printf("[%Rf] e_xy %.10Re %.10Re r_xy %.10Re %.10Re diff %.10Re %.10Re\n", t, e_x, e_y,
				r_x, r_y, d_x, d_y);
			mpfr_printf("      e_vxy %.10Re %.10Re r_vxy %.10Re %.10Re diff %.10Re %.10Re\n", e_vx, e_vy,
				r_vx, r_vy, d_vx, d_vy);
			mpfr_printf("      e_axy %.10Re %.10Re r_axy %.10Re %.10Re diff %.10Re %.10Re\n", e_ax, e_ay,
				r_ax, r_ay, d_ax, d_ay);
			mpfr_printf("      e_jxy %.10Re %.10Re r_jxy %.10Re %.10Re diff %.10Re %.10Re\n", e_jx, e_jy,
				r_jx, r_jy, d_jx, d_jy);

			mpfr_t dist_x, dist_v, dist_a, dist_j;
			mpfr_inits(dist_x, dist_v, dist_a, dist_j, NULL);

			mpfr_sub(dist_x, e_x, r_x, rnd);
			mpfr_mul(dist_x, dist_x, dist_x, rnd);
			mpfr_sub(tmp, e_y, r_y, rnd);
			mpfr_mul(tmp, tmp, tmp, rnd);
			mpfr_add(dist_x, dist_x, tmp, rnd);
			mpfr_sqrt(dist_x, dist_x, rnd);

			mpfr_sub(dist_v, e_vx, r_vx, rnd);
			mpfr_mul(dist_v, dist_v, dist_v, rnd);
			mpfr_sub(tmp, e_vy, r_vy, rnd);
			mpfr_mul(tmp, tmp, tmp, rnd);
			mpfr_add(dist_v, dist_v, tmp, rnd);
			mpfr_sqrt(dist_v, dist_v, rnd);

			mpfr_sub(dist_a, e_ax, r_ax, rnd);
			mpfr_mul(dist_a, dist_a, dist_a, rnd);
			mpfr_sub(tmp, e_ay, r_ay, rnd);
			mpfr_mul(tmp, tmp, tmp, rnd);
			mpfr_add(dist_a, dist_a, tmp, rnd);
			mpfr_sqrt(dist_a, dist_a, rnd);

			mpfr_sub(dist_j, e_jx, r_jx, rnd);
			mpfr_mul(dist_j, dist_j, dist_j, rnd);
			mpfr_sub(tmp, e_jy, r_jy, rnd);
			mpfr_mul(tmp, tmp, tmp, rnd);
			mpfr_add(dist_j, dist_j, tmp, rnd);
			mpfr_sqrt(dist_j, dist_j, rnd);

			mpfr_printf("      point distance %.10Re\n", dist_x);
			mpfr_printf("      speed distance %.10Re\n", dist_v);
			mpfr_printf("      accel distance %.10Re\n", dist_a);
			mpfr_printf("      jerk  distance %.10Re\n", dist_j);

			mpfr_t d_tmp;
			mpfr_init(d_tmp);

			mpfr_set_z(tmp, f_cx, rnd);
			mpfr_div(tmp, tmp, scale, rnd);
			mpfr_printf("f_cx %Zd == %.10Re\n", f_cx, tmp);
			mpfr_set_z(tmp, f_sx, rnd);
			mpfr_div(tmp, tmp, scale, rnd);
			mpfr_printf("f_sx %Zd == %.10Re\n", f_sx, tmp);
			mpfr_set_z(tmp, f_jx, rnd);
			mpfr_div(tmp, tmp, scale, rnd);
			mpfr_sub(d_tmp, r_jx, tmp, rnd);
			mpfr_printf("f_jx %Zd == %.10Re diff %.10Re\n", f_jx, tmp, d_tmp);
			mpfr_set_z(tmp, f_ax, rnd);
			mpfr_div(tmp, tmp, scale, rnd);
			mpfr_sub(d_tmp, r_ax, tmp, rnd);
			mpfr_printf("f_ax %Zd == %.10Re diff %.10Re\n", f_ax, tmp, d_tmp);
			mpfr_set_z(tmp, f_vx, rnd);
			mpfr_div(tmp, tmp, scale, rnd);
			mpfr_sub(d_tmp, r_vx, tmp, rnd);
			mpfr_printf("f_vx %Zd == %.10Re diff %.10Re\n", f_vx, tmp, d_tmp);
			mpfr_set_z(tmp, f_x, rnd);
			mpfr_div(tmp, tmp, scale, rnd);
			mpfr_sub(d_tmp, r_x, tmp, rnd);
			mpfr_printf("f_x %Zd == %.10Re diff %.10Re\n", f_x, tmp, d_tmp);

			mpfr_set_z(tmp, f_cy, rnd);
			mpfr_div(tmp, tmp, scale, rnd);
			mpfr_printf("f_cy %Zd == %.10Re\n", f_cy, tmp);
			mpfr_set_z(tmp, f_sy, rnd);
			mpfr_div(tmp, tmp, scale, rnd);
			mpfr_printf("f_sy %Zd == %.10Re\n", f_sy, tmp);
			mpfr_set_z(tmp, f_jy, rnd);
			mpfr_div(tmp, tmp, scale, rnd);
			mpfr_printf("f_jy %Zd == %.10Re\n", f_jy, tmp);
			mpfr_set_z(tmp, f_ay, rnd);
			mpfr_div(tmp, tmp, scale, rnd);
			mpfr_printf("f_ay %Zd == %.10Re\n", f_ay, tmp);
			mpfr_set_z(tmp, f_vy, rnd);
			mpfr_div(tmp, tmp, scale, rnd);
			mpfr_printf("f_vy %Zd == %.10Re\n", f_vy, tmp);
			mpfr_set_z(tmp, f_y, rnd);
			mpfr_div(tmp, tmp, scale, rnd);
			mpfr_printf("f_y %Zd == %.10Re\n", f_y, tmp);

			step_dev(stepsize, f_x, f_vx, f_ax, f_jx, f_sx, f_cx,
				f_y, f_vy, f_ay, f_jy, f_sy, f_cy);
		}
		printf("==================================\n");
	}

	return 0;
}
