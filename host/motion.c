#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <mpfr.h>

#include "conan.h"

/*
 * operations:
 * line:
 * starting speed
 * end speed
 * starting curvature
 * ending curvature
 */

/*
 * calculate parameters for 5th oder polynomial
 */
static void
fit5th(mpfr_t u,
	mpfr_t x0, mpfr_t v0, mpfr_t a0, mpfr_t x1, mpfr_t v1, mpfr_t a1,
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
	mpfr_printf("%.10Re, %.10Re, %.10Re, %.10Re, %.10Re, %.10Re\n",
		x, v, a, j, s, c);
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

int
motion_init(motion_t *mp, kinematics_t kinematics, const char *coeff_file,
	const char *path_file, int simsteps)
{
	memset(mp, 0, sizeof(*mp));

	mpz_init_set_ui(mp->m_x, 0);
	mpfr_init_set_ui(mp->m_vx, 0, rnd);
	mpfr_init_set_ui(mp->m_ax, 0, rnd);
	mpz_init_set_ui(mp->m_y, 0);
	mpfr_init_set_ui(mp->m_vy, 0, rnd);
	mpfr_init_set_ui(mp->m_ay, 0, rnd);

	mp->m_coeff = fopen(coeff_file, "w");
	if (mp->m_coeff == NULL) {
		printf("failed to open coeff output file %s\n", coeff_file);
		return -1;
	}
	mp->m_path = fopen(path_file, "w");
	if (mp->m_path == NULL) {
		printf("failed to open output file %s\n", path_file);
		return -1;
	}

	/*
	 * init controller/printer parameters
	 */
	/* scale microsteps to mm, 256 * 400 microsteps make 32mm */
	mpfr_init_set_ui(mp->m_scale_f, 256 * 400, rnd);
	mpfr_div_ui(mp->m_scale_f, mp->m_scale_f, 32, rnd);

	/* device clock */
	mp->m_hz = 20000000;
	mpz_init_set_ui(mp->m_hz_z, mp->m_hz);
	mpfr_init_set_ui(mp->m_freq, 1, rnd);
	mpfr_div_ui(mp->m_freq, mp->m_freq, mp->m_hz, rnd);

	/* precision in device */
	mp->m_bits = 208;

	/* pi */
	mpfr_init_set_si(mp->m_pi, -1, rnd);
	mpfr_acos(mp->m_pi, mp->m_pi, rnd);
mpfr_printf("pi is %Re\n", mp->m_pi);
	/* scaling factor to fixed point, 2^bits, as mpfr_t */
	mpfr_init_set_ui(mp->m_scale, 1, rnd);
	mpfr_mul_2exp(mp->m_scale, mp->m_scale, mp->m_bits, rnd);

	/* maximum and minimum values representable with these bits, as mpz_t */
	mpz_init(mp->m_minval);
	mpz_init_set_ui(mp->m_maxval, 1);
	mpz_mul_2exp(mp->m_maxval, mp->m_maxval, mp->m_bits);
	mpz_init_set(mp->m_scale_z, mp->m_maxval);
	mpz_neg(mp->m_minval, mp->m_maxval);
	mpz_sub_ui(mp->m_maxval, mp->m_maxval, 1);

	/*
	 * tmp variables
	 */
	mpfr_init(mp->m_tmp);

	mp->m_simsteps = simsteps;
	mp->m_kinematics = kinematics;

	return 0;
}

void
motion_fini(motion_t *mp)
{
	mpfr_clears(mp->m_scale_f, mp->m_scale, mp->m_tmp, mp->m_freq, NULL);
	mpz_clears(mp->m_scale_z, mp->m_maxval, mp->m_minval, mp->m_hz_z, mp->m_x, mp->m_y, NULL);
	mpfr_clears(mp->m_vx, mp->m_ax, mp->m_vy, mp->m_ay, mp->m_pi, NULL);

	fclose(mp->m_coeff);
	fclose(mp->m_path);
}

static void
check_limits(motion_t *mp, mpz_t z, char *name)
{
	if (mpz_cmp(z, mp->m_maxval) > 0) {
		printf("%s exceeded maxval\n", name);
		exit(1);
	}
	if (mpz_cmp(z, mp->m_minval) < 0) {
		printf("%s exceeded minval\n", name);
		exit(1);
	}
}

/*
 * scale one parameter into dev range. The transformation depends of the
 * derivate number
 */
static void
scale_to_dev(motion_t *mp, mpfr_t in, mpz_t out)
{
	/* mm -> microsteps */
	mpfr_mul(mp->m_tmp, in, mp->m_scale_f, rnd);

	/* to fixed point */
	mpfr_mul(mp->m_tmp, mp->m_tmp, mp->m_scale, rnd);

	mpfr_get_z(out, mp->m_tmp, rnd);
}

/* XXX TODO replace div by mul 1/x if needed for performance */
static void
scale_from_dev(motion_t *mp, mpz_t in, mpfr_t out)
{
	mpfr_set_z(out, in, rnd);

	mpfr_div(out, out, mp->m_scale, rnd);
	mpfr_div(out, out, mp->m_scale_f, rnd);
}

/*
 * the device calculates the function
 *
 *   f'(n) = x'
 *         + v'    *n
 *         + a'/2  *n*(n-1)
 *         + j'/6  *n*(n-1)*(n-2)
 *         + s'/24 *n*(n-1)*(n-2)*(n-3)
 *         + c'/120*n*(n-1)*(n-2)*(n-3)*(n-4)
 *
 * while we have our parameters in the form
 *
 *   f(n) = x + v*n + a/2*n^2 + j/6*n^3 + s/24*n^4 + c/120*n^5
 *
 * which is equivalent to
 *   f(n) = x
 *        + (v + a/2 + j/6 + s/24 + c/120)    * n
 *        +          (a + j + 7/12*s +c/4)/2  * n*(n-1)
 *        +            (j + 3/2*s + 5/4*c)/6  * n*(n-1)*(n-2)
 *        +                        (s+2*c)/24 * n*(n-1)*(n-2)*(n-3)
 *        +                             c /120* n*(n-1)*(n-2)*(n-3)*(n-4)
 *
 * so we have
 *   x' = x
 *   v' = v + a/2 + j/6 + s/24 + c/120
 *   a' = a + j + 7/12*s +c/4
 *   j' = j + 3/2*s + 5/4*c
 *   s' = s + 2*c
 *   c' = c
 *
 * and for the other direction
 *   x = x'
 *   v = v' - a'/2 + j'/3 - s'/4 + c'/5
 *   a = a' - j' + 11/12*s' - 5/6*c'
 *   j = j' - 3/2*s' + 7/4*c'
 *   s = s' - 2*c'
 *   c = c'
 */
static void
transform_to_dev(motion_t *mp, mpfr_t v, mpfr_t a, mpfr_t j, mpfr_t s, mpfr_t c,
	mpz_t d_v, mpz_t d_a, mpz_t d_j, mpz_t d_s, mpz_t d_c)
{
	mpfr_t f_v, f_a, f_j, f_s, f_c, tmp;
	mpfr_init(tmp);

	mpfr_init_set(f_v, v, rnd);
	mpfr_init_set(f_a, a, rnd);
	mpfr_init_set(f_j, j, rnd);
	mpfr_init_set(f_s, s, rnd);
	mpfr_init_set(f_c, c, rnd);

	/*
	 * transform to device domain
	 * for each derivate, devide by time scale once (* 1/hz)
	 */
	mpfr_mul(f_v, f_v, mp->m_freq, rnd);
	mpfr_mul(f_a, f_a, mp->m_freq, rnd);
	mpfr_mul(f_a, f_a, mp->m_freq, rnd);
	mpfr_mul(f_j, f_j, mp->m_freq, rnd);
	mpfr_mul(f_j, f_j, mp->m_freq, rnd);
	mpfr_mul(f_j, f_j, mp->m_freq, rnd);
	mpfr_mul(f_s, f_s, mp->m_freq, rnd);
	mpfr_mul(f_s, f_s, mp->m_freq, rnd);
	mpfr_mul(f_s, f_s, mp->m_freq, rnd);
	mpfr_mul(f_s, f_s, mp->m_freq, rnd);
	mpfr_mul(f_c, f_c, mp->m_freq, rnd);
	mpfr_mul(f_c, f_c, mp->m_freq, rnd);
	mpfr_mul(f_c, f_c, mp->m_freq, rnd);
	mpfr_mul(f_c, f_c, mp->m_freq, rnd);
	mpfr_mul(f_c, f_c, mp->m_freq, rnd);

	/* v' = v + a/2 + j/6 + s/24 + c/120 */
	mpfr_div_ui(tmp, f_a, 2, rnd);
	mpfr_add(f_v, f_v, tmp, rnd);
	mpfr_div_ui(tmp, f_j, 6, rnd);
	mpfr_add(f_v, f_v, tmp, rnd);
	mpfr_div_ui(tmp, f_s, 24, rnd);
	mpfr_add(f_v, f_v, tmp, rnd);
	mpfr_div_ui(tmp, f_c, 120, rnd);
	mpfr_add(f_v, f_v, tmp, rnd);

	/* a' = a + j + 7/12*s +c/4 */
	mpfr_add(f_a, f_a, f_j, rnd);
	mpfr_mul_ui(tmp, f_s, 7, rnd);
	mpfr_div_ui(tmp, tmp, 12, rnd);
	mpfr_add(f_a, f_a, tmp, rnd);
	mpfr_div_ui(tmp, f_c, 4, rnd);
	mpfr_add(f_a, f_a, tmp, rnd);

	/* j' = j + 3/2*s + 5/4*c */
	mpfr_mul_ui(tmp, f_s, 3, rnd);
	mpfr_div_ui(tmp, tmp, 2, rnd);
	mpfr_add(f_j, f_j, tmp, rnd);
	mpfr_mul_ui(tmp, f_c, 5, rnd);
	mpfr_div_ui(tmp, tmp, 4, rnd);
	mpfr_add(f_j, f_j, tmp, rnd);

	/* s' = s + 2*c */
	mpfr_mul_ui(tmp, f_c, 2, rnd);
	mpfr_add(f_s, f_s, tmp, rnd);

mpfr_printf("trans %.20Re %.20Re %.20Re %.20Re %.20Re\n", f_v, f_a, f_j, f_s, f_c);
	scale_to_dev(mp, f_v, d_v);
	scale_to_dev(mp, f_a, d_a);
	scale_to_dev(mp, f_j, d_j);
	scale_to_dev(mp, f_s, d_s);
	scale_to_dev(mp, f_c, d_c);

	mpfr_clears(f_v, f_a, f_j, f_s, f_c, tmp, NULL);
}

static void
transform_from_dev(motion_t *mp, mpz_t d_v, mpz_t d_a, mpz_t d_j, mpz_t d_s, mpz_t d_c,
	mpfr_t v, mpfr_t a, mpfr_t j, mpfr_t s, mpfr_t c)
{
	mpfr_t f_v, f_a, f_j, f_s, f_c, tmp;
	mpfr_inits(f_v, f_a, f_j, f_s, f_c, tmp, NULL);

	scale_from_dev(mp, d_v, f_v);
	scale_from_dev(mp, d_a, f_a);
	scale_from_dev(mp, d_j, f_j);
	scale_from_dev(mp, d_s, f_s);
	scale_from_dev(mp, d_c, f_c);

	/* v = v' - a'/2 + j'/3 - s'/4 + c'/5 */
	if (v != NULL) {
		mpfr_div_ui(tmp, f_a, 2, rnd);
		mpfr_sub(f_v, f_v, tmp, rnd);
		mpfr_div_ui(tmp, f_j, 3, rnd);
		mpfr_add(f_v, f_v, tmp, rnd);
		mpfr_div_ui(tmp, f_s, 4, rnd);
		mpfr_sub(f_v, f_v, tmp, rnd);
		mpfr_div_ui(tmp, f_c, 5, rnd);
		mpfr_add(f_v, f_v, tmp, rnd);
		mpfr_mul_z(f_v, f_v, mp->m_hz_z, rnd);
		mpfr_set(v, f_v, rnd);
	}

	/* a = a' - j' + 11/12*s' - 5/6*c' */
	if (a != NULL) {
		mpfr_sub(f_a, f_a, f_j, rnd);
		mpfr_mul_ui(tmp, f_s, 11, rnd);
		mpfr_div_ui(tmp, tmp, 12, rnd);
		mpfr_add(f_a, f_a, tmp, rnd);
		mpfr_mul_ui(tmp, f_c, 5, rnd);
		mpfr_div_ui(tmp, tmp, 6, rnd);
		mpfr_sub(f_a, f_a, tmp, rnd);
		mpfr_mul_z(f_a, f_a, mp->m_hz_z, rnd);
		mpfr_mul_z(f_a, f_a, mp->m_hz_z, rnd);
		mpfr_set(a, f_a, rnd);
	}

	/* j = j' - 3/2*s' + 7/4*c' */
	if (j != NULL) {
		mpfr_mul_ui(tmp, f_s, 3, rnd);
		mpfr_div_ui(tmp, tmp, 2, rnd);
		mpfr_sub(f_j, f_j, tmp, rnd);
		mpfr_mul_ui(tmp, f_c, 7, rnd);
		mpfr_div_ui(tmp, tmp, 4, rnd);
		mpfr_add(f_j, f_j, tmp, rnd);
		mpfr_mul_z(f_j, f_j, mp->m_hz_z, rnd);
		mpfr_mul_z(f_j, f_j, mp->m_hz_z, rnd);
		mpfr_mul_z(f_j, f_j, mp->m_hz_z, rnd);
		mpfr_set(j, f_j, rnd);
	}

	if (s != NULL) {
		mpfr_mul_ui(tmp, f_c, 2, rnd);
		mpfr_sub(f_s, f_s, tmp, rnd);
		mpfr_mul_z(f_s, f_s, mp->m_hz_z, rnd);
		mpfr_mul_z(f_s, f_s, mp->m_hz_z, rnd);
		mpfr_mul_z(f_s, f_s, mp->m_hz_z, rnd);
		mpfr_mul_z(f_s, f_s, mp->m_hz_z, rnd);
		mpfr_set(s, f_s, rnd);
	}

	if (c != NULL) {
		mpfr_mul_z(f_c, f_c, mp->m_hz_z, rnd);
		mpfr_mul_z(f_c, f_c, mp->m_hz_z, rnd);
		mpfr_mul_z(f_c, f_c, mp->m_hz_z, rnd);
		mpfr_mul_z(f_c, f_c, mp->m_hz_z, rnd);
		mpfr_mul_z(f_c, f_c, mp->m_hz_z, rnd);
		mpfr_set(c, f_c, rnd);
	}

	mpfr_clears(f_v, f_a, f_j, f_s, f_c, tmp, NULL);
}

static void
dev_step(motion_t *mp, mpz_t x, mpz_t v, mpz_t a, mpz_t j, mpz_t s, mpz_t c,
	int n)
{
	mpz_t _x, _v, _a, _j, _s, tmp;
	mpz_t _n0, _n1, _n2, _n3, _n4, _ntmp;
	mpz_inits(_x, _v, _a, _j, _s, tmp, NULL);
	mpz_inits(_n0, _n1, _n2, _n3, _n4, _ntmp, NULL);

	mpz_set_ui(_n0, n);
	mpz_sub_ui(_ntmp, _n0, 1);
	mpz_mul(_n1, _n0, _ntmp);
	mpz_sub_ui(_ntmp, _ntmp, 1);
	mpz_mul(_n2, _n1, _ntmp);
	mpz_sub_ui(_ntmp, _ntmp, 1);
	mpz_mul(_n3, _n2, _ntmp);
	mpz_sub_ui(_ntmp, _ntmp, 1);
	mpz_mul(_n4, _n3, _ntmp);
	mpz_div_ui(_n1, _n1, 2);
	mpz_div_ui(_n2, _n2, 6);
	mpz_div_ui(_n3, _n3, 24);
	mpz_div_ui(_n4, _n4, 120);

	/* sn = s0 + n * c0 */
	mpz_init_set(_s, s);
	mpz_mul(tmp, c, _n0);
	mpz_add(_s, _s, tmp);

	/* jn = j0 + n * s0  + 1/2 * n*(n-1) * c0 */
	mpz_init_set(_j, j);
	mpz_mul(tmp, s, _n0);
	mpz_add(_j, _j, tmp);
	mpz_mul(tmp, c, _n1);
	mpz_add(_j, _j, tmp);

	mpz_init_set(_a, a);
	mpz_mul(tmp, j, _n0);
	mpz_add(_a, _a, tmp);
	mpz_mul(tmp, s, _n1);
	mpz_add(_a, _a, tmp);
	mpz_mul(tmp, c, _n2);
	mpz_add(_a, _a, tmp);

	mpz_init_set(_v, v);
	mpz_mul(tmp, a, _n0);
	mpz_add(_v, _v, tmp);
	mpz_mul(tmp, j, _n1);
	mpz_add(_v, _v, tmp);
	mpz_mul(tmp, s, _n2);
	mpz_add(_v, _v, tmp);
	mpz_mul(tmp, c, _n3);
	mpz_add(_v, _v, tmp);

	mpz_init_set(_x, x);
	mpz_mul(tmp, v, _n0);
	mpz_add(_x, _x, tmp);
	mpz_mul(tmp, a, _n1);
	mpz_add(_x, _x, tmp);
	mpz_mul(tmp, j, _n2);
	mpz_add(_x, _x, tmp);
	mpz_mul(tmp, s, _n3);
	mpz_add(_x, _x, tmp);
	mpz_mul(tmp, c, _n4);
	mpz_add(_x, _x, tmp);

#if 0
	int i;
	for (i = 0; i < n; ++i) {
		mpz_add(x, x, v);
		mpz_add(v, v, a);
		mpz_add(a, a, j);
		mpz_add(j, j, s);
		mpz_add(s, s, c);
	}

	if (mpz_cmp(x, _x) != 0 ||
	    mpz_cmp(v, _v) != 0 ||
	    mpz_cmp(a, _a) != 0 ||
	    mpz_cmp(j, _j) != 0 ||
	    mpz_cmp(s, _s) != 0) {
		mpfr_printf("values differ:\n");
		mpfr_printf("  x: %Zd %Zd\n", x, _x);
		mpfr_printf("  v: %Zd %Zd\n", v, _v);
		mpfr_printf("  a: %Zd %Zd\n", a, _a);
		mpfr_printf("  j: %Zd %Zd\n", j, _j);
		mpfr_printf("  s: %Zd %Zd\n", s, _s);
		exit(1);
	}
#else
	mpz_set(x, _x);
	mpz_set(v, _v);
	mpz_set(a, _a);
	mpz_set(j, _j);
	mpz_set(s, _s);
#endif
}

static void
to_kinematics(motion_t *mp, mpfr_t x, mpfr_t y)
{
	switch (mp->m_kinematics) {
	case KIN_CARTESIAN:
		/* nothing to do */
		break;
	case KIN_COREXY:
		mpfr_neg(mp->m_tmp, x, rnd);
		mpfr_add(x, x, y, rnd);
		mpfr_add(y, mp->m_tmp, y, rnd);
		break;
	default:
		abort();
	}
}

static void
from_kinematics(motion_t *mp, mpfr_t x, mpfr_t y)
{
	switch (mp->m_kinematics) {
	case KIN_CARTESIAN:
		/* nothing to do */
		break;
	case KIN_COREXY:
		mpfr_neg(mp->m_tmp, y, rnd);
		mpfr_add(y, x, y, rnd);
		mpfr_add(x, mp->m_tmp, x, rnd);
		mpfr_div_ui(x, x, 2, rnd);
		mpfr_div_ui(y, y, 2, rnd);
		break;
	default:
		abort();
	}
}

void
motion_move_origin(motion_t *mp, mpfr_t e0_x, mpfr_t e0_y)
{
	to_kinematics(mp, e0_x, e0_y);
	scale_to_dev(mp, e0_x, mp->m_x);
	scale_to_dev(mp, e0_y, mp->m_y);
}

/*
 * move from the current position to the given position in given time t
 */
void
motion_move(motion_t *mp, mpfr_t t,
	mpfr_t e1_x, mpfr_t e1_vx, mpfr_t e1_ax,
	mpfr_t e1_y, mpfr_t e1_vy, mpfr_t e1_ay,
	func_cb_t func, void *func_ctx)
{
	mpfr_t e0_x, e0_y, e0_vx, e0_vy, e0_ax, e0_ay;
	mpfr_t x, vx, ax, jx, sx, cx;
	mpfr_t y, vy, ay, jy, sy, cy;
	mpz_t d_vx, d_ax, d_jx, d_sx, d_cx;
	mpz_t d_vy, d_ay, d_jy, d_sy, d_cy;
	mpfr_t u;
	uint64_t cyc;
	mpfr_t r_t, r_x, r_y, r_xy, r_vx, r_vy, r_v, r_ax, r_ay, r_a;

	mpfr_inits(r_t, r_x, r_y, r_xy, r_vx, r_vy, r_v, r_ax, r_ay, r_a, NULL);
	mpfr_inits(e0_x, e0_y, e0_vx, e0_vy, e0_ax, e0_ay, NULL);
	mpfr_inits(x, vx, ax, jx, sx, cx, NULL);
	mpfr_inits(y, vy, ay, jy, sy, cy, NULL);
	mpz_inits(d_vx, d_ax, d_jx, d_sx, d_cx, NULL);
	mpz_inits(d_vy, d_ay, d_jy, d_sy, d_cy, NULL);

	/*
	 * transform the input parameters to target kinematics
	 */
	to_kinematics(mp, e1_x, e1_y);
	to_kinematics(mp, e1_vx, e1_vy);
	to_kinematics(mp, e1_ax, e1_ay);

	/*
	 * get the starting point of the move. It is encoded in dev-coordinates
	 * in mp, so first transform it to host units
	 */
	scale_from_dev(mp, mp->m_x, e0_x);
	mpfr_set(e0_vx, mp->m_vx, rnd);
	mpfr_set(e0_ax, mp->m_ax, rnd);
	scale_from_dev(mp, mp->m_y, e0_y);
	mpfr_set(e0_vy, mp->m_vy, rnd);
	mpfr_set(e0_ay, mp->m_ay, rnd);

	/*
	 * do the fitting, calculate the coefficients
	 */
	fit5th(t, e0_x, e0_vx, e0_ax, e1_x, e1_vx, e1_ax,
		x, vx, ax, jx, sx, cx);
	fit5th(t, e0_y, e0_vy, e0_ay, e1_y, e1_vy, e1_ay,
		y, vy, ay, jy, sy, cy);

mpfr_printf("t %.10Re e0_x %.10Re v %.10Re a %.10Re\n e1_x %.10Re v %.10Re a %.10Re\nfit parameters: %.10Re %.10Re %.10Re %.10Re %.10Re %.10Re\n",
	t, e0_x, e0_vx, e0_ax, e1_x, e1_vx, e1_ax, x, vx, ax, jx, sx, cx);
	/*
	 * transform coefficients to dev
	 */
	transform_to_dev(mp, vx, ax, jx, sx, cx, d_vx, d_ax, d_jx, d_sx, d_cx);
	transform_to_dev(mp, vy, ay, jy, sy, cy, d_vy, d_ay, d_jy, d_sy, d_cy);

	/*
	 * transform time to dev clocks
	 */
	mpfr_init_set(u, t, rnd);
	mpfr_mul_z(u, u, mp->m_hz_z, rnd);
	cyc = mpfr_get_ui(u, rnd);

	/*
	 * send parameters to dev
	 */
	if (mp->m_coeff != NULL)
		mpfr_fprintf(mp->m_coeff, "%d %Zx %Zx %Zx %Zx %Zx %Zx %Zx %Zx %Zx %Zx\n", cyc,
			d_vx, d_ax, d_jx, d_sx, d_cx, d_vy, d_ay, d_jy, d_sy, d_cy);

	/*
	 * if function callback or path output is given, simulate in small
	 * steps and compare with func / write to path output
	 */
	if (func != NULL || mp->m_path != NULL) {
		mpz_t sim_x, sim_vx, sim_ax, sim_jx, sim_sx, sim_cx;
		mpz_t sim_y, sim_vy, sim_ay, sim_jy, sim_sy, sim_cy;
		mpfr_t path_x, path_vx, path_ax, path_jx;
		mpfr_t path_y, path_vy, path_ay, path_jy;
		mpfr_t path_v, path_a, path_j;
		uint64_t step = mp->m_absstep;
		uint64_t steps = cyc;
		uint64_t n;
		uint64_t relstep = 0;

		mpz_init_set(sim_x, mp->m_x);
		mpz_init_set(sim_vx, d_vx);
		mpz_init_set(sim_ax, d_ax);
		mpz_init_set(sim_jx, d_jx);
		mpz_init_set(sim_sx, d_sx);
		mpz_init_set(sim_cx, d_cx);
		mpz_init_set(sim_y, mp->m_y);
		mpz_init_set(sim_vy, d_vy);
		mpz_init_set(sim_ay, d_ay);
		mpz_init_set(sim_jy, d_jy);
		mpz_init_set(sim_sy, d_sy);
		mpz_init_set(sim_cy, d_cy);

		mpfr_inits(path_x, path_vx, path_ax, path_jx,
			path_y, path_vy, path_ay, path_jy,
			path_v, path_a, path_j, NULL);

		while (steps) {
			scale_from_dev(mp, sim_x, path_x);
			transform_from_dev(mp, sim_vx, sim_ax, sim_jx, sim_sx, sim_cx,
				path_vx, path_ax, path_jx, NULL, NULL);
			scale_from_dev(mp, sim_y, path_y);
			transform_from_dev(mp, sim_vy, sim_ay, sim_jy, sim_sy, sim_cy,
				path_vy, path_ay, path_jy, NULL, NULL);
mpfr_printf("x/y from sim %.10Re %.10Re\n", path_x, path_y);
			from_kinematics(mp, path_x, path_y);
			from_kinematics(mp, path_vx, path_vy);
			from_kinematics(mp, path_ax, path_ay);
			norm(path_v, path_vx, path_vy);
			norm(path_a, path_ax, path_ay);
			norm(path_j, path_jx, path_jy);

			if (func != NULL) {
				/* get comparison values from func */
				/* step to time */
				mpfr_set_ui(r_t, relstep, rnd);
				mpfr_mul(r_t, r_t, mp->m_freq, rnd);
				func(func_ctx, r_t, r_x, r_y, r_vx, r_vy, r_ax, r_ay);
				/*
				 * calculate deviation from real path
				 */
				mpfr_sub(r_x, r_x, path_x, rnd);
				mpfr_sub(r_y, r_y, path_y, rnd);
				norm(r_xy, r_x, r_y);
				norm(r_v, r_vx, r_vy);
				mpfr_sub(r_v, r_v, path_v, rnd);
				norm(r_a, r_ax, r_ay);
				mpfr_sub(r_a, r_a, path_a, rnd);
			}

			if (mp->m_path != NULL) {
				mpfr_fprintf(mp->m_path,
					"%.5f %.3Rf %.3Rf %.3Rf %.3Rf %.3Rf | %.3Re %.3Re %.3Re\n",
					(double)step / mp->m_hz,
					path_x, path_y, path_v, path_a, path_j,
					r_xy, r_v, r_a);
			}
			n = steps > mp->m_simsteps ? mp->m_simsteps : steps;

			dev_step(mp, sim_x, sim_vx, sim_ax, sim_jx,
				sim_sx, sim_cx, n);
			dev_step(mp, sim_y, sim_vy, sim_ay, sim_jy,
				sim_sy, sim_cy, n);

			check_limits(mp, sim_sx, "sim_sx");
			check_limits(mp, sim_jx, "sim_jx");
			check_limits(mp, sim_ax, "sim_ax");
			check_limits(mp, sim_vx, "sim_vx");
			check_limits(mp, sim_sy, "sim_sy");
			check_limits(mp, sim_jy, "sim_jy");
			check_limits(mp, sim_ay, "sim_ay");
			check_limits(mp, sim_vy, "sim_vy");

			step += n;
			steps -= n;
			relstep += n;
		}

		mpz_clears(sim_x, sim_vx, sim_ax, sim_jx, sim_sx, sim_cx,
			sim_y, sim_vy, sim_ay, sim_jy, sim_sy, sim_cy, NULL);
		mpfr_clears(path_x, path_vx, path_ax, path_jx,
			path_y, path_vy, path_ay, path_jy,
			path_v, path_a, path_j, NULL);
	}

	dev_step(mp, mp->m_x, d_vx, d_ax, d_jx, d_sx, d_cx, cyc);
	dev_step(mp, mp->m_y, d_vy, d_ay, d_jy, d_sy, d_cy, cyc);

{
mpfr_t px, py; mpfr_inits(px, py, NULL);
scale_from_dev(mp, mp->m_x, px);
scale_from_dev(mp, mp->m_y, py);
mpfr_printf("x/y after sim %.20Re %.20Re\n", px, py);
mpfr_clears(px, py, NULL);
}
	transform_from_dev(mp, d_vx, d_ax, d_jx, d_sx, d_cx, vx, ax, NULL, NULL, NULL);
	transform_from_dev(mp, d_vy, d_ay, d_jy, d_sy, d_cy, vy, ay, NULL, NULL, NULL);

	/* save endpoint as next start */
	mpfr_set(mp->m_vx, vx, rnd);
	mpfr_set(mp->m_ax, ax, rnd);
	mpfr_set(mp->m_vy, vy, rnd);
	mpfr_set(mp->m_ay, ay, rnd);

	mp->m_absstep += cyc;

	mpfr_clears(e0_x, e0_y, e0_vx, e0_vy, e0_ax, e0_ay, NULL);
	mpfr_clears(x, vx, ax, jx, sx, cx, NULL);
	mpfr_clears(y, vy, ay, jy, sy, cy, NULL);
	mpz_clears(d_jx, d_sx, d_cx, d_jy, d_sy, d_cy, NULL);
	mpfr_clear(u);
	mpfr_clears(r_t, r_x, r_y, r_xy, r_vx, r_vy, r_v, r_ax, r_ay, r_a, NULL);
}

#if 0
int
main(int argc, char **argv)
{
	int ret;
	motion_t m;

	mpfr_set_default_prec(300);

	ret = motion_init(&m, KIN_COREXY, "coeff.txt", "path.csv", 10000);
	if (ret < 0)
		exit(1);


	mpfr_t omega, t, r, delta;
	mpfr_t e0_x, e0_vx, e0_ax, e0_y, e0_vy, e0_ay;
	mpfr_inits(omega, t, r, delta, NULL);
	mpfr_inits(e0_x, e0_vx, e0_ax, e0_y, e0_vy, e0_ay, NULL);


	exit(0);
#if 0
#if 1
	mpfr_set_d(e0_x,  -5.0, rnd);
	mpfr_set_d(e0_y,  75.0, rnd);
	motion_move_origin(&m, e0_x, e0_y);
#endif
#if 0
	/*
	 * a = 100mm/s^2
	 * j = 1000mm/s^3
	 */
mpfr_t _a, _j;
mpfr_inits(_a, _j, NULL);
mpfr_set_ui(_a, 100, rnd);
mpfr_set_ui(_j, 1000, rnd);
	mpfr_div(delta, _a, _j, rnd);
	mpfr_set_ui(e0_x, 1, rnd);
	mpfr_div_ui(e0_x, e0_x, 6, rnd);
	mpfr_mul(e0_x, e0_x, _j, rnd);
	mpfr_mul(e0_x, e0_x, delta, rnd);
	mpfr_mul(e0_x, e0_x, delta, rnd);
	mpfr_mul(e0_x, e0_x, delta, rnd);
	mpfr_set_ui(e0_y, 0, rnd);
	mpfr_set_ui(e0_vx, 1, rnd);
	mpfr_div_ui(e0_vx, e0_vx, 2, rnd);
	mpfr_mul(e0_vx, e0_vx, _j, rnd);
	mpfr_mul(e0_vx, e0_vx, delta, rnd);
	mpfr_mul(e0_vx, e0_vx, delta, rnd);
	mpfr_set_ui(e0_vy, 0, rnd);
	mpfr_set(e0_ax, _a, rnd);
	mpfr_set_ui(e0_ay, 0, rnd);
	motion_move(&m, delta, e0_x, e0_vx, e0_ax, e0_y, e0_vy, e0_ay, NULL, NULL);
#endif

#if 1
	mpfr_set_d(omega, 1, rnd);
	mpfr_set_d(t, 0, rnd);
	mpfr_set_d(r, 80, rnd);

	int div = 40;
	int i;

	mpfr_add(delta, m.m_pi, m.m_pi, rnd); /* 2PI */
	mpfr_div_ui(delta, delta, div, rnd);

	mpfr_mul_z(delta, delta, m.m_hz_z, rnd);
	mpfr_floor(delta, delta);
	mpfr_div_z(delta, delta, m.m_hz_z, rnd);

	circle_ctx_t cc;
	mpfr_init_set(cc.omega, omega, rnd);
	mpfr_init_set(cc.r, r, rnd);
	mpfr_init(cc.start_t);
	mpfr_init(cc.t);

	for (i = 0; i <= div; ++i) {
		calccircle(omega, t, r, e0_x, e0_vx, e0_ax, NULL,
			e0_y, e0_vy, e0_ay, NULL);

		mpfr_printf("circle at %.10Re %.10Re %.10Re %.10Re %.10Re %.10Re %.10Re\n",
			t, e0_x, e0_vx, e0_ax, e0_y, e0_vy, e0_ay);

		motion_move(&m, delta, e0_x, e0_vx, e0_ax, e0_y, e0_vy, e0_ay, i == 0 ? NULL : circle_cb, &cc);
		mpfr_set(cc.start_t, t, rnd);

		mpfr_add(t, t, delta, rnd);
	}

	mpfr_clears(cc.omega, cc.r, cc.t, cc.start_t, NULL);

	mpfr_set_d(e0_x,   5.0, rnd);
	mpfr_set_d(e0_y,  75.0, rnd);
	mpfr_set_d(e0_vx,    0, rnd);
	mpfr_set_d(e0_vy,    0, rnd);
	mpfr_set_d(e0_ax,    0, rnd);
	mpfr_set_d(e0_ay,    0, rnd);

	motion_move(&m, delta, e0_x, e0_vx, e0_ax, e0_y, e0_vy, e0_ay, NULL, NULL);

#if 0
	enter_circle(&params, 30);
	approx_circle(&params, 80, 30, 40);
#endif
#endif

#endif
	motion_fini(&m);

	return 0;
}
#endif
