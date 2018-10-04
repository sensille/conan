#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <mpfr.h>

#include "conan.h"

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

void
circle_cb(void *ctx, mpfr_t t, mpfr_t x, mpfr_t y,
	mpfr_t vx, mpfr_t vy, mpfr_t ax, mpfr_t ay)
{
	circle_ctx_t *cp = ctx;

	mpfr_add(cp->t, t, cp->start_t, rnd);
	calccircle(cp->omega, cp->t, cp->r, x, vx, ax, NULL, y, vy, ay, NULL);
}
