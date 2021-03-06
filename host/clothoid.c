#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdint.h>
#include <string.h>
#include <assert.h>
#include <mpfr.h>

#include "conan.h"

#undef FRES_DEBUG

/*
 * calculate Fresnel-integrals. See
 * Smith, D. M. (2011).
 * Algorithm 911. ACM Transactions on Mathematical Software, 37(4), 1-16.
 * doi:10.1145/1916461.1916470
 */
static void
calcfresnel(mpfr_t C, mpfr_t S, mpfr_t x)
{
	mpfr_t tmp, t1, t2, v, half, pi;
	mpfr_inits(tmp, t1, t2, v, half, NULL);
	int i;
	int bits = mpfr_get_prec(C);

	if (mpfr_cmp_ui(x, 0) == 0) {
		mpfr_set_ui(C, 0, rnd);
		mpfr_set_ui(S, 0, rnd);
		return;
	}

	if (bits < mpfr_get_prec(S))
		bits = mpfr_get_prec(S);

	mpfr_init2(pi, bits*2);
	mpfr_const_pi(pi, rnd);

#ifdef FRES_DEBUG
	int rounds;
#endif
	int add = 0;

	mpfr_set_ui(half, 1, rnd);
	mpfr_div_ui(half, half, 2, rnd);

	/*
	 * decide which algorithm to use
	 */
	/* v = (pi*x^2 - 1)/2 */
	mpfr_mul(v, x, x, rnd);
	mpfr_mul(v, v, pi, rnd);
	mpfr_sub_ui(v, v, 1, rnd);
	mpfr_div_ui(v, v, 2, rnd);

	/* (2v + 1.5)ln(2v+2) - (2v+1) -vln(2) - (v+0.5)ln(v+1)+v-v*ln(2v+1) */
	mpfr_add(t1, v, v, rnd);
	mpfr_add_ui(t1, t1, 2, rnd);
	mpfr_log(t1, t1, rnd);
	mpfr_add(t2, v, v, rnd);
	mpfr_add(t2, t2, half, rnd);
	mpfr_add_ui(t2, t2, 1, rnd);
	mpfr_mul(tmp, t1, t2, rnd);
	mpfr_add(t1, v, v, rnd);
	mpfr_add_ui(t1, t1, 1, rnd);
	mpfr_sub(tmp, tmp, t1, rnd);
	mpfr_set_ui(t1, 2, rnd);
	mpfr_log(t1, t1, rnd);
	mpfr_mul(t1, t1, v, rnd);
	mpfr_sub(tmp, tmp, t1, rnd);
	mpfr_add(t1, v, half, rnd);
	mpfr_add_ui(t2, v, 1, rnd);
	mpfr_log(t2, t2, rnd);
	mpfr_mul(t1, t1, t2, rnd);
	mpfr_sub(tmp, tmp, t1, rnd);
	mpfr_add(tmp, tmp, v, rnd);
	mpfr_add(t1, v, v, rnd);
	mpfr_add_ui(t1, t1, 1, rnd);
	mpfr_log(t1, t1, rnd);
	mpfr_mul(t1, t1, v, rnd);
	mpfr_sub(tmp, tmp, t1, rnd);

	mpfr_set_ui(t1, 2, rnd);
	mpfr_log(t1, t1, rnd);
	mpfr_mul_ui(t1, t1, bits + 1, rnd);
	mpfr_neg(t1, t1, rnd);
#ifdef FRES_DEBUG
	mpfr_printf("x %.10Rf decide: %.10Rf threshold %.10Rf\n", x, tmp, t1);
#endif

	if (mpfr_cmp(tmp, t1) < 0) {
		mpfr_t f, g, px2, px22, f_prev, g_prev, tmp;
		mpfr_inits(f, g, px2, px22, tmp, f_prev, g_prev, NULL);

#ifdef FRES_DEBUG
		printf("using asymptotic series\n");
#endif
		/* asymptotic series */
		mpfr_set_ui(f, 1, rnd);
		mpfr_set_ui(g, 1, rnd);
		mpfr_set_ui(tmp, 1, rnd);
		mpfr_mul(px2, pi, x, rnd);
		mpfr_mul(px2, px2, x, rnd);
		mpfr_mul(px22, px2, px2, rnd);
		mpfr_div(px22, tmp, px22, rnd); /* 1/(pi*x^2)^2 */

		for (i = 3; i; i += 4) {
			mpfr_set(f_prev, f, rnd);
			mpfr_set(g_prev, g, rnd);
			mpfr_mul(tmp, tmp, px22, rnd);
			mpfr_mul_ui(tmp, tmp, i, rnd);
			if (add)
				mpfr_add(f, f, tmp, rnd);
			else
				mpfr_sub(f, f, tmp, rnd);
			mpfr_mul_ui(tmp, tmp, i + 2, rnd);
			if (add)
				mpfr_sub(g, g, tmp, rnd);
			else
				mpfr_add(g, g, tmp, rnd);
			add = !add;
#ifdef FRES_DEBUG
			mpfr_printf("[%d] tmp is %.10Re px22 %.10Re f "
				"%.10Re g %.10Re\n", i, tmp, px22, f, g);
#endif
			if (mpfr_cmp(f, f_prev) == 0 &&
			    mpfr_cmp(g, g_prev) == 0)
				break;
		}
#ifdef FRES_DEBUG
		rounds = (i - 3) / 4;
#endif
		mpfr_div(f, f, pi, rnd);
		mpfr_div(f, f, x, rnd);
		mpfr_div(g, g, px2, rnd);
		mpfr_div(g, g, pi, rnd);
		mpfr_div(g, g, x, rnd);
		mpfr_div_ui(t1, px2, 2, rnd);
		mpfr_sin(t2, t1, rnd);	/* t2 = sin(...) */
		mpfr_cos(t1, t1, rnd);	/* t1 = cos(...) */
		mpfr_set(C, half, rnd);
		mpfr_mul(tmp, f, t2, rnd);
		mpfr_add(C, C, tmp, rnd);
		mpfr_mul(tmp, g, t1, rnd);
		mpfr_sub(C, C, tmp, rnd);
		mpfr_set(S, half, rnd);
		mpfr_mul(tmp, f, t1, rnd);
		mpfr_sub(S, S, tmp, rnd);
		mpfr_mul(tmp, g, t2, rnd);
		mpfr_sub(S, S, tmp, rnd);
		mpfr_clears(f, g, px2, px22, tmp, f_prev, g_prev, NULL);
	} else {
		mpfr_t tmp, _C, _S, px2_2, fac, acc, C_prev, S_prev, prec;
		/* XXX find out why BITS are not enough, change the algorithm */
		mpfr_inits2(bits * 2, tmp, _C, _S, px2_2, fac, acc, C_prev, S_prev, prec, NULL);

		/* convergent series */
		mpfr_set_ui(prec, 1, rnd);
		mpfr_div_2exp(prec, prec, bits, rnd);

		/* pi*x^2 / 2 */
		mpfr_mul(px2_2, x, x, rnd);
		mpfr_mul(px2_2, px2_2, pi, rnd);
		mpfr_div_ui(px2_2, px2_2, 2, rnd);

		mpfr_set_ui(fac, 1, rnd);
		mpfr_set_ui(_C, 0, rnd);
		mpfr_set_ui(_S, 0, rnd);
		mpfr_set_ui(acc, 1, rnd);
#ifdef FRES_DEBUG
		printf("using convergent series\n");
#endif
		add = 1;
		for (i = 0; i >= 0; i += 1) {
			mpfr_set(C_prev, _C, rnd);
			mpfr_set(S_prev, _S, rnd);
			mpfr_div(tmp, acc, fac, rnd);
			mpfr_div_ui(tmp, tmp, 4 * i + 1, rnd);
			if (add)
				mpfr_add(_C, _C, tmp, rnd);
			else
				mpfr_sub(_C, _C, tmp, rnd);
			mpfr_mul_ui(fac, fac, 2 * i + 1, rnd);
			mpfr_mul(acc, acc, px2_2, rnd);
			mpfr_div(tmp, acc, fac, rnd);
			mpfr_div_ui(tmp, tmp, 4 * i + 3, rnd);
			if (add)
				mpfr_add(_S, _S, tmp, rnd);
			else
				mpfr_sub(_S, _S, tmp, rnd);
			add = !add;
			mpfr_mul_ui(fac, fac, 2 * i + 2, rnd);
			mpfr_mul(acc, acc, px2_2, rnd);
#ifdef FRES_DEBUG
			mpfr_printf("[%d] tmp is %.10Re fac %.10Re acc %.10Re "
				"C %.10Re S %.10Re\n", i, tmp, fac, acc, _C, _S);
#endif
			mpfr_sub(tmp, _C, C_prev, rnd);
			mpfr_abs(tmp, tmp, rnd);
			if (mpfr_cmp(tmp, prec) > 0)
				continue;
			mpfr_sub(tmp, _S, S_prev, rnd);
			mpfr_abs(tmp, tmp, rnd);
			if (mpfr_cmp(tmp, prec) > 0)
				continue;
			break;
		}
#ifdef FRES_DEBUG
		rounds = i;
#endif
		mpfr_mul(C, _C, x, rnd);
		mpfr_mul(S, _S, x, rnd);
		mpfr_clears(px2_2, fac, acc, C_prev, S_prev, _C, _S, tmp, NULL);
	}
#ifdef FRES_DEBUG
	mpfr_printf("x %.10Re => C %.10Re S %.10Re in %d rounds\n",
		x, C, S, rounds);
#endif
	mpfr_clears(tmp, t1, t2, v, half, NULL);
}

static void
_calcclothoid(mpfr_t x, mpfr_t y, mpfr_t a, mpfr_t b, mpfr_t t)
{
	mpfr_t tmp;
	mpfr_inits(tmp, NULL);

	mpfr_mul(tmp, b, t, rnd);
	calcfresnel(x, y, tmp);
	mpfr_mul(x, x, a, rnd);
	mpfr_mul(y, y, a, rnd);
}

void
cloth_reset(path_elem_t *pe)
{
	clothoid_t *c = pe->clothoid;
	mpfr_t tmp;

	mpfr_inits(tmp, NULL);

	/*
	 * simple approach for now: no segment longer than 10mm, but have
	 * at least 4 segments. equidistant.
	 */
	mpfr_div_ui(tmp, c->L, 10, rnd);
	c->nsegments = mpfr_get_ui(tmp, rnd);
	if (c->nsegments < 2)
		c->nsegments = 2;
	c->ptr = 1;
}

int
cloth_next(path_elem_t *pe, mpfr_t t)
{
	clothoid_t *c = pe->clothoid;

	if (c->ptr > c->nsegments * 2)
		return 0;

	mpfr_div_ui(t, c->t, c->nsegments, rnd);
	mpfr_mul_ui(t, t, c->ptr, rnd);
	++c->ptr;

	return 1;
}

void
cloth_calc(void *ctx, mpfr_t t, mpfr_t x, mpfr_t vx,
	mpfr_t ax, mpfr_t y, mpfr_t vy, mpfr_t ay)
{
	path_elem_t *pe = ctx;
	clothoid_t *c = pe->clothoid;

	calcclothoid(c, t, x, y, vx, vy, ax, ay);
}

/*
 * - all clothoids we need start from angle 0 or end there
 * - the largest angle we need is 90 degree
 * - the calcfresnel calculates integral({sin,cos}(pi/2*t^2)dt) from 0 to x
 * - with C'(x) = 0 => x = pi/2 as first occurence, so we only need the range
 *   for 0 to 1
 * - the curve still needs to be stretched and rotated
 * - as a first step, we only fit line-to-line
 * - we fit with 2 clothoids
 * - we calculate the angle between the lines and find the x where the angle
 *   is half
 *   the angle between the lines
 * - at this point, we reach maximum perpendicular accelaration a = v^2/r
 * - scale the clothoid so that the target acceleration is reached
 * - find starting and end point
 * - done fitting
 * - the lines are given by their intersection point and 2 two vectors
 * - output is starting point, the point where the two clothoids meet, A of the
 *   the clothoids and the transformation matrix
 */
path_elem_t *
cloth_line(mpfr_t inter_x, mpfr_t inter_y, mpfr_t _vx1, mpfr_t _vy1,
	mpfr_t _vx2, mpfr_t _vy2, mpfr_t v, mpfr_t acc)
{
	clothoid_t *c;
	path_elem_t *pe;
	mpfr_t phi, tmp, pi, x, dx, dy, vx1, vy1, vx2, vy2, rx, ry;
	mpfr_t mx, my, ml, vl;
	int left = 0;

	pe = calloc(sizeof(*pe) + sizeof(*c), 1);
	assert(pe != NULL);
	c = (clothoid_t *)(pe + 1);
	pe->type = PT_CLOTHOID;
	pe->clothoid = c;
	pe->reset = cloth_reset;
	pe->next = cloth_next;
	pe->calc = cloth_calc;

	mpfr_inits(c->a, c->b, c->t, c->vl, c->L, NULL);
	mpfr_inits(c->px1, c->py1, c->px2, c->py2, NULL);
	mpfr_inits(c->r1_11, c->r1_12, c->r1_21, c->r1_22, NULL);
	mpfr_inits(c->r2_11, c->r2_12, c->r2_21, c->r2_22, NULL);

	mpfr_inits(phi, tmp, pi, x, dx, dy, vx1, vy1, vx2, vy2, rx, ry, NULL);
	mpfr_inits(mx, my, ml, vl, NULL);

	mpfr_printf("inter (%.4Rf/%.4Rf) v1 (%.4Rf/%.4Rf) v2 (%.4Rf/%.4Rf) v %.5Rf acc %0.5Rf\n",
		inter_x, inter_y, _vx1, _vy1, _vx2, _vy2, v, acc);

	mpfr_const_pi(pi, rnd);

	/*
	 * normate v1 and v2
	 */
	mpfr_hypot(tmp, _vx1, _vy1, rnd);
	mpfr_div(vx1, _vx1, tmp, rnd);
	mpfr_div(vy1, _vy1, tmp, rnd);
	mpfr_hypot(tmp, _vx2, _vy2, rnd);
	mpfr_div(vx2, _vx2, tmp, rnd);
	mpfr_div(vy2, _vy2, tmp, rnd);

	mpfr_printf("after normation: %.10Rf %.10Rf %.10Rf %.10Rf\n",
		vx1, vy1, vx2, vy2);
	/*
	 * calculate the angle between the lines
	 */
	/* dot product for 0-180 deg */
	mpfr_mul(phi, vx1, vx2, rnd);
	mpfr_fma(phi, vy1, vy2, phi, rnd);
	mpfr_printf("phi is %.20Rf\n", phi);
	mpfr_acos(phi, phi, rnd);
	/* cross product to determine left/right turn */
	mpfr_mul(tmp, vx1, vy2, rnd);
	mpfr_fms(tmp, vx2, vy1, tmp, rnd);
	if (mpfr_cmp_ui(tmp, 0) < 0) {
		mpfr_printf("left turn detected\n");
		left = 1;
		/* swap v1 and v2 */
		mpfr_set(tmp, vx1, rnd);
		mpfr_neg(vx1, vx2, rnd);
		mpfr_neg(vx2, tmp, rnd);
		mpfr_set(tmp, vy1, rnd);
		mpfr_neg(vy1, vy2, rnd);
		mpfr_neg(vy2, tmp, rnd);
	}

	mpfr_div(tmp, phi, pi, rnd);
	mpfr_mul_ui(tmp, tmp, 180, rnd);
	mpfr_printf("full angle is %.5Rf\n", tmp);

	/*
	 * 90 deg - half the angle is the endpoint of the clothoid
	 */
	mpfr_div_ui(phi, phi, 2, rnd);

	mpfr_div(tmp, phi, pi, rnd);
	mpfr_mul_ui(tmp, tmp, 180, rnd);
	mpfr_printf("angle is %.5Rf\n", tmp);

	/*
	 * calculate parameter from phi, v and acc
	 */
	mpfr_set_ui(c->t, 2, rnd);
	mpfr_div(c->t, c->t, acc, rnd);
	mpfr_mul(c->t, c->t, phi, rnd);
	mpfr_mul(c->t, c->t, v, rnd);

	mpfr_mul(c->b, phi, pi, rnd);
	mpfr_mul_ui(c->b, c->b, 2, rnd);
	mpfr_sqrt(c->b, c->b, rnd);
	mpfr_mul(c->b, c->b, v, rnd);
	mpfr_div(c->b, acc, c->b, rnd);

	mpfr_div(c->a, v, c->b, rnd);

	/*
	 * distance traveled by the clothoid
	 */
	_calcclothoid(dx, dy, c->a, c->b, c->t);

	mpfr_printf("a %.5Rf b %.5Rf t %.5Rf dx %.5Rf dy %.5Rf\n",
		c->a, c->b, c->t, dx, dy);

	/*
	 * middle line
	 */
	mpfr_sub(mx, vx2, vx1, rnd);
	mpfr_sub(my, vy2, vy1, rnd);

	mpfr_printf("mx %.5Rf my %.5Rf\n", mx, my);
	/*
	 * rotate middle line so that v1 is vertical
	 * |vy -vx|
	 * |vx  vy|
	 */
	mpfr_mul(rx, my, vx1, rnd);
	mpfr_fms(rx, mx, vy1, rx, rnd);
	mpfr_mul(ry, my, vy1, rnd);
	mpfr_fma(ry, mx, vx1, ry, rnd);

	mpfr_printf("rx %.5Rf ry %.5Rf\n", rx, ry);

	/*
	 * match clothoid to it
	 */
	mpfr_div(ml, dy, rx, rnd);
	mpfr_mul(vl, ml, ry, rnd);
	mpfr_sub(vl, dx, vl, rnd);

	mpfr_printf("ml %.5Rf vl %.5Rf\n", ml, vl);

	/*
	 * calculate start and end points
	 */
	mpfr_mul(c->px1, vl, vx1, rnd);
	mpfr_sub(c->px1, inter_x, c->px1, rnd);
	mpfr_mul(c->py1, vl, vy1, rnd);
	mpfr_sub(c->py1, inter_y, c->py1, rnd);
	mpfr_mul(c->px2, vl, vx2, rnd);
	mpfr_add(c->px2, inter_x, c->px2, rnd);
	mpfr_mul(c->py2, vl, vy2, rnd);
	mpfr_add(c->py2, inter_y, c->py2, rnd);

	/*
	 * store rotation matrix for first clothoid
	 * mirror on x, rotate to v1
	 * |vx1 -vy1| \/ |1  0|  _ | vx1  vy1|
	 * |vy1  vx1| /\ |0 -1|  - | vy1 -vx1|
	 */
	mpfr_set(c->r1_11, vx1, rnd);
	mpfr_set(c->r1_12, vy1, rnd);
	mpfr_set(c->r1_21, vy1, rnd);
	mpfr_neg(c->r1_22, vx1, rnd);

	/*
	 * store rotation matrix for second clothoid
	 * mirror in O, rotate to v2
	 * |vx2 -vy2| \/ |-1 0|  _ |-vx2  vy2|
	 * |vy2  vx2| /\ |0 -1|  - |-vy2 -vx2|
	 */
	mpfr_neg(c->r2_11, vx2, rnd);
	mpfr_set(c->r2_12, vy2, rnd);
	mpfr_neg(c->r2_21, vy2, rnd);
	mpfr_neg(c->r2_22, vx2, rnd);

	mpfr_set(c->vl, vl, rnd);

	/*
	 * set total path length
	 */
	mpfr_mul(c->L, c->a, c->b, rnd);
	mpfr_mul(c->L, c->L, c->t, rnd);

	if (left) {
		/* swap results back */
		mpfr_set(tmp, c->px1, rnd);
		mpfr_set(c->px1, c->px2, rnd);
		mpfr_set(c->px2, tmp, rnd);
		mpfr_set(tmp, c->py1, rnd);
		mpfr_set(c->py1, c->py2, rnd);
		mpfr_set(c->py2, tmp, rnd);

		mpfr_set(tmp, c->r1_11, rnd);
		mpfr_set(c->r1_11, c->r2_11, rnd);
		mpfr_set(c->r2_11, tmp, rnd);
		mpfr_set(tmp, c->r1_12, rnd);
		mpfr_set(c->r1_12, c->r2_12, rnd);
		mpfr_set(c->r2_12, tmp, rnd);
		mpfr_set(tmp, c->r1_21, rnd);
		mpfr_set(c->r1_21, c->r2_21, rnd);
		mpfr_set(c->r2_21, tmp, rnd);
		mpfr_set(tmp, c->r1_22, rnd);
		mpfr_set(c->r1_22, c->r2_22, rnd);
		mpfr_set(c->r2_22, tmp, rnd);
	}

	mpfr_clears(phi, tmp, pi, x, dx, dy, vx1, vy1, vx2, vy2, mx, my, ml,
		vl, rx, ry, NULL);

	return pe;
}

void
free_clothoid(path_elem_t *pe)
{
	clothoid_t *c = C(pe);

	mpfr_clears(c->a, c->b, c->t, c->vl, c->L, NULL);
	mpfr_clears(c->px1, c->py1, c->px2, c->py2, NULL);
	mpfr_clears(c->r1_11, c->r1_12, c->r1_21, c->r1_22, NULL);
	mpfr_clears(c->r2_11, c->r2_12, c->r2_21, c->r2_22, NULL);
	free(pe);
}

void
print_clothoid(clothoid_t *c)
{
	mpfr_printf(
		"a %.5Re b %.5Re t %.5Re start %.5Rf/%.5Rf end %.5Rf/%.5Rf\n",
		c->a, c->b, c->t, c->px1, c->py1, c->px2, c->py2);
	mpfr_printf("r1 %.5Rf %.5Rf %.5Rf %.5Rf\n",
		c->r1_11, c->r1_12, c->r1_21, c->r1_22);
	mpfr_printf("r2 %.5Rf %.5Rf %.5Rf %.5Rf\n",
		c->r2_11, c->r2_12, c->r2_21, c->r2_22);
}

/*
 * calculate a point on the clothoid c at time t, return x/y
 */
void
calcclothoid(clothoid_t *c, mpfr_t t, mpfr_t x, mpfr_t y,
	mpfr_t vx, mpfr_t vy, mpfr_t ax, mpfr_t ay)
{
	mpfr_t tmp, pab3, t_sin, t_cos, _t;
	mpfr_inits(tmp, pab3, t_sin, t_cos, _t, NULL);

	if (mpfr_cmp(t, c->t) <= 0) {
		mpfr_set(_t, t, rnd);
		_calcclothoid(x, y, c->a, c->b, _t);
		/* rotate + translate point */
		mpfr_mul(tmp, x, c->r1_11, rnd);
		mpfr_fma(tmp, y, c->r1_12, tmp, rnd);
		mpfr_mul(y, y, c->r1_22, rnd);
		mpfr_fma(y, x, c->r1_21, y, rnd);
		mpfr_set(x, tmp, rnd);

		mpfr_add(x, x, c->px1, rnd);
		mpfr_add(y, y, c->py1, rnd);
	} else {
		mpfr_sub(tmp, c->t, t, rnd);
		mpfr_add(_t, tmp, c->t, rnd);
		_calcclothoid(x, y, c->a, c->b, _t);

		/* rotate + translate point */
		mpfr_mul(tmp, x, c->r2_11, rnd);
		mpfr_fma(tmp, y, c->r2_12, tmp, rnd);
		mpfr_mul(y, y, c->r2_22, rnd);
		mpfr_fma(y, x, c->r2_21, y, rnd);
		mpfr_set(x, tmp, rnd);

		mpfr_add(x, x, c->px2, rnd);
		mpfr_add(y, y, c->py2, rnd);
	}
	/* x', y' */
	mpfr_const_pi(tmp, rnd);
	mpfr_div_ui(tmp, tmp, 2, rnd);
	mpfr_mul(tmp, tmp, _t, rnd);
	mpfr_mul(tmp, tmp, _t, rnd);
	mpfr_mul(tmp, tmp, c->b, rnd);
	mpfr_mul(tmp, tmp, c->b, rnd);
	mpfr_sin(t_sin, tmp, rnd);
	mpfr_cos(t_cos, tmp, rnd);
	mpfr_mul(vx, t_cos, c->a, rnd);
	mpfr_mul(vx, vx, c->b, rnd);
	mpfr_mul(vy, t_sin, c->a, rnd);
	mpfr_mul(vy, vy, c->b, rnd);

	/* x", y" */
	mpfr_const_pi(ax, rnd);
	mpfr_mul(ax, ax, vy, rnd);
	mpfr_mul(ax, ax, _t, rnd);
	mpfr_mul(ax, ax, c->b, rnd);
	mpfr_mul(ax, ax, c->b, rnd);
	mpfr_neg(ax, ax, rnd);
	mpfr_const_pi(ay, rnd);
	mpfr_mul(ay, ay, vx, rnd);
	mpfr_mul(ay, ay, _t, rnd);
	mpfr_mul(ay, ay, c->b, rnd);
	mpfr_mul(ay, ay, c->b, rnd);

#if 0
mpfr_printf("before reversal: x %.05Rf y %.05Rf vx %.05Rf vy %.05Rf ax %.05Rf ay %.05Rf\n", x, y, vx, vy, ax, ay);
#endif
	if (mpfr_cmp(t, c->t) > 0) {
		/*
		 * as we're going backward along the second
		 * clothoid, reverse the direction of v
		 */
		mpfr_neg(vx, vx, rnd);
		mpfr_neg(vy, vy, rnd);

#if 0
		/*
		 * Reversing a is a bit more involved.
		 * By differentiating C(t0-t) two times, we
		 * get a corrective term in addition to the
		 * negation
		 */
		mpfr_const_pi(pab3, rnd);
		mpfr_mul(pab3, pab3, c->a, rnd);
		mpfr_mul(pab3, pab3, c->b, rnd);
		mpfr_mul(pab3, pab3, c->b, rnd);
		mpfr_mul(pab3, pab3, c->b, rnd);
		mpfr_mul(pab3, pab3, c->t, rnd);
		mpfr_neg(tmp, pab3, rnd);
		mpfr_mul(tmp, tmp, t_sin, rnd);
		mpfr_sub(ax, tmp, ax, rnd);
		mpfr_mul(tmp, pab3, t_cos, rnd);
		mpfr_sub(ay, tmp, ay, rnd);
#endif
	}

#if 0
mpfr_printf("before rotation: x %.05Rf y %.05Rf vx %.05Rf vy %.05Rf ax %.05Rf ay %.05Rf\n", x, y, vx, vy, ax, ay);
#endif

	/* rotate */
	if (mpfr_cmp(t, c->t) <= 0) {
		mpfr_mul(tmp, vx, c->r1_11, rnd);
		mpfr_fma(tmp, vy, c->r1_12, tmp, rnd);
		mpfr_mul(vy, vy, c->r1_22, rnd);
		mpfr_fma(vy, vx, c->r1_21, vy, rnd);
		mpfr_set(vx, tmp, rnd);

		mpfr_mul(tmp, ax, c->r1_11, rnd);
		mpfr_fma(tmp, ay, c->r1_12, tmp, rnd);
		mpfr_mul(ay, ay, c->r1_22, rnd);
		mpfr_fma(ay, ax, c->r1_21, ay, rnd);
		mpfr_set(ax, tmp, rnd);
	} else {
		mpfr_mul(tmp, vx, c->r2_11, rnd);
		mpfr_fma(tmp, vy, c->r2_12, tmp, rnd);
		mpfr_mul(vy, vy, c->r2_22, rnd);
		mpfr_fma(vy, vx, c->r2_21, vy, rnd);
		mpfr_set(vx, tmp, rnd);

		mpfr_mul(tmp, ax, c->r2_11, rnd);
		mpfr_fma(tmp, ay, c->r2_12, tmp, rnd);
		mpfr_mul(ay, ay, c->r2_22, rnd);
		mpfr_fma(ay, ax, c->r2_21, ay, rnd);
		mpfr_set(ax, tmp, rnd);
	}
	mpfr_clears(tmp, pab3, t_sin, t_cos, _t, NULL);
}

#ifdef TEST
int
main(int argc, char **argv)
{
	int i;

	mpfr_set_default_prec(300);

	mpfr_t ix, iy, vx1, vy1, vx2, vy2, v, acc;
	mpfr_inits(ix, iy, vx1, vy1, vx2, vy2, v, acc, NULL);

	mpfr_set_ui(ix,   30, rnd);
	mpfr_set_ui(iy,  100, rnd);
	mpfr_set_si(vx1,  20, rnd);
	mpfr_set_si(vy1,  90, rnd);
	mpfr_set_si(vx2, 100, rnd);
	mpfr_set_si(vy2, -10, rnd);
	mpfr_set_ui(v,   100, rnd);
	mpfr_set_ui(acc, 500, rnd);

	clothoid_t *c;

	c = cloth_line(ix, iy, vx1, vy1, vx2, vy2, v, acc);

	print_clothoid(c);

	mpfr_t t, dt, x, y, vx, vy, ax, ay;
	mpfr_inits(t, dt, x, y, vx, vy, ax, ay, NULL);

	mpfr_set_ui(t, 0, rnd);
	mpfr_div_ui(dt, c->t, 400, rnd);

	for (i = 0; i < 800; ++i) {
		calcclothoid(c, t, x, y, vx, vy, ax, ay);
#if 1
		mpfr_printf("[%d] [%.5Rf] %.5Rf %.5Rf %.5Rf %.5Rf %.5Rf %.5Rf\n", i, t, x, y, vx, vy, ax, ay);
#endif
		mpfr_add(t, t, dt, rnd);
	}
	free_clothoid(c);

	return 0;
}
#endif
