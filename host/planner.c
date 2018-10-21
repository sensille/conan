#include <stdlib.h>
#include <assert.h>
#include <string.h>

#include "conan.h"

/*
motion planner v1
input: series of lines
output: list of geom_t pointers to play
two passes
first pass: determine max speed for each 3 lines, step to next line
second pass: render at min[max_speeds]
*/

/*
 * calculate the maximum speed across two junctions.
 * Do this by calculating the total line length that would be needed at max
 * speed. If this exceeds the actual line length, approximatively find the
 * max possible speed.
 * lx1/ly1: first junction
 j lx2/ly2: second junction
 * vx1/vy1: line direction and length into the first junction (optional)
 * vx2/vy2: line direction and length out of the second junction (optional)
 * The line under consideration is the line between lx1/ly1 and lx2/ly2.
 * The function calculates both the clothoid to get into this line and the one
 * the get out of it. It might be that at max speed the line is not long enough
 * to accomodate both clothoids. In that case lower the speed until it fits.
 * If v1 is not given, it is assumed that the second clothoid may take the
 * full line. If v2 is not given, the second clothoid may take the full line.
 */
static int
get_max_speed(mpfr_t max_acc, mpfr_t max_v,
	mpfr_t lx1, mpfr_t ly1, mpfr_t lx2, mpfr_t ly2,
	mpfr_t vx1, mpfr_t vy1, mpfr_t vx2, mpfr_t vy2, mpfr_t v)
{
	path_elem_t *c;
	int i;

mpfr_t nan; mpfr_init(nan);
mpfr_printf("xx: lx1 %.5Rf ly1 %.5Rf lx2 %.5Rf ly2 %.5Rf vx1 %.5Rf vy1 %.5Rf vx2 %.5Rf vy2 %.5Rf v %.5Rf\n",
	lx1, ly1, lx2, ly2, vx1 ? vx1 : nan, vy1 ? vy1 : nan, vx2 ? vx2 : nan, vy2 ? vy2 : nan, max_v);

	mpfr_t dx, dy, dl, l, low_v, hi_v;
	mpfr_inits(dx, dy, dl, l, low_v, hi_v, NULL);
	mpfr_set_ui(l, 0, rnd);

	mpfr_sub(dx, lx2, lx1, rnd);
	mpfr_sub(dy, ly2, ly1, rnd);
	mpfr_hypot(dl, dx, dy, rnd);

	if (vx1 != NULL) {
		c = cloth_line(lx1, ly1, vx1, vy1, dx, dy, max_v, max_acc);
		if (c == NULL)
			return -1;
mpfr_printf("xx: from c1: l %.5Rf\n", C(c)->vl);
		mpfr_add(l, l, C(c)->vl, rnd);
		free_path_elem(c);
	}
	if (vx2 != NULL) {
		c = cloth_line(lx2, ly2, dx, dy, vx2, vy2, max_v, max_acc);
		if (c == NULL)
			return -1;
mpfr_printf("xx: from c2: l %.5Rf\n", C(c)->vl);
		mpfr_add(l, l, C(c)->vl, rnd);
		free_path_elem(c);
	}

mpfr_printf("xx: l %.5Rf dl %.5Rf\n", l, dl);
	if (mpfr_cmp(l, dl) < 0) {
		/*
		 * we can run with max_v
		 */
		mpfr_set(v, max_v, rnd);
		return 0;
	}

printf("xx: approximation needed\n");
	/*
	 * approximate a max v to use
	 */
	mpfr_set_ui(low_v, 0, rnd);
	mpfr_set(hi_v, max_v, rnd);

	/* XXX TODO use better approximation */
	for (i = 0; i < 20; ++i) {
		mpfr_add(v, low_v, hi_v, rnd);
		mpfr_div_ui(v, v, 2, rnd);
		mpfr_set_ui(l, 0, rnd);

		if (vx1 != NULL) {
			c = cloth_line(lx1, ly1, vx1, vy1, dx, dy, v, max_acc);
			if (c == NULL)
				return -1;
mpfr_printf("xx: a. from c1: l %.5Rf start %.5Rf/%.5Rf end %.5Rf/%.5Rf\n", C(c)->vl, C(c)->px1, C(c)->py1, C(c)->px2, C(c)->py2);
			mpfr_add(l, l, C(c)->vl, rnd);
			free_path_elem(c);
		}
		if (vx2 != NULL) {
			c = cloth_line(lx2, ly2, dx, dy, vx2, vy2, v, max_acc);
			if (c == NULL)
				return -1;
mpfr_printf("xx: a. from c2: l %.5Rf start %.5Rf/%.5Rf end %.5Rf/%.5Rf\n", C(c)->vl, C(c)->px1, C(c)->py1, C(c)->px2, C(c)->py2);
			mpfr_add(l, l, C(c)->vl, rnd);
			free_path_elem(c);
		}

mpfr_printf("xx: v %.5Rf low_v %.5Rf high_v %.5Rf l %.5Rf dl %.5Rf\n", v, low_v, hi_v, l, dl);
		if (mpfr_cmp(l, dl) < 0)
			mpfr_set(low_v, v, rnd);
		else
			mpfr_set(hi_v, v, rnd);
			
	
	}
	mpfr_set(v, low_v, rnd);
mpfr_printf("xx: approximation done: v %.5Rf\n", v);

	mpfr_clears(dx, dy, dl, l, low_v, hi_v, NULL);

	return 0;
}

int
plan(mpfr_t max_acc, mpfr_t max_v, point_t *_points, int _npoints,
	path_t **path)
{
	int i;
	int ret;
	mpfr_t vx1, vy1, vx2, vy2, v, target_v, tmp, thresh;
	path_elem_t *pe;
	path_t *pp;
	point_t prev_endpoint;
	path_elem_t *c;
	path_elem_t *l;
	int eix = 0;
	int max_elem;
	point_t *points;
	int npoints = _npoints;

	mpfr_inits(vx1, vy1, vx2, vy2, v, target_v, tmp, thresh, NULL);

printf("--- planner start ---\n");
	/*
	 * copy points for preprocessing
	 */
	points = calloc(sizeof(*points), npoints);
	for (i = 0; i < npoints; ++i) {
		mpfr_inits(points[i].x, points[i].y, NULL);
		mpfr_set(points[i].x, _points[i].x, rnd);
		mpfr_set(points[i].y, _points[i].y, rnd);
	}

	/*
	 * preprocessing: remove duplicate points
	 */
	mpfr_set_d(thresh, 0.001, rnd);
	for (i = 0; npoints > 1 && i < npoints - 1; ++i) {
		mpfr_sub(tmp, points[i].x, points[i+1].x, rnd);
		mpfr_abs(tmp, tmp, rnd);
		if (mpfr_cmp(tmp, thresh) > 0)
			continue;
		mpfr_sub(tmp, points[i].y, points[i+1].y, rnd);
		mpfr_abs(tmp, tmp, rnd);
		if (mpfr_cmp(tmp, thresh) > 0)
			continue;
		printf("remove point %d\n", i);
		memmove(points + i, points + i + 1,
			sizeof(*points) * (npoints - i));
		--npoints;
	}

	/*
	 * preprocessing: remove points in a line
	 */
	mpfr_set_d(thresh, 0.00001, rnd);
	for (i = 0; npoints > 2 && i < npoints - 2; ++i) {
		mpfr_sub(vx1, points[i].x, points[i+1].x, rnd);
		mpfr_sub(vy1, points[i].y, points[i+1].y, rnd);
		mpfr_sub(vx2, points[i+1].x, points[i+2].x, rnd);
		mpfr_sub(vy2, points[i+1].y, points[i+2].y, rnd);
#if 0
		mpfr_hypot(tmp, vx1, vy1, rnd);
		mpfr_div(vx1, vx1, tmp, rnd);
		mpfr_div(vy1, vy1, tmp, rnd);
		mpfr_hypot(tmp, vx2, vy2, rnd);
		mpfr_div(vx2, vx2, tmp, rnd);
		mpfr_div(vy2, vy2, tmp, rnd);
#endif
		mpfr_mul(tmp, vx1, vy2, rnd);
		mpfr_fms(tmp, vx2, vy1, tmp, rnd);
		mpfr_abs(tmp, tmp, rnd);
mpfr_printf("point %d v1 %.5Re/%.5Re v2 %.5Re/%.5Re cp %.5Re\n", i, vx1, vy1, vx2, vy2, tmp);
		if (mpfr_cmp(tmp, thresh) > 0)
			continue;
		mpfr_printf("remove point %d from line, cross product is %.5Re\n",
			i, tmp);
		memmove(points + i + 1, points + i + 2,
			sizeof(*points) * (npoints - i - 1));
		--npoints;
	}

	pp = calloc(sizeof(*pp), 1);
	mpfr_init(pp->target_v);
	max_elem = 2 * npoints - 3;
	pp->p = calloc(sizeof(*pe), max_elem);

	mpfr_set(target_v, max_v, rnd);

	/* pass one: get max speed */
	for (i = 0; i < npoints - 1; ++i) {
		if (i > 0) {
			mpfr_sub(vx1, points[i].x, points[i-1].x, rnd);
			mpfr_sub(vy1, points[i].y, points[i-1].y, rnd);
		}
		if (i < npoints - 2) {
			mpfr_sub(vx2, points[i+2].x, points[i+1].x, rnd);
			mpfr_sub(vy2, points[i+2].y, points[i+1].y, rnd);
		}
		ret = get_max_speed(max_acc, target_v,
			points[i].x, points[i].y,
			points[i+1].x, points[i+1].y,
			i > 0 ? vx1 : NULL, i > 0 ? vy1 : NULL,
			i < npoints-2 ? vx2 : NULL, i < npoints-2 ? vy2 : NULL,
			v);
		if (ret)
			return ret;
		if (mpfr_cmp(v, target_v) < 0)
			mpfr_set(target_v, v, rnd);
	}

mpfr_printf("xx: max speed is %.5Rf\n", target_v);
	prev_endpoint = points[0];
	/* pass two: render path */
	for (i = 1; i < npoints - 1; ++i) {
		point_t *p1 = &points[i];
		point_t *p2 = &points[i + 1];

		mpfr_sub(vx1, p1->x, prev_endpoint.x, rnd);
		mpfr_sub(vy1, p1->y, prev_endpoint.y, rnd);
		mpfr_sub(vx2, p2->x, p1->x, rnd);
		mpfr_sub(vy2, p2->y, p1->y, rnd);

		/*
		 * calc the clothoid, so we know the start point of it, which
		 * is the endpoint of the previous line
		 */
		c = cloth_line(p1->x, p1->y, vx1, vy1, vx2, vy2, target_v, max_acc);
		assert(c != NULL);
mpfr_printf("xx: c: l %.5Rf start %.5Rf/%.5Rf end %.5Rf/%.5Rf\n", C(c)->vl, C(c)->px1, C(c)->py1, C(c)->px2, C(c)->py2);

		/*
		 * draw the line connecting to the clothoid
		 */
		l = const_line(target_v, prev_endpoint.x, prev_endpoint.y,
			C(c)->px1, C(c)->py1);
		
		if (l != NULL)
			pp->p[eix++] = l;
		pp->p[eix++] = c;

		mpfr_set(prev_endpoint.x, C(c)->px2, rnd);
		mpfr_set(prev_endpoint.y, C(c)->py2, rnd);
	}
	/*
	 * finish the path by a line from the endpoint of the clothoid to the
	 * final point
	 */
	l = const_line(target_v, prev_endpoint.x, prev_endpoint.y,
		points[npoints - 1].x, points[npoints - 1].y);
	if (l != NULL)
		pp->p[eix++] = l;

	assert(eix <= max_elem);

	pp->nelem = eix;

	mpfr_set(pp->target_v, target_v, rnd);

	mpfr_clears(vx1, vy1, vx2, vy2, v, target_v, NULL);

printf("--- planner end ---\n");
	*path = pp;

	return 0;
}

