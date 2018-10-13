#include <stdio.h>
#include <stdlib.h>

#include "conan.h"

mpfr_rnd_t rnd = MPFR_RNDN;

void
free_path_elem(path_elem_t *pe)
{
	switch (pe->type) {
	case PT_LINE:
		free_line(pe);
		break;
	case PT_CLOTHOID:
		free_clothoid(pe);
		break;
	default:
		printf("unknown type %d\n", pe->type);
		abort();
	}
}

void
print_path_elem(int i, path_elem_t *pe)
{
	printf("[%d] type %d\n", i, pe->type);
	switch (pe->type) {
	case PT_LINE:
		print_line(pe->line);
		break;
	case PT_CLOTHOID:
		print_clothoid(pe->clothoid);
		break;
	default:
		printf("unknown type %d\n", pe->type);
		abort();
	}
}

void
print_path(path_t *p)
{
	int i;

	printf("path has %d elements\n", p->nelem);

	for (i = 0; i < p->nelem; ++i) {
		path_elem_t *pe = p->p[i];

		print_path_elem(i, pe);
	}
}

int
main(int argc, char **argv)
{
	path_t *path;
	int ret;
	int npoints;
	point_t *points;
	int i;

	mpfr_set_default_prec(300);

	mpfr_t max_acc, max_v;
	mpfr_inits(max_acc, max_v, NULL);

	mpfr_set_ui(max_acc, 1000, rnd);
	mpfr_set_d(max_v, 100, rnd);

	npoints = 5;
	points = calloc(sizeof(*points), npoints);
	
	mpfr_init_set_ui(points[0].x, 20, rnd);
	mpfr_init_set_ui(points[0].y, 20, rnd);
	mpfr_init_set_ui(points[1].x, 30, rnd);
	mpfr_init_set_ui(points[1].y, 80, rnd);
	mpfr_init_set_ui(points[2].x, 70, rnd);
	mpfr_init_set_ui(points[2].y, 90, rnd);
	mpfr_init_set_ui(points[3].x, 120, rnd);
	mpfr_init_set_ui(points[3].y, 50, rnd);
	mpfr_init_set_ui(points[4].x, 180, rnd);
	mpfr_init_set_ui(points[4].y, 70, rnd);
	ret = plan(max_acc, max_v, points, npoints, &path);
	printf("plan retourned %d\n", ret);

	print_path(path);

	motion_t m;

	ret = motion_init(&m, KIN_COREXY, "coeff.txt", "path.csv", 10000);
	if (ret < 0)
		exit(1);


	mpfr_t x, y, vx, vy, ax, ay, t, delta, toff;
	mpfr_inits(x, y, vx, vy, ax, ay, t, delta, toff, NULL);

	/* for now: initialize motion for first point */
	path_elem_t *pe = path->p[0];
	pe->reset(pe);
	mpfr_set_ui(t, 0, rnd);
	pe->calc(pe, t, x, vx, ax, y, vy, ay);
	motion_set_xy_v_a(&m, x, vx, ax, y, vy, ay);

	for (i = 0; i < path->nelem; ++i) {
		path_elem_t *pe = path->p[i];

print_path_elem(i, pe);
		pe->reset(pe);
		mpfr_set_ui(t, 0, rnd);
		while (1) {
			mpfr_set(delta, t, rnd);
			mpfr_set(toff, t, rnd);
			ret = pe->next(pe, t);
			if (ret == 0)
				break;
			mpfr_sub(delta, t, delta, rnd);
			pe->calc(pe, t, x, vx, ax, y, vy, ay);
mpfr_printf("calc: t %.5Rf x %.5Rf y %.5Rf vx %.5Rf vy %.5Rf ax %.5Rf ay %.5Rf \n",
		t, x, y, vx, vy, ax, ay);
#if 1
			motion_move(&m, delta, x, vx, ax, y, vy, ay,
				pe->calc, pe, toff);
#endif
		}
	}

	return 0;
}
