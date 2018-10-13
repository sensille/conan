#include <stdlib.h>
#include <assert.h>

#include "conan.h"

void
line_reset(path_elem_t *pe)
{
	line_t *l = pe->line;

	l->ptr = 1;
}

int
line_next(path_elem_t *pe, mpfr_t t)
{
	line_t *l = pe->line;

	assert(l->nsegments == 1);

	if (l->ptr > 1)
		return 0;

	mpfr_set(t, l->segments[0].time, rnd);
	++l->ptr;

	return 1;
}

void
line_calc(void *ctx, mpfr_t t, mpfr_t x, mpfr_t vx,
	mpfr_t ax, mpfr_t y, mpfr_t vy, mpfr_t ay)
{
	path_elem_t *pe = ctx;
	line_t *l = pe->line;
	mpfr_t tmp;
	mpfr_init(tmp);

	mpfr_sub(x, l->end_x, l->start_x, rnd);
	mpfr_sub(y, l->end_y, l->start_y, rnd);
	mpfr_set(vx, x, rnd);
	mpfr_set(vy, y, rnd);
	mpfr_div(tmp, t, l->segments[0].time, rnd);
	mpfr_mul(x, x, tmp, rnd);
	mpfr_mul(y, y, tmp, rnd);
	mpfr_add(x, x, l->start_x, rnd);
	mpfr_add(y, y, l->start_y, rnd);
	mpfr_hypot(tmp, vx, vy, rnd);
	mpfr_div(vx, vx, tmp, rnd);
	mpfr_div(vy, vy, tmp, rnd);
	mpfr_mul(vx, vx, l->start_v, rnd);
	mpfr_mul(vy, vy, l->start_v, rnd);
	mpfr_set_ui(ax, 0, rnd);
	mpfr_set_ui(ay, 0, rnd);

	mpfr_clear(tmp);
}

path_elem_t *
const_line(mpfr_t v, mpfr_t x1, mpfr_t y1, mpfr_t x2, mpfr_t y2)
{
	line_t *l;
	path_elem_t *pe = calloc(sizeof(*pe) + sizeof(*l) +
		sizeof(line_segment_t), 1);
	
	assert(pe);

	l = (line_t *)(pe + 1);
	pe->line = l;
	pe->type = PT_LINE;
	pe->reset = line_reset;
	pe->next = line_next;
	pe->calc = line_calc;
	
	mpfr_t dx, dy, time;
	mpfr_inits(dx, dy, time, NULL);
	mpfr_sub(dx, x2, x1, rnd);
	mpfr_sub(dy, y2, y1, rnd);
	mpfr_hypot(time, dx, dy, rnd);
	mpfr_div(time, time, v, rnd);

	mpfr_init_set(l->start_x, x1, rnd);
	mpfr_init_set(l->start_y, y1, rnd);
	mpfr_init_set(l->start_v, v, rnd);
	mpfr_init_set(l->end_x, x2, rnd);
	mpfr_init_set(l->end_y, y2, rnd);
	mpfr_init_set(l->end_v, v, rnd);

	l->nsegments = 1;
	mpfr_init_set_ui(l->segments[0].jerk_x, 0, rnd);
	mpfr_init_set_ui(l->segments[0].jerk_y, 0, rnd);
	mpfr_init_set(l->segments[0].time, time, rnd);

	mpfr_clears(dx, dy, time, NULL);

	return pe;
}

void
free_line(path_elem_t *pe)
{
	int i;
	line_t *l = L(pe);

	for (i = 0; i < l->nsegments; ++i) {
		mpfr_clear(l->segments[i].jerk_x);
		mpfr_clear(l->segments[i].jerk_y);
		mpfr_clear(l->segments[i].time);
	}
	mpfr_clear(l->start_x);
	mpfr_clear(l->start_y);
	mpfr_clear(l->start_v);
	mpfr_clear(l->end_x);
	mpfr_clear(l->end_y);
	mpfr_clear(l->end_v);
	free(pe);
}

void
print_line(line_t *l)
{
	mpfr_printf("start %.5Rf/%.5Rf v %.5Rf end %.5Rf/%.Rf v %.Rf\n",
		l->start_x, l->start_y, l->start_v,
		l->end_x, l->end_y, l->end_v);
}
