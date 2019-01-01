#ifndef __CONAN_H__
#define __CONAN_H__

#include <stdint.h>
#include <stdio.h>
#include <assert.h>
#include <mpfr.h>

extern mpfr_rnd_t rnd;

struct _path_elem;

#ifndef min
#define min(x, y) ((x)<(y)?(x):(y))
#endif
#ifndef max
#define max(x, y) ((x)>(y)?(x):(y))
#endif

/*
 * motion.c
 */

typedef enum {
	KIN_CARTESIAN	= 1,
	KIN_COREXY	= 2,
} kinematics_t;

typedef void (*func_cb_t)(void *ctx, mpfr_t t, mpfr_t x, mpfr_t y,
	mpfr_t vx, mpfr_t vy, mpfr_t ax, mpfr_t ay);

typedef struct _motion {
	kinematics_t	m_kinematics;

	/*
	 * scaling info
	 */
	mpfr_t		m_scale_f;
	mpfr_t		m_scale;
	mpz_t		m_scale_z;
	mpz_t		m_maxval;
	mpz_t		m_minval;
	uint64_t	m_hz;
	mpz_t		m_hz_z;
	mpfr_t		m_freq;
	int		m_bits;
	mpfr_t		m_pi;

	/*
	 * device state
	 */
	uint64_t	m_absstep;
	mpz_t		m_x;
	mpz_t		m_y;
	/* already transformed values, for next starting point */
	mpfr_t		m_vx;
	mpfr_t		m_ax;
	mpfr_t		m_vy;
	mpfr_t		m_ay;
	/*
	 * output
	 */
	FILE		*m_coeff;
	FILE		*m_path;
	int		m_first_line;

	/*
	 * scratchpad, we keep them around here to save inits
	 */
	mpfr_t		m_tmp;

	/*
	 * simulation info
	 * callback function to get reference value, for error checking
	 */
	uint64_t	m_simsteps;

	/*
	 * limits
	 */
	mpfr_t		m_max_acc;
	mpfr_t		m_max_v;
} motion_t;

int
motion_init(motion_t *mp, mpfr_t max_acc, mpfr_t max_v,
	kinematics_t kinematics, const char *coeff_file,
	const char *path_file, int simsteps);
void motion_move_origin(motion_t *mp, mpfr_t e0_x, mpfr_t e0_y);
void
motion_set_xy_v_a(motion_t *mp,
	mpfr_t x, mpfr_t vx, mpfr_t ax,
	mpfr_t y, mpfr_t vy, mpfr_t ay);
void
motion_move(motion_t *mp, mpfr_t t,
	mpfr_t e1_x, mpfr_t e1_vx, mpfr_t e1_ax,
	mpfr_t e1_y, mpfr_t e1_vy, mpfr_t e1_ay,
	func_cb_t func, void *func_ctx, mpfr_t func_toff);

/*
 * clothoid.c
 */
typedef struct _clothoid {
	int	ptr;
	int	nsegments;
	mpfr_t	a;
	mpfr_t	b;
	mpfr_t	t;
	mpfr_t	vl;
	mpfr_t	L;		/* total path length */
	mpfr_t	px1, py1;	/* start point */
	mpfr_t	px2, py2;	/* end point */
	mpfr_t	r1_11, r1_12, r1_21, r1_22; /* rotation matrix for 1st half */
	mpfr_t	r2_11, r2_12, r2_21, r2_22; /* rotation matrix for 2nd half */
} clothoid_t;

void
calcclothoid(clothoid_t *c, mpfr_t t, mpfr_t x, mpfr_t y,
	mpfr_t vx, mpfr_t vy, mpfr_t ax, mpfr_t ay);
void print_clothoid(clothoid_t *c);
void free_clothoid(struct _path_elem *pe);
struct _path_elem *
cloth_line(mpfr_t inter_x, mpfr_t inter_y, mpfr_t _vx1, mpfr_t _vy1,
	mpfr_t _vx2, mpfr_t _vy2, mpfr_t v, mpfr_t acc);

/*
 * circle.c
 */

typedef struct _circle_ctx {
	mpfr_t	omega;
	mpfr_t	r;
	mpfr_t	start_t;
	mpfr_t	t;
} circle_ctx_t;
void
calccircle(mpfr_t omega, mpfr_t t, mpfr_t r,
	mpfr_t x, mpfr_t vx, mpfr_t ax, mpfr_t jx,
	mpfr_t y, mpfr_t vy, mpfr_t ay, mpfr_t jy);

/*
 * line.c
 */
typedef struct _line {
	int	nsegments;
	int	ptr;
	mpfr_t	start_x;
	mpfr_t	start_y;
	mpfr_t	start_v;
	mpfr_t	end_x;
	mpfr_t	end_y;
	mpfr_t	end_v;
	struct _line_segment {
		mpfr_t	jerk_x;
		mpfr_t	jerk_y;
		mpfr_t	time;
	} segments[];
} line_t;
typedef struct _line_segment line_segment_t;

struct _path_elem *
const_line(mpfr_t v, mpfr_t x1, mpfr_t y1, mpfr_t x2, mpfr_t y2);
void print_line(line_t *l);
void free_line(struct _path_elem *pe);

/*
 * planner.c
 */

typedef struct _point {
	mpfr_t	x;
	mpfr_t	y;
} point_t;

#define PT_CLOTHOID	1
#define PT_LINE		2
typedef struct _path_elem {
	union {
		clothoid_t	*clothoid;
		line_t		*line;
	};
	int type;
	void (*reset)(struct _path_elem *);
	int (*next)(struct _path_elem *, mpfr_t t);
	void (*calc)(void *, mpfr_t t, mpfr_t x, mpfr_t y,
		mpfr_t vx, mpfr_t vy, mpfr_t ax, mpfr_t ay);
} path_elem_t;

typedef struct _path {
	int		nelem;
	path_elem_t	**p;
	mpfr_t		target_v;
} path_t;

void free_path_elem(path_elem_t *pe);

static inline clothoid_t *C(path_elem_t *pe) {
	assert(pe->type == PT_CLOTHOID);
	return pe->clothoid;
}

static inline line_t *L(path_elem_t *pe) {
	assert(pe->type == PT_LINE);
	return pe->line;
}

int
plan(mpfr_t max_acc, mpfr_t max_v, point_t *points, int npoints, path_t **path);

/*
 * gcode.c
 */
int read_gcode(motion_t *m, const char *file);

/*
 * conan.c
 */
int plan_and_execute(motion_t *m, point_t *points, int npoints);

/*
 * b-approx.c
 */
void match_path(motion_t *m, point_t *points, int npoints);

#endif
