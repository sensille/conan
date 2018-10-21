#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "conan.h"

#define MAXLINE	1000
#define MAXCODES	26
#define P_CHUNK	10

typedef struct _gccmd {
	char	codes[MAXCODES];
	double	values[MAXCODES];
} gccmd_t;

static int
parse_line(const char *line, gccmd_t *gc)
{
	const char *l = line;
	char val[100];
	int vp = 0;
	int cp = 0;
	int current = -1;

	val[vp] = 0;

	/* XXX valid with any double implementation? */
	memset(gc, 0, sizeof(*gc));

	while (*l) {
		/* skip all whitespace */
		if (strchr(" \t\r\n", *l)) {
			++l;
			continue;
		}
		if (*l == ';')	/* comment, done */
			break;

		/* one letter code */
		if ((*l >= 'A' && *l <= 'Z') || (*l >= 'a' && *l <= 'z')) {
			if (current >= 0)
				gc->values[current] = atof(val);
			val[0] = 0;
			vp = 0;
			gc->codes[cp++] = toupper(*l);
			current = toupper(*l) - 'A';
			if (cp == MAXCODES) {
				printf("too many codes in line %s\n", line);
				return -1;
			}
			++l;
		} else if (strchr("0123456789+-.", *l)) {
			val[vp++] = *l++;
			val[vp] = 0;
			if (vp == sizeof(val)) {
				printf("val too long in %s\n", line);
				return -1;
			}
		} else {
			printf("failed to parse cmd/param %c\n", *l);
			return -1;
		}
	}
	if (current >= 0)
		gc->values[current] = atof(val);

	return 0;
}

void
dump_gccmd(gccmd_t *gc)
{
	int i;

	for (i = 0; i < MAXCODES; ++i) {
		if (gc->codes[i] == 0)
			break;
		if (i != 0)
			printf("  ");
		printf("%c %f\n", gc->codes[i], gc->values[gc->codes[i] - 'A']);
	}
}

static int
check_params(gccmd_t *gc, const char *mandatory, const char *optional)
{
	const char *p;

	/* see if all mandatory parameters are present */
	for (p = mandatory; *p; ++p) {
		if (strchr(gc->codes + 1, *p) == NULL) {
			printf("missing mandatory parameter %c\n", *p);
			return -1;
		}
	}
	/* see if all parameters are either mandatory or optional */
	for (p = gc->codes + 1; *p; ++p) {
		if (strchr(mandatory, *p) == NULL &&
		    strchr(optional, *p) == NULL) {
			printf("invalid parameter %c given\n", *p);
			return -1;
		}
	}

	return 0;
}

static void
set_if_set(gccmd_t *gc, char p, double *v)
{
	if (strchr(gc->codes + 1, p) == NULL)
		return;
	*v = gc->values[p - 'A'];
}

static void
push_points(motion_t *m, point_t *points, int npoints)
{
	int i;

	for (i = 0; i < npoints; ++i) {
		mpfr_printf("point: %.5Rf %.5Rf\n", points[i].x, points[i].y);
	}
	if (npoints <= 1) {
		printf("WARNING: not enough points for plan, fix planner!\n");
		return;
	}
	plan_and_execute(m, points, npoints);
}

int
read_gcode(motion_t *m, const char *file)
{
	FILE *fp;
	char buf[MAXLINE];
	char *p;
	int ret;
	gccmd_t gc;
	int i;

	int	extruder_mode_abs = 0;
	double	X = 0;
	double	Y = 0;
	double	Z = 0;
	double	E = 0;
	double	F = 0;
	point_t	*points = NULL;
	int	npoints = 0;
	int	alloc_points = 0;

	fp = fopen(file, "r");
	if (fp == NULL) {
		printf("read_gcode: failed to open %s\n", file);
		return -1;
	}

	while (1) {
		double v;
		int vi;

		p = fgets(buf, sizeof(buf), fp);
		if (p == NULL) {
			printf("reading gcode finished\n");
			break;
		}
		if (strncmp(buf, ";LAYER:1", 8) == 0) {
			printf("layer1 finished\n");
			break;
		}
		ret = parse_line(buf, &gc);
		if (ret < 0) {
			printf("parsing line failed\n");
			goto err;
		}
		dump_gccmd(&gc);
		if (gc.codes[0] == 0)
			continue;
		v = gc.values[gc.codes[0] - 'A'];
		vi = (int)v;
		if (v != vi) {
			printf("non-integer codes not supported: %c%d\n",
				gc.codes[0], vi);
			goto err;
		}
		if (gc.codes[0] == 'M') {
			switch(vi) {
			case 82:
				ret = check_params(&gc, "", "");
				if (ret)
					goto err;
				printf("setting extruder to abs\n");
				extruder_mode_abs = 1;
				break;
			case 84:
				ret = check_params(&gc, "", "");
				if (ret)
					goto err;
				printf("stopping idle holds, power down "
					"steppers\n");
				break;
			case 104:
				ret = check_params(&gc, "S", "");
				if (ret)
					goto err;
				printf("setting temp to %f\n",
					gc.values['S'-'A']);
				break;
			case 105:
				ret = check_params(&gc, "", "");
				if (ret)
					goto err;
				printf("sending extruder temp\n");
				break;
			case 106:
				ret = check_params(&gc, "", "SP");
				if (ret)
					goto err;
				printf("turning fan on\n");
				break;
			case 107:
				ret = check_params(&gc, "", "");
				if (ret)
					goto err;
				printf("turning fan off\n");
				break;
			case 109:
				ret = check_params(&gc, "S", "");
				if (ret)
					goto err;
				printf("setting temp to %f and wait\n",
					gc.values['S'-'A']);
				break;
			case 140:
				ret = check_params(&gc, "S", "H");
				if (ret)
					goto err;
				printf("setting bed temp to %f\n",
					gc.values['S'-'A']);
				break;
			default:
				printf("unknown command M%d\n", vi);
				goto err;
			}
		} else if (gc.codes[0] == 'G') {
			switch(vi) {
			case 0:
				ret = check_params(&gc, "", "XYZEF");
				if (ret)
					goto err;
				set_if_set(&gc, 'X', &X);
				set_if_set(&gc, 'Y', &Y);
				set_if_set(&gc, 'Z', &Z);
				set_if_set(&gc, 'E', &E);
				set_if_set(&gc, 'F', &F);
				printf("rapid move\n");
				if (npoints)
					push_points(m, points, npoints);
				npoints = 0;
				/*
				 * ignore it
				 */
				break;
			case 1:
				ret = check_params(&gc, "", "XYZEF");
				if (ret)
					goto err;
				set_if_set(&gc, 'X', &X);
				set_if_set(&gc, 'Y', &Y);
				set_if_set(&gc, 'Z', &Z);
				set_if_set(&gc, 'E', &E);
				set_if_set(&gc, 'F', &F);
				if (strchr(gc.codes, 'Z') &&
				    strpbrk(gc.codes, "XYE")) {
					printf("linear move in Z combined with "
						"X/Y/E not supported\n");
					exit(1);
				}
				printf("linear move to %f/%f/%f/%f at %f\n",
					X, Y, Z, E, F);
				if (strpbrk(gc.codes, "XY") == NULL)
					break;
				if (npoints == alloc_points) {
					points = realloc(points,
						sizeof(*points) *
						 (alloc_points + P_CHUNK));
					for (i = 0; i < P_CHUNK; ++i) {
						mpfr_inits(
						   points[i + alloc_points].x,
						   points[i + alloc_points].y,
						   NULL);
					}
					alloc_points += P_CHUNK;
				}
				mpfr_set_d(points[npoints].x, X, rnd);
				mpfr_set_d(points[npoints].y, Y, rnd);
				++npoints;
				break;
			case 28:
				ret = check_params(&gc, "", "XYZ");
				if (ret)
					goto err;
				printf("homing\n");
				break;
			case 92:
				ret = check_params(&gc, "", "XYZE");
				if (ret)
					goto err;
				printf("set position\n");
				break;
			default:
				printf("unknown command G%d\n", vi);
				goto err;
			}
		} else {
			printf("unknown command %c\n", gc.codes[0]);
			goto err;
		}
	}
	if (npoints)
		push_points(m, points, npoints);

	/* XXX TODO free points */
	return 0;

err:
	printf("%d", extruder_mode_abs);/* keep compiler happy */
	return -1;
}

