#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <sys/stat.h>
#include <sys/types.h>

#include "conan.h"

#define MAXLINE	1000
#define MAXCODES	26
#define P_CHUNK	10

typedef struct _gccmd {
	char	codes[MAXCODES];
	double	values[MAXCODES];
} gccmd_t;

static int prepnum = 0;
FILE *gp;
FILE *html;
static void
init_prep(void)
{
	prepnum = 0;
	mkdir("gc", 0777);
	gp = fopen("gc/data.gp", "w");
	assert(gp);
	html = fopen("gc/index.html", "w");
	assert(html);
	fprintf(html, "<html>\n<head>\n<title>3d printing from scratch</title>\n");
	fprintf(html, "</head>\n<body>\n<table>\n<tr>\n");
}

static void
fini_prep(void)
{
	fprintf(html, "</tr>\n</body>\n</html>\n");
	fclose(html);
	fclose(gp);
}

static double
check_middle_deviation(point_t *p1, point_t *p2, point_t *p3)
{
	mpfr_t d, dx21, dy21, dx10, dy10, tmp;
	double rv;

	mpfr_inits(d, dx21, dy21, dx10, dy10, tmp, NULL);

	mpfr_sub(dx21, p3->x, p2->x, rnd);
	mpfr_sub(dy21, p3->y, p2->y, rnd);
	mpfr_sub(dx10, p2->x, p1->x, rnd);
	mpfr_sub(dy10, p2->y, p1->y, rnd);

	mpfr_mul(tmp, dx10, dy21, rnd);
	mpfr_fms(tmp, dx21, dy10, tmp, rnd);
	mpfr_abs(d, tmp, rnd);
	mpfr_hypot(tmp, dx21, dy21, rnd);
	mpfr_div(d, d, tmp, rnd);
	rv = mpfr_get_d(d, rnd);

	mpfr_printf("check_middle_deviation: p1 (%.3Rf/%.3Rf) p2 (%.3Rf/%.3Rf) p3 (%.3Rf/%.3Rf) d %.5Rf\n", p1->x, p1->y, p2->x, p2->y, p3->x, p3->y, d);

	mpfr_clears(d, dx21, dy21, dx10, dy10, tmp, NULL);

	return rv;
}

static void
remove_redundant_points(motion_t *m, point_t *points, int *npoints)
{
	int i;
	int ri;
	double min_dev = 0.01;
	double d;
	double d2;

	for (i = 0; i < *npoints - 2; ++i) {
		printf("check point %d/%d\n", i, *npoints);
		d = check_middle_deviation(points + i, points + i+1, points + i+2);
		if (d > min_dev)
			continue;
		ri = i + 1;
		if (i < *npoints - 3) {
			d2 = check_middle_deviation(points + i+1, points + i+2,
				points + i+3);
			if (d2 < d)
				ri = i + 2;
		}
		printf("remove point %d\n", ri);
		mpfr_clears(points[ri].x, points[ri].y, NULL);
		memmove(points + ri, points + ri + 1,
			sizeof(*points) * (*npoints - ri - 1));
		mpfr_inits(points[*npoints-1].x, points[*npoints-1].y, NULL);
		--*npoints;
		--i; /* check from current point again */
	}
}

static void
output_path(motion_t *m, point_t *points, int npoints)
{
	int i;
	mpfr_t dx, dy, tmp, t;
	char fn[50];
	char csv[50];
	FILE *fp;

	sprintf(fn, "path%04d", prepnum++);
	sprintf(csv, "gc/%s.csv", fn);

	fp = fopen(csv, "w");
	assert(fp);

	mpfr_inits(dx, dy, tmp, t, NULL);

	mpfr_set_ui(t, 0, rnd);
	mpfr_set_ui(tmp, 0, rnd);

	for (i = 0; i < npoints; ++i) {
		if (i > 0) {
			mpfr_sub(dx, points[i-1].x, points[i].x, rnd);
			mpfr_sub(dy, points[i-1].y, points[i].y, rnd);
			mpfr_hypot(tmp, dx, dy, rnd);
			mpfr_add(t, t, tmp, rnd);
		}
		mpfr_fprintf(fp, "%.7Rf %.5Rf %.5Rf %.5Rf\n",
			t, points[i].x, points[i].y, tmp);
	}
	fclose(fp);
	mpfr_clears(dx, dy, tmp, NULL);

	fprintf(gp, "set term png size 1000,1000\n");
	fprintf(gp, "set output \"gc/%s-path.png\"\n", fn);
	fprintf(gp, "plot \"%s\" using 2:3 title \"Path\" with lines lt rgb \"red\", \\\n", csv);
	fprintf(gp, "\"%s\" using 2:3 title \"Path\" with points lt rgb \"blue\"\n", csv);
	fprintf(gp, "set output \"gc/%s-x.png\"\n", fn);
	fprintf(gp, "plot \"%s\" using 1:2 title \"X\" with lines lt rgb \"red\", \\\n", csv);
	fprintf(gp, "\"%s\" using 1:2 title \"X\" with points lt rgb \"blue\"\n", csv);
	fprintf(gp, "set output \"gc/%s-y.png\"\n", fn);
	fprintf(gp, "plot \"%s\" using 1:3 title \"Y\" with lines lt rgb \"red\", \\\n", csv);
	fprintf(gp, "\"%s\" using 1:3 title \"Y\" with points lt rgb \"blue\"\n", csv);
	fprintf(html, "<td>\n\t<table>\n\t\t<tr>\n");
	fprintf(html, "\t\t\t<img src=\"%s-path.png\" alt=\"Path\">\n", fn);
	fprintf(html, "\t\t</tr>\n");
	fprintf(html, "\t\t<tr>\n");
	fprintf(html, "\t\t\t<img src=\"%s-x.png\" alt=\"X\">\n", fn);
	fprintf(html, "\t\t</tr>\n");
	fprintf(html, "\t\t<tr>\n");
	fprintf(html, "\t\t\t<img src=\"%s-y.png\" alt=\"Y\">\n", fn);
	fprintf(html, "\t\t</tr>\n");
	fprintf(html, "\t</table>\n");
	fprintf(html, "</td>\n");
}

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

	printf("push_points\n");

	for (i = 0; i < npoints; ++i) {
		mpfr_printf("pp: [%d] %.5Rf %.5Rf\n", i, points[i].x, points[i].y);
	}
	if (npoints <= 1) {
		printf("WARNING: not enough points for plan, fix planner!\n");
		return;
	}

	remove_redundant_points(m, points, &npoints);
	output_path(m, points, npoints);
#if 0
	match_path(m, points, npoints);
#endif
#if 1
	plan_and_execute(m, points, npoints);
#endif
	printf("pp: -------\n");
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
	double	last_g0_X = 0;
	double	last_g0_Y = 0;
	int	first_z = 1;

	init_prep();

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
			case 83:
				ret = check_params(&gc, "", "");
				if (ret)
					goto err;
				printf("setting extruder to rel\n");
				extruder_mode_abs = 0;
				break;
			case 84:
				ret = check_params(&gc, "", "");
				if (ret)
					goto err;
				printf("stopping idle holds, power down "
					"steppers\n");
				break;
			case 104:
				ret = check_params(&gc, "S", "T");
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
				ret = check_params(&gc, "S", "T");
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
			case 190:
				ret = check_params(&gc, "", "S");
				if (ret)
					goto err;
				printf("wait for bed to reach temp\n");
				break;
			case 204:
				ret = check_params(&gc, "", "PNS");
				if (ret)
					goto err;
				printf("set default acceleration\n");
				break;
			case 205:
				ret = check_params(&gc, "", "BEJSTXYZ");
				if (ret)
					goto err;
				printf("set advanced settings\n");
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
				 * save last point as start for g1 series
				 */
				if (strpbrk(gc.codes, "XY") == NULL)
					break;
				last_g0_X = X;
				last_g0_Y = Y;
				break;
			case 1: {
				double old_e = E;
				int is_extruding = 0;

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
				if (strchr(gc.codes, 'Z')) {
					if (first_z) {
						first_z = 0;
						break;
					}
					goto done;
				}
				if (strpbrk(gc.codes, "XY") == NULL) {
					printf("pure extracting move\n");
					break;
				}
				printf("linear move to %f/%f/%f/%f at %f\n",
					X, Y, Z, E, F);
				if (strchr(gc.codes, 'E')) {
					if (extruder_mode_abs) {
						if (E > old_e)
							is_extruding = 1;
					} else {
						if (E > 0)
							is_extruding = 1;
					}
				}
				if (!is_extruding) {
					if (npoints)
						push_points(m, points, npoints);
					npoints = 0;
					last_g0_X = X;
					last_g0_Y = Y;
					break;
				}
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
				if (npoints == 0) {
					printf("add %.3f/%.3f as start point\n",
						last_g0_X, last_g0_Y);
					mpfr_set_d(points[0].x, last_g0_X, rnd);
					mpfr_set_d(points[0].y, last_g0_Y, rnd);
					++npoints;
				}
				printf("add %.3f/%.3f as point %d\n", X, Y, npoints);
				mpfr_set_d(points[npoints].x, X, rnd);
				mpfr_set_d(points[npoints].y, Y, rnd);
				++npoints;
				break;
			}
			case 21:
				ret = check_params(&gc, "", "");
				if (ret)
					goto err;
				printf("set units to mm\n");
				break;
			case 28:
				ret = check_params(&gc, "", "XYZW");
				if (ret)
					goto err;
				printf("homing\n");
				break;
			case 80:
				ret = check_params(&gc, "", "");
				if (ret)
					goto err;
				printf("stop cycle\n");
				break;
			case 90:
				ret = check_params(&gc, "", "XYZE");
				if (ret)
					goto err;
				printf("move tool\n");
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
		} else if (gc.codes[0] == 'T') {
			ret = check_params(&gc, "", "");
			if (ret)
				goto err;
			printf("set tool\n");
		} else {
			printf("unknown command %c\n", gc.codes[0]);
			goto err;
		}
	}
done:
	if (npoints)
		push_points(m, points, npoints);

	fini_prep();

	/* XXX TODO free points */
	return 0;

err:
	printf("%d", extruder_mode_abs);/* keep compiler happy */
	return -1;
}

