#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <cblas.h>
#include <lapacke.h>

#include "conan.h"

typedef struct _p {
	double	x;
	double	y;
	double	t;
} p_t;

void
match_path(motion_t *m, point_t *points, int npoints)
{
	int i;
	p_t p[npoints];
	double cl;
#if 0
	double epsilon = 0.005;	/* noise level: data is approx good to 0.001 */
#endif

if (npoints <= 2) return;

	/* copy to local array */
	for (i = 0; i < npoints; ++i) {
		p[i].x = mpfr_get_d(points[i].x, rnd);
		p[i].y = mpfr_get_d(points[i].y, rnd);
	}

	/* accumulate chord lengths */
	cl = 0;
	for (i = 0; i < npoints - 1; ++i) {
		double dx = p[i + 1].x - p[i].x;
		double dy = p[i + 1].y - p[i].y;
		double l = sqrt(dx * dx + dy * dy);
printf("l %.5f\n", l);
		p[i].t = cl;
		cl += l;
	}

	/* parametrize in chord length */
	for (i = 0; i < npoints - 1; ++i)
		p[i].t /= cl;
	p[i].t = 1;

	for (i = 0; i < npoints; ++i)
		printf("ba: %.4f %.3f %.3f\n", p[i].t, p[i].x, p[i].y);
	printf("ba: done\n");

#if 0
	serial_bisect(p, npoints);
#endif

exit(1);
}

#if 0

serial_bisect(p_t p, int np, int degree, double epsilon)
{
	int StartIdx = 1;
	int EndIdx = n;
	int LeftIdx;
	int RightIdx;
	int Success;

	while (EndIdx > StartIdx) {
		LeftIdx = StartIdx;
		RightIdx = EndIdx;
		while (RightIdx - LeftIdx) <= 1) {
			if (StartIdx + degree) >= n
				break;
			LocalBsplineFitting(&Emax);
			Success = 0;
			if (Emax < epsilon) {
				Success = 1;
				Save temporary local spline
				Save temporary
			}
			if  (Success == 1)
				LeftIdx = EndIdx;
			else
				RightIdx = EndIdx;
			EndIdx = floor((LeftIdx + RightIdx) / 2);
		}
		StartIdx = LeftIdx + 1;
		EndIdx = n;
		Save local spline
		Save knot
	}
	endwhile
}

LocalBsplineFitting()
{
    % Computing b-spline basis function
     ;
    ;
    % identify control points for the last segment by solving equation 3.3

     % Computing fitting error

     ;   % element-wise multiplication
    ;
}
#endif

void
print_matrix(const char *name, int rows, int cols, double *m, int lda)
{
	int i;
	int j;

	printf("%s: (%dx%d)\n", name, rows, cols);
	for (i = 0; i < rows; ++i) {
		for (j = 0; j < cols; ++j) {
			printf("%.05f ", m[i * lda + j]);
		}
		printf("\n");
	}
	printf("----\n");
}

void NmatCal_1P(double Knotin[2], int Degree, double *Nmat);
/*
 * file name: SerialBisection.m
 * Description: This function employes parallel method to divide input data
 * into small single-b-spline pieces having the maximum fitted error being
 * smaller than a certain input value
 * Prototype:
 * VectorUX = ParallelBisection(DataIn,Degree,CtrlError,NumOfPieceStart)
 * Input parameters:
 * - DataIn: A matrix contains input data, having at least 2 columns. 1st
 * column is parametric, 2nd column is X, 3rd column is Y and so on.
 * - Degree: Degree of the fitted B-spline
 * - CtrlError: Maximum fitted error of a single piece b-spline
 * Output parameters:
 * - VectorUX: a matrix contains information of each single piece b-spline.
 * each column stores the information of each single piece b-spline.
 * 1st row: start index, 2nd: end index, 3 to 2+Order^2: Coefficient vector,
 * 3 + Order^2 to end - 1: B-spline control points, end: maximum fitted
 * error.
 * Version: 1.0
 * Date: 30-June-2016
 * Author: Dvthan
 */
//function VectorUX = SerialBisection(DataIn,Degree,CtrlError)

typedef struct _VectorUX {
	int	StartIdx;
	int	LeftIdx;
	double	*Nmatrix;
	double	*CtrlPointY;
	double	ErrorMax;
} VectorUX_t;

void
SerialBisection(double *DataIn, int rows, int cols, int Degree, double CtrlError,
	VectorUX_t **VectorUX, int *VectorUXlen)
{
	int i;
	int j;
	int ret;

// %% Public variables
// % prepare Amat
	/* Public variables */
	/* Amat */

// dataT = DataIn(:,1);
	double dataT[rows];
	for (i = 0; i < rows; ++i)
		dataT[i] = DataIn[i * cols];

// [datasize, S] = size(DataIn);
	int datasize = rows;
	int S = cols;

// % Amat = [1,x1,...,x1^N ;1,x2,...,x2^N ;...];
// Amat = zeros(datasize,Degree+1);
	/* Amat = [1,x1,...,x1^N ;1,x2,...,x2^N ;...]; */
// Amat = zeros(datasize,Degree+1);
	double Amat[datasize][Degree + 1];
	for (i = 0; i < datasize; ++i)
		for (j = 0; j < Degree + 1; ++j)
			Amat[i][j] = 0;

// Amat(:,1) = 1;
	for (i = 0; i < datasize; ++i)
		Amat[i][0] = 1;

// for ii = 1:(Degree)
//     Amat(:,1+ii) = Amat(:,ii).*dataT;
// end
	for (i = 1; i <= Degree; ++i)
		for (j = 0; j < datasize; ++j)
			Amat[j][i] = Amat[j][i - 1] * dataT[j];

// Ymat = DataIn(:,2:end);
	double Ymat[datasize][S - 1];
	for (i = 0; i < datasize; ++i)
		for (j = 0; j < S - 1; ++j)
			Ymat[i][j] = DataIn[i * cols + j + 1];

//
// %% Bisecting
// StartIdx = 1;
// EndIdx = datasize;
// ptr = 1;
	/* Bisecting */
	int StartIdx = 0;
	int EndIdx = datasize - 1;
	*VectorUXlen = 0;

// % one piece B-spline fitting
// % computing Nmatrix
// Order = Degree + 1;
	/* one piece B-spline fitting */
	/* computing Nmatrix */
	int Order = Degree + 1;

// Order2= Order*Order;
	int Order2 = Order * Order;

#if 0
// VectorUX1 = zeros(3+Order*(S-1)+Order*Order, fix(datasize/Order));
	int ux1 = 3 + Order * (S - 1) + Order2;
	int ux2 = datasize / Order;
	double VectorUX1[ux1][ux2];
	for (i = 0; i < ux1; ++i)
		for (j = 0; j < ux2; ++j)
			VectorUX1[i][j] = 0;
#endif

	double CtrlPointYsave[Order][S - 1];
	double Nmatrixsave[Order2];
	double ErrorMaxsave;

// % Calculation of Ni,0 format [a0,a1,...,an-1,an]
//
// while EndIdx > StartIdx
	/* Calculation of Ni,0 format [a0,a1,...,an-1,an] */
	while (EndIdx > StartIdx) {

//
//     LeftIdx = StartIdx;
//     RightIdx = EndIdx;

		int LeftIdx = StartIdx;
		int RightIdx = EndIdx;

//     Knotin(1)= dataT(StartIdx);

		double Knotin[2];
		Knotin[0] = dataT[StartIdx];

//     while(1)

		while(1) {

//         if (StartIdx + Degree)>=datasize
//             CtrlPointYsave = zeros(Order,(S-1));
//             Nmatrixsave = zeros(Order2,1);
//             ErrorMaxsave  = 0;
//             break;
//         end

			if ((StartIdx + Degree) >= datasize) {
				for (i = 0; i < Order; ++i)
					for (j = 0; j < S - 1; ++j)
						CtrlPointYsave[i][j] = 0;
				for (i = 0; i < Order2; ++i)
					Nmatrixsave[i] = 0;
				ErrorMaxsave = 0;
				break;
			}

			Knotin[1]= dataT[EndIdx];
//         % Basis function calculation
//         SPNmat = NmatCal_1P(Knotin,Degree);

			/* Basis function calculation */
			double SPNmat[Order2];
			NmatCal_1P(Knotin, Degree, SPNmat);

//        %Spline fitting
//         clear Nmatrix  Ymatrix
//         Nmatrix = Amat(StartIdx:EndIdx,:)*reshape(SPNmat,[Order,Order]);
//         Ymatrix = Ymat(StartIdx:EndIdx,:);
//         CtrlP = Nmatrix\Ymatrix;

			/* Spline fitting */
			int nm1 = EndIdx - StartIdx + 1;
			double Nmatrix[nm1][Order];
			double r_SPNmat[Order][Order];
			/* reshape */
			/* XXX TODO reshape kann wahrscheinlich auch direkt
			 * ueber den aufruf an dgemm gemacht werden */
			for (i = 0; i < Order; ++i)
				for (j = 0; j < Order; ++j)
					r_SPNmat[i][j] = SPNmat[i + j * Order];

print_matrix("Amat", nm1, Order, Amat[StartIdx], Order);
print_matrix("r_SPNmat", Order, Order, *r_SPNmat, Order);

//	double Amat[datasize][Degree + 1];
//	double r_SPNmat[Order][Order];
//	double Nmatrix[nm1][Order];

			/* mult */
			cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
				nm1, Order, Order, 1.0,
				Amat[StartIdx], Order, *r_SPNmat, Order,
				0, *Nmatrix, Order);

#if 0
void cblas_dgemm (const CBLAS_LAYOUT Layout, const CBLAS_TRANSPOSE transa, const CBLAS_TRANSPOSE transb, const MKL_INT m, const MKL_INT n, const MKL_INT k, const double alpha, const double *a, const MKL_INT lda, const double *b, const MKL_INT ldb, const double beta, double *c, const MKL_INT ldc);
#endif

print_matrix("Nmatrix", nm1, Order, *Nmatrix, Order);
print_matrix("Ymatrix", nm1, S - 1, Ymat[StartIdx], S - 1);

//         Ymatrix = Ymat(StartIdx:EndIdx,:);
//         CtrlP = Nmatrix\Ymatrix;

			double Ymatrix[nm1][S - 1];
			for (i = 0; i < nm1; ++i)
				for (j = 0; j < S - 1; ++j)
					Ymatrix[i][j] = Ymat[StartIdx + i][j];

			/* Nmatrix gets overwritten, pass in a copy */
			double Nm_cpy[nm1][Order];
			for (i = 0; i < nm1; ++i)
				for (j = 0; j < Order; ++j)
					Nm_cpy[i][j] = Nmatrix[i][j];

print_matrix("Ymatrix", nm1, S - 1, *Ymatrix, S - 1);
			ret = LAPACKE_dgels(LAPACK_ROW_MAJOR, 'N',
				nm1, Order, S - 1, *Nm_cpy, Order, *Ymatrix, S - 1);

			double CtrlP[Order][S - 1];
			for (i = 0; i < Order; ++i)
				for (j = 0; j < S - 1; ++j)
					CtrlP[i][j] = Ymatrix[i][j];

			printf("dgels returned %d\n", ret);
print_matrix("Nmatrix", nm1, Order, *Nmatrix, Order);
print_matrix("CtrlP", Order, S - 1, *CtrlP, S - 1);

//         R = Ymatrix-Nmatrix*CtrlP;
//         R = R.*R;
			/*
			 * calculate residual
			 */
			double R[nm1][S - 1];
			cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
				nm1, S - 1, Order, 1.0,
				*Nmatrix, Order, *CtrlP, S - 1,
				0, *R, S - 1);

print_matrix("R (pre, first 10)", 10, S - 1, *R, S - 1);
			for (i = 0; i < nm1; ++i) {
				for (j = 0; j < S - 1; ++j) {
					R[i][j] = Ymat[i+StartIdx][j] - R[i][j];
					R[i][j] *= R[i][j];
				}
			}
print_matrix("R (post)", nm1, S - 1, *R, S - 1);
//         Errorcal = max(sum(R,2));
//         Errorcal = sqrt(Errorcal);
			double Errorcal = 0;
			for (i = 0; i < nm1; ++i) {
				double sum = 0;
				for (j = 0; j < S - 1; ++j)
					sum += R[i][j];
				if (sum > Errorcal)
					Errorcal = sum;
			}
			Errorcal = sqrt(Errorcal);

printf("Errorcal: %.5f\n", Errorcal);
//
//         % Fitting error evaluation
			/* Fitting error evaluation */
//
//         success  = 0;
//         if Errorcal <= CtrlError
//             success = 1;
//             CtrlPointYsave = CtrlP;
//             Nmatrixsave = SPNmat;
//             ErrorMaxsave  = Errorcal;
//         end
//         if success ==1
//             LeftIdx = EndIdx;
//         else
//             RightIdx = EndIdx;
//         end
//         if RightIdx - LeftIdx <=1
//             break;
//         end
//         EndIdx = floor((LeftIdx+RightIdx)/2);
			int success = 0;
printf("round: left %d right %d error %f\n", LeftIdx, RightIdx, Errorcal);
			if (Errorcal <= CtrlError) {
printf("saving Errorcal %f\n", Errorcal);
				success = 1;
				for (i = 0; i < Order; ++i)
					for (j = 0; j < S - 1; ++j)
						CtrlPointYsave[i][j] = CtrlP[i][j];
				for (i = 0; i < Order2; ++i)
					Nmatrixsave[i] = SPNmat[i];
				ErrorMaxsave = Errorcal;
			}
			if (success == 1)
				LeftIdx = EndIdx;
			else
				RightIdx = EndIdx;
			if (RightIdx - LeftIdx <= 1)
				break;
			EndIdx = (LeftIdx + RightIdx) / 2;
printf("new EndIdx(+1): %d\n", EndIdx + 1);
		}
printf("ErrorMaxSave: %f\n", ErrorMaxsave);
//     % fill vector Ux
			/* fill vector Ux */
//     VectorUX1(1, ptr) = StartIdx;
//     VectorUX1(2, ptr) = LeftIdx;
//     VectorUX1(3:(2+Order*Order), ptr) = Nmatrixsave';
//     VectorUX1(3+Order2:end-1,ptr) = reshape(CtrlPointYsave,[Order*(S-1),1]);
//     VectorUX1(end, ptr) = ErrorMaxsave;
//     ptr = ptr + 1;
//     StartIdx = LeftIdx+1;
//     EndIdx = datasize;
		VectorUX_t *v = calloc(sizeof(**VectorUX), 1);
		assert(v);
		v->StartIdx = StartIdx;
		v->LeftIdx = LeftIdx;
		v->Nmatrix = malloc(sizeof(double) * Order2);
		memcpy(v->Nmatrix, Nmatrixsave, sizeof(double) * Order2);
		v->CtrlPointY = malloc(sizeof(double) * Order * (S - 1));
		memcpy(v->CtrlPointY, CtrlPointYsave, sizeof(double) * Order * (S - 1));
		v->ErrorMax = ErrorMaxsave;
		VectorUX[(*VectorUXlen)++] = v;
		StartIdx = LeftIdx + 1;
		EndIdx = datasize - 1;
	}
// end
// VectorUX = VectorUX1(:,1:(ptr-1));
// end
}

/*
 * return result in Nmat of size (Degree + 1) ^ 2
 */
void
NmatCal_1P(double Knotin[2], int Degree, double *Nmat)
{
	int i;
	int jj;
	int kk;

printf("NmatCal_1P: knotin %f/%f, degree %d\n", Knotin[0], Knotin[1], Degree);
	int Order = Degree + 1;
	int rownumbers = Order * (Order + 1) * (Order / 2.0);

	double CoeNmat[rownumbers];
	for (i = 0; i < rownumbers; ++i)
		CoeNmat[i] = 0;

	CoeNmat[0] = 1;

	double Knot[2 * Order];

	for (i = 0; i < Order; ++i)
		Knot[i] = Knotin[0];

	for (i = Order; i < 2 * Order; ++i)
		Knot[i] = Knotin[1];

	for (jj = 1; jj <= Degree; ++jj) {
		int id0 = Order * (jj - 1) / 2.0 * jj;
		int id1 = Order * jj / 2.0 * (jj + 1);

		for (kk = 0; kk <= jj; ++kk) {
			/*
			 * effective internal we focus,eg when ii=4,jj=1,then
			 * id 2=3,4
			 */
			int id2 = Order - jj + kk - 1;

			/* effective knot num 1 */
			int id2Knot00 = id2 + jj;
			int id2Knot01 = id2Knot00 + 1;

			if ((id2 >= 0) && (id2Knot01 < 2 * Order)) {
				/* Access previous data Ni-1,j-1 Ni,j-1 and
				 * Ni+1,j-1
				 */
				int id00 = id0 + (kk - 1) * Order;
				int id01 = id0 + kk * Order;

				double N0[Order];
				double N1[Order];

				if (kk == 0) {	/* first box of matrix */
					for (i = 0; i < Order; ++i) {
						N0[i] = 0;
						N1[i] = CoeNmat[id01 + i];
					}
				} else if (kk == jj) {
					for (i = 0; i < Order; ++i) {
						N0[i] = CoeNmat[id00 + i];
						N1[i] = 0;
					}
				} else {
					for (i = 0; i < Order; ++i) {
						N0[i] = CoeNmat[id00 + i];
						N1[i] = CoeNmat[id01 + i];
					}
				}

				/* calculate a1x+a0, */
				double aden = Knot[id2Knot00] - Knot[id2];
				double bden = Knot[id2Knot01] - Knot[id2 + 1];

				double a0 = 0;
				double a1 = 0;
				if (aden != 0) {
					a1 = 1 / aden;
					a0 = -Knot[id2] / aden;
				}

				/* calculate b1x+b0 */
				double b0 = 0;
				double b1 = 0;
				if (bden != 0) {
					b1 = -1 / bden;
					b0 = Knot[id2Knot01] / bden;
				}

				/* Multiplication, */
				double Acoef[Order];
				double N00[Order];
				double N01[Order];
				double N10[Order];
				double N11[Order];
				for (i = 0; i < Order; ++i) {
					Acoef[i] = 0;
					N00[i] = a0 * N0[i];
					N01[i] = a1 * N0[i];
					N10[i] = b0 * N1[i];
					N11[i] = b1 * N1[i];
				}
				Acoef[0] = N00[0] + N10[0];
				for (i = 1; i <= Order; ++i)
					Acoef[i] = N00[i] + N10[i] +
					           N01[i - 1] + N11[i - 1];
				int id11 = id1 + kk * Order;
				for (i = 0; i <= Degree; ++i)
					CoeNmat[id11 + i] = Acoef[i];
			}
		}
	}

	int id10 = Order * (Degree / 2.0) * (Degree + 1);

	for (i = id10; i < rownumbers; ++i)
		Nmat[i - id10] = CoeNmat[i];

printf("id10: %d\n", id10);
	for (i = id10; i < rownumbers; ++i)
		printf("Nmat[%d] = %0.5f\n", i - id10, CoeNmat[i]);
}

#ifdef BATEST
mpfr_rnd_t rnd = MPFR_RNDN;

double demo[94][3] = {
	{ 0.0000, 94.550, 112.880 },
	{ 0.0049, 95.116, 112.436 },
	{ 0.0083, 95.508, 112.157 },
	{ 0.0137, 96.174, 111.732 },
	{ 0.0191, 96.879, 111.377 },
	{ 0.0246, 97.626, 111.092 },
	{ 0.0283, 98.135, 110.927 },
	{ 0.0337, 98.897, 110.722 },
	{ 0.0418, 100.063, 110.562 },
	{ 0.2077, 124.151, 109.086 },
	{ 0.2129, 124.916, 109.076 },
	{ 0.2264, 126.867, 109.148 },
	{ 0.2331, 127.845, 109.246 },
	{ 0.2456, 129.644, 109.542 },
	{ 0.2507, 130.361, 109.695 },
	{ 0.2573, 131.289, 109.939 },
	{ 0.2621, 131.963, 110.149 },
	{ 0.2742, 133.611, 110.746 },
	{ 0.2788, 134.238, 111.005 },
	{ 0.2854, 135.102, 111.406 },
	{ 0.2903, 135.739, 111.738 },
	{ 0.3020, 137.206, 112.591 },
	{ 0.3096, 138.118, 113.212 },
	{ 0.3273, 140.137, 114.806 },
	{ 0.3327, 140.732, 115.326 },
	{ 0.3355, 141.020, 115.621 },
	{ 0.3493, 142.381, 117.087 },
	{ 0.3555, 142.955, 117.781 },
	{ 0.3656, 143.824, 118.966 },
	{ 0.3723, 144.353, 119.793 },
	{ 0.3808, 144.951, 120.870 },
	{ 0.3862, 145.299, 121.579 },
	{ 0.3916, 145.574, 122.319 },
	{ 0.3971, 145.772, 123.083 },
	{ 0.4016, 145.877, 123.730 },
	{ 0.4043, 145.924, 124.123 },
	{ 0.4097, 145.978, 124.911 },
	{ 0.4151, 145.952, 125.700 },
	{ 0.4202, 145.856, 126.427 },
	{ 0.4226, 145.793, 126.778 },
	{ 0.4281, 145.615, 127.547 },
	{ 0.4335, 145.359, 128.294 },
	{ 0.4413, 144.864, 129.313 },
	{ 0.4504, 144.196, 130.465 },
	{ 0.4574, 143.630, 131.308 },
	{ 0.4697, 142.536, 132.730 },
	{ 0.4772, 141.814, 133.545 },
	{ 0.4941, 140.058, 135.266 },
	{ 0.4996, 139.466, 135.789 },
	{ 0.5029, 139.079, 136.079 },
	{ 0.5192, 137.137, 137.442 },
	{ 0.5246, 136.468, 137.863 },
	{ 0.5277, 136.069, 138.072 },
	{ 0.5445, 133.872, 139.141 },
	{ 0.5499, 133.145, 139.450 },
	{ 0.5534, 132.668, 139.608 },
	{ 0.5702, 130.325, 140.300 },
	{ 0.5756, 129.558, 140.485 },
	{ 0.5791, 129.057, 140.562 },
	{ 0.5965, 126.539, 140.868 },
	{ 0.6035, 125.534, 140.925 },
	{ 0.6126, 124.203, 140.914 },
	{ 0.6154, 123.790, 140.899 },
	{ 0.7782, 100.149, 139.450 },
	{ 0.7837, 99.364, 139.362 },
	{ 0.7891, 98.592, 139.195 },
	{ 0.7942, 97.888, 138.969 },
	{ 0.8000, 97.093, 138.671 },
	{ 0.8055, 96.369, 138.356 },
	{ 0.8109, 95.680, 137.971 },
	{ 0.8163, 95.033, 137.517 },
	{ 0.8217, 94.436, 137.000 },
	{ 0.8252, 94.084, 136.641 },
	{ 0.8301, 93.604, 136.118 },
	{ 0.8355, 93.101, 135.510 },
	{ 0.8409, 92.661, 134.854 },
	{ 0.8464, 92.290, 134.157 },
	{ 0.8511, 92.024, 133.518 },
	{ 0.8554, 91.812, 132.937 },
	{ 0.8608, 91.579, 132.182 },
	{ 0.8662, 91.423, 131.408 },
	{ 0.8739, 91.339, 130.291 },
	{ 0.9437, 91.302, 120.141 },
	{ 0.9516, 91.382, 119.000 },
	{ 0.9553, 91.460, 118.459 },
	{ 0.9607, 91.612, 117.684 },
	{ 0.9688, 91.978, 116.577 },
	{ 0.9718, 92.146, 116.173 },
	{ 0.9772, 92.486, 115.460 },
	{ 0.9826, 92.896, 114.785 },
	{ 0.9859, 93.173, 114.402 },
	{ 0.9891, 93.464, 114.025 },
	{ 0.9946, 93.978, 113.425 },
	{ 1.0000, 94.550, 112.880 }
};

void
nmat_test()
{
	double Knotin[2];
	int Degree = 3;
	double Nmat[(Degree + 1) * (Degree + 1)];

	Knotin[0] = 0;
	Knotin[1] = 1;

	NmatCal_1P(Knotin, Degree, Nmat);
}

int
main(int argc, char **argv)
{
	int len = 0;
	VectorUX_t *VectorUX[94];
	int i;

	SerialBisection(*demo, 94, 3, 3, 0.05, VectorUX, &len);

	for (i = 0; i < len; ++i) {
		printf("[%d]: %d-%d ErrorMax %f\n", i,
			VectorUX[i]->StartIdx,
			VectorUX[i]->LeftIdx,
			VectorUX[i]->ErrorMax);
	}

	return 0;
}
#endif
