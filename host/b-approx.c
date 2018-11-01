#include <math.h>
#include <stdlib.h>

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

void
SerialBisection(double *DataIn, int rows, int cols, int Degree, double CtrlError)
{
	int i;
	int j;

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
	for (i = 1; i < datasize; ++i)
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
	int ptr = 0;

// % one piece B-spline fitting
// % computing Nmatrix
// Order = Degree + 1;
	/* one piece B-spline fitting */
	/* computing Nmatrix */
	int Order = Degree + 1;

// Order2= Order*Order;
	int Order2= Order * Order;

// VectorUX1 = zeros(3+Order*(S-1)+Order*Order, fix(datasize/Order));
	int ux1 = 3 + Order * (S - 1) + Order2;
	int ux2 = datasize / Order;
	double VectorUX1[ux1][ux2];
	for (i = 0; i < ux1; ++i)
		for (j = 0; j < ux2; ++j)
			VectorUX1[i][j] = 0;

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

		LeftIdx = StartIdx;
		RightIdx = EndIdx;

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

			if (StartIdx + Degree) > datasize) {
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
			double Nmatrix[EndIdx-StartIdx+1][
//         Nmatrix = Amat(StartIdx:EndIdx,:)*reshape(SPNmat,[Order,Order]);
//         Ymatrix = Ymat(StartIdx:EndIdx,:);
//         CtrlP = Nmatrix\Ymatrix;
//
#if 0
//         R = Ymatrix-Nmatrix*CtrlP;
//         R = R.*R;
//         Errorcal = max(sum(R,2));
//         Errorcal = sqrt(Errorcal);
//         % Fitting error evaluation
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
//     end
//     % fill vector Ux
//     VectorUX1(1, ptr) = StartIdx;
//     VectorUX1(2, ptr) = LeftIdx;
//     VectorUX1(3:(2+Order*Order), ptr) = Nmatrixsave';
//     VectorUX1(3+Order2:end-1,ptr) = reshape(CtrlPointYsave,[Order*(S-1),1]);
//     VectorUX1(end, ptr) = ErrorMaxsave;
//     ptr = ptr + 1;
//     StartIdx = LeftIdx+1;
//     EndIdx = datasize;
// end
// VectorUX = VectorUX1(:,1:(ptr-1));
// end
#endif
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
int
main(int argc, char **argv)
{
	double Knotin[2];
	int Degree = 3;
	double Nmat[(Degree + 1) * (Degree + 1)];

	Knotin[0] = 0;
	Knotin[1] = 1;

	NmatCal_1P(Knotin, Degree, Nmat);

	return 0;
}
#endif
