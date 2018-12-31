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
void NewNmatrix(double *Knot, int nknots, int Degree, double *Nmat);

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

	/* Public variables */
	/* Amat */

	double dataT[rows];
	for (i = 0; i < rows; ++i)
		dataT[i] = DataIn[i * cols];

	int datasize = rows;
	int S = cols;

	/* Amat = [1,x1,...,x1^N ;1,x2,...,x2^N ;...]; */
	double Amat[datasize][Degree + 1];
	for (i = 0; i < datasize; ++i)
		for (j = 0; j < Degree + 1; ++j)
			Amat[i][j] = 0;

	for (i = 0; i < datasize; ++i)
		Amat[i][0] = 1;

	for (i = 1; i <= Degree; ++i)
		for (j = 0; j < datasize; ++j)
			Amat[j][i] = Amat[j][i - 1] * dataT[j];

	double Ymat[datasize][S - 1];
	for (i = 0; i < datasize; ++i)
		for (j = 0; j < S - 1; ++j)
			Ymat[i][j] = DataIn[i * cols + j + 1];

	/* Bisecting */
	int StartIdx = 0;
	int EndIdx = datasize - 1;
	*VectorUXlen = 0;

	/* one piece B-spline fitting */
	/* computing Nmatrix */
	int Order = Degree + 1;

	int Order2 = Order * Order;

	double CtrlPointYsave[Order][S - 1];
	double Nmatrixsave[Order2];
	double ErrorMaxsave;

	/* Calculation of Ni,0 format [a0,a1,...,an-1,an] */
	while (EndIdx > StartIdx) {
		int LeftIdx = StartIdx;
		int RightIdx = EndIdx;

		double Knotin[2];
		Knotin[0] = dataT[StartIdx];

		while(1) {
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

			/* Basis function calculation */
			double SPNmat[Order2];
			NmatCal_1P(Knotin, Degree, SPNmat);

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

#if 0
print_matrix("Amat", nm1, Order, Amat[StartIdx], Order);
print_matrix("r_SPNmat", Order, Order, *r_SPNmat, Order);
#endif

			/* mult */
			cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
				nm1, Order, Order, 1.0,
				Amat[StartIdx], Order, *r_SPNmat, Order,
				0, *Nmatrix, Order);

#if 0
print_matrix("Nmatrix", nm1, Order, *Nmatrix, Order);
print_matrix("Ymatrix", nm1, S - 1, Ymat[StartIdx], S - 1);
#endif

			double Ymatrix[nm1][S - 1];
			for (i = 0; i < nm1; ++i)
				for (j = 0; j < S - 1; ++j)
					Ymatrix[i][j] = Ymat[StartIdx + i][j];

			/* Nmatrix gets overwritten, pass in a copy */
			double Nm_cpy[nm1][Order];
			for (i = 0; i < nm1; ++i)
				for (j = 0; j < Order; ++j)
					Nm_cpy[i][j] = Nmatrix[i][j];

#if 0
print_matrix("Ymatrix", nm1, S - 1, *Ymatrix, S - 1);
#endif
			ret = LAPACKE_dgels(LAPACK_ROW_MAJOR, 'N',
				nm1, Order, S - 1, *Nm_cpy, Order, *Ymatrix, S - 1);

			double CtrlP[Order][S - 1];
			for (i = 0; i < Order; ++i)
				for (j = 0; j < S - 1; ++j)
					CtrlP[i][j] = Ymatrix[i][j];

			printf("dgels returned %d\n", ret);
			assert(ret == 0);
#if 0
print_matrix("Nmatrix", nm1, Order, *Nmatrix, Order);
print_matrix("CtrlP", Order, S - 1, *CtrlP, S - 1);
#endif

			/*
			 * calculate residual
			 */
			double R[nm1][S - 1];
			cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
				nm1, S - 1, Order, 1.0,
				*Nmatrix, Order, *CtrlP, S - 1,
				0, *R, S - 1);

			for (i = 0; i < nm1; ++i) {
				for (j = 0; j < S - 1; ++j) {
					R[i][j] = Ymat[i+StartIdx][j] - R[i][j];
					R[i][j] *= R[i][j];
				}
			}
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
			/* Fitting error evaluation */
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
		}
		/* fill vector Ux */
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

#if 0

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% file name: TwoPieceOptimalKnotSolver.m
% Description: This function finds optimal knot and its continuity. The
% input data must fully define each single pieces
% Prototype:
% [OptimalKnotOut,MultipleOut]= TwoPieceOptimalKnotSolver(DataIn,Degree,...
%	   SearchRange,MinAngle,N0ScanForDiscontinuosCase,GaussNewtonLoopTime)
% Input parameters:
% - DataIn: A matrix contains input data, having at least 2 columns. 1st
% column is parametric, 2nd column is X, 3rd column is Y and so on.
% - Degree: Degree of the fitted B-spline
% - SearchRange: [a,b] is a range to find optimal knot
% - MinAngle: Minimum joining angle at knot that we can accept the finding
% optimal knot.
% - N0ScanForDiscontinuosCase: Number of evaluation times in calculating
% discontinuity case.
% - GaussNewtonLoopTime: Maximum number of loops in Gauss-Newton solving.
% Output parameters:
% - OptimalKnotOut: Optimum knot location
% - MultipleOut: Multiple knot at the optimal knot join. If there is no
% optimal knot found, this variable will return 0 (Eliminate the join)
% error.
% Version: 1.1
% Date: 19-July-2016
% Author: Dvthan
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [OptimalKnotOut,MultipleOut] = TwoPieceOptimalKnotSolver1(DataIn,DataKnot,...
	Degree,MinAngle,MaxSmooth,N0ScanForDiscontinuosCase,GaussNewtonLoopTime,ScanKnot)

	Order = Degree + 1;
	[m,n] = size(DataIn);

	T1 = DataIn(:,1);
	T = zeros(numel(T1),Order);
	T(:,1) = 1;
	for cnt = 2:Order
		T(:,cnt) = T(:,cnt-1).*T1;
	end
	Ymat = DataIn(:,2:end);

	idx01 = DataKnot;
	idx10 = DataKnot + 1;
	idx00 = 1;
	idx11 = m;
	DP1 = idx01 - idx00 - Degree;
	DP2 = idx11 - idx10 - Degree;
	LeftSearch = zeros(1,Order);
	RightSearch = zeros(1,Order);
	for cnt = 1:Order
		LeftSearch(1,cnt) = min(DP1,3)+Degree-cnt;
		RightSearch(1,cnt) = min(DP2,3)+Degree-cnt;
	end
	% if (DP1>0)&&(DP2>0)
	%	 LeftSearch(1,Order) = LeftSearch(1,Order) + 1;
	%	 RightSearch(1,Order) = RightSearch(1,Order) + 1;
	% elseif (DP1<=0)&&(DP2>0)
	%	 RightSearch(1,Order) = RightSearch(1,Order) + 1;
	% elseif (DP2<=0)&&(DP1>0)
	%	 LeftSearch(1,Order) = LeftSearch(1,Order) + 1;
	% end
	%% Uniform scanning the lowest error
	OptimalKnotSave = zeros(1,Order);
	ErrorSave = zeros(1,Order);
	AngleSave = zeros(1,Order);
	SearchRangeOut = zeros(2,Order);
	MultipleMax = Order;
	for ii = Order:-1:1
	%	 LeftSearch1 = max(LeftSearch(ii),0);
	%	 RightSearch1 = max(RightSearch(ii),0);
		LeftSearch1 = LeftSearch(ii);
		RightSearch1 = RightSearch(ii);
		SearchRange = [DataIn(max(idx01-LeftSearch1,1),1),DataIn(min(idx10+RightSearch1,idx11),1)];
		if (SearchRange(2)>SearchRange(1))
	%%keyboard();
			ExpandRange = max(ceil(0.5*(N0ScanForDiscontinuosCase+1)/(LeftSearch1+RightSearch1)),1);
			Multiple = ii;
			KnotLocation = SearchRange(1):((SearchRange(2)-SearchRange(1))/N0ScanForDiscontinuosCase):SearchRange(2);
			DisErrorSave = zeros(1,numel(KnotLocation));
			DisAnglePsave = zeros(1,numel(KnotLocation));
			for cnt = 1:numel(KnotLocation)
				StartPoint = KnotLocation(cnt);
				[~,DisErrorSave(cnt),DisAnglePsave(cnt)] = TwoPieceBspineKnotEval1(T,Ymat,Degree,Multiple,StartPoint);
			end
			% compute fourfold knot position
			KnotLocation = KnotLocation(~isnan(DisAnglePsave));
			DisErrorSave = DisErrorSave(~isnan(DisAnglePsave));
			DisAnglePsave = DisAnglePsave(~isnan(DisAnglePsave));
			%% Select region
			if ii == Order
				CtrlError1 = min(DisErrorSave)+1e-10;
				disIdx = find(DisErrorSave<CtrlError1);
				MidPosition = round(0.5*(min(disIdx)+max(disIdx)));
				%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
				SearchRangeOut(:,ii) = [KnotLocation(max(disIdx(1)-1,1));KnotLocation(min(disIdx(end)+1,numel(KnotLocation)))];
				OptimalKnotSave(ii) = KnotLocation(MidPosition);
				AngleSave(ii) = DisAnglePsave(MidPosition);
				ErrorSave(ii) = DisErrorSave(MidPosition);
				if ~ScanKnot
					OptimalKnotSave(1:Degree) = OptimalKnotSave(ii);
					for cnt = 1:Degree
						SearchRangeOut(:,cnt) = [KnotLocation(max(disIdx(1)-1,1));KnotLocation(min(disIdx(end)+1,numel(KnotLocation)))];
					end
					break;
				end
			else
				[MinError,disIdx] = min(DisErrorSave);
				SearchRangeOut(:,ii) = [KnotLocation(max(disIdx-ExpandRange,1));KnotLocation(min(disIdx+ExpandRange,numel(KnotLocation)))];
				OptimalKnotSave(ii) = KnotLocation(disIdx);
				AngleSave(ii) = DisAnglePsave(disIdx);
				ErrorSave(ii) = MinError;
			end
		else
			MultipleMax = MultipleMax - 1;
			ScanKnot = 1;
		end
	end
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	for cnt = 1:(min(MultipleMax,Degree))
		StartPoint = OptimalKnotSave(cnt);
		SearchRange = SearchRangeOut(:,cnt)';			%'
		Multiple = cnt;
		[OptimalKnotO,ErrorO,~,~,Angle] = GNKnotSolver1(T,Ymat,...
		Degree,Multiple,SearchRange,StartPoint, GaussNewtonLoopTime);
		if (ScanKnot)&&(MultipleMax==Order)
			StartPoint = OptimalKnotSave(Order);
			SearchRange = SearchRangeOut(:,Order)';		%'
			[OptimalKnotOs,ErrorOs,~,~,Angles] = GNKnotSolver1(T,Ymat,...
			Degree,Multiple,SearchRange,StartPoint, GaussNewtonLoopTime);
			if ErrorOs<ErrorO
				OptimalKnotO = OptimalKnotOs;
				ErrorO = ErrorOs;
				Angle = Angles;
			end
		end
		OptimalKnotSave(cnt) = OptimalKnotO;
		ErrorSave(cnt) = ErrorO;
		AngleSave(cnt) = Angle;
	end
	%% decide multiple knot
	SmoothOutput = Degree - MaxSmooth;
	SmoothOutput = min(SmoothOutput,MultipleMax);
	AngleOut = AngleSave(1:SmoothOutput);
	idx1 = find(AngleOut>MinAngle);

	if numel(idx1)
		ErrorOut = ErrorSave(idx1);
		[~, minidx] = min(ErrorOut);
		OptimalIdx = idx1(minidx);

		OptimalKnotOut = OptimalKnotSave(OptimalIdx);
		MultipleOut = OptimalIdx;
	else
		% eliminate the knot
		[~,idx] = max(AngleOut);
		OptimalKnotOut = OptimalKnotSave(idx);
		MultipleOut = idx;
	end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Subfunctions
%% Find error of a knot
function [OptimalKnot,Error,JoinAngle] = TwoPieceBspineKnotEval1(T,Ymat,Degree,Multiple,StartPoint)
	OptimalKnot = StartPoint;
	Order = Degree + 1;
	T1 = T(:,2);
	[datalength,datasize12] = size(Ymat);
%	 datasize12 = datasize12 + 1;
	% Derivative
	Tcoef = ones(Order,Order);
	for ii = Degree:-1:1
		for cnt = 1:Order
			mul = (cnt - (Order -ii));
			if mul<0
				mul = 0;
			end
			Tcoef(ii,cnt) = Tcoef(ii+1,cnt)*mul;
		end
	end
	Tcal = ones(1,Order);
	for cnt = 2:Order
		Tcal(cnt) = Tcal(cnt-1)*OptimalKnot;
	end
	Tcal = Tcal.*Tcoef(Multiple,:);
	% B-spline evaluation
	Knot = zeros(1,2*Order+Multiple);
	Knot(1:Order) = T1(1);
	Knot((Order+Multiple + 1):end) = T1(end);
	Knot((Order+1):(Order+Multiple))= OptimalKnot;
	left = 1;
	right = datalength;
	middle = fix(0.5*(left+right));
	loopknot = ceil(log2(datalength));
	for cnt = 1:loopknot
		if OptimalKnot < T1(middle)
			right = middle;
		else
			left = middle;
		end
		middle = fix(0.5*(left+right));
	end
	% update Nmat
	Nmat = NewNmatrix(Knot,Degree);
	% generate the knot1 for G'mat		'
	% Compute Basis function base an KnotL and Knot1L

	N = zeros(datalength,Order+Multiple);
	Ncal = reshape(Nmat(:,Order),[Order,Order]);
	Sleft = Tcal*Ncal;
	N(1:middle,1:Order) = T(1:middle,:)*Ncal;
	Ncal = reshape(Nmat(:,Order + Multiple),[Order,Order]);
	Sright = Tcal*Ncal;
	N(middle+1:end,1+Multiple:(Multiple+Order)) = T(middle+1:end,:)*Ncal;
	if datasize12 == 1
		Pctrl = N(:,1:(Multiple+Order))\Ymat;
		Gmat = abs(Ymat-N(:,1:(Multiple+Order))*Pctrl);
	else
		Pctrl = N(:,1:(Multiple+Order))\Ymat;
		G = (Ymat-N(:,1:(Multiple+Order))*Pctrl);
		G = G.*G;
		Gmat = sum(G,2);
		Gmat = sqrt(Gmat);
	end
	Error = max(abs(Gmat));
	%% Calculating join angle
	Sleft1 = Sleft*Pctrl(1:Order,:);
	Sright1 = Sright*Pctrl(Multiple + 1:Multiple+Order,:);
	if datasize12 == 1
		JoinAngle = abs(atan(Sright1/OptimalKnot)-atan(Sleft1/OptimalKnot));
	else
		CosPhi = dot(Sleft1,Sright1)/(norm(Sleft1)*norm(Sright1));
		JoinAngle = acos(CosPhi);
	end
	JoinAngle = JoinAngle*180/pi();
end

#endif

#undef GNKNOTSOLVER1_DEBUG
/*
 * Gauss Newton solver
 */
void
GNKnotSolver1(
	double *T, int T_rows, int T_cols,
	double *Ymat, int ym_cols,
	int Degree, int Multiple, double SearchRange[2],
	double StartPoint, int GaussNewtonLoopTime,
	/* return values */
	double *pOptimalKnotO, double *pErrorO, double *pOptimalKnot,
	double *pError, double *pAngle)
{
	int i;
	int j;
	int cnt;
	double Stepsize = sqrt(2.2204e-16);	/* resolution of double */
	int Order = Degree + 1;

	int datalength = T_rows;
	int datasize12 = ym_cols;
	double KnotLeft = SearchRange[0];
	double KnotRight = SearchRange[1];
	int loopknot = ceil(log2(datalength));

	double T1[datalength];
	for (i = 0; i < datalength; ++i)
		T1[i] = T[i * T_cols + 1];

	double OptimalKnot = StartPoint;
	int Knotsize = 2 * Order + Multiple;
	double Knot[Knotsize];
	for (i = 0; i < Knotsize; ++i)
		Knot[i] = 0;
	for (i = 0; i < Order; ++i)
		Knot[i] = T1[0];

	for (i = 0; i < Order; ++i)
		Knot[Order + Multiple + i] = T1[datalength - 1];

	double Knot1[Knotsize];
	for (i = 0; i < Knotsize; ++i)
		Knot1[i] = Knot[i];

	/* Derivative */
	double Tcoef[Order][Order];
	for (i = 0; i < Order; ++i)
		for (j = 0; j < Order; ++j)
			Tcoef[i][j] = 1;

	int ii;
	for (ii = Degree - 1; ii >= 0; --ii) {
		for (cnt = 0; cnt < Order; ++cnt) {
			int mul = (cnt + 1) - (Order - (ii + 1));
			if (mul < 0)
				mul = 0;
			Tcoef[ii][cnt] = Tcoef[ii + 1][cnt] * mul;
		}
	}

#ifdef GNKNOTSOLVER1_DEBUG
print_matrix("Tcoef", Order, Order, *Tcoef, Order);
#endif

	int iterationstep;
	double deltaXLOld = 0;
	double OptimalKnotLOld = 0;
	double OptimalKnotO = 0;
	double ErrorO = 0;
	double Error = 0;
	double Angle = 0;
	int checkflag = 0;
	for (iterationstep = 1; iterationstep <= GaussNewtonLoopTime; ++iterationstep) {
		/* generate the knot based on the multiple type */
		for (i = Order; i < Order + Multiple; ++i)
			Knot[i] = OptimalKnot;

		double OptimalKnot1 = OptimalKnot + Stepsize;

		for (i = Order; i < Order + Multiple; ++i)
			Knot1[i] = OptimalKnot1;

		/* find knot location */
		int left = 1;
		int right = datalength;
		int left1 = left;
		int right1 = right;
		int middle = (left + right) / 2;
		int middle1 = (left1 + right1) / 2;

#ifdef GNKNOTSOLVER1_DEBUG
printf("middles: %d/%d\n", middle, middle1);
#endif
		for (cnt = 1; cnt <= loopknot; ++cnt) {
			if (OptimalKnot < T1[middle - 1])
				right = middle;
			else
				left = middle;
			if (OptimalKnot1 < T1[middle1 - 1])
				right1 = middle1;
			else
				left1 = middle1;
			middle = (left + right) / 2;
			middle1 = (left1 + right1) / 2;
		}
#ifdef GNKNOTSOLVER1_DEBUG
printf("middles: %d/%d\n", middle, middle1);
#endif

		/* update Nmat */
		int nmsz = (Degree + 1) * (Degree + 1) * (Knotsize - 1);
		double Nmat[nmsz];
		double Nmat1[nmsz];
		NewNmatrix(Knot, Knotsize, Degree, Nmat);
		NewNmatrix(Knot1, Knotsize, Degree, Nmat1);

#ifdef GNKNOTSOLVER1_DEBUG
print_matrix("Nmat", Order * Order, Knotsize - 1, Nmat, Knotsize - 1);
print_matrix("Nmat1", Order * Order, Knotsize - 1, Nmat1, Knotsize - 1);
#endif
		/* generate the knot1 for G'mat	*/
		double Tcal[Order];
		for (i = 0; i < Order; ++i)
			Tcal[i] = 1;

		for (cnt = 1; cnt < Order; ++cnt)
			Tcal[cnt] = Tcal[cnt - 1] * OptimalKnot;

#ifdef GNKNOTSOLVER1_DEBUG
print_matrix("Tcal(1)", 1, Order, Tcal, 1);
#endif
		for (i =0; i < Order; ++i)
			Tcal[i] *= Tcoef[Multiple - 1][i];
#ifdef GNKNOTSOLVER1_DEBUG
print_matrix("Tcal(2)", 1, Order, Tcal, 1);
#endif

		/* Compute Basis function base an KnotL and Knot1L */
		double N[datalength][Order + Multiple];
		for (i = 0; i < datalength; ++i)
			for (j = 0; j < Order + Multiple; ++j)
				N[i][j] = 0;

#ifdef GNKNOTSOLVER1_DEBUG
print_matrix("N", datalength, Order + Multiple, *N, Order + Multiple);
#endif
		double N1[datalength][Order + Multiple];
		for (i = 0; i < datalength; ++i)
			for (j = 0; j < Order + Multiple; ++j)
				N1[i][j] = N[i][j];

		double Ncal[Order][Order];
		for (i = 0; i < Order; ++i)
			for (j = 0; j < Order; ++j)
				Ncal[i][j] = Nmat[(j * Order + i) * (Knotsize - 1) + (Order - 1)];

#ifdef GNKNOTSOLVER1_DEBUG
print_matrix("Ncal", Order, Order, *Ncal, Order);
#endif
		double Sleft[Order];
		cblas_dgemv(CblasRowMajor, CblasTrans, Order, Order, 1.0, *Ncal, Order, Tcal, 1, 0, Sleft, 1);

#ifdef GNKNOTSOLVER1_DEBUG
print_matrix("Sleft", 1, Order, Sleft, 1);
print_matrix("T", datalength, Order, T, Order);
#endif
		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
			middle, Order, Order, 1.0,
			T, Order, *Ncal, Order,
			0, *N, Order + Multiple);
#ifdef GNKNOTSOLVER1_DEBUG
print_matrix("N", datalength, Order + Multiple, *N, Order + Multiple);
#endif

		for (i = 0; i < Order; ++i)
			for (j = 0; j < Order; ++j)
				Ncal[i][j] = Nmat[(j * Order + i) * (Knotsize - 1) + (Order + Multiple - 1)];
#ifdef GNKNOTSOLVER1_DEBUG
print_matrix("Ncal", Order, Order, *Ncal, Order);
#endif
		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
			datalength - middle, Order, Order, 1.0,
			&T[middle * Order], Order, *Ncal, Order,
			0, &N[middle][Multiple], Order + Multiple);
#ifdef GNKNOTSOLVER1_DEBUG
print_matrix("N", datalength, Order + Multiple, *N, Order + Multiple);
#endif
		double Sright[Order];
		cblas_dgemv(CblasRowMajor, CblasTrans, Order, Order, 1.0, *Ncal, Order, Tcal, 1, 0, Sright, 1);
#ifdef GNKNOTSOLVER1_DEBUG
print_matrix("Sright", 1, Order, Sright, 1);
#endif

		for (i = 0; i < Order; ++i)
			for (j = 0; j < Order; ++j)
				Ncal[i][j] = Nmat1[(j * Order + i) * (Knotsize - 1) + (Order - 1)];
		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
			middle1, Order, Order, 1.0,
			T, Order, *Ncal, Order,
			0, *N1, Order + Multiple);

		for (i = 0; i < Order; ++i)
			for (j = 0; j < Order; ++j)
				Ncal[i][j] = Nmat1[(j * Order + i) * (Knotsize - 1) + (Order + Multiple - 1)];

		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
			datalength - middle, Order, Order, 1.0,
			&T[middle * Order], Order, *Ncal, Order,
			0, &N1[middle][Multiple], Order + Multiple);
#ifdef GNKNOTSOLVER1_DEBUG
print_matrix("N1", datalength, Order + Multiple, *N, Order + Multiple);
#endif

		double Pctrl[datalength][datasize12];
		double Pctrl1[datalength][datasize12];
		double G[datalength][datasize12];
		double G1[datalength][datasize12];
		for (i = 0; i < datalength; ++i) {
			for (j = 0; j < datasize12; ++j) {
				Pctrl[i][j] = Ymat[i * datasize12 + j];
				Pctrl1[i][j] = Ymat[i * datasize12 + j];
				G[i][j] = Ymat[i * datasize12 + j];
				G1[i][j] = Ymat[i * datasize12 + j];
			}
		}

		double N_copy[datalength][Order + Multiple];
		double N1_copy[datalength][Order + Multiple];
		for (i = 0; i < datalength; ++i) {
			for (j = 0; j < Order + Multiple; ++j) {
				N_copy[i][j] = N[i][j];
				N1_copy[i][j] = N1[i][j];
			}
		}
		int ret = LAPACKE_dgels(LAPACK_ROW_MAJOR, 'N',
			datalength, Order + Multiple, datasize12, *N_copy, Order + Multiple, *Pctrl, datasize12);
#ifdef GNKNOTSOLVER1_DEBUG
printf("dgels returns %d\n", ret);
#endif
		assert(ret == 0);
#ifdef GNKNOTSOLVER1_DEBUG
print_matrix("Pctrl", Order + Multiple, datasize12, *Pctrl, datasize12);
#endif

		ret = LAPACKE_dgels(LAPACK_ROW_MAJOR, 'N',
			datalength, Order + Multiple, datasize12, *N1_copy, Order + Multiple, *Pctrl1, datasize12);
#ifdef GNKNOTSOLVER1_DEBUG
printf("dgels returns %d\n", ret);
#endif
		assert(ret == 0);
#ifdef GNKNOTSOLVER1_DEBUG
print_matrix("Pctrl1", Order + Multiple, datasize12, *Pctrl1, datasize12);
#endif

		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
			datalength, datasize12, Order + Multiple, -1.0,
			*N, Order + Multiple, *Pctrl, datasize12,
			1, *G, datasize12);
#ifdef GNKNOTSOLVER1_DEBUG
print_matrix("G", datalength, datasize12, *G, datasize12);
#endif

		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
			datalength, datasize12, Order + Multiple, -1.0,
			*N1, Order + Multiple, *Pctrl1, datasize12,
			1, *G1, datasize12);
#ifdef GNKNOTSOLVER1_DEBUG
print_matrix("G1", datalength, datasize12, *G1, datasize12);
#endif

		double Gmat[datalength];
		double Gmat1[datalength];
		if (datasize12 == 1) {
			for (i = 0; i < datalength; ++i) {
				Gmat[i] = abs(G[i][0]);
				Gmat1[i] = abs(G1[i][0]);
			}
		} else {
			for (i = 0; i < datalength; ++i) {
				for (j = 0; j < datasize12; ++j) {
					G[i][j] *= G[i][j];
					G1[i][j] *= G1[i][j];
				}
			}
			for (i = 0; i < datalength; ++i) {
				Gmat[i] = 0;
				Gmat1[i] = 0;
				for (j = 0; j < datasize12; ++j) {
					Gmat[i] += G[i][j];
					Gmat1[i] += G1[i][j];
				}
				Gmat[i] = sqrt(Gmat[i]);
				Gmat1[i] = sqrt(Gmat1[i]);
			}
		}
#ifdef GNKNOTSOLVER1_DEBUG
print_matrix("Gmat", 1, datalength, Gmat, datalength);
print_matrix("Gmat1", 1, datalength, Gmat1, datalength);
#endif

		/* Compute Jacobian matrix G'mat */
		double Jmat[datalength];
		for (i = 0; i < datalength; ++i)
			Jmat[i] = (Gmat1[i] - Gmat[i]) / Stepsize;

#ifdef GNKNOTSOLVER1_DEBUG
print_matrix("Jmat", 1, datalength, Jmat, datalength);
#endif

		double d1 = 0;
		double d2 = 0;
		for (i = 0; i < datalength; ++i) {
			d1 += Jmat[i] * Jmat[i];
			d2 += Jmat[i] * Gmat[i];
		}
		double deltaX = d2 / d1;
#ifdef GNKNOTSOLVER1_DEBUG
printf("deltaX %f\n", deltaX);
#endif

		if (iterationstep == 1) {
			deltaXLOld = deltaX;
		} else {
			if (deltaXLOld * deltaX < 0)
				deltaX = 0.5 * deltaX;
			deltaXLOld = deltaX;
		}
		OptimalKnot = OptimalKnot - deltaX;
		/* calculating angle */
		double Sleft1[datasize12];
		cblas_dgemv(CblasRowMajor, CblasTrans, Order, datasize12, 1.0, *Pctrl, datasize12, Sleft, 1, 0, Sleft1, 1);
#ifdef GNKNOTSOLVER1_DEBUG
print_matrix("Sleft1", 1, datasize12, Sleft1, datasize12);
#endif
		double Sright1[datasize12];
		cblas_dgemv(CblasRowMajor, CblasTrans, Order, datasize12, 1.0, &Pctrl[Multiple][0], datasize12, Sright, 1, 0, Sright1, 1);
#ifdef GNKNOTSOLVER1_DEBUG
print_matrix("Sright1", 1, datasize12, Sright1, datasize12);
#endif
		double JoinAngle;
		if (datasize12 == 1) {
			JoinAngle = abs(atan(Sright1[0]/OptimalKnot) - atan(Sleft1[0]/OptimalKnot));
		} else {
			double dot = 0;
			double norml = 0;
			double normr = 0;
			for (i = 0; i < datasize12; ++i) {
				dot += Sleft1[i] * Sright1[i];
				norml += Sleft1[i] * Sleft1[i];
				normr += Sright1[i] * Sright1[i];
			}
			double CosPhi = dot / sqrt(norml * normr);
#ifdef GNKNOTSOLVER1_DEBUG
printf("dot %f norml %f normr %f CosPhi %f\n", dot, norml, normr, CosPhi);
#endif
			JoinAngle = acos(CosPhi);
		}
		JoinAngle = JoinAngle * 180.0 / M_PI;
#ifdef GNKNOTSOLVER1_DEBUG
printf("JoinAngle %f\n", JoinAngle);
#endif
		/* saturation Optimal knot */
		if (OptimalKnot < KnotLeft)
			OptimalKnot = KnotLeft;
		if (OptimalKnot > KnotRight)
			OptimalKnot = KnotRight;
		Error = 0;
		for (i = 0; i < datalength; ++i) {
			double a = fabs(Gmat[i]);
			if (a > Error)
				Error = a;
		}
#ifdef GNKNOTSOLVER1_DEBUG
printf("Error: %f\n", Error);
#endif
		/* check for stop */
		if (iterationstep > 1) {
			if (Error < ErrorO) {
				OptimalKnotO = OptimalKnot;
				ErrorO = Error;
				Angle = JoinAngle;
			}
			if (fabs(OptimalKnotLOld - OptimalKnot) < 1e-12) {
				if (checkflag == 0)
					checkflag = 1;
				else
					break;
			} else {
				OptimalKnotLOld = OptimalKnot;
				checkflag = 0;
			}
		} else {
			OptimalKnotLOld = OptimalKnot;
			OptimalKnotO = OptimalKnot;
			ErrorO = Error;
			Angle = JoinAngle;
			checkflag = 0;
		}
#ifdef GNKNOTSOLVER1_DEBUG
printf("iterationstep %d OptimalKnotO %f ErrorO %f OptimalKnot %f Error %f Angle %f\n", iterationstep, OptimalKnotO, ErrorO, OptimalKnot, Error, Angle);
#endif
	}
printf("OptimalKnotO %f ErrorO %f OptimalKnot %f Error %f Angle %f\n", OptimalKnotO, ErrorO, OptimalKnot, Error, Angle);
	*pOptimalKnotO = OptimalKnotO;
	*pErrorO = ErrorO;
	*pOptimalKnot = OptimalKnot;
	*pError = Error;
	*pAngle = Angle;
}

#undef NEWNMATRIX_DEBUG
/*
 * return result in Nmat of size (Degree+1)^2 * (nknots - 1)
 */
void
NewNmatrix(double *Knot, int nknots, int Degree, double *Nmat)
{
	int i;
	int j;
	int ii;
	int jj;
	int kk;

	int Order = Degree + 1;
	int rownumbers = Order * (Order + 1) * (Order / 2.0);
	int Intervals = nknots - 1;
	int activeKnot[Intervals];

	double CoeNmat[rownumbers][Intervals];
	for (i = 0; i < rownumbers; ++i)
		for (j = 0; j < Intervals; ++j)
			CoeNmat[i][j] = 0;

	for (ii = 0; ii < Intervals; ++ii) {
		activeKnot[ii] = 0;
		if (Knot[ii + 1] - Knot[ii] != 0) {
			CoeNmat[0][ii] = 1;
			activeKnot[ii] = 1;
		}
	}

#ifdef NEWNMATRIX_DEBUG
for (i=0; i < Intervals; ++i) printf("activeKnot[%d]=%d\n", i, activeKnot[i]);
print_matrix("CoeNmat", rownumbers, Intervals, *CoeNmat, Intervals);
#endif
	for (ii = 0; ii < Intervals; ++ii) {
		if (activeKnot[ii] == 0)
			continue;
		for (jj = 1; jj <= Degree; ++jj) {
			int id0 = Order * (jj - 1) / 2.0 * jj;
			int id1 = Order * jj / 2.0 * (jj + 1);

			for (kk = 0; kk <= jj; ++kk) {
				/*
				 * effective internal we focus,eg when ii=4,jj=1,then
				 * id 2=3,4
				 */
				int id2 = ii - jj + kk;

				/* effective knot num 1 */
				int id2Knot00 = id2 + jj;
				int id2Knot01 = id2Knot00 + 1;

				if ((id2 >= 0) && (id2Knot01 < nknots)) {
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
							N1[i] = CoeNmat[id01 + i][ii];
						}
					} else if (kk == jj) {
						for (i = 0; i < Order; ++i) {
							N0[i] = CoeNmat[id00 + i][ii];
							N1[i] = 0;
						}
					} else {
						for (i = 0; i < Order; ++i) {
							N0[i] = CoeNmat[id00 + i][ii];
							N1[i] = CoeNmat[id01 + i][ii];
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
#ifdef NEWNMATRIX_DEBUG
printf("aden %f a0 %f a1 %f bden %f b0 %f b1 %f\n", aden, a0, a1, bden, b0, b1);
printf("id2Knot00 %d id2Knot01 %d\n", id2Knot00, id2Knot01);
#endif

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
					for (i = 1; i < Order; ++i)
						Acoef[i] = N00[i] + N10[i] +
							   N01[i - 1] + N11[i - 1];
					int id11 = id1 + kk * Order;
					for (i = 0; i < Order; ++i)
						CoeNmat[id11 + i][ii] = Acoef[i];
#ifdef NEWNMATRIX_DEBUG
print_matrix("N0", 1, Order, N0, 1);
print_matrix("N1", 1, Order, N1, 1);
print_matrix("N00", 1, Order, N00, 1);
print_matrix("N01", 1, Order, N01, 1);
print_matrix("N10", 1, Order, N10, 1);
print_matrix("N11", 1, Order, N11, 1);
printf("ii %d jj %d kk %d id0 %d id1 %d id2 %d id11 %d\n", ii, jj, kk, id0, id1, id2, id11);
print_matrix("CoeNmat", rownumbers, Intervals, *CoeNmat, Intervals);
#endif
				}
			}
		}
	}

	int id10 = Order * (Degree / 2.0) * (Degree + 1);

#ifdef NEWNMATRIX_DEBUG
print_matrix("CoeNmat", rownumbers, Intervals, *CoeNmat, Intervals);
printf("id10: %d\n", id10);
#endif

	for (i = id10; i < rownumbers; ++i)
		for (j = 0; j < Intervals; ++j)
			Nmat[(i - id10) * Intervals + j] = CoeNmat[i][j];

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

void
gnknotsolver1_test()
{
	double T[33 * 4] = {
		1.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000,
		1.000000000000000, 0.004900000000000, 0.000024010000000, 0.000000117649000,
		1.000000000000000, 0.008300000000000, 0.000068890000000, 0.000000571787000,
		1.000000000000000, 0.013700000000000, 0.000187690000000, 0.000002571353000,
		1.000000000000000, 0.019100000000000, 0.000364810000000, 0.000006967871000,
		1.000000000000000, 0.024600000000000, 0.000605160000000, 0.000014886936000,
		1.000000000000000, 0.028300000000000, 0.000800890000000, 0.000022665187000,
		1.000000000000000, 0.033700000000000, 0.001135690000000, 0.000038272753000,
		1.000000000000000, 0.041800000000000, 0.001747240000000, 0.000073034632000,
		1.000000000000000, 0.207700000000000, 0.043139290000000, 0.008960030533000,
		1.000000000000000, 0.212900000000000, 0.045326410000000, 0.009649992689000,
		1.000000000000000, 0.226400000000000, 0.051256960000000, 0.011604575744000,
		1.000000000000000, 0.233100000000000, 0.054335610000000, 0.012665630691000,
		1.000000000000000, 0.245600000000000, 0.060319360000000, 0.014814434816000,
		1.000000000000000, 0.250700000000000, 0.062850490000000, 0.015756617843000,
		1.000000000000000, 0.257300000000000, 0.066203290000000, 0.017034106517000,
		1.000000000000000, 0.262100000000000, 0.068696410000000, 0.018005329061000,
		1.000000000000000, 0.274200000000000, 0.075185640000000, 0.020615902488000,
		1.000000000000000, 0.278800000000000, 0.077729440000000, 0.021670967872000,
		1.000000000000000, 0.285400000000000, 0.081453160000000, 0.023246731864000,
		1.000000000000000, 0.290300000000000, 0.084274090000000, 0.024464768327000,
		1.000000000000000, 0.302000000000000, 0.091204000000000, 0.027543608000000,
		1.000000000000000, 0.309600000000000, 0.095852160000000, 0.029675828736000,
		1.000000000000000, 0.327300000000000, 0.107125290000000, 0.035062107417000,
		1.000000000000000, 0.332700000000000, 0.110689290000000, 0.036826326783000,
		1.000000000000000, 0.335500000000000, 0.112560250000000, 0.037763963875000,
		1.000000000000000, 0.349300000000000, 0.122010490000000, 0.042618264157000,
		1.000000000000000, 0.355500000000000, 0.126380250000000, 0.044928178875000,
		1.000000000000000, 0.365600000000000, 0.133663360000000, 0.048867324416000,
		1.000000000000000, 0.372300000000000, 0.138607290000000, 0.051603494067000,
		1.000000000000000, 0.380800000000000, 0.145008640000000, 0.055219290112000,
		1.000000000000000, 0.386200000000000, 0.149150440000000, 0.057601899928000,
		1.000000000000000, 0.391600000000000, 0.153350560000000, 0.060052079296000,
	};
	double Ymat[33 * 2] = {
		94.550, 112.880,
		95.116, 112.436,
		95.508, 112.157,
		96.174, 111.732,
		96.879, 111.377,
		97.626, 111.092,
		98.135, 110.927,
		98.897, 110.722,
		100.063, 110.562,
		124.151, 109.086,
		124.916, 109.076,
		126.867, 109.148,
		127.845, 109.246,
		129.644, 109.542,
		130.361, 109.695,
		131.289, 109.939,
		131.963, 110.149,
		133.611, 110.746,
		134.238, 111.005,
		135.102, 111.406,
		135.739, 111.738,
		137.206, 112.591,
		138.118, 113.212,
		140.137, 114.806,
		140.732, 115.326,
		141.020, 115.621,
		142.381, 117.087,
		142.955, 117.781,
		143.824, 118.966,
		144.353, 119.793,
		144.951, 120.870,
		145.299, 121.579,
		145.574, 122.319,
	};
	double SearchRange[2] = { 0.033700, 0.213160 };
	double OptimalKnotO;
	double ErrorO;
	double OptimalKnot;
	double Error;
	double Angle;
	GNKnotSolver1(T, 33, 4, Ymat, 2, 3, 1, SearchRange, 0.13340, 10,
		&OptimalKnotO, &ErrorO, &OptimalKnot, &Error, &Angle);
	printf("OptimalKnotO %f ErrorO %f OptimalKnot %f Error %f Angle %f\n", OptimalKnotO, ErrorO, OptimalKnot, Error, Angle);
}

// int
// function [BSpline,MaxError,FittedData] = BSplineCurveFittingSerialBisection...
//   (DataIn,Degree,CtrlError,MinAngle,MaxSmooth,N0ScanForDiscontinuosCase,...
//   GaussNewtonLoopTime,ScanKnot)
// % data seperation
// tic
// VectorUX1 = SerialBisection(DataIn,Degree,CtrlError);
// keyboard();
// toc
// tic
// %% find optimal knot
// [~,n]= size(VectorUX1);
// OptimalKnotOut = zeros(1,n-1);
// MultipleOut = zeros(1,n-1);
// %%keyboard();
// for ii = 1:(n-1)
//     idx00 =  VectorUX1(1,ii);
//     DataKnot =  VectorUX1(1,ii+1)-VectorUX1(1,ii);
//     idx11 =  VectorUX1(2,ii+1);
//     DataInKnot = DataIn(idx00:idx11,:);
//     [OptimalKnotOut(ii),MultipleOut(ii)] = TwoPieceOptimalKnotSolver1(DataInKnot,DataKnot,...
//     Degree,MinAngle,MaxSmooth,N0ScanForDiscontinuosCase,GaussNewtonLoopTime,ScanKnot);
// end
// [BSpline,MaxError,FittedData] = BsplineFitting(OptimalKnotOut,MultipleOut,DataIn,Degree);
// toc
// end

int
main(int argc, char **argv)
{
	int len = 0;
	VectorUX_t *VectorUX[94];
	int i;
	int j;
	int k;
	int dim = 2;
	int degree = 3;
	int order = degree + 1;
	double max_error = 0.05;

gnknotsolver1_test();
exit(1);
	SerialBisection(*demo, 94, dim + 1, degree, max_error, VectorUX, &len);

	for (i = 0; i < len; ++i) {
		printf("[%d]: %d-%d ErrorMax %f ", i,
			VectorUX[i]->StartIdx,
			VectorUX[i]->LeftIdx,
			VectorUX[i]->ErrorMax);
		for (j = 0; j < order; ++j) {
			printf("(");
			for (k = 0; k < dim; ++k) {
				if (k > 0)
					printf("/");
				printf("%.3f", VectorUX[i]->CtrlPointY[j * dim + k]);
			}
			printf(") ");
		}
		printf("knots ");
		for (j = 0; j < order * order; ++j) {
			printf("%.3f ", VectorUX[i]->Nmatrix[j]);
		}
		printf("\n");
	}

	return 0;
}
#endif
