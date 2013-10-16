/*
 **      Author: Tapas Kanungo, kanungo@cfar.umd.edu
 **      Date:   15 December 1997
 **      File:   baumwelch.c
 **      Purpose: Baum-Welch algorithm for estimating the parameters
 **              of a HMM model, given an observation sequence.
 **      Organization: University of Maryland
 **
 **	Update:
 **	Author: Tapas Kanungo
 **	Date:	19 April 1999
 **	Purpose: Changed the convergence criterion from ratio
 **		to absolute value.
 **
 **      $Id: baumwelch.c,v 1.6 1999/04/24 15:58:43 kanungo Exp kanungo $
 */


#include "hmmTree.h"
#include <stdio.h>
#include "nrutil.h"
#include <math.h>

#include "specialfunctions/digamma_function.c"
#include "specialfunctions/asa121.c"

static char rcsid[] = "$Id: baumwelch.c,v 1.6 1999/04/24 15:58:43 kanungo Exp kanungo $";

#define DELTA 0.001




void BaumWelchTree(HMMT *phmm, int T, int *O, int *P, double **logalpha, double **logalpha2, double **logbeta,
		double **gamma, int *pniter, BaumConfig *baumConf, int maxiter)
{
	int	i, j, k;
	int	t, l = 0;



	double	logprobf, logprobb,  threshold;
	double	 denominatorA;
	double	 denominatorB;

	double ***xi, *scale;
	double **beta, **alpha2;
	double delta, deltaprev, logprobprev;

	double sum;

	deltaprev = 10e-70;


	FindSiblings(baumConf->backConf->bro, P, baumConf->numLeaf, T);
	xi = AllocXi(T, phmm->N);
	scale = dvector(1, T);

	beta = ExpMatrix(logalpha, T, phmm->N);
	alpha2 = ExpMatrix(logbeta, T, phmm->N * phmm->N);
	ForwardTree(phmm, T, O, baumConf->numLeaf, logalpha, logalpha2, &logprobf, baumConf->forwardConf);
	*baumConf->plogprobinit = logprobf; /* log P(O |initial model) */
	BackwardTree(phmm, T, O, baumConf->numLeaf, beta, baumConf->forwardConf->phi, baumConf->backConf);

	ComputeGamma(phmm, T, alpha2, logbeta, gamma);
	ComputeXi(phmm, T, O, baumConf->numLeaf, alpha2, beta, xi);

	logprobprev = logprobf;

	do  {

		/* reestimate frequency of state i in time t=1 */
		for (i = 1; i <= phmm->N; i++)
			for (j = 1; j <= baumConf->numLeaf; j++)
				phmm->pi[i] = gamma[j][i];


		/* R Code Translated MStep*/
		for (i = 1; i < phmm->N*phmm->N; i++) {
			sum = 0.0;
			for (j = 1; j <= phmm->N; j++) {
				for (t = 1; t <= T - baumConf->numLeaf - 1; t++) {
					baumConf->F[i][j] += xi[t][i][j];
				}
				sum += baumConf->F[i][j];
			}

			for (j = 1; j <= phmm->N; j++) {
				phmm->AF[i][j] = baumConf->F[i][j] / sum;
			}
		}

		/* Beta Maximization Step */
		MstepBeta(phmm, T, baumConf, gamma, O, 200);

		MakeSymmetric(phmm->AF, phmm->AB, phmm->N, phmm->N);
		ForwardTree(phmm, T, O, baumConf->numLeaf, logalpha, logalpha2, &logprobf, baumConf->forwardConf);
		BackwardTree(phmm, T, O, baumConf->numLeaf, beta, baumConf->forwardConf->phi, baumConf->backConf);
		ComputeGamma(phmm, T, logalpha, logbeta, gamma);
		Exp_Matrices(alpha2, beta, phmm->N, T);
		ComputeXi(phmm, T, O, baumConf->numLeaf, alpha2, beta, xi);

		/* compute difference between log probability of
		   two iterations */
		delta = logprobf - logprobprev;
		logprobprev = logprobf;
		l++;

	}
	while (delta > DELTA); /* if log probability does not
                                  change much, exit */
	FreeMatrix(alpha2);
	FreeMatrix(beta);
	*pniter = l;
	*baumConf->plogprobfinal = logprobf; /* log P(O|estimated model) */
	FreeXi(xi, T, phmm->N);
	free_dvector(scale, 1, T);
}

void MstepBeta(HMMT *phmm, int T, BaumConfig baumConf, double **gamma, int *O, maxiter) {
	int i, j, t, iter;
	double sum1, sum2;
	double shape1J, shape2J;
	int i_fault;
	double dLL;
	double a, b, c, d, determinant; /* top left, top right, bottom left, and bottom right entries */
	double incr1, incr2;

	for (i = 1; i <= phmm->N; i++) {
		for (t = 1; t <= T; t++) {
			baumConf->betaDenom[i] += gamma[t][i];
		}
	}

	for (i = 1; i <= phmm->N; i++) {
		sum1 = 0.0;
		sum2 = 0.0;
		for (t = 1; t <= T; t++) {
			sum1 += gamma[t][i] * log(O[t]);
			sum2 += gamma[t][i] * log(1-O[t]);
		}
		baumConf->betaY1[i] = sum1/baumConf->betaDenom[i];
		baumConf->betaY2[i] = sum2/baumConf->betaDenom[i];

	}

	for (j = 1; j <= phmm->N; j++) {

		for (iter = 1; iter <= maxiter; i++) {
			sum1 = 0.0;
			shape1J = phmm->pmshape1[j] + phmm->pmshape2[j] - DiGamma_Function(phmm->pmshape1[j]) + baumConf->betaY1[j]; //digamma is a placeholder
			shape2J = phmm->pmshape1[j] + phmm->pmshape2[j] - DiGamma_Function(phmm->pmshape2[j]) + baumConf->betaY2[j];
			sum1 = trigamma(phmm->pmshape1[j], &i_fault) + trigamma(phmm->pmshape2[j], &i_fault);
			/* general inverse of 2x2 matrix with pmshape1J and pmshape2J as diagonals */
			a = -trigamma(phmm->pmshape1[j], &i_fault) + sum1;
			b = sum1;
			c = sum1;
			d = -trigamma(phmm->pmshape2[j], &i_fault) + sum1;
			determinant = 1/((a * d)-(b * c));
			a = d/determinant;
			b = -b/determinant;
			c = -c/determinant;
			d = a/determinant;
			incr1 = a * shape1J + b * shape2J;
			incr2 = c * shape1J + d * shape2J;
			if (phmm->pmshape1[j] - incr1 <= 0) { /* new estimates of shape1 are negative */
				phmm->pmshape1[j] = .01;
			} else {
				phmm->pmshape1[j] -= incr1;
			}
			if (phmm->pmshape2[j] - incr2 <= 0) {
				phmm->pmshape2[j] = .01;
			} else {
				phmm->pmshape2[j] -= incr2;
			}
			if (phmm->pmshape1[j] <= 0 || phmm->pmshape2[j] <= 0) {
				/* invalid parameters */
			}
			if (abs(incr1) < .00001 || abs(incr2) < .00001)
				break;

		}
	}
}


void ComputeGamma(HMMT *phmm, int T, double **logalpha, double **logbeta, int numLeaf, double **gamma, double LL) {
	int i,t;
	for (i = 1; i <= phmm->N; i++) {
		for (t = 1; t <= T; t++) {
			exp(gamma[t][i] = logalpha[t][i] + logbeta[t][i] - LL);
		}
	}
}


/* Xi dimensions: (T-numLeaf) X N^2 X N */
void ComputeXi(HMMT* phmm, int T, int *O, int numLeaf, double **alpha2, double **beta,
		double ***xi)
{
	int i, j;
	int t;
	double sum;

	for (t = 1; t <= T - numLeaf-1; t++) {
		sum = 0.0;
		for (i = 1; i <= phmm->N*phmm->N; i++)
			for (j = 1; j <= phmm->N; j++) {
				xi[t][i][j] = alpha2[t+numLeaf][i]*beta[t+numLeaf+1][j] //Fix this for N^2 Cols
				             *(phmm->AF[i][j])
				             *(phmm->B[j][O[t+numLeaf+1]]);
				sum += xi[t][i][j];
			}

		for (i = 1; i <= phmm->N*phmm->N; i++)
			for (j = 1; j <= phmm->N; j++)
				xi[t][i][j]  /= sum;
	}
}

double *** AllocXi(int T, int N)
{
	int t;
	double ***xi;

	xi = (double ***) malloc(T*sizeof(double **));

	xi --;

	for (t = 1; t <= T; t++)
		xi[t] = dmatrix(1, N*N, 1, N*N); //Change this to N^2 Cols
	return xi;

}

void FreeXi(double *** xi, int T, int N)
{
	int t;

	for (t = 1; t <= T; t++)
		free_dmatrix(xi[t], 1, N, 1, N);

	xi ++;
	free(xi);

}

void FindSiblings(int *sib, int *P, int numleaf, int T) {
	int i, j, flag = 1, s = -1;
	for (i = numleaf + 1; i <= T; i++) {
		for (j = 1; j <= T; j++) {
			if (P[j] == i && flag) {
				s = j;
				flag = 0;
			} else if (j == s) {
				sib[j] = s;
				sib[s] = j;
			}
		}
		flag = 1;
	}
}

double ** ExpMatrix(double **mat, int row, int col) {
	int i, j;
	double **newMat = dmatrix(1,row,1,col);

	for (i = 1; i <= row; i++)
		for (j = 1;  j <= col; j++)
			newMat[i][j] = exp(mat[i][j]);

	return newMat;
}

void ExpMatrices(double **mat1, double **mat2, int row, int col) {
	int i,j;
	for (i = 1; i <= row; i++) {
		for (j = 1; j <= col; j++) {
			mat1[i][j] = exp(mat1[i][j]);
			mat2[i][j] = exp(mat2[i][j]);
		}
	}
}

void FreeMatrix(double **mat, int row, int col) {
	int i;
	for (i = 1; i <= row; i++) {
		free(mat[i]);
	}
	free(mat);
}

/* This function is probably way inefficient */
void MakeSymmetric(double **asym, double **sym, int row, int col) {
	int i, j, l;
	int x1, y1, x2, y2;
	for (j = 1; j <= row; j++) {
		y1 = 0;
		for (i = 1; i <= col; i++) {
			if (j % col == 1)
				y1++;
			x1 = (j-1) % row;
			y2 = 0;
			if (x1 != y1)
				for (l = 1; l <= col; l++) {
					if (l % row == 1)
						y2++;
					x2 = (l-1) % row;
					if (x1 == y2 && x2 == y1) {
						sym[i][j] = (asym[i][j] + asym[i][l])/2;
						sym[i][l] = (asym[i][j] + asym[i][l])/2;
					}
				}
		}
	}
}

void ThreeToTwo(double **threeD, double **twoD, int row, int col) {
	int i,j,k,l;
	for (i = 1; i <= row; i++) {
		for (j = 1; j <= col; j++) {
			k = 1;
			l = 1 + (col*(j-1));
			while (k <= col*col) {
				twoD[i][j] += threeD[k][j];
				twoD[i][j] += threeD[l][j];
				k += col;
				l += 1;
			}
			twoD[i][j] /= 6;
		}
	}
}
