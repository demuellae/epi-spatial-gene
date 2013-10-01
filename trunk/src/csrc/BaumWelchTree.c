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

static char rcsid[] = "$Id: baumwelch.c,v 1.6 1999/04/24 15:58:43 kanungo Exp kanungo $";

#define DELTA 0.001




void BaumWelch(HMMT *phmm, int T, int *O, int *P, double **logalpha, double **logalpha2, double **logbeta,
	double **gamma, int *pniter, ForwardConfig *fConf, BackwardConfig *bConf, int numLeaf,
	double *plogprobinit, double *plogprobfinal)
{
	int	i, j, k;
	int	t, l = 0;

	double	logprobf, logprobb,  threshold;
	double	numeratorA, denominatorA;
	double	numeratorB, denominatorB;

	double ***xi, *scale;
	double **beta, **alpha2;
	double delta, deltaprev, logprobprev;

	deltaprev = 10e-70;



	FindSiblings(bConf->bro, P, numLeaf, T);
	xi = AllocXi(T, phmm->N);
	scale = dvector(1, T);

	ForwardTree(phmm, T, O, numLeaf, logalpha, logalpha2, &logprobf, fConf);
	*plogprobinit = logprobf; /* log P(O |initial model) */
	BackwardTree(phmm, T, O, numLeaf, beta, fConf->phi, bConf);

	alpha2 = ExpMatrix(logalpha2, T, phmm->N * phmm->N);
	beta = ExpMatrix(logbeta, T, phmm->N);

	ComputeGamma(phmm, T, alpha2, logbeta, gamma);
	ComputeXi(phmm, T, O, alpha2, beta, xi);
	FreeMatrix(alpha2);
	FreeMatrix(beta);
	logprobprev = logprobf;

	do  {

		/* reestimate frequency of state i in time t=1 */
		for (i = 1; i <= phmm->N; i++)
			phmm->pi[i] = .001 + .999*gamma[1][i];

		/* reestimate transition matrix  and symbol prob in
		   each state */
		for (i = 1; i <= phmm->N; i++) {
			denominatorA = 0.0;
			for (t = 1; t <= T - 1; t++)
				denominatorA += gamma[t][i];

			for (j = 1; j <= phmm->N; j++) {
				numeratorA = 0.0;
				for (t = 1; t <= T - 1; t++)
					numeratorA += xi[t][i][j];
				phmm->AF[i][j] = .001 +
						.999*numeratorA/denominatorA;
			}

			denominatorB = denominatorA + gamma[T][i];
			for (k = 1; k <= phmm->M; k++) {
				numeratorB = 0.0;
				for (t = 1; t <= T; t++) {
					if (O[t] == k)
						numeratorB += gamma[t][i];
				}

				phmm->B[i][k] = .001 +
						.999*numeratorB/denominatorB;
			}
		}
		MakeSymmetric(phmm->AF, phmm->N * phmm->N, phmm->N);

		ForwardTree(phmm, T, O, numLeaf, logalpha, logalpha2, &logprobf, fConf);
		BackwardTree(phmm, T, O, numLeaf, beta, fConf->phi, bConf);
		ComputeGamma(phmm, T, logalpha, logbeta, gamma);
		ComputeXi(phmm, T, O, alpha2, beta, xi);


		/* compute difference between log probability of
		   two iterations */
		delta = logprobf - logprobprev;
		logprobprev = logprobf;
		l++;

	}
	while (delta > DELTA); /* if log probability does not
                                  change much, exit */

	*pniter = l;
	*plogprobfinal = logprobf; /* log P(O|estimated model) */
	FreeXi(xi, T, phmm->N);
	free_dvector(scale, 1, T);
}

void ComputeGamma(HMM *phmm, int T, double **alpha, double **beta,
	double **gamma)
{

	int 	i, j;
	int	t;
	double	denominator;

	for (t = 1; t <= T; t++) {
		denominator = 0.0;
		for (j = 1; j <= phmm->N; j++) {
			gamma[t][j] = alpha[t][j]*beta[t][j];
			denominator += gamma[t][j];
		}

		for (i = 1; i <= phmm->N; i++)
			gamma[t][i] = gamma[t][i]/denominator;
	}
}

/* Xi dimensions: (T-numLeaf) X N^2 X N */
void ComputeXi(HMMT* phmm, int T, int *O, int numLeaf, double **alpha2, double **beta,
	double ***xi)
{
	int i, j;
	int t;
	double sum;



	for (t = 1; t <= T - numLeaf; t++) {
		sum = 0.0;
		for (i = 1; i <= phmm->N*phmm->N; i++)
			for (j = 1; j <= phmm->N; j++) {
				xi[t][i][j] = alpha2[t+numLeaf][i]*beta[t+numLeaf+1][j] //Fix this for N^2 Cols
					*(phmm->AF[i][j])
					*(phmm->B[j][O[t+numLeaf+1]]);
				sum += xi[t][i][j];
			}

		for (i = 1; i <= phmm->N; i++)
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
		xi[t] = dmatrix(1, N, 1, N); //Change this to N^2 Cols
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
			}
			if (j == s) {
				sib[j] = s;
				sib[s] = j;
			}
		}
		flag = 1;
	}
}

double ** ExpMatrix(double **mat, int row, int col) {
	int i, j;
	double **newMat = malloc(row * sizeof(double *));
	for (i = 1; i <= row; i++) {
		newMat[i] = malloc(col * sizeof(double));
	}

	for (i = 1; i <= row; i++)
		for (j = 1;  j <= col; j++)
			newMat[i][j] = exp(mat[i][j]);

	return newMat;
}

void FreeMatrix(double **mat, int row, int col) {
	int i;
	for (i = 1; i <= row; i++) {
		free(mat[i]);
	}
	free(mat);
}

void MakeSymmetric(double **mat, int row, int col) {

}
