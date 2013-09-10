/*
 **      Author: Tapas Kanungo, kanungo@cfar.umd.edu
 **      Date:   15 December 1997
 **      File:   forward.c
 **      Purpose: Foward algorithm for computing the probabilty
 **		of observing a sequence given a HMM model parameter.
 **      Organization: University of Maryland
 **
 **      $Id: forward.c,v 1.2 1998/02/19 12:42:31 kanungo Exp kanungo $
 */
#include <stdio.h>
#include "hmmTree.h"
static char rcsid[] = "$Id: forward.c,v 1.2 1998/02/19 12:42:31 kanungo Exp kanungo $";

void Forward(HMMT *phmm, int T, int *O, double **alpha, double *pprob)
{
	int     i, j;   /* state indices */
	int     t;      /* time index */

	double sum;     /* partial sum */

	/* 1. Initialization */

	for (i = 1; i <= phmm->N; i++)
		alpha[1][i] = phmm->pi[i]* phmm->B[i][O[1]];

	/* 2. Induction */

	for (t = 1; t < T; t++) {
		for (j = 1; j <= phmm->N; j++) {
			sum = 0.0;
			for (i = 1; i <= phmm->N; i++)
				sum += alpha[t][i]* (phmm->A[i][j]);

			alpha[t+1][j] = sum*(phmm->B[j][O[t+1]]);
		}
	}

	/* 3. Termination */
	*pprob = 0.0;
	for (i = 1; i <= phmm->N; i++)
		*pprob += alpha[T][i];

}

void ForwardWithScale(HMMT *phmm, int T, int *O, int numLeaf, double **alpha,
		double *scale1, double *scale2 double *pprob)
{
	int	i, j; 	/* state indices */
	int	t;	/* time index */

	double sum;	/* partial sum */


	/* 1. Initialization */

	/*
	scale[1] = 0.0;
	for (i = 1; i <= numLeaf; i++) {
	  scale[i+1] = 0.0;
	  sum = 0.0;


	  alpha[i][j] /= sum

	}


	for (i = 1; i <= phmm->N; i++)
	alpha[1][i] /= scale[1]; */

	/* 2. Induction */

	for (t = 1; t <= T - 1; t++) {
		scale[t+1] = 0.0;
		if (t < numLeaf) {

			for (j = 1; j <= phmm->N; j++) {
				alpha[i][j] = phmm->pi[j]* (phmm->B[j][O[i]]);
				scale[i] += alpha[i][j];
				sum += alpha[i][j];
			}

		}
		else {

			for (j = 1; j <= phmm->N; j++) {
				sum = 0.0;
				for (i = 1; i <= phmm->N; i++)
					sum += alpha[t][i]* (phmm->A[i][j]);

				alpha[t+1][j] = sum*(phmm->B[j][O[t+1]]);
			}

			for (j = 1; j <= phmm->N; j++)
				alpha[t+1][j] /= scale[t+1];
		}
	}

	/* 3. Termination */
	*pprob = 0.0;

	for (t = 1; t <= T; t++)
		*pprob += log(scale[t]);

}
