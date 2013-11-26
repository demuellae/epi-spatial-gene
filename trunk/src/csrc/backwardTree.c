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



void BackwardTree(HMMT *phmm, int T, double *O, int numLeaf, double **beta, double **phi, double *scale1, BackwardConfig *conf)
{
	int	i, j, k, idx; 	/* state indices */
	int	t;	/* time index */

	double sum;	/* partial sum */
	double sumTheta;

	/* Initialize theta */
	for (i = 1; i <= phmm->N; i++) {
		conf->theta[T][i] = 1.0/phmm->N;
	}

	conf->scaleB[T] = log(phmm->N);
	for (t = 1; t < T; t++) {
		conf->scaleB[t] = 0.0;
	}

	for (t = T-1; t >= 1; t--) {
		idx = 1;
		k = 0;
		for (i = 1; i <= phmm->N*phmm->N; i++) {
			sum = 0.0;
			for (j = 1; j <= phmm->N; j++) {
				sum += conf->theta[conf->P[t]][j] * phmm->B[conf->P[t]][j] * (phmm->AF[i][j]);
			}
			/* Result of matrix multiplication should be in N * N matrix rather than 1 * N^2 */
			conf->thetaT[(k%phmm->N)+1][idx] = sum;
			k++;
			if (k % phmm->N == 0)
				idx++;
		}

		/* phi X thetaT (1 * N times N * N) */
		for (i = 1; i <= phmm->N; i++) {
			sum = 0.0;
			for (j = 1; j <= phmm->N; j++) {
				sum += phi[conf->bro[t]][j] * conf->thetaT[j][i];
			}
			conf->thetaTRow[i] = sum;
		}

		sumTheta = 0.0;
		for (i = 1; i <= phmm->N; i++) {
			sumTheta += conf->thetaTRow[i];
		}

		for (i = 1; i <= phmm->N; i++) {
			if (sumTheta != 0.0)
				conf->theta[t][i] = conf->thetaTRow[i]/sumTheta;
			else
				conf->theta[t][i] = 0.0;
		}

		conf->scaleB[t] = conf->scaleB[conf->P[t]] + scale1[conf->bro[t]] + log(sumTheta);

		for (i = 1; i <= phmm->N; i++) {
			beta[t][i] = log(conf->theta[t][i]) + conf->scaleB[t];
		}

	}


	/* 3. Termination */


}

