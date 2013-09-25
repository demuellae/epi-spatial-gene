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
#include <math.h>
static char rcsid[] = "$Id: forward.c,v 1.2 1998/02/19 12:42:31 kanungo Exp kanungo $";



void BackwardWithScale(HMMT *phmm, int T, int *O, int numLeaf, double **beta, double **phi, BackwardConfig **conf)
{
	int	i, j, k; 	/* state indices */
	int	t;	/* time index */

	double sum;	/* partial sum */
	double sumTheta;

	/* Initialize theta */
	for (i = 1; i < phmm->N; i++) {
		conf->theta[phmm->N][i] = 1.0/phmm->N;
	}

	/* Initialize scale */
	for (i = 1; i <= T-1; i++) {
		conf->scale[i] = 0.0;
	}
	conf->scale[T] = log(phmm->N);

	for (t = T; t >= 1; t--) {

		k = 0;
		for (i = 1; i <= phmm->N*phmm->N; i++) {
			sum = 0.0;
			for (j = 1; j <= phmm->N; j++) {
				sum += conf->thetaT[t][j] * phmm->B[conf->P[t]][j] * (phmm->AF[j][i]);
			}
			/* Result of matrix multiplication should be in N * N matrix rather than 1 * N^2 */
			conf->thetaT[k%phmm->N][i] = sum;
			k++;
		}

		/* phi X thetaT (1 * N times N * N) */
		for (i = 1; i <= phmm->N; i++) {
			sum = 0.0;
			for (j = 1; j <= phmm->N; j++) {
				sum += phi[t][j] * conf->thetaT[j][i];
			}
			conf->thetaTRow[i] = sum;
		}

		sumTheta = 0.0;
		for (i = 1; i <= phmm->N; i++) {
			sumTheta += conf->thetaTRow[i];
		}

		for (i = 1; i <= phmm->N; i++) {
			conf->theta[t][i] = conf->thetaTRow[i]/sumTheta;
		}

		conf->scale[t] = conf->scale[conf->P[t]] + conf->scale[conf->bro[t]] + log(sumTheta);

		for (i = 1; i <= T; i++) {
			beta[t][i] = log(conf->theta[t][i]) + conf->scale[t];
		}

	}


	/* 3. Termination */


}

