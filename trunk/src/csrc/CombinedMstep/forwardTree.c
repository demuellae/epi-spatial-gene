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
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "hmmTree.h" /* All HMM declarations */
#include "specfunc.h" /* All distribution calculations and special function declarations */
//static char rcsid[] = "$Id: forward.c,v 1.2 1998/02/19 12:42:31 kanungo Exp kanungo $";


void ForwardTree(HMMT *phmm, int T, double *O, int numLeaf, double **logalpha, double **logalpha2, double *pmshape1, double *pmshape2, double **B, int g, BaumConfig *baumConf)
{
	int	i, j, k; 	/* state indices */
	int	t;	/* time index */

	double sum;	/* partial sum */
	double sumPhi;
	ForwardConfig *conf = baumConf->forwardConf
	CalcObsProb(phmm, O, T);

	/* Initialize phi2 to -1.0 */
	for (t = 1; t <= T; t++) {
		conf->phi2[t][1] = -1.0;
		conf->scale1[t] = 0.0;
		conf->scale2[t] = 0.0;
		for (j = 2; j <= phmm->N * phmm->N; j++) {
			conf->phi2[t][j] = 0.0;
		}
	}

	for (i = 1; i <= phmm->N; i++) {
		conf->phiT[i] = 0.0;
	}

	/* Loop over sequence */
	for (t = 1; t <= T; t++) {

		/* Initialization */
		if (t <= numLeaf) {
			for (i = 1; i <= phmm->N; i++) {
				conf->phiT[i] = phmm->pi[t][i] * B[t][i];
			}
		} else {
			/* Matrix multiplication of transition prob and row t of phi2
			 * Only occurs if O[t] is not a leaf */
			for (j = 1; j <= phmm->N; j++) {
				sum = 0.0;
				for (i = 1; i <= phmm->N*phmm->N; i++)
					sum += conf->phi2[t][i] * (phmm->AF[i][j]);

				conf->phiT[j] = sum*(B[t][j]);
			}
		}

		sumPhi = 0.0;
		/* sum(Phi(t)) */
		for (i = 1; i <= phmm->N; i++) {
			sumPhi += conf->phiT[i];
		}
		/* set row t of phi to phiT/sumPhi */
		for (i = 1; i <= phmm->N; i++) {
			if (sumPhi != 0)
				conf->phi[t][i] = conf->phiT[i]/sumPhi;
			else
				conf->phi[t][i] = 0.0;
		}
		conf->scale1[t] = conf->scale2[t] + log(sumPhi);

		for (i = 1; i <= phmm->N; i++) {
			logalpha[t][i] = log(conf->phi[t][i]) + conf->scale1[t];
		}

		if (t < T) {

			/*Create a 1-D copy of phi2[T] for the outer product in case we need it
			 * also check if phi2 contains negative values in range 1:N
			 */
			/* Don't compute outer product if phi2[t] from 1:N has negative values */
			if (conf->phi2[conf->P[t]][1] < -.1) {
				for (i = 1; i <= phmm->N; i++) {
					conf->phi2[conf->P[t]][i] = conf->phi[t][i];
				}
			} else {
				/* Store a temporary copy of phi2T[i] for the computation of the outerproduct */
				for (i = 1; i <= phmm->N; i++) {
					conf->phi2Temp[i] = conf->phi2[conf->P[t]][i];
				}
				/* outer product of row t of phi2 and row t of phi */
				/* store result as a row vector in phi2 instead of a matrix */
				k = 1;
				for (i = 1; i <= phmm->N; i++) {
					for (j = 1; j <= phmm->N; j++) {
						conf->phi2[conf->P[t]][k] = conf->phi2Temp[j] * conf->phi[t][i];
						k++;
					}
				}
			}
			conf->scale2[conf->P[t]] = conf->scale1[t] + conf->scale2[conf->P[t]];
			for (i = 1; i <= phmm->N*phmm->N; i++) {
				logalpha2[conf->P[t]][i] = log(conf->phi2[conf->P[t]][i]) + conf->scale2[conf->P[t]];
			}
		}
	}
	baumConf->LL[g] = conf->scale1[T];
	//printf("Log-Likelihood %f\n", *LL);
}

/* 0 = Beta
 * 1 = Binomial
 */
void CalcObsProb(HMMT *phmm, double *O, int T) {
	int i, j;
	/* iterate through all observations */
	if (phmm->dist == 0) {
		for (i = 1; i <= T; i++) {
			/* Beta distribution pdf */
			for (j = 1; j <= phmm->N; j++) {
				if (O[i] <= 0 || O[i] >= 1)
					B[i][j] = 0.0;
				else
					B[i][j] = (powl(O[i], pmshape1[j] - 1.0) * powl(1.0 - O[i],pmshape2[j] - 1.0))
					/ Beta_Function(pmshape1[j],pmshape2[j]);
			}
		}
	} else {
		for (i = 1; i <= T; i++) {
			/* Binom distribution pdf */
			for (j = 1; j <= phmm->N; j++) {
				if (O[i] == 1) {
					B[i][j] = pmshape1[j];
				} else {
					B[i][j] = 1.0-pmshape1[j];
				}
			}
		}
	}
}




