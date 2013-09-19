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
#include "hmm.h"
#include <math.h>
static char rcsid[] = "$Id: forward.c,v 1.2 1998/02/19 12:42:31 kanungo Exp kanungo $";



void ForwardWithScale(HMMT *phmm, int T, int *O, int numLeaf, double **alpha, double **alpha2, double *LL,
		      Config *conf)
{
	int	i, j, k; 	/* state indices */
	int	t;	/* time index */

	double sum;	/* partial sum */
	double sumPhi;



	/* Forward Section */
	for (t = 1; t <= T - 1; t++) {

		scale1[t+1] = 0.0;
		/* Initialization */
		if (t < numLeaf) {
			for (i = 1; i <= phmm->N*phmm->N; i++) {
				conf->phiT[i] = phmm->pi[t][i];
				conf->phiT[i] *= phmm->B[t][i];
			}
		} else {
			/* Matrix multiplication of transition prob and row t of phi2
			 * Only occurs if O[t] is a leaf */
			for (j = 1; j <= phmm->N; j++) {
				sum = 0.0;
				for (i = 1; i <= phmm->N*phmm->N; i++)
					sum += conf->phi2[t][i]* (phmm->AF[i][j]);

				conf->phiT[j] = sum*(phmm->B[t][j]);
			}
		}

		sumPhi = 0.0;
		/* sum(Phi(t)) */
		for (i = 1; i <= phmm->N; i++) {
			sumPhi += conf->phiT[i];
		}
		/* set row t of phi to phiT/sumPhi */
		for (i = 1; i <= phmm->N; i++) {
			conf->phi[t][i] = conf->phiT[i]/sumPhi;
		}
		conf->scale1[t] = conf->scale2[t] + log(sumPhi);

		for (i = 1; i < phmm->N; i++) {
			alpha[t][i] = conf->phi[t][i] + conf->scale1[i];
		}

		if (t < T) {
			/* If t is one, set phi2 at row t from 1 to N to phi at row t from 1 to n */
			if (t == 1) {
				for (i = 1; i < phmm->N; i++) {
					conf->phi2[t][i] = conf->phi[t][i];
				}
			} else {

				for (i = 1; i < phmm->N; i++) {
					conf->phi2Temp[i] = phi2[t][i];
				}
				/* outer product of row t of phi2 and row t of phi */
				/* store as a vector instead of a matrix */
				k = 1;
				for (i = 1; i < phmm->N; i++) {
					for (j = 1; j < phmm->N; j++) {
						conf->phi2[t][k] = conf->phi2Temp[j] * conf->phi[t][i];
						k++;
					}
				}
			}
			conf->scale2[P[i]] = conf->scale1[i] + conf->scale2[P[i]];
			for (i = 1; i < phmm->N*phmm->N; i++) {
				alpha2[P[t]][i] = log(conf->phi2[P[t]][i]) + scale1[P[i]];
			}
		}
	}
	LL = conf->scale1[T];


	/* 3. Termination */


}

