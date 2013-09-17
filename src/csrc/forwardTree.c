/*
 **      Author: Tapas Kanungo, kanungo@cfar.umd.edu
 **      Date:   15 December 1997
 **      File:   forward.c
 **      Purpose: Foward algorithm for computing the probabilty
 **		of observing a sequence given a HMMT model parameter.
 **      Organization: University of Maryland
 **
 **      $Id: forward.c,v 1.2 1998/02/19 12:42:31 kanungo Exp kanungo $
 */
#include <stdio.h>
#include "hmm.h"
#include <math.h>
static char rcsid[] = "$Id: forward.c,v 1.2 1998/02/19 12:42:31 kanungo Exp kanungo $";



void ForwardWithScale(HMMT *phmm, int T, int *O, int numLeaf, double **alpha, int *P,
		double *scale1, double *scale2, double **phi, double **phi2)
{
	int	i, j; 	/* state indices */
	int	t;	/* time index */

	double sum;	/* partial sum */
	double sumPhi;
	double phiT[];

	/* 1. Initialization */
	/*
	scale1[1] = 0.0;
	for (i = 1; i <= phmm->N; i++) {
		for (j = 1; j <= numLeaf; j++) {
		alpha[1][i] = phmm->pi[i]* (phmm->B[i][O[1]]);
		scale1[1] += alpha[1][i];
		}
	}

	for (i = 1; i <= phmm->N; i++)
		alpha[1][i] /= scale1[1];
	 */


	/* Forward Section */
	for (t = 1; t <= T - 1; t++) {

		scale1[t+1] = 0.0;
		/* Initialization */
		if (t < numLeaf) {
			for (i = 1; i <= phmm->N*phmm->N; i++) {
				phiT[i] = phmm->pi[t][i];
				phiT[i] *= phmm->B[i][O[t+1]];
			}
		} else {

			for (j = 1; j <= phmm->N*phmm->N; j++) {
				sum = 0.0;
				for (i = 1; i <= phmm->N*phmm->N; i++)
					sum += phi2[t][i]* (phmm->AF[i][j]);

				phiT[j] = sum*(phmm->B[j][O[t+1]]);
			}
		}

		sumPhi = 0.0;
		for (i = 1; i <= phmm->N; i++) {
			sumPhi += phiT[i];
		}

		for (i = 1; i <= phmm->N; i++) {
			phi[t][i] = phiT[i]/sumPhi;
		}
		scale1[t] = scale2[t] + log(sumPhi);

		for (i = 1; i < phmm->N; i++) {
			alpha[t][i] = phi[t][i] + scale1[i];
		}

		if (i < T) {
			for (i = 1; i < phmm->N; i++) {
				sum = 0.0;
				for (j = 1; j < phmm->N; j++) {
					sum += phi2[P[t]][j]*phi[i][j];
				}
				phi2[P[t]][i] = sum;
			}
			scale2[P[t]] = scale2[P[t]] + scale1[t];

		}
	}

	/* 3. Termination */


}

