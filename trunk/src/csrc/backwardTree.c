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



void ForwardWithScale(HMMT *phmm, int T, int *O, int numLeaf, double **alpha, double **alpha2, int *P,
		double *scale1, double *scale2, double **phi, double **phi2)
{
	int	i, j, k; 	/* state indices */
	int	t;	/* time index */

	double LL;
	double sum;	/* partial sum */
	double sumPhi;
	double *phiT = malloc(sizeof(double)*phmm->N);
	double *phi2Temp = malloc(sizeof(double)*phmm->N);



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
			/* Matrix multiplication of transition prob and row t of phi2
			 * Only occurs if O[t] is a leaf */
			for (j = 1; j <= phmm->N; j++) {
				sum = 0.0;
				for (i = 1; i <= phmm->N*phmm->N; i++)
					sum += phi2[t][i]* (phmm->AF[i][j]);

				phiT[j] = sum*(phmm->B[t][j]);
			}
		}

		sumPhi = 0.0;
		/* sum(Phi(t)) */
		for (i = 1; i <= phmm->N; i++) {
			sumPhi += phiT[i];
		}
		/* set row t of phi to phiT/sumPhi */
		for (i = 1; i <= phmm->N; i++) {
			phi[t][i] = phiT[i]/sumPhi;
		}
		scale1[t] = scale2[t] + log(sumPhi);

		for (i = 1; i < phmm->N; i++) {
			alpha[t][i] = phi[t][i] + scale1[i];
		}

		if (t < T) {
			/* If t is one, set phi2 at row t from 1 to N to phi at row t from 1 to n */
			if (t == 1) {
				for (i = 1; i < phmm->N; i++) {
					phi2[t][i] = phi[t][i];
				}
			} else {

				for (i = 1; i < phmm->N; i++) {
					phi2Temp[i] = phi2[t][i];
				}
				/* outer product of row t of phi2 and row t of phi */
				/* store as a vector instead of a matrix */
				k = 1;
				for (i = 1; i < phmm->N; i++) {
					for (j = 1; j < phmm->N; j++) {
						phi2[t][k] = phi2Temp[i] * phi[t][j];
						k++;
					}
				}
			}
			free(phi2Temp);
			scale2[P[i]] = scale1[i] + scale2[P[i]];
			for (i = 1; i < phmm->N*phmm->N; i++) {
				alpha2[P[t]][i] = log(phi2[P[t]][i]) + scale1[P[i]];
			}
		}
	}
	LL = scale1[T];


	/* 3. Termination */


}

