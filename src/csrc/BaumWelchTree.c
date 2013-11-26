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
#include <stdlib.h>
#include "nrutil.h"
#include <math.h>
#include "specfunc.h"

static char rcsid[] = "$Id: baumwelch.c,v 1.6 1999/04/24 15:58:43 kanungo Exp kanungo $";



void BaumWelchTree(HMMT *phmm, int T, double *O, int *P, double **logalpha, double **logalpha2, double **logbeta,
		double **gamma, int *pniter, BaumConfig *baumConf, int maxiter)
{
	int	i, j;
	int	t, l = 0; /* l is number of iterations */

	double	logprobf;

	double ***xi;
	double **beta, **alpha2;
	double delta, logprobprev;

	double sum;


	double **temp = dmatrix(1, phmm->N * phmm->N, 1, phmm->N);


	FindSiblings(baumConf->backConf->bro, P, baumConf->numLeaf, T);
	xi = AllocXi(T, phmm->N);



	/* Intial forward backward call */
	ForwardTree(phmm, T, O, baumConf->numLeaf, logalpha, logalpha2, &logprobf, baumConf->forwardConf);
	*baumConf->plogprobinit = logprobf; /* log P(O |initial model) */
	BackwardTree(phmm, T, O, baumConf->numLeaf, logbeta, baumConf->forwardConf->phi, baumConf->forwardConf->scale1, baumConf->backConf);

	/* Get exponentiated versions of beta and alpha2 */
	beta = (double **) ExpMatrix(logbeta, T, phmm->N);
	alpha2 = (double **) ExpMatrix(logalpha2, T, phmm->N * phmm->N);

	ComputeGamma(phmm, T, logalpha, logbeta, baumConf->numLeaf, gamma, logprobf);
	ComputeXi(phmm, T, O, baumConf->numLeaf, alpha2, beta, logprobf, xi);

	logprobprev = logprobf;

	do  {

		/* reestimate frequency of state i in time t=1 */
		for (i = 1; i <= phmm->N; i++)
			for (j = 1; j <= baumConf->numLeaf; j++)
				phmm->pi[j][i] = gamma[j][i];


		/* R Code Translated EStep*/
		for (i = 1; i <= phmm->N*phmm->N; i++) {
			sum = 0.0;
			for (j = 1; j <= phmm->N; j++) {
				baumConf->F[i][j] = 0.0;
				for (t = 1; t <= T - baumConf->numLeaf; t++) {
					baumConf->F[i][j] += xi[t][i][j];
				}
				sum += baumConf->F[i][j];
			}

			for (j = 1; j <= phmm->N; j++) {
				if (sum == 0.0)
					phmm->AF[i][j] = 0.0;
				else
					phmm->AF[i][j] = baumConf->F[i][j] / sum;
			}
		}

		/* Beta Maximization Step */
		MstepBeta(phmm, T, baumConf, gamma, O, 200);

		MakeSymmetric(phmm->AF, temp, phmm->N*phmm->N, phmm->N);
		ForwardTree(phmm, T, O, baumConf->numLeaf, logalpha, logalpha2, &logprobf, baumConf->forwardConf);
		BackwardTree(phmm, T, O, baumConf->numLeaf, logbeta, baumConf->forwardConf->phi, baumConf->forwardConf->scale1, baumConf->backConf);
		ComputeGamma(phmm, T, logalpha, logbeta, baumConf->numLeaf, gamma, logprobf);
		ExpMatrices(beta, alpha2, logbeta, logalpha2, T, phmm->N);
		ComputeXi(phmm, T, O, baumConf->numLeaf, alpha2, beta, logprobf, xi);

		/* compute difference between log probability of
		   two iterations */
		delta = logprobf - logprobprev;
		logprobprev = logprobf;
		if (l > MAXITER)
			break;
		l++;

	}
	while (delta > DELTA); /* if log probability does not
                                  change much, exit */
	printf("%f\n",logprobf);
	free_dmatrix(alpha2, 1, T, 1, phmm->N * phmm->N);
	free_dmatrix(beta, 1, T, 1, phmm->N);
	free_dmatrix(temp, 1, phmm->N * phmm->N, 1, phmm->N);
	*pniter = l;
	*baumConf->plogprobfinal = logprobf; /* log P(O|estimated model) */
	//FreeXi(xi, T, phmm->N);
}

/* Compute Maximization step for emission probabilities */
void MstepBeta(HMMT *phmm, int T, BaumConfig *baumConf, double **gamma, double *O, int maxiter) {
	int i, j, t, iter;
	double sum1, sum2;
	double trigammaSum;
	double shape1J, shape2J;
	int i_fault;
	double a, b, c, d, determinant; /* top left, top right, bottom left, and bottom right entries */
	double aNew, bNew, cNew, dNew;
	double incr1, incr2;

	for (i = 1; i <= phmm->N; i++) {
		baumConf->betaDenom[i] = 0.0;
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
		for (iter = 1; iter <= maxiter; iter++) {
			trigammaSum = trigamma(phmm->pmshape1[j] + phmm->pmshape2[j], &i_fault);
			shape1J = DiGamma_Function(phmm->pmshape1[j] + phmm->pmshape2[j]) - DiGamma_Function(phmm->pmshape1[j]) + baumConf->betaY1[j];
			shape2J = DiGamma_Function(phmm->pmshape1[j] + phmm->pmshape2[j]) - DiGamma_Function(phmm->pmshape2[j]) + baumConf->betaY2[j];
			/* general inverse of 2x2 matrix with pmshape1J and pmshape2J as diagonals */
			a = -trigamma(phmm->pmshape1[j], &i_fault) + trigammaSum; // [ a  b ]
			b = trigammaSum;										   // [ c  d ]
			c = trigammaSum;
			d = -trigamma(phmm->pmshape2[j], &i_fault) + trigammaSum;
			determinant = (a * d)-(b * c);
			if (!((determinant/(fabs(shape1J) + fabs(shape2J))/2) < 0.001)) {
				aNew = d/determinant;
				bNew = -b/determinant;
				cNew = -c/determinant;
				dNew = a/determinant;
				/* solution to system of 2 equations */
				incr1 = aNew * shape1J + cNew * shape2J;             // new matrix: [ a  b ][ shape1J ]
				incr2 = bNew * shape1J + dNew * shape2J;			   //             [ c  d ][ shape2J ]
			} else {
				incr1 = 0;
				incr2 = shape1J/b;
			}
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
			if ((incr1 < .00001 && incr1 > -.00001)  && (incr2 < .00001 && incr2 > -.00001))
				break;

		}
	}
}


void ComputeGamma(HMMT *phmm, int T, double **logalpha, double **logbeta, int numLeaf, double **gamma, double LL) {
	int i,t;
	for (i = 1; i <= phmm->N; i++) {
		for (t = 1; t <= T; t++) {
			if (logalpha[t][i] + logbeta[t][i] - LL == -1.0/0.0) {
				gamma[t][i] = 0.0;
			} else {
				gamma[t][i] = exp(logalpha[t][i] + logbeta[t][i] - LL);
			}
		}
	}
}


/* Xi dimensions: (T-numLeaf) X N^2 X N */
void ComputeXi(HMMT* phmm, int T, double *O, int numLeaf, double **alpha2, double **beta, double LL,
		double ***xi)
{
	int i, j;
	int t;

	/* Note: Xi is indexed from 1 to T-numLeaf-1 but represents values of N from numleaf to T-1 */
	for (t = 1; t <= T - numLeaf; t++) {
		for (i = 1; i <= phmm->N*phmm->N; i++)
			for (j = 1; j <= phmm->N; j++) {
				xi[t][i][j] = alpha2[t+numLeaf][i]*beta[t+numLeaf][j]
				                                                   *(phmm->AF[i][j])
				                                                   *(phmm->B[t+numLeaf][j])
				                                                   /exp(LL);
			}


	}
}

double *** AllocXi(int T, int N)
{
	int t;
	double ***xi;

	xi = (double ***) malloc(T*sizeof(double **));

	xi --;

	for (t = 1; t <= T; t++)
		xi[t] = dmatrix(1, N*N, 1, N);
	return xi;

}

void FreeXi(double *** xi, int T, int N)
{
	int t;

	for (t = 1; t <= T; t++)
		free_dmatrix(xi[t], 1, N*N, 1, N);

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
			} else if (P[j] == i) {
				sib[j] = s;
				sib[s] = j;
			}
		}
		flag = 1;
	}
}

/* allocated and return an exponentiated version of matrix */
double ** ExpMatrix(double **mat, int row, int col) {
	int i, j;
	double **newMat = (double **) dmatrix(1,row,1,col);

	for (i = 1; i <= row; i++) {
		for (j = 1;  j <= col; j++) {
			if (mat[i][j] == -1.0/0.0)
				newMat[i][j] = 0.0;
			else
				newMat[i][j] = exp(mat[i][j]);
		}
	}

	return newMat;
}

/* sets values of res1 and res2 to exponentiated versions of mat1 and mat2 */
void ExpMatrices(double ** res1, double **res2, double **mat1, double **mat2, int row, int col) {
	int i,j;
	for (i = 1; i <= row; i++) {
		for (j = 1; j <= col; j++) {
			if (mat1[i][j] == -1.0/0.0)
				res1[i][j] = 0.0;
			else
				res1[i][j] = exp(mat1[i][j]);
		}
	}

	for (i = 1; i <= row; i++) {
		for (j = 1; j <= col * col; j++) {
			if (mat2[i][j] == -1.0/0.0)
				res2[i][j] = 0.0;
			else
				res2[i][j] = exp(mat2[i][j]);
		}
	}


}

void MakeSymmetric(double **asym, double **temp, int row, int col) {
	int i, j;
	int x;
	for (j = 1; j <= col; j++) {
		x = 1;
		for (i = 1; i <= row; i++) {
			if ((i-1) % col == 0 && i != 1)
				x++;
			x = ((x-1) % row) + 1;
			temp[i][j] = (asym[x][j] + asym[i][j])/2;
			x += col;
		}
	}

	for (i = 1; i <= row; i++)
		for (j = 1; j <= col; j++)
			asym[i][j] = temp[i][j];

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
