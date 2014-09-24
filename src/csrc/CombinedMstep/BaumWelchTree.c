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



void BaumWelchTree(HMMT *phmm, int T, double **O, int *P, double ***logalpha, double ***logalpha2, double ***logbeta,
		double **gamma, int *pniter, BaumConfig **baumConf, TreeConfig *treeConf, int numGenes, int maxiter)
{
	int	i, j, g;
	int	t, l = 0; /* l is number of iterations */

	double ****xi = (double ****) malloc(sizeof(double ****));
	double delta, logprobprev, logprobf;

	double sum, denom;


	double **temp = dmatrix(1, phmm->N * phmm->N, 1, phmm->N);


	FindSiblings(treeConf->bro, P, treeConf->numLeaf, T);



	/* Intial forward backward call */
	for (g = 1; g <= numGenes; g++) {
		xi[g] = AllocXi(T, phmm->N);
		ForwardTree(phmm, T, O[g], treeConf->numLeaf, logalpha[g], logalpha2[g], baumConf[g], treeConf);
		//*baumConf->plogprobinit = logprobf; /* log P(O |initial model) */
		BackwardTree(phmm, T, O[g], treeConf->numLeaf, logbeta[g], baumConf[g]->forwardConf->phi,
				baumConf[g]->forwardConf->scale1, baumConf[g]->backConf);
		ComputeXi(phmm, T, O[g], treeConf->numLeaf, logalpha2[g], logbeta[g], baumConf[g], xi[g]);
	}
	ComputeGamma(phmm, T, numGenes, logalpha, logbeta, logalpha2, treeConf->numLeaf,
			gamma, baumConf);
	logprobf = MaxLL(baumConf, numGenes);

	logprobprev = logprobf;

	do  {

		/* reestimate frequency of state i in time t=1 */
		for (i = 1; i <= phmm->N; i++)
			for (j = 1; j <= treeConf->numLeaf; j++)
				phmm->pi[j][i] = gamma[j][i];


		/* R Code Translated EStep*/
		for (i = 1; i <= phmm->N*phmm->N; i++) {
			denom = 0.0;
			for (j = 1; j <= phmm->N; j++) {
				//baumConf->F[i][j] = 0.0;
				sum = 0.0;
				for (g = 1; g <= numGenes; g++) {
					for (t = 1; t <= T - treeConf->numLeaf; t++) {
						sum += xi[g][t][i][j]; //sum = F[i, j, g]
					}
				denom += sum;
				}
				if (denom == 0.0)
					phmm->AF[i][j] = 0.0;
				else
					phmm->AF[i][j] = sum / denom;

			}

		}

		/* Beta Maximization Step */
		MakeSymmetric(phmm->AF, temp, phmm->N*phmm->N, phmm->N);
		// MultiThreading here
		for (g = 1; g <= numGenes; g++) {
			Mstep(phmm, T, baumConf[g], gamma, O[g]);
			ForwardTree(phmm, T, O[g], treeConf->numLeaf, logalpha[g], logalpha2[g], baumConf[g], treeConf);
			BackwardTree(phmm, T, O[g], treeConf->numLeaf, logbeta[g], baumConf[g]->forwardConf->phi,
					baumConf[g]->forwardConf->scale1, baumConf[g]->backConf);
			ComputeXi(phmm, T, O[g], treeConf->numLeaf, logalpha2[g], logbeta[g], baumConf[g], xi[g]);

		}
		logprobf = MaxLL(baumConf, numGenes);
		ComputeGamma(phmm, T, numGenes, logalpha, logbeta, logalpha2, treeConf->numLeaf,
				gamma, baumConf);
		//ComputeGamma(phmm, T, logalpha, logbeta, baumConf->numLeaf, gamma);

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
	//printf("%f\n", );
	free_dmatrix(temp, 1, phmm->N * phmm->N, 1, phmm->N);
	*pniter = l;
	//*baumConf->plogprobfinal = logprobf; /* log P(O|estimated model) */
	//FreeXi(xi, T, phmm->N);
}

void Mstep(HMMT *phmm, int T, BaumConfig *baumConf, double **gamma, double *O) {
	if (phmm->dist == 0) {
		MstepBeta(phmm, T, baumConf, gamma, O, 200);
	} else {
		MstepBinom(phmm, T, baumConf, gamma, O);
	}
}

// void MstepBinom(HMMT *phmm, int T, BaumConfig *baumConf, double ***gamma, double **O) {
void MstepBinom(HMMT *phmm, int T, BaumConfig *baumConf, double **gamma, double *O) {
	int i, j, t;
	for (i = 1; i <= phmm->N; i++) {
		baumConf->pmshape1[i] = 0.0;
		baumConf->betaDenom[i] = 0.0;
	}
	for (t = 1; t <= T; t++) {
		for (i = 1; i <= phmm->N; i++) {
			baumConf->pmshape1[i] += gamma[t][i] * (double) O[t];
			baumConf->betaDenom[i] += gamma[t][i];
		}
	}
	for (i = 1; i <= phmm->N; i++) {
		baumConf->pmshape1[i] /= baumConf->betaDenom[i];
	}

}

double MaxLL(BaumConfig **baumConf, int G) {
	int g;
	double max = -1.0/0;
	for (g = 1; g <= G; g++) {
		if (baumConf[g]->LL > max) {
			max = baumConf[g]->LL;
		}
	}
	return max;
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
			trigammaSum = trigamma(baumConf->pmshape1[j] + baumConf->pmshape2[j], &i_fault);
			shape1J = DiGamma_Function(baumConf->pmshape1[j] + baumConf->pmshape2[j]) - DiGamma_Function(baumConf->pmshape1[j]) + baumConf->betaY1[j];
			shape2J = DiGamma_Function(baumConf->pmshape1[j] + baumConf->pmshape2[j]) - DiGamma_Function(baumConf->pmshape2[j]) + baumConf->betaY2[j];
			/* general inverse of 2x2 matrix with pmshape1J and pmshape2J as diagonals */
			a = -trigamma(baumConf->pmshape1[j], &i_fault) + trigammaSum; // [ a  b ]
			b = trigammaSum;										   // [ c  d ]
			c = trigammaSum;
			d = -trigamma(baumConf->pmshape2[j], &i_fault) + trigammaSum;
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
			if (baumConf->pmshape1[j] - incr1 <= 0) { /* new estimates of shape1 are negative */
				baumConf->pmshape1[j] = .01;
			} else {
				baumConf->pmshape1[j] -= incr1;
			}
			if (baumConf->pmshape2[j] - incr2 <= 0) {
				baumConf->pmshape2[j] = .01;
			} else {
				baumConf->pmshape2[j] -= incr2;
			}
			if (baumConf->pmshape1[j] <= 0 || baumConf->pmshape2[j] <= 0) {
				/* invalid parameters */
			}
			if ((incr1 < .00001 && incr1 > -.00001)  && (incr2 < .00001 && incr2 > -.00001))
				break;

		}
	}
}


void ComputeGamma(HMMT *phmm, int T, int numGenes, double ***logalpha, double ***logbeta, double ***logalpha2, int numLeaf,
		double **gamma, BaumConfig **baumConf) {
	int i,j,t,g;
	for (i = 1; i <= phmm->N; i++) {
		for (t = 1; t <= T; t++) {
			gamma[t][i] = 0.0;
		}

	}

	for (g = 1; g <= numGenes; g++) {
		for (i = 1; i <= phmm->N; i++) {
			for (t = 1; t <= T; t++) {
				if (logalpha[g][t][i] + logbeta[g][t][i] - baumConf[g]->LL == -1.0/0.0) {
					gamma[t][i] += 0.0;
				} else {
					gamma[t][i] += exp(logalpha[g][t][i] + logbeta[g][t][i] - baumConf[g]->LL);
				}
			}
		}
	}
	for (i = 1; i <= phmm->N; i++) {
		for (t = 1; t <= T; t++) {
			gamma[t][i] /= numGenes;
		}
	}

}


/* Xi dimensions: (T-numLeaf) X N^2 X N */


void ComputeXi(HMMT *phmm, int T, double *O, int numLeaf, double **logalpha2, double **logbeta, BaumConfig *baumConf,
		double ***xi)
{
	int i, j;
	int t;
	double g = log(phmm->AF[1][1]);

	/* Note: Xi is indexed from 1 to T-numLeaf-1 but represents values of N from numleaf to T-1 */
	for (t = 1; t <= T - numLeaf; t++) {
		for (i = 1; i <= phmm->N*phmm->N; i++)
			for (j = 1; j <= phmm->N; j++) {
				xi[t][i][j] = logalpha2[t+numLeaf][i] + logbeta[t+numLeaf][j] + log(phmm->AF[i][j]) + log(phmm->B[t+numLeaf][j]) - baumConf->LL;
			}
		for (i = 1; i <= phmm->N*phmm->N; i++) {
			for (j = 1; j <= phmm->N; j++) {
				xi[t][i][j] = exp(xi[t][i][j]);
			}
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