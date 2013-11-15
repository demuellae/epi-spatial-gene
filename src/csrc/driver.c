#include <stdio.h>
#include <stdlib.h>
#include "hmmTree.h"
#include "nrutil.h"
#include "specfunc.h"
#include <math.h>

int main() {
	int i,j;
	int T = 5;
	int N = 3;
	int NbyN = N*N;
	int iter;
	int numLeaf = ((T+1)/2);

	HMMT hmm;
	BaumConfig *baumConf = (BaumConfig *) malloc(sizeof(BaumConfig));
	double *pmshape1 = (double *) dvector(1, N);
	double *pmshape2 = (double *) dvector(1, N);
	double *O = (double *) dvector(1, T);

	pmshape1[1] = 1;
	pmshape1[2] = .5;
	pmshape1[3] = 3;
	pmshape2[1] = 3;
	pmshape2[2] = .5;
	pmshape2[3] = 1;
	/*
	for (i = 1; i <= N; i++) {
			pmshape1[i] = .5;
			pmshape2[i] = .5;
	} */

	/* Need to create a function to generate O and pmshape1/pmshape2 */

	/* Allocates configuration structures and HMM using nrutil */
	AllocateConfigs(baumConf, T, N, numLeaf);
	AllocateHMM(&hmm, T, N, numLeaf, pmshape1, pmshape2);

	for (i = 1; i <= NbyN; i++)
		for (j = 1; j <= N; j++)
			hmm.AF[i][j] = 1.0/3.0;

	/*for (i = 1; i <= numLeaf; i++) {
		hmm.pi[i][1] = 1;
	}*/
	hmm.pi[1][3] = 1;
	hmm.pi[2][3] = 1;
	hmm.pi[3][1] = 1;

	/*
	for (i = 1; i <= T-1; i++) {
		baumConf->forwardConf->P[i] = i+1;
	}
	*/
	baumConf->forwardConf->P[1] = 4;
	baumConf->forwardConf->P[2] = 4;
	baumConf->forwardConf->P[3] = 5;
	baumConf->forwardConf->P[4] = 5;
	baumConf->forwardConf->P[5] = 0;


	/*for (i = 1; i <= T; i++) {
		O[i] = .5;
	} */

	O[1] = .95;
	O[2] = .95;
	O[3] = .05;
	O[4] = .95;
	O[5] = .5;

	baumConf->backConf->bro[1] = 2;
	baumConf->backConf->bro[2] = 1;
	baumConf->backConf->bro[3] = 4;
	baumConf->backConf->bro[4] = 3;
	baumConf->backConf->bro[5] = 0;

	baumConf->backConf->P = baumConf->forwardConf->P;


	double ** logalpha = (double **) dmatrix(1,T,1,N);
	double ** logalpha2 = (double **) dmatrix(1,T,1,NbyN);
	double ** logbeta = (double **) dmatrix(1,T,1,N);
	double ** gamma = (double **) dmatrix(1,T,1,N);
	double LL;

	//ForwardTree(&hmm, T, O, numLeaf, logalpha, logalpha2, &LL, baumConf->forwardConf);
	//printf("Log Likelihood: %f\n", LL);
	//printf("%f %f %f\n%f %f %f\n %f %f %f\n", logalpha[1][1], logalpha[1][2], logalpha[1][3], logalpha[2][1], logalpha[2][2], logalpha[2][3],
	//		logalpha[3][1], logalpha[3][2], logalpha[3][3]);
	//BackwardTree(&hmm, T, O, numLeaf, logbeta, baumConf->forwardConf->phi, baumConf->forwardConf->scale1, baumConf->backConf);
	//printf("logbeta[1][1]: %f\n", logbeta[1][1]);

	BaumWelchTree(&hmm, T, O, baumConf->forwardConf->P, logalpha, logalpha2, logbeta, gamma, &iter, baumConf, 200);


	return 0;
}

void AllocateHMM(HMMT *phmm, int T, int N, int numLeaf, double *pmshape1, double *pmshape2) {
	int i,j;
	phmm->AF = (double **) dmatrix(1,N*N,1,N);
	phmm->B = (double **) dmatrix(1,T,1,N);
	phmm->pi = (double **) dmatrix(1,numLeaf,1,N);
	phmm->N = N;
	phmm->pmshape1 = pmshape1;
	phmm->pmshape2 = pmshape2;

}

void AllocateConfigs(BaumConfig *baumConf, int T, int N, int numLeaf) {

	baumConf->forwardConf = malloc(sizeof(ForwardConfig));
	baumConf->backConf = malloc(sizeof(BackwardConfig));

	baumConf->forwardConf->phi = (double **) dmatrix(1,T,1,N);
	baumConf->forwardConf->phi2 = (double **) dmatrix(1,T,1,N*N);
	baumConf->forwardConf->phiT = (double *) dvector(1,N);
	baumConf->forwardConf->phi2Temp = (double *) dvector(1,N*N);
	baumConf->forwardConf->P = (int *) ivector(1,T);
	baumConf->forwardConf->scale1 = (double *) dvector(1, T-1);
	baumConf->forwardConf->scale2 = (double *) dvector(1, T-1);

	baumConf->backConf->P = (int *) ivector(1,T);
	baumConf->backConf->bro = (int *) ivector(1,T);
	baumConf->backConf->scaleB = (double *) dvector(1, T);
	baumConf->backConf->thetaT = (double **) dmatrix(1,T,1,N);
	baumConf->backConf->theta = (double **) dmatrix(1,T,1,N*N);
	baumConf->backConf->thetaTRow = (double *) dvector(1, N);

	baumConf->F = (double **) dmatrix(1, N*N, 1, N);
	baumConf->numLeaf = numLeaf;
	baumConf->betaDenom = (double *) dvector(1, N);
	baumConf->betaY1 = (double *) dvector(1, N);
	baumConf->betaY2 = (double *) dvector(1, N);
	baumConf->plogprobinit = malloc(sizeof(double));
	baumConf->plogprobfinal = malloc(sizeof(double));
}

