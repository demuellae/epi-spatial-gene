#include <stdio.h>
#include <stdlib.h>
#include "hmmTree.h"
#include "nrutil.h"

int main() {
	int i,j;
	int T = 10;
	int N = 3;
	int M = 3;
	int numLeaf = (T+1/2);

	HMMT *hmm = (HMMT *) malloc(sizeof(HMMT));
	BaumConfig *baumConf = (BaumConfig *) malloc(sizeof(BaumConfig));
	double *pmshape1 = (double *) dvector(1, N);
	double *pmshape2 = (double *) dvector(1, N);
	double *O = (double *) dvector(1, N);

	/* Need to create a function to generate O and pmshape1/pmshape2 */

	/* Allocates configuration structures and HMM using nrutil */
	AllocateConfigs(baumConf, T, N, numLeaf);
	AllocateHMM(hmm, T, N, M, pmshape1, pmshape2);


	for (i = 1; i <= N; i++)
		for (j = 1; j <= N*N; j++)
			hmm->AF[i][j] = 1.0/3.0;

	double ** logalpha = (double **) dmatrix(1,T,1,N);
	double ** logalpha2 = (double **) dmatrix(1,T,1,N*N);
	double ** logbeta = (double **) dmatrix(1,T,1,N);
	double *LL = (double *) malloc(sizeof(double));

	ForwardTree(hmm, T, O, numLeaf, logalpha, logalpha2, LL, baumConf->forwardConf);

	return 0;
}

void AllocateHMM(HMMT *phmm, int T, int N, int M, double *pmshape1, double *pmshape2) {
	phmm->AF = (double **) dmatrix(1,N,1,N);
	phmm->B = (double **) dmatrix(1,T,1,N);
	phmm->M = M;
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
	baumConf->forwardConf->P = (int *) ivector(1,T-1);
	baumConf->forwardConf->scale1 = (double **) dvector(1, T-1);
	baumConf->forwardConf->scale2 = (double **) dvector(1, T-1);

	baumConf->backConf->P = (int *) ivector(1,T-1);
	baumConf->backConf->bro = (int *) ivector(1,T-1);
	baumConf->backConf->scale = (double *) dvector(1, T-1);
	baumConf->backConf->thetaT = (double **) dmatrix(1,T,1,N);
	baumConf->backConf->theta = (double **) dmatrix(1,T,1,N*N);

	baumConf->F = (double **) dmatrix(1, N*N, 1, N);
	baumConf->numLeaf = numLeaf;
	baumConf->betaDenom = (double *) dvector(1, N);
	baumConf->betaY1 = (double *) dvector(1, N);
	baumConf->betaY2 = (double *) dvector(1, N);
	baumConf->plogprobinit = malloc(sizeof(double));
	baumConf->plogprobfinal = malloc(sizeof(double));
}

