#include <stdio.h>
#include <stdlib.h>
#include "hmmTree.h"
#include "nrutil.h"
#include "specfunc.h"
#include <math.h>

void AllocateHMM(HMMT *phmm, int T, int N, int numLeaf) {
	phmm->AF = (double **) dmatrix(1,N*N,1,N);
	phmm->B = (double **) dmatrix(1,T,1,N);
	phmm->pi = (double **) dmatrix(1,numLeaf,1,N);
	phmm->N = N;
	phmm->pmshape1 = (double *) dvector(1, N);
	phmm->pmshape2 = (double *) dvector(1, N);

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

void FreeConfigs(BaumConfig *baumConf, int T, int N, int numLeaf) {

	free_dmatrix(baumConf->forwardConf->phi,1,T,1,N);
	free_dmatrix(baumConf->forwardConf->phi2,1,T,1,N*N);
	free_dvector(baumConf->forwardConf->phiT,1,N);
	free_dvector(baumConf->forwardConf->phi2Temp,1,N*N);
	free_ivector(baumConf->forwardConf->P,1,T);
	free_dvector(baumConf->forwardConf->scale1, 1, T-1);
	free_dvector(baumConf->forwardConf->scale2, 1, T-1);

	free_ivector(baumConf->backConf->P,1,T);
	free_ivector(baumConf->backConf->bro,1,T);
	free_dvector(baumConf->backConf->scaleB,1, T);
	free_dmatrix(baumConf->backConf->thetaT,1,T,1,N);
	free_dmatrix(baumConf->backConf->theta,1,T,1,N*N);
	free_dvector(baumConf->backConf->thetaTRow, 1, N);

	free_dmatrix(baumConf->F,1, N*N, 1, N);
	free_dvector(baumConf->betaDenom,1, N);
	free_dvector(baumConf->betaY1,1, N);
	free_dvector(baumConf->betaY2,1, N);
	free(baumConf->plogprobinit);
	free(baumConf->plogprobfinal);

	free(baumConf->forwardConf);
	free(baumConf->backConf);

}

void FreeHMM(HMMT *phmm, int T, int N, int numLeaf) {
	free_dmatrix(phmm->AF, 1,N*N,1,N);
	free_dmatrix(phmm->B, 1,T,1,N);
	free_dmatrix(phmm->pi,1,numLeaf,1,N);
	free_dvector(phmm->pmshape1, 1, N);
	free_dvector(phmm->pmshape2, 1, N);
}

void InitTrans(double **trans, int N) {
	int i, j;
	for (i = 1; i <= N * N; i++) {
		for (j = 1; j <= N; j++) {
			trans[i][j] = 1.0/3;
		}
	}
}

void InitDelta(double **delta, double *O, int numLeaf) {
	int i, j;
	for (i = 1; i <= numLeaf; i++) {
		if (O[i] > .5) {
			delta[i][3] = 1.0;
		} else {
			delta[i][1] = 1.0;
		}
	}
}