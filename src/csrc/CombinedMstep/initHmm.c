#include <stdio.h>
#include <stdlib.h>
#include "hmmTree.h"
#include "nrutil.h"
#include "specfunc.h"
#include <math.h>

HMMT InitHMM(double *O, int T, int N, int *P, int numLeaf, BaumConfig *baumConf) {
	HMMT hmm;
	AllocateConfigs(baumConf, T, N, numLeaf);
	baumConf->forwardConf->P = P;
	baumConf->backConf->P = P;
	AllocateHMM(&hmm, T, N, numLeaf);
	InitTrans(hmm.AF, N);
	InitDelta(hmm.pi, O, numLeaf);
	return hmm;
}

void AllocateHMM(HMMT *phmm, int T, int N, int numLeaf) {
	phmm->AF = (double **) dmatrix(1,N*N,1,N);
	phmm->B = (double **) dmatrix(1,T,1,N);
	phmm->pi = (double **) dmatrix(1,numLeaf,1,N);
	phmm->N = N;
	phmm->pmshape1 = (double *) dvector(1, N);
	phmm->pmshape2 = (double *) dvector(1, N);

}

void AllocateConfigs(BaumConfig **baumConf, int T, int G, int N, int numLeaf) {
	int g;
	for (g = 1; g <= G; g++) {
		baumConf[g] = (BaumConfig *) malloc(sizeof(BaumConfig));
	}

	for (g = 1; g <= G; g++) {
		baumConf[g]->forwardConf = malloc(sizeof(ForwardConfig));
		baumConf[g]->backConf = malloc(sizeof(BackwardConfig));

		baumConf[g]->forwardConf->phi = (double **) dmatrix(1,T,1,N);
		baumConf[g]->forwardConf->phi2 = (double **) dmatrix(1,T,1,N*N);
		baumConf[g]->forwardConf->phiT = (double *) dvector(1,N);
		baumConf[g]->forwardConf->phi2Temp = (double *) dvector(1,N*N);
		baumConf[g]->forwardConf->P = (int *) ivector(1,T);
		baumConf[g]->forwardConf->scale1 = (double *) dvector(1, T);
		baumConf[g]->forwardConf->scale2 = (double *) dvector(1, T);

		baumConf[g]->backConf->P = (int *) ivector(1,T);
		baumConf[g]->backConf->bro = (int *) ivector(1,T);
		baumConf[g]->backConf->scaleB = (double *) dvector(1, T);
		baumConf[g]->backConf->thetaT = (double **) dmatrix(1,N,1,N);
		baumConf[g]->backConf->theta = (double **) dmatrix(1,T,1,N*N);
		baumConf[g]->backConf->thetaTRow = (double *) dvector(1, N);

		baumConf[g]->F = (double **) dmatrix(1, N*N, 1, N);
		baumConf[g]->numLeaf = numLeaf;
		baumConf[g]->betaDenom = (double *) dvector(1, N);
		baumConf[g]->betaY1 = (double *) dvector(1, N);
		baumConf[g]->betaY2 = (double *) dvector(1, N);
		baumConf[g]->plogprobinit = malloc(sizeof(double));
		baumConf[g]->plogprobfinal = malloc(sizeof(double));
	}
}

void FreeConfigs(BaumConfig *baumConf, int T, int N, int numLeaf) {

	free_dmatrix(baumConf->forwardConf->phi,1,T,1,N);
	free_dmatrix(baumConf->forwardConf->phi2,1,T,1,N*N);
	free_dvector(baumConf->forwardConf->phiT,1,N);
	free_dvector(baumConf->forwardConf->phi2Temp,1,N*N);
	free_ivector(baumConf->forwardConf->P,1,T);
	free_dvector(baumConf->forwardConf->scale1, 1, T);
	free_dvector(baumConf->forwardConf->scale2, 1, T);

	free_ivector(baumConf->backConf->P,1,T);
	free_ivector(baumConf->backConf->bro,1,T);
	free_dvector(baumConf->backConf->scaleB,1, T);
	free_dmatrix(baumConf->backConf->thetaT,1,N,1,N);
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

void TestInitTrans(double **trans, int N) {
	int i, j;

	for (i = 1; i <= N*N; i++) {
		for (j = 1; j <= N; j++) {
			trans[i][j] = 0.002;
			if (j == 2)
				trans[i][j] = .996;
		}
	}
	trans[1][1] = .996;
	trans[N*N][N] = .996;
	trans[1][2] = .002;
	trans[9][2] = .002;

}

void GenerateZ(double **gamma, int T, int *Z) {
	int t;
	double max;
	for (t = 1; t <= T; t++) {
		max = 0.0;
		if (gamma[t][1] > gamma[t][2]) {
			max = gamma[t][1];
			Z[t] = -1;
		} else {
			max = gamma[t][2];
			Z[t] = 0;
		}

		if (gamma[t][3] > max) {
			max = gamma[t][3];
			Z[t] = 1;
		}
	}
}

void ReadCommaSequence(char *buffer, double *O) {
	int i = 0;
	int j = 0;
	int oIdx = 0;
	double sum = 0.0;
	while (buffer[i] != '\0') {
		while (buffer[i] != ',') {
			if (buffer[i] == '\0')
				break;
			i++;
		}
		buffer[i] = '\0';
		O[oIdx] = atof(&buffer[j]);
		i++;
		j = i;
		oIdx++;
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
