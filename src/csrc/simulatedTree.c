#include <stdio.h>
#include <stdlib.h>
#include "hmmTree.h"
#include "nrutil.h"
#include "specfunc.h"
#include <math.h>

int main() {
	HMMT hmm;
	BaumConfig *baumConf = (BaumConfig *) malloc(sizeof(BaumConfig));

	int i, iter;
	int T, numLeaf, N = 3;

	int x = abs(3);

	FILE *obsFile = fopen("Y", "r");
	FILE *parFile = fopen("P", "r");

	double *O = ReadInputD(obsFile, &T);
	int *P = ReadInputI(parFile, &T);

	numLeaf = (T+1)/2;

	AllocateConfigs(baumConf, T, N, numLeaf);
	baumConf->forwardConf->P = P;
	baumConf->backConf->P = P;
	AllocateHMM(&hmm, T, N, numLeaf);
	InitTrans(hmm.AF, N);
	InitDelta(hmm.pi, O, numLeaf);

	double ** logalpha = (double **) dmatrix(1,T,1,N);
	double ** logalpha2 = (double **) dmatrix(1,T,1,N * N);
	double ** logbeta = (double **) dmatrix(1,T,1,N);
	double ** gamma = (double **) dmatrix(1,T,1,N);
	double LL;

	hmm.pmshape1[1] = 1;
	hmm.pmshape1[2] = 2;
	hmm.pmshape1[3] = 3;
	hmm.pmshape2[1] = 3;
	hmm.pmshape2[2] = 5;
	hmm.pmshape2[3] = 1;


	BaumWelchTree(&hmm, T, O, baumConf->forwardConf->P, logalpha, logalpha2, logbeta, gamma, &iter, baumConf, 200);

	FreeHMM(&hmm, T, N, numLeaf);
	free_dmatrix(logalpha, 1, T, 1, N);
	free_dmatrix(logalpha2, 1, T, 1, N*N);
	free_dmatrix(logbeta, 1, T, 1, N);
	free_dmatrix(gamma, 1, T, 1, N);
	//FreeConfigs(baumConf, T, N, numLeaf);

	return 0;
}
