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
	double *O = (double *) dvector(1, T);

	/*
	for (i = 1; i <= N; i++) {
			pmshape1[i] = .5;
			pmshape2[i] = .5;
	} */

	/* Need to create a function to generate O and pmshape1/pmshape2 */

	/* Allocates configuration structures and HMM using nrutil */
	AllocateConfigs(baumConf, T, N, numLeaf);
	AllocateHMM(&hmm, T, N, numLeaf);

	hmm.pmshape1[1] = 1;
	hmm.pmshape1[2] = .5;
	hmm.pmshape1[3] = 3;
	hmm.pmshape2[1] = 3;
	hmm.pmshape2[2] = .5;
	hmm.pmshape2[3] = 1;


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

	BaumWelchTree(&hmm, T, O, baumConf->forwardConf->P, logalpha, logalpha2, logbeta, gamma, &iter, baumConf, 200);


	return 0;
}


