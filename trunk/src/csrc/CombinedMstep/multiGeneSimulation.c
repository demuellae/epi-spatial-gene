#include <stdio.h>
#include <stdlib.h>
#include "hmmTree.h"
#include "nrutil.h"
#include "specfunc.h"
#include <math.h>

int main() {
	BaumConfig *baumConf = (BaumConfig **) malloc(sizeof(BaumConfig **));
	TreeConfig *treeConf = (TreeConfig *) malloc(sizeof(TreeConfig *));

	int i, j, t, iter;
	int T = 19, numLeaf, numGenes = 100, N = 3;
	double error = 0;

	FILE *obsFile = fopen("Y", "r");
	FILE *parFile = fopen("P", "r");
	FILE *zFile= fopen("z", "r");


	double **O = ReadInputMatD(obsFile, numGenes, T);
	/*
	int *P = ReadInputI(parFile, numGenes, T);
	int **correctZ = ReadInputMatI(zFile, numGenes, T);

	numLeaf = (T+1)/2;

	HMMT hmm = InitHmm(&hmm, O, T, N, P, numLeaf, baumConf);

	double *** logalpha = (double **) d3array(numGenes, T, N);
	double *** logalpha2 = (double **) d3array(numGenes, T, N);
	double *** logbeta = (double **) d3array(numGenes, T, N);
	double ** gamma = (double **) dmatrix(1,T,1,N);
	double LL;


	int *Z = (int *) ivector(1,T);
	/* 0 -- Beta
	   1 -- Binomial */
	/*
	hmm.dist = 1;

	hmm.pmshape1[1] = .5;
	hmm.pmshape1[2] = .5;
	hmm.pmshape1[3] = .5;
	/*hmm.pmshape2[1] = 3;
	hmm.pmshape2[2] = 5;
	hmm.pmshape2[3] = .7;*/
/*

	BaumWelchTree(&hmm, T, O, baumConf->forwardConf->P, logalpha, logalpha2, logbeta, gamma, &iter, baumConf, 600);
	/*
	GenerateZ(gamma, T, Z);



	for (t = 1; t <= T; t++) {
	  printf("%d\n", Z[t]);
	  if (Z[t] != correctZ[t]) {
	    error+= 1.0;
	  } 
	}
	printf("Percent Error: %f\n", error/(T-numLeaf));
	printf("Pi: \n");
	for (i = 1; i <= N*N; i++) {
	  for (j = 1; j <= N; j++) {
	    printf("%f ", hmm.AF[i][j]);
	  }
	  printf("\n");
	}

	printf("hmm.shape1: \n");
	for (i = 1; i <= N; i++) {
	  printf("%f ", hmm.pmshape1[i]);
	}
	printf("\n");

	
	printf("hmm.shape2: \n");
	for (i = 1; i <= N; i++) {
	  printf("%f ", hmm.pmshape2[i]);
	}

	printf("\n");

	FreeHMM(&hmm, T, N, numLeaf);
	free_dmatrix(logalpha, 1, T, 1, N);
	free_dmatrix(logalpha2, 1, T, 1, N*N);
	free_dmatrix(logbeta, 1, T, 1, N);
	free_dmatrix(gamma, 1, T, 1, N);
	//FreeConfigs(baumConf, T, N, numLeaf);
*/


	return 0;
}
