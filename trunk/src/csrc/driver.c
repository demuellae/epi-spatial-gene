#include "hmmTree.h"

int main() {
	int seqLength = 10;
	int numStates = 3;
	int numSymbols = 3;
	int numLeaf = 4;

	/* indices for malloc of 2-d arrays */
	int i;

	HMMT *hmm = malloc(sizeof(HMMT));

	double **alpha = (double **) malloc(sizeof(double *) * seqLength);
	double **alpha2 = (double **) malloc(sizeof(double *) * (seqLength));
	hmm->AF = (double **) malloc(sizeof(double *) * numStates);

	for (i = 1; i < numStates; i++) {
		hmm->AF[i] = (double *) malloc(sizeof(double) * (numStates * numStates));
	}

	for (i = 1; i < seqLength; i++) {
		alpha[i] = (double *) malloc(sizeof(double) * numStates);
		alpha2[i] = (double *) malloc(sizeof(double) * (numStates * numStates));

	}

	hmm->N = numStates;
	hmm->M = numSymbols;
	/* O = obsSeq */


	ForwardWithScale(hmm, seqLength, O, numLeaf, alpha, double **alpha2, int *P,
			double *scale1, double *scale2, double **phi, double **phi2);
}
