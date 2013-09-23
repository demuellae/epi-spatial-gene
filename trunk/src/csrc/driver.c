#include "hmmTree.h"

int main() {
	int seqLength = 10;
	int numStates = 3;
	int numSymbols = 3;
	int numLeaf = 4;

	/* indices for malloc of 2-d arrays */
	int i, j;

	HMMT *hmm = malloc(sizeof(HMMT));
	ForwardConfig *conf = malloc(sizeof(ForwardConfig));


	/* Initial memory allocation */
	double **alpha = (double **) malloc(sizeof(double *) * seqLength);
	double **alpha2 = (double **) malloc(sizeof(double *) * (seqLength));
	hmm->AF = (double **) malloc(sizeof(double *) * numStates);
	hmm->N = numStates;
	hmm->M = seqLength;
	hmm->pi = malloc(sizeof(double) * numStates);
	hmm->B = malloc(sizeof(double) * seqLength);
	conf->phi = (double **) malloc(sizeof(double *) * seqLength);
	conf->phi2 = (double **) malloc(sizeof(double *) * seqLength);


	/* Initialize 2D Arrays with as many rows as number of states */
	for (i = 1; i < numStates; i++) {
		hmm->AF[i] = (double *) malloc(sizeof(double) * (numStates * numStates));
	}

	/* Initialize 2D Arrays with as many rows as sequence length */
	for (i = 1; i < seqLength; i++) {

		/* Arrays with as many columns as number of states */
		alpha[i] = (double *) malloc(sizeof(double) * numStates);
		conf->phi = (double *) malloc(sizeof(double) * (numStates));

		/*Arrays with as many columns as number of states squared */
		alpha2[i] = (double *) malloc(sizeof(double) * (numStates * numStates));
		hmm->B = (double *) malloc(sizeof(double) * (numStates * numStates));
		conf->phi2 = (double *) malloc(sizeof(double) * (numStates * numStates));



	}

	hmm->N = numStates;
	hmm->M = numSymbols;



	//ForwardWithScale(hmm, seqLength, O, numLeaf, alpha,**alpha2, conf);
}
