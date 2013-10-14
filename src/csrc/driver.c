#include "hmmTree.h"
#include "nrutil.h"

int main() {
	int T = 10;
	int N = 3;
	int M = 3;
	int numLeaf = 4;

	HMMT *hmm = (HMMT *) malloc(sizeof(HMMT));
	ForwardConfig *fConf = (ForwardConfig *) malloc(sizeof(ForwardConfig));
	BackwardConfig *bConf = (BackwardConfig *) malloc(sizeof(BackwardConfig));

	double ** logalpha = dmatrix(1,T,1,N);
	double ** logalpha2 = dmatrix(1,T,1,N*N);
	double ** logbeta = dmatrix(1,T,1,N);
	double *LL = (double *) malloc(sizeof(double));

	//ForwardWithScale(hmm, T, O, numLeaf, alpha,**alpha2, fConf);
}

void AllocateHMM(HMMT *phmm, int T, int N, int M) {
	phmm->AF = dmatrix(1,N,1,N);
	phmm->B = dmatrix(1,N,1,M);
	phmm->M = M;
	phmm->N = N;

}

void AllocateConfigs(ForwardConfig *fConf, BackwardConfig *bConf, int T, int N) {
	fConf->phi = dmatrix(1,T,1,N);
	fConf->phi2 = dmatrix(1,T,1,N*N);
	fConf->phiT = dvector(1,N);
	fConf->phi2Temp = dvector(1,N*N);
	fConf->P = ivector(1,T-1);
	fConf->scale1 = dvector(1, T-1);
	fConf->scale2 = dvector(1, T-1);

	bConf->P = ivector(1,T-1);
	bConf->bro = ivector(1,T-1);
	bConf->scale = dvector(1, T-1);
	bConf->thetaT = dmatrix(1,T,1,N);
	bConf->theta = dmatrix(1,T,1,N*N);
}

