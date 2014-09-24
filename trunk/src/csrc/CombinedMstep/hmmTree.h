/*
 **      Author: Tapas Kanungo, kanungo@cfar.umd.edu
 **      Date:   15 December 1997
 **      File:   hmm.h
 **      Purpose: datastructures used for HMM.
 **      Organization: University of Maryland
 **
 **	Update:
 **	Author: Tapas Kanungo
 **	Purpose: include <math.h>. Not including this was
 **		creating a problem with forward.c
 **      $Id: hmm.h,v 1.9 1999/05/02 18:38:11 kanungo Exp kanungo $
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define DELTA 0.001
#define MAXITER 200

typedef struct {
	int N;		/* number of states;  Q={1,2,...,N} */
	int M; 		/* number of observation symbols; V={1,2,...,M}*/
	int dist;       /* probability distribution for Mstep (0 or 1) */
	double	**AF;	/* A[1..N*N][1..N]. a[i][j] is the transition prob
			   of going from state i at time t to state j at time t+1 */
	double **AB;

	double	**B;	/* B[1..N][1..M]. b[j][k] is the probability of
			   of observing symbol k in state j */
	double	**pi;	/* pi[1..N] pi[i] is the initial state distribution.
	 */
	//double *pmshape1;
	//double *pmshape2; /* parameter of beta distribution */
} HMMT;

typedef struct {
	int *bro;
	int numLeaf;
	int *P;

} TreeConfig;

typedef struct {
	double **phi;
	double *phiT;
	double **phi2;
	double *phi2Temp;
	double *scale1;
	double *scale2;
} ForwardConfig;

typedef struct {
	double **theta;
	double **thetaT;
	double *thetaTRow;
	double *scale;
	double *scaleB;
} BackwardConfig;

//Each gene has its own BaumConfig
typedef struct {
	ForwardConfig *forwardConf;
	BackwardConfig *backConf;
	double LL;
	double *plogprobinit;
	double *plogprobfinal;
	double **F;
	double *betaDenom;
	double *betaY2;
	double *betaY1;
	double **B;
	double *pmshape1;
	double *pmshape2;
} BaumConfig;


void ReadHMM(FILE *fp, HMMT *phmm);
void PrintHMM(FILE *fp, HMMT *phmm);
void CopyHMM(HMMT *phmm1, HMMT *phmm2);

double *ReadInputD(FILE *file, int *L);
double *ReadInputMatD(FILE *file, int numRows, int numcols);
int *ReadInputI(FILE *file, int *L);
void InitDelta(double **delta, double *O, int numLeaf);

HMMT InitHMM(double *O, int T, int N, int *P, int numLeaf, BaumConfig *baumConf, TreeConfig *treeConf);

void PrintSequence(FILE *fp, int T, double *O);
void GenSequenceArray(HMMT *phmm, int seed, int T, double *O, int *q);
int GenInitalState(HMMT *phmm);
int GenNextState(HMMT *phmm, int q_t);
int GenSymbol(HMMT *phmm, int q_t);

void FindSiblings(int *B, int *P, int numLeaf, int T);
void ForwardTree(HMMT *phmm, int T, double *O, int numLeaf, double **logalpha, double **logalpha2, BaumConfig *baumConf, TreeConfig *treeConf);
void BackwardTree(HMMT *phmm, int T, double *O, int numLeaf, double **logbeta, double **phi, double *scale, BackwardConfig *conf);
void BaumWelchTree(HMMT *phmm, int T, double **O, int *P, double ***logalpha, double ***logalpha2, double ***logbeta,
		double **gamma, int *pniter, BaumConfig **baumConf, TreeConfig *treeConf, int numGenes, int maxiter);

double *** AllocXi(int T, int N);
void FreeXi(double *** xi, int T, int N);
void ComputeGamma(HMMT *phmm, int T, int numGenes, double ***logalpha, double ***logbeta, double ***logalpha2, int numLeaf,
		double **gamma, BaumConfig **baumConf);
void ComputeXi(HMMT *phmm, int T, double *O, int numLeaf, double **logalpha2, double **logbeta, BaumConfig *baumConf,
		double ***xi);
void AllocateHMM(HMMT *phmm, int T, int N, int M);
void FreeHMM(HMMT *phmm, int T, int N, int M);
void FreeConfigs(BaumConfig *baumConf, int T, int N, int numLeaf);
void AllocateConfigs(BaumConfig *baumConf, TreeConfig *treeConf, int T, int N, int numLeaf, int *P);
void MakeSymmetric(double **three, double ** temp, int row, int col);
void MstepBeta(HMMT *phmm, int T, BaumConfig *baumConf, double **gamma, double *O, int maxiter);
void Mstep(HMMT *phmm, int T, BaumConfig *baumConf, double **gamma, double *O);
void MstepBinom(HMMT *phmm, int T, BaumConfig *baumConf, double **gamma, double *O);
void CombinedEstep(HMMT *phmm, int T, int numGenes, double ***logalpha, double ***logbeta, double ***logalpha2, int numLeaf,
		double **gamma, double ***xi, double LL);
double MaxLL(BaumConfig **baumConf, int G);



void CalcObsProb(HMMT *phmm, double *O, int T, BaumConfig *baumConf);
double ** ExpMatrix(double **mat, int row, int col);
void ExpMatrices(double ** res1, double **res2, double **mat1, double **mat2, int row, int col);

/* random number generator related functions*/


#define MAX(x,y)        ((x) > (y) ? (x) : (y))
#define MIN(x,y)        ((x) < (y) ? (x) : (y))