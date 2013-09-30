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

typedef struct {
	int N;		/* number of states;  Q={1,2,...,N} */
	int M; 		/* number of observation symbols; V={1,2,...,M}*/
	double	**AF;	/* A[1..N][1..N]. a[i][j] is the transition prob
			   of going from state i at time t to state j at time t+1 */
	double **AB;

	double	**B;	/* B[1..N][1..M]. b[j][k] is the probability of
			   of observing symbol k in state j */
	double	**pi;	/* pi[1..N] pi[i] is the initial state distribution.
	 */
} HMMT;

typedef struct {
  double **phi;
  double *phiT;
  double **phi2;
  double *phi2Temp;
  double *scale1;
  double *scale2;
  int *P;
} ForwardConfig;

typedef struct {
  double **theta;
  double **thetaT;
  double *thetaTRow;
  double *bro;
  double *scale;
  int *P;
} BackwardConfig;

void InitConfForward(ForwardConfig *f);
void InitConfBackward(BackwardConfig *b);

void ReadHMM(FILE *fp, HMMT *phmm);
void PrintHMM(FILE *fp, HMMT *phmm);
void InitHMM(HMMT *phmm, int N, int M, int seed);
void CopyHMM(HMMT *phmm1, HMMT *phmm2);
void FreeHMM(HMMT *phmm);

void ReadSequence(FILE *fp, int *pT, int **pO);
void PrintSequence(FILE *fp, int T, int *O);
void GenSequenceArray(HMMT *phmm, int seed, int T, int *O, int *q);
int GenInitalState(HMMT *phmm);
int GenNextState(HMMT *phmm, int q_t);
int GenSymbol(HMMT *phmm, int q_t);


void FindSiblings(int *B, int *P, int numLeaf);
void ForwardTree(HMMT *phmm, int T, int *O, int numLeaf, double **logalpha, double **logalpha2, double *LL,
	      ForwardConfig *conf);
void BackwardTree(HMMT *phmm, int T, int *O, int numLeaf, double **logbeta, double **phi, BackwardConfig **conf);
void Backward(HMMT *phmm, int T, int *O, double **beta, double *pprob);
void BackwardWithScale(HMMT *phmm, int T, int *O, double **beta,
        double *scale, double *pprob);
void BaumWelch(HMMT *phmm, int T, int *O, double **alpha, double **beta,
        double **gamma, int *niter,
	double *plogprobinit, double *plogprobfinal);

double *** AllocXi(int T, int N);
void FreeXi(double *** xi, int T, int N);
void ComputeGamma(HMMT *phmm, int T, double **alpha, double **beta,
        double **gamma);
void ComputeXi(HMMT* phmm, int T, int *O, double **alpha, double **beta,
        double ***xi);
void Viterbi(HMMT *phmm, int T, int *O, double **delta, int **psi,
        int *q, double *pprob);
void ViterbiLog(HMMT *phmm, int T, int *O, double **delta, int **psi,
        int *q, double *pprob);

/* random number generator related functions*/

int hmmgetseed(void);
void hmmsetseed(int seed);
double hmmgetrand(void);

#define MAX(x,y)        ((x) > (y) ? (x) : (y))
#define MIN(x,y)        ((x) < (y) ? (x) : (y))
