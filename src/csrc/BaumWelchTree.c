#include "hmmTree.h"
#include <stdio.h>

void BaumWelch(HMM *phmm, int T, int *O, double **alpha, double **beta,
	double **gamma, int *pniter,
	double *plogprobinit, double *plogprobfinal)
{
	int	i, j, k;
	int	t, l = 0;

	double	logprobf, logprobb,  threshold;
	double	numeratorA, denominatorA;
	double	numeratorB, denominatorB;

	double ***xi, *scale;
	double delta, deltaprev, logprobprev;

	deltaprev = 10e-70;

	xi = AllocXi(T, phmm->N);
	scale = dvector(1, T);

	ForwardWithScale(phmm, T, O, alpha, scale, &logprobf);
	*plogprobinit = logprobf; /* log P(O |intial model) */
	BackwardWithScale(phmm, T, O, beta, scale, &logprobb);
	ComputeGamma(phmm, T, alpha, beta, gamma);
	ComputeXi(phmm, T, O, alpha, beta, xi);
	logprobprev = logprobf;

	do  {

		/* reestimate frequency of state i in time t=1 */
		for (i = 1; i <= phmm->N; i++)
			phmm->pi[i] = .001 + .999*gamma[1][i];

		/* reestimate transition matrix  and symbol prob in
		   each state */
		for (i = 1; i <= phmm->N; i++) {
			denominatorA = 0.0;
			for (t = 1; t <= T - 1; t++)
				denominatorA += gamma[t][i];

			for (j = 1; j <= phmm->N; j++) {
				numeratorA = 0.0;
				for (t = 1; t <= T - 1; t++)
					numeratorA += xi[t][i][j];
				phmm->A[i][j] = .001 +
						.999*numeratorA/denominatorA;
			}

			denominatorB = denominatorA + gamma[T][i];
			for (k = 1; k <= phmm->M; k++) {
				numeratorB = 0.0;
				for (t = 1; t <= T; t++) {
					if (O[t] == k)
						numeratorB += gamma[t][i];
				}

				phmm->B[i][k] = .001 +
						.999*numeratorB/denominatorB;
			}
		}
