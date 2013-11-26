
#ifndef SPECFUNC_H_
#define SPECFUNC_H_

/* Beta_Function.c */
double Beta_Function(double a, double b);
long double xBeta_Function(long double a, long double b);

void timestamp ( void );
double trigamma ( double x, int *ifault );
void trigamma_values ( int *n_data, double *x, double *fx );

/* beta_distribution.c */
double Beta_Distribution(double x, double a, double b);

static long double Beta_Continued_Fraction( long double x, long double a,
                                                               long double b);
static long double xBeta_Distribution(double x, double a, double b);

/* digamma_function.c */
double DiGamma_Function( double x );
long double xDiGamma_Function( long double x );

static long double xDiGamma(long double x);
static long double xDiGamma_Asymptotic_Expansion( long double x );

/* ln_gamma_function .c */
double Ln_Gamma_Function(double x);
long double xLn_Gamma_Function(long double x);

static long double xLnGamma_Asymptotic_Expansion( long double x );

/* gamma_function.c */
double Gamma_Function(double x);
long double xGamma_Function(long double x);
double Gamma_Function_Max_Arg( void );
long double xGamma_Function_Max_Arg( void );

static long double xGamma(long double x);
static long double Duplication_Formula( long double two_x );





#endif /* SPECFUNC_H_ */
