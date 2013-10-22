/*
 * specfunc.h
 *
 *  Created on: Oct 15, 2013
 *      Author: Dan
 */

#ifndef SPECFUNC_H_
#define SPECFUNC_H_

/* Beta_Function.c */
double Beta_Function(double a, double b);
long double xBeta_Function(long double a, long double b);

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

static long double const e =  2.71828182845904523536028747L;
static long double const pi = 3.14159265358979323846264338L;
static long double const g =  9.65657815377331589457187L;
static long double const exp_g_o_sqrt_2pi = +6.23316569877722552586386e+3L;
static double max_double_arg = 171.0;
static long double max_long_double_arg = 1755.5L;

static const long double ln_LDBL_MAX =  1.13565234062941435e+4L;


static long double const a[] = {
                                 +1.14400529453851095667309e+4L,
                                 -3.23988020152318335053598e+4L,
                                 +3.50514523505571666566083e+4L,
                                 -1.81641309541260702610647e+4L,
                                 +4.63232990536666818409138e+3L,
                                 -5.36976777703356780555748e+2L,
                                 +2.28754473395181007645155e+1L,
                                 -2.17925748738865115560082e-1L,
                                 +1.08314836272589368860689e-4L
                              };

static const long double B[] = {   1.0L / (long double)(6 * 2 * 1),
                                  -1.0L / (long double)(30 * 4 * 3),
                                   1.0L / (long double)(42 * 6 * 5),
                                  -1.0L / (long double)(30 * 8 * 7),
                                   5.0L / (long double)(66 * 10 * 9),
                                -691.0L / (long double)(2730 * 12 * 11),
                                   7.0L / (long double)(6 * 14 * 13),
                               -3617.0L / (long double)(510 * 16 * 15),
                               43867.0L / (long double)(796 * 18 * 17),
                             -174611.0L / (long double)(330 * 20 * 19)
                           };

static const int n = sizeof(B) / sizeof(long double);



#endif /* SPECFUNC_H_ */
