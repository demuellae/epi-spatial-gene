#include <stdio.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include "specfunc.h"


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

/******************************************************************************/

void timestamp ( void )

/******************************************************************************/
/*
  Purpose:

    TIMESTAMP prints the current YMDHMS date as a time stamp.

  Example:

    31 May 2001 09:45:54 AM

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    24 September 2003

  Author:

    John Burkardt

  Parameters:

    None
*/
{
# define TIME_SIZE 40

  static char time_buffer[TIME_SIZE];
  const struct tm *tm;
  size_t len;
  time_t now;

  now = time ( NULL );
  tm = localtime ( &now );

  len = strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

  printf ( "%s\n", time_buffer );

  return;
# undef TIME_SIZE
}
/******************************************************************************/

double trigamma ( double x, int *ifault )

/******************************************************************************/
/*
  Purpose:

    TRIGAMMA calculates trigamma(x) = d^2 log(gamma(x)) / dx^2

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 November 2010

  Author:

    Original FORTRAN77 version by BE Schneider.
    C version by John Burkardt.

  Reference:

    BE Schneider,
    Algorithm AS 121:
    Trigamma Function,
    Applied Statistics,
    Volume 27, Number 1, pages 97-99, 1978.

  Parameters:

    Input, double X, the argument of the trigamma function.
    0 < X.

    Output, int *IFAULT, error flag.
    0, no error.
    1, X <= 0.

    Output, double TRIGAMMA, the value of the trigamma function at X.
*/
{
  double a = 0.0001;
  double b = 5.0;
  double b2 =  0.1666666667;
  double b4 = -0.03333333333;
  double b6 =  0.02380952381;
  double b8 = -0.03333333333;
  double value;
  double y;
  double z;
/*
  Check the input.
*/
  if ( x <= 0.0 )
  {
    *ifault = 1;
    value = 0.0;
    return value;
  }

  *ifault = 0;
  z = x;
/*
  Use small value approximation if X <= A.
*/
  if ( x <= a )
  {
    value = 1.0 / x / x;
    return value;
  }
/*
  Increase argument to ( X + I ) >= B.
*/
  value = 0.0;

  while ( z < b )
  {
    value = value + 1.0 / z / z;
    z = z + 1.0;
  }
/*
  Apply asymptotic formula if argument is B or greater.
*/
  y = 1.0 / z / z;

  value = value + 0.5 *
      y + ( 1.0
    + y * ( b2
    + y * ( b4
    + y * ( b6
    + y *   b8 )))) / z;

  return value;
}
/******************************************************************************/

void trigamma_values ( int *n_data, double *x, double *fx )

/******************************************************************************/
/*
  Purpose:

    TRIGAMMA_VALUES returns some values of the TriGamma function.

  Discussion:

    In Mathematica, the function can be evaluated by:

      PolyGamma[1,x]

    TriGamma(X) = d^2 ln ( Gamma ( X ) ) / d X^2

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 September 2004

  Author:

    John Burkardt

  Reference:

    Milton Abramowitz, Irene Stegun,
    Handbook of Mathematical Functions,
    National Bureau of Standards, 1964,
    ISBN: 0-486-61272-4,
    LC: QA47.A34.

    Stephen Wolfram,
    The Mathematica Book,
    Fourth Edition,
    Cambridge University Press, 1999,
    ISBN: 0-521-64314-7,
    LC: QA76.95.W65.

  Parameters:

    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
    first call.  On each call, the routine increments N_DATA by 1, and
    returns the corresponding data; when there is no more data, the
    output value of N_DATA will be 0 again.

    Output, double *X, the argument of the function.

    Output, double *FX, the value of the function.
*/
{
# define N_MAX 11

  double fx_vec[N_MAX] = {
    0.1644934066848226E+01,
    0.1433299150792759E+01,
    0.1267377205423779E+01,
    0.1134253434996619E+01,
    0.1025356590529597E+01,
    0.9348022005446793E+00,
    0.8584318931245799E+00,
    0.7932328301639984E+00,
    0.7369741375017002E+00,
    0.6879720582426356E+00,
    0.6449340668482264E+00 };

  double x_vec[N_MAX] = {
     1.0E+00,
     1.1E+00,
     1.2E+00,
     1.3E+00,
     1.4E+00,
     1.5E+00,
     1.6E+00,
     1.7E+00,
     1.8E+00,
     1.9E+00,
     2.0E+00 };

  if ( *n_data < 0 )
  {
    *n_data = 0;
  }

  *n_data = *n_data + 1;

  if ( N_MAX < *n_data )
  {
    *n_data = 0;
    *x = 0.0;
    *fx = 0.0;
  }
  else
  {
    *x = x_vec[*n_data-1];
    *fx = fx_vec[*n_data-1];
  }

  return;
# undef N_MAX
}

static double cutoff = 171.0;

////////////////////////////////////////////////////////////////////////////////
// double DiGamma_Function( double x )                                        //
//                                                                            //
//  Description:                                                              //
//     This function uses the derivative of the log of Lanczos's expression   //
//     for the Gamma function to calculate the DiGamma function for x > 1 and //
//     x <= cutoff = 171.  An asymptotic expression for the DiGamma function  //
//     for x > cutoff.  The reflection formula                                //
//                      DiGamma(x) = DiGamma(1+x) - 1/x                       //
//     for 0 < x < 1. and the reflection formula                              //
//                DiGamma(x) = DiGamma(1-x) - pi * cot(pi*x)                  //
//     for x < 0.                                                             //
//     The DiGamma function has singularities at the nonpositive integers.    //
//     At a singularity, DBL_MAX is returned.                                 //
//                                                                            //
//  Arguments:                                                                //
//     double x   Argument of the DiGamma function.                           //
//                                                                            //
//  Return Values:                                                            //
//     If x is a nonpositive integer then DBL_MAX is returned otherwise       //
//     DiGamma(x) is returned.                                                //
//                                                                            //
//  Example:                                                                  //
//     double x, psi;                                                         //
//                                                                            //
//     psi = DiGamma_Function( x );                                           //
////////////////////////////////////////////////////////////////////////////////
double DiGamma_Function(double x)
{
   long double psi = xDiGamma_Function((long double) x);

   if (fabsl(psi) < DBL_MAX) return (double) psi;
   return (psi < 0.0L) ? -DBL_MAX : DBL_MAX;
}


////////////////////////////////////////////////////////////////////////////////
// long double xDiGamma_Function( long double x )                             //
//                                                                            //
//  Description:                                                              //
//     This function uses the derivative of the log of Lanczos's expression   //
//     for the Gamma function to calculate the DiGamma function for x > 1 and //
//     x <= cutoff = 171.  An asymptotic expression for the DiGamma function  //
//     for x > cutoff.  The reflection formula                                //
//                      DiGamma(x) = DiGamma(1+x) - 1/x                       //
//     for 0 < x < 1. and the reflection formula                              //
//                DiGamma(x) = DiGamma(1-x) - pi * cot(pi*x)                  //
//     for x < 0.                                                             //
//     The DiGamma function has singularities at the nonpositive integers.    //
//     At a singularity, LDBL_MAX is returned.                                //
//                                                                            //
//  Arguments:                                                                //
//     long double x   Argument of the DiGamma function.                      //
//                                                                            //
//  Return Values:                                                            //
//     If x is a nonpositive integer then LDBL_MAX is returned otherwise      //
//     DiGamma(x) is returned.                                                //
//                                                                            //
//  Example:                                                                  //
//     long double x, psi;                                                    //
//                                                                            //
//     psi = xDiGamma_Function( x );                                          //
////////////////////////////////////////////////////////////////////////////////
long double xDiGamma_Function(long double x)
{
   long double sin_x, cos_x;
   long int ix;

             // For a positive argument (x > 0)                 //
             //    if x <= cutoff return Lanzcos approximation  //
             //    otherwise return Asymptotic approximation.   //

   if ( x > 0.0L )
      if (x <= cutoff)
         if ( x >= 1.0L) return xDiGamma(x);
         else return xDiGamma( x + 1.0L ) - (1.0L / x);
      else return xDiGamma_Asymptotic_Expansion(x);

                  // For a nonpositive argument (x <= 0). //
               // If x is a singularity then return LDBL_MAX. //

   if ( x > -(long double)LONG_MAX) {
      ix = (long int) x;
      if ( x == (long double) ix) return LDBL_MAX;
   }

   sin_x = sinl(pi * x);
   if (sin_x == 0.0L) return LDBL_MAX;
   cos_x = cosl(pi * x);
   if (fabsl(cos_x) == 1.0L) return LDBL_MAX;

            // If x is not a singularity then return DiGamma(x). //

   return xDiGamma(1.0L - x) - pi * cos_x / sin_x;
}


////////////////////////////////////////////////////////////////////////////////
// static long double xDiGamma( long double x )                               //
//                                                                            //
//  Description:                                                              //
//     This function uses the derivative of the log of Lanczos's expression   //
//     for the Gamma function to calculate the DiGamma function for x > 1 and //
//     x <= cutoff = 171.                                                     //
//                                                                            //
//  Arguments:                                                                //
//     double x   Argument of the Gamma function.                             //
//                                                                            //
//  Return Values:                                                            //
//     DiGamma(x)                                                             //
//                                                                            //
//  Example:                                                                  //
//     long double x;                                                         //
//     long double psi;                                                       //
//                                                                            //
//     g = xDiGamma( x );                                                     //
////////////////////////////////////////////////////////////////////////////////
static long double xDiGamma(long double x) {

   long double lnarg = (g + x - 0.5L);
   long double temp;
   long double term;
   long double numerator = 0.0L;
   long double denominator = 0.0L;
   int const k = sizeof(a) / sizeof(long double);
   int i;

   for (i = k-1; i >= 0; i--) {
      temp = x + (long double) i;
      term = a[i] / temp;
      denominator += term;
      numerator += term / temp;
   }
   denominator += 1.0L;

   return logl(lnarg) - (g / lnarg) - numerator / denominator;
}


////////////////////////////////////////////////////////////////////////////////
// static long double xDiGamma_Asymptotic_Expansion( long double x )          //
//                                                                            //
//  Description:                                                              //
//     This function estimates DiGamma(x) by evaluating the asymptotic        //
//     expression:                                                            //
//         DiGamma(x) ~ ln(x) - (1/2) x +                                     //
//                        Sum B[2j] / [ 2j * x^(2j) ], summed over            //
//     j from 1 to 3 and where B[2j] is the 2j-th Bernoulli number.           //
//                                                                            //
//  Arguments:                                                                //
//     long double x   Argument of the DiGamma function. The argument x must  //
//                     be positive.                                           //
//                                                                            //
//  Return Values:                                                            //
//     DiGamma(x) where x > cutoff.                                           //
//                                                                            //
//  Example:                                                                  //
//     long double x;                                                         //
//     long double g;                                                         //
//                                                                            //
//     g = xDiGamma_Asymptotic_Expansion( x );                                //
////////////////////////////////////////////////////////////////////////////////

// Bernoulli numbers B(2j) / 2j: B(2)/2,B(4)/4,B(6)/6,...,B(20)/20.  Only     //
//  B(2)/2,..., B(6)/6 are currently used.                                    //



static long double xDiGamma_Asymptotic_Expansion(long double x ) {
   const int  m = 3;
   long double term[3];
   long double sum = 0.0L;
   long double xx = x * x;
   long double xj = x;
   long double digamma = logl(xj) - 1.0L / (xj + xj);
   int i;

   xj = xx;
   for (i = 0; i < m; i++) { term[i] = B[i] / xj; xj *= xx; }
   for (i = m - 1; i >= 0; i--) sum += term[i];
   return digamma - sum;
}

////////////////////////////////////////////////////////////////////////////////
// double Gamma_Function( double x )                                          //
//                                                                            //
//  Description:                                                              //
//     This function uses Lanczos' expression to calculate Gamma(x) for real  //
//     x, where -(max_double_arg - 1) < x < max_double_arg.                   //
//     Note the Gamma function is meromorphic in the complex plane and has    //
//     poles at the nonpositive integers.                                     //
//     Tests for x a positive integer or a half positive integer give a       //
//     maximum absolute relative error of about 1.9e-16.                      //
//                                                                            //
//     If x > max_double_arg, then one should either use xGamma_Function(x)   //
//     or calculate lnGamma(x).                                               //
//     Note that for x < 0, ln (Gamma(x)) may be a complex number.            //
//                                                                            //
//  Arguments:                                                                //
//     double x   Argument of the Gamma function.                             //
//                                                                            //
//  Return Values:                                                            //
//     If x is positive and is less than max_double_arg then Gamma(x) is      //
//     returned and if x > max_double_arg then DBL_MAX is returned.  If x is  //
//     a nonpositive integer i.e. x is a pole, then DBL_MAX is returned       //
//     ( note that Gamma(x) changes sign on each side of the pole).  If x is  //
//     nonpositive nonintegral, then if Gamma(x) > DBL_MAX, then DBL_MAX is   //
//     returned and if Gamma(x) < -DBL_MAX, then -DBL_MAX is returned.        //
//                                                                            //
//  Example:                                                                  //
//     double x, g;                                                           //
//                                                                            //
//     g = Gamma_Function( x );                                               //
////////////////////////////////////////////////////////////////////////////////
double Gamma_Function(double x)
{
   long double g;

   if ( x > max_double_arg ) return DBL_MAX;
   g = xGamma_Function( (long double) x);
   if (fabsl(g) < DBL_MAX) return (double) g;
   return (g < 0.0L) ? -DBL_MAX : DBL_MAX;

}


////////////////////////////////////////////////////////////////////////////////
// long double xGamma_Function( long double x )                               //
//                                                                            //
//  Description:                                                              //
//     This function uses Lanczos' expression to calculate Gamma(x) for real  //
//     x, where -(max_long_double_arg - 1) < x < max_long_double_arg.         //
//     Note the Gamma function is meromorphic in the complex plane and has    //
//     poles at the nonpositive integers.                                     //
//     Tests for x a positive integer or a half positive integer give a       //
//     maximum absolute relative error of about 3.5e-16.                      //
//                                                                            //
//     If x > max_long_double_arg, then one should use lnGamma(x).            //
//     Note that for x < 0, ln (Gamma(x)) may be a complex number.            //
//                                                                            //
//  Arguments:                                                                //
//     long double x   Argument of the Gamma function.                        //
//                                                                            //
//  Return Values:                                                            //
//     If x is positive and is less than max_long_double_arg, then Gamma(x)   //
//     is returned and if x > max_long_double_arg, then LDBL_MAX is returned. //
//     If x is a nonpositive integer i.e. x is a pole, then LDBL_MAX is       //
//     returned ( note that Gamma(x) changes sign on each side of the pole).  //
//     If x is nonpositive nonintegral, then if x > -(max_long_double_arg + 1)//
//     then Gamma(x) is returned otherwise 0.0 is returned.                   //
//                                                                            //
//  Example:                                                                  //
//     long double x, g;                                                      //
//                                                                            //
//     g = xGamma_Function( x );                                              //
////////////////////////////////////////////////////////////////////////////////
long double xGamma_Function(long double x)
{
   long double sin_x;
   long double rg;
   long int ix;

             // For a positive argument (x > 0)                 //
             //    if x <= max_long_double_arg return Gamma(x)  //
             //    otherwise      return LDBL_MAX.              //

   if ( x > 0.0L )
      if (x <= max_long_double_arg) return xGamma(x);
      else return LDBL_MAX;

                   // For a nonpositive argument (x <= 0) //
                   //    if x is a pole return LDBL_MAX   //

   if ( x > -(long double)LONG_MAX) {
      ix = (long int) x;
      if ( x == (long double) ix) return LDBL_MAX;
   }
   sin_x = sinl(pi * x);
   if ( sin_x == 0.0L ) return LDBL_MAX;

          // if x is not a pole and x < -(max_long_double_arg - 1) //
          //                                     then return 0.0L  //

   if ( x < -(max_long_double_arg - 1.0L) ) return 0.0L;

            // if x is not a pole and x >= -(max_long_double - 1) //
            //                               then return Gamma(x) //

   rg = xGamma(1.0L - x) * sin_x / pi;
   if ( rg != 0.0L ) return (1.0L / rg);
   return LDBL_MAX;
}


////////////////////////////////////////////////////////////////////////////////
// static long double xGamma( long double x )                                 //
//                                                                            //
//  Description:                                                              //
//     This function uses Lanczos' expression to calculate Gamma(x) for real  //
//     x, where 0 < x <= 900. For 900 < x < 1755.5, the duplication formula   //
//     is used.                                                               //
//     The major source of relative error is in the use of the c library      //
//     function powl().  The results have a relative error of about 10^-16.   //
//     except near x = 0.                                                     //
//                                                                            //
//     If x > 1755.5, then one should calculate lnGamma(x).                   //
//                                                                            //
//  Arguments:                                                                //
//     long double x   Argument of the Gamma function.                        //
//                                                                            //
//  Return Values:                                                            //
//     If x is positive and is less than 1755.5 then Gamma(x) is returned and //
//     if x > 1755.5 then LDBL_MAX is returned.                               //
//                                                                            //
//  Example:                                                                  //
//     long double x;                                                         //
//     long double g;                                                         //
//                                                                            //
//     g = xGamma_Function( x );                                              //
////////////////////////////////////////////////////////////////////////////////
static long double xGamma(long double x)
{

   long double xx = (x < 1.0L) ? x + 1.0L : x;
   long double temp;
   int const n = sizeof(a) / sizeof(long double);
   int i;

   if (x > 1755.5L) return LDBL_MAX;

   if (x > 900.0L) return Duplication_Formula(x);

   temp = 0.0L;
   for (i = n-1; i >= 0; i--) {
      temp += ( a[i] / (xx + (long double) i) );
   }
   temp += 1.0L;
   temp *= ( powl((g + xx - 0.5L) / e, xx - 0.5L) / exp_g_o_sqrt_2pi );
   return (x < 1.0L) ?  temp / x : temp;
}


////////////////////////////////////////////////////////////////////////////////
// static long double Duplication_Formula(long double two_x)                  //
//                                                                            //
//  Description:                                                              //
//     This function returns the Gamma(two_x) using the duplication formula   //
//     Gamma(2x) = (2^(2x-1) / sqrt(pi)) Gamma(x) Gamma(x+1/2).               //
//                                                                            //
//  Arguments:                                                                //
//     none                                                                   //
//                                                                            //
//  Return Values:                                                            //
//     Gamma(two_x)                                                           //
//                                                                            //
//  Example:                                                                  //
//     long double two_x, g;                                                  //
//                                                                            //
//     g = Duplication_Formula(two_x);                                        //
////////////////////////////////////////////////////////////////////////////////
static long double Duplication_Formula( long double two_x )
{
   long double x = 0.5L * two_x;
   long double g;
   double two_n = 1.0;
   int n = (int) two_x - 1;

   g = powl(2.0L, two_x - 1.0L - (long double) n);
   g = ldexpl(g,n);
   g /= sqrt(pi);
   g *= xGamma_Function(x);
   g *= xGamma_Function(x + 0.5L);

   return g;
}


////////////////////////////////////////////////////////////////////////////////
// double Gamma_Function_Max_Arg( void )                                      //
//                                                                            //
//  Description:                                                              //
//     This function returns the maximum argument of Gamma_Function for which //
//     a number < DBL_MAX is returned, for arguments greater than 1.          //
//                                                                            //
//  Arguments:                                                                //
//     none                                                                   //
//                                                                            //
//  Return Values:                                                            //
//     max_double_arg (171.0).                                                //
//                                                                            //
//  Example:                                                                  //
//     double x;                                                              //
//                                                                            //
//     x = Gamma_Function_Max_Arg();                                          //
////////////////////////////////////////////////////////////////////////////////
double Gamma_Function_Max_Arg( void ) { return max_double_arg; }


////////////////////////////////////////////////////////////////////////////////
// long double xGamma_Function_Max_Arg( void )                                //
//                                                                            //
//  Description:                                                              //
//     This function returns the maximum argument of Gamma_Function for which //
//     a number < LDBL_MAX is returned, for arguments greater than 1.         //
//                                                                            //
//  Arguments:                                                                //
//     none                                                                   //
//                                                                            //
//  Return Values:                                                            //
//     max_long_double_arg (1755.5).                                          //
//                                                                            //
//  Example:                                                                  //
//     long double x;                                                         //
//                                                                            //
//     x = xGamma_Function_Max_Arg();                                         //
////////////////////////////////////////////////////////////////////////////////
long double xGamma_Function_Max_Arg( void ) { return max_long_double_arg; }

////////////////////////////////////////////////////////////////////////////////
// double Ln_Gamma_Function( double x )                                       //
//                                                                            //
//  Description:                                                              //
//     This function calculates the natural log of Gamma(x) for positive real //
//     x.                                                                     //
//     Assuming that Gamma_Function_Max_Arg() = 171, then                     //
//     If 0 < x <= 171, then ln(gamma(x)) is calculated by taking the natural //
//     log of the result from Gamma_Function(x).  If x > 171, then            //
//     ln(gamma(x)) is calculated using the asymptotic expansion              //
//         ln(gamma(x)) ~ ln(2sqrt(2pi)) - x + (x - 1/2) ln x +               //
//                        Sum B[2j] / [ 2j * (2j-1) * x^(2j-1) ], summed over //
//     j from 1 to 3 and where B[2j] is the 2j-th Bernoulli number.           //
//                                                                            //
//  Arguments:                                                                //
//     double x   Argument of the ln Gamma function. The argument x must be   //
//                positive.                                                   //
//                                                                            //
//  Return Values:                                                            //
//     ln(Gamma(x)) where x > 0.                                              //
//                                                                            //
//  Example:                                                                  //
//     double x, g;                                                           //
//                                                                            //
//     g = Ln_Gamma_Function( x );                                            //
////////////////////////////////////////////////////////////////////////////////

double Ln_Gamma_Function(double x)
{

       // For a positive argument, 0 < x <= Gamma_Function_Max_Arg() //
       // then  return log Gamma(x).                                 //

   if (x <= Gamma_Function_Max_Arg()) return log(Gamma_Function(x));

    // otherwise return result from asymptotic expansion of ln Gamma(x). //

   return (double) xLnGamma_Asymptotic_Expansion( (long double) x );
}


////////////////////////////////////////////////////////////////////////////////
// long double xLn_Gamma_Function( long double x )                            //
//                                                                            //
//  Description:                                                              //
//     This function calculates the natural log of Gamma(x) for positive real //
//     x.                                                                     //
//     Assuming that Gamma_Function_Max_Arg() = 171, then                     //
//     If 0 < x <= 171, then ln(gamma(x)) is calculated by taking the natural //
//     log of the result from Gamma_Function(x).  If x > 171, then            //
//     ln(gamma(x)) is calculated using the asymptotic expansion              //
//         ln(gamma(x)) ~ ln(2sqrt(2pi)) - x + (x - 1/2) ln x +               //
//                        Sum B[2j] / [ 2j * (2j-1) * x^(2j-1) ], summed over //
//     j from 1 to 3 and where B[2j] is the 2j-th Bernoulli number.           //
//                                                                            //
//  Arguments:                                                                //
//     long double x   Argument of the ln Gamma function. The argument x must //
//                     be positive.                                           //
//                                                                            //
//  Return Values:                                                            //
//     ln(Gamma(x)) where x > 0.                                              //
//                                                                            //
//  Example:                                                                  //
//     double x;                                                              //
//     long double g;                                                         //
//                                                                            //
//     g = xLn_Gamma_Function( x );                                           //
////////////////////////////////////////////////////////////////////////////////
long double xLn_Gamma_Function(long double x)
{

       // For a positive argument, 0 < x <= Gamma_Function_Max_Arg() //
       // then  return log Gamma(x).                                 //

   if (x <= Gamma_Function_Max_Arg()) return logl(xGamma_Function(x));

    // otherwise return result from asymptotic expansion of ln Gamma(x). //

   return xLnGamma_Asymptotic_Expansion( x );
}


////////////////////////////////////////////////////////////////////////////////
// static long double xLnGamma_Asymptotic_Expansion( long double x )          //
//                                                                            //
//  Description:                                                              //
//     This function estimates log(gamma(x)) by evaluating the asymptotic     //
//     expression:                                                            //
//         ln(Gamma(x)) ~ ln(2sqrt(2pi)) - x + (x - 1/2) ln x +               //
//                        Sum B[2j] / [ 2j * (2j-1) * x^(2j-1) ], summed over //
//     j from 1 to 3 and where B[2j] is the 2j-th Bernoulli number.           //
//                                                                            //
//  Arguments:                                                                //
//     long double x   Argument of the ln Gamma function. The argument x must //
//                     be  positive.                                          //
//                                                                            //
//  Return Values:                                                            //
//     ln(Gamma(x)) where x > Gamma_Function_Max_Arg()                        //
//                                                                            //
//  Example:                                                                  //
//     double x;                                                              //
//     long double g;                                                         //
//                                                                            //
//     g = xlnGamma_Asymptotic_Expansion( x );                                //
////////////////////////////////////////////////////////////////////////////////

static const long double log_sqrt_2pi = 9.18938533204672741780329736e-1L;

// Bernoulli numbers B(2),B(4),B(6),...,B(20).  Only B(2),...,B(6) currently //
// used.                                                                     //




static long double xLnGamma_Asymptotic_Expansion( long double x ) {
   const int  m = 3;
   long double term[3];
   long double sum = 0.0L;
   long double xx = x * x;
   long double xj = x;
   long double lngamma = log_sqrt_2pi - xj + (xj - 0.5L) * logl(xj);
   int i;

   for (i = 0; i < m; i++) { term[i] = B[i] / xj; xj *= xx; }
   for (i = m - 1; i >= 0; i--) sum += term[i];
   return lngamma + sum;
}

////////////////////////////////////////////////////////////////////////////////
// double Beta_Function( double a, double b)                                  //
//                                                                            //
//  Description:                                                              //
//     This function returns beta(a,b) = gamma(a) * gamma(b) / gamma(a+b),    //
//     where a > 0, b > 0.                                                    //
//                                                                            //
//  Arguments:                                                                //
//     double a   Argument of the Beta function, a must be positive.          //
//     double b   Argument of the Beta function, b must be positive.          //
//                                                                            //
//  Return Values:                                                            //
//     If beta(a,b) exceeds DBL_MAX then DBL_MAX is returned otherwise        //
//     beta(a,b) is returned.                                                 //
//                                                                            //
//  Example:                                                                  //
//     double a, b, beta;                                                     //
//                                                                            //
//     beta = Beta_Function( a, b );                                          //
////////////////////////////////////////////////////////////////////////////////
double Beta_Function(double a, double b)
{
   long double beta = xBeta_Function( (long double) a, (long double) b);
   return (beta < DBL_MAX) ? (double) beta : DBL_MAX;
}


////////////////////////////////////////////////////////////////////////////////
// long double xBeta_Function( long double a, long double b)                  //
//                                                                            //
//  Description:                                                              //
//     This function returns beta(a,b) = gamma(a) * gamma(b) / gamma(a+b),    //
//     where a > 0, b > 0.                                                    //
//                                                                            //
//  Arguments:                                                                //
//     long double a   Argument of the Beta function, a must be positive.     //
//     long double b   Argument of the Beta function, b must be positive.     //
//                                                                            //
//  Return Values:                                                            //
//     If beta(a,b) exceeds LDBL_MAX then LDBL_MAX is returned otherwise      //
//     beta(a,b) is returned.                                                 //
//                                                                            //
//  Example:                                                                  //
//     long double a, b;                                                      //
//     long double beta;                                                      //
//                                                                            //
//     beta = xBeta_Function( a, b );                                         //
////////////////////////////////////////////////////////////////////////////////
long double xBeta_Function(long double a, long double b)
{
   long double lnbeta;

     // If (a + b) <= Gamma_Function_Max_Arg() then simply return //
     //  gamma(a)*gamma(b) / gamma(a+b).                          //

   if ( (a + b) <= Gamma_Function_Max_Arg() )
      return xGamma_Function(a) / (xGamma_Function(a + b) / xGamma_Function(b));

     // If (a + b) > Gamma_Function_Max_Arg() then simply return //
     //  exp(lngamma(a) + lngamma(b) - lngamma(a+b) ).           //

   lnbeta = xLn_Gamma_Function(a) + xLn_Gamma_Function(b)
                                                 - xLn_Gamma_Function(a + b);
   return (lnbeta > ln_LDBL_MAX) ? (long double) LDBL_MAX : expl(lnbeta);
}

