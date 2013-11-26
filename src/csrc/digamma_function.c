////////////////////////////////////////////////////////////////////////////////
// File: digamma_function.c                                                   //
// Routine(s):                                                                //
//    DiGamma_Function                                                        //
//    xDiGamma_Function                                                       //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//  Description:                                                              //
//     The digamma function, also called the psi function, evaluated at a     //
//     x is the derivative of the log of the gamma function evaluated at x,   //
//     i.e.  psi(x) = d ln gamma(x) / dx = (1 / gamma(x)) d gamma(x) / dx.    //
////////////////////////////////////////////////////////////////////////////////
#include <float.h>                 // required for DBL_MAX and LDBL_MAX.
#include <limits.h>     // required for LONG_MAX
#include "specfunc.h"
#include <math.h>


//                         Internally Defined Routines                        //




//                         Internally Defined Constants                       //

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
