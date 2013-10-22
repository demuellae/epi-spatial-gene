////////////////////////////////////////////////////////////////////////////////
// File: ln_gamma_function.c                                                  //
// Routine(s):                                                                //
//    Ln_Gamma_Function                                                       //
//    xLn_Gamma_Function                                                      //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//  Description:                                                              //
//     These functions, Ln_Gamma_Function(x) and xLn_Gamma_Function(x),       //
//      calculate the natural log of Gamma(x) for positive real x.            //
////////////////////////////////////////////////////////////////////////////////

#include <math.h>                  // required for log(), logl() and sqrtl()
#include "specfunc.h"

//                         Externally Defined Routines                        //

//                         Internally Defined Routines                        //



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
