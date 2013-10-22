////////////////////////////////////////////////////////////////////////////////
// File: beta_density.c                                                       //
// Routine(s):                                                                //
//    Beta_Density                                                            //
////////////////////////////////////////////////////////////////////////////////

#include <math.h>
#include "specfunc.h"// required for pow()

//                         Externally Defined Routines                        //

extern double Beta_Function(double a, double b);

////////////////////////////////////////////////////////////////////////////////
// double Beta_Density( double x, double a, double b)                         //
//                                                                            //
//  Description:                                                              //
//     The density of the beta distribution is                                //
//                               0                if x < 0,                   //
//                 x^(a-1) (1-x)^(b-1) / B(a,b)   if 0 <= x <= 1,             //
//                               0                if x > 1,                   //
//     where a > 0, b > 0, and B(a,b) is the (complete) beta function.        //
//                                                                            //
//  Arguments:                                                                //
//     double x                                                               //
//                Argument of the density function.  If x <= 0, the result    //
//                is 0 and if x >= 1, the result is 0, otherwise the result   //
//                is the function above except if x = 0 when a < 1 or x = 1   //
//                when b < 1 in which case DBL_MAX is returned.               //
//     double a                                                               //
//                A positive shape parameter of the beta distriubtion,        //
//                a - 1 is the exponent of the factor x in the integrand.     //
//     double b                                                               //
//                A positive shape parameter of the beta distribution,        //
//                b - 1 is the exponent of the factor (1-x) in the integrand. //
//                                                                            //
//  Return Values:                                                            //
//     If x <= 0 or x >= 1, then 0 is returned.  If 0 < x < 1 then            //
//     x^(a-1) (1-x)^(b-1) / B(a,b) is returned.                              //
//                                                                            //
//  Example:                                                                  //
//     double a, b, p, x;                                                     //
//                                                                            //
//     p = Beta_Density(x, a, b);                                             //
////////////////////////////////////////////////////////////////////////////////

double Beta_Density(double x, double a, double b)
{
   if ( x <= 0.0 ) return 0.0;
   if ( x >= 1.0 ) return 0.0;

   return pow(x, a - 1.0) * pow(1.0 - x, b - 1.0) / Beta_Function(a,b);
}
