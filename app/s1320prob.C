//===========================================================================
//                                                                           
// File: s1320prob.C                                                         
//                                                                           
// Created: Mon Oct  9 11:56:37 2000                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: s1320prob.C,v 1.1 2001-09-05 14:30:11 afr Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================

#include "sisl.h"
#include <cmath>
#include <stdlib.h>
#include <iostream>
 

using namespace std;


int main()
{
  unsigned n;
  const unsigned dim = 2;
  const unsigned NoPoints = 3;
  const unsigned Order = 3;
  int state, NoIntPoints, NoIntCurves;
  double* IntPointsOnCurve;
  double* IntPointsOnLine;
  SISLIntcurve** IntCurves;
  SISLCurve* Line = NULL;
  SISLCurve* Curve = NULL;
  double endPar;
  double knots[NoPoints+Order];
  double coef[(dim+1)*NoPoints];
  // define line
  double a[] = { 0.0, 4.0};
  double b[] = { 4.0, 4.0};
 
  // generate control polygon

   
  coef[ 0] =  0.0;
  coef[ 1] =  0.0;
  coef[ 2] =  1.0;
 
  coef[ 3] =  2.0 * sqrt(2.0)/2.0;
  coef[ 4] =  2.0 * sqrt(2.0)/2.0;
  coef[ 5] =  1.0 * sqrt(2.0)/2.0;
 
  coef[ 6] =  4.0;
  coef[ 7] =  0.0;
  coef[ 8] =  1.0;


  /*
  coef[ 0] =  0.0;
  coef[ 1] =  0.0;
  coef[ 2] =  0.0;
  coef[ 3] =  1.0;
 
  coef[ 4] =  2.0 * sqrt(2.0)/2.0;
  coef[ 5] =  2.0 * sqrt(2.0)/2.0;
  coef[ 6] =  0.0;
  coef[ 7] =  1.0 * sqrt(2.0)/2.0;
 
  coef[ 8] =  4.0;
  coef[ 9] =  0.0;
  coef[ 10] =  0.0;
 coef[ 11] =  1.0;
  */

  /*
  coef[0] = 0.0;
  coef[1] =  0.0;
  coef[2] =  2.0;
  coef[3] =  2.0;
  coef[4] =  4.0;
  coef[5] =  0.0;
  */

  // generate knots
  knots[0] = 0.0;
  knots[1] = 0.0;
  knots[2] = 0.0;
  knots[3] = 1.0;
  knots[4] = 1.0;
  knots[5] = 1.0;
 
  // ask SISL for a curve.
  Curve = newCurve(NoPoints, Order, knots, coef, 4, dim, 1);
 
  // ask SISL for a line.
  s1602(a, b, Order, dim, 0.0, &endPar, &Line, &state);
 
  if (state!=0)
    return 1;
 
  // least distance between lines
  // Denne linja krasjer
  s1955(Curve, Line, 1E-3, 1E-3, &NoIntPoints, &IntPointsOnCurve,
        &IntPointsOnLine, &NoIntCurves, &IntCurves, &state);
 
  cout << NoIntPoints << ' ' << IntPointsOnCurve[0] << ' ' 
       << IntPointsOnLine[0] << endl;

  // cleanup
  freeCurve(Curve);
  freeCurve(Line);
  free(IntPointsOnCurve);
  free(IntPointsOnLine);
  freeIntcrvlist(IntCurves, NoIntCurves);
}







