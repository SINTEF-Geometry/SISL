/*
 * Copyright (C) 1998, 2000-2007, 2010, 2011, 2012, 2013 SINTEF ICT,
 * Applied Mathematics, Norway.
 *
 * Contact information: E-mail: tor.dokken@sintef.no                      
 * SINTEF ICT, Department of Applied Mathematics,                         
 * P.O. Box 124 Blindern,                                                 
 * 0314 Oslo, Norway.                                                     
 *
 * This file is part of SISL.
 *
 * SISL is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version. 
 *
 * SISL is distributed in the hope that it will be useful,        
 * but WITHOUT ANY WARRANTY; without even the implied warranty of         
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public
 * License along with SISL. If not, see
 * <http://www.gnu.org/licenses/>.
 *
 * In accordance with Section 7(b) of the GNU Affero General Public
 * License, a covered work must retain the producer line in every data
 * file that is created or manipulated using SISL.
 *
 * Other Usage
 * You can be released from the requirements of the license by purchasing
 * a commercial license. Buying such a license is mandatory as soon as you
 * develop commercial activities involving the SISL library without
 * disclosing the source code of your own applications.
 *
 * This file may be used in accordance with the terms contained in a
 * written agreement between you and SINTEF ICT. 
 */

#include "sisl.h"
#include <cmath>
#include <stdlib.h>
#include <iostream>
 

using namespace std;


int main()
{
  // unsigned n;
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







