/*
 * Copyright (C) 1998, 2000-2007, 2010, 2011, 2012, 2013 SINTEF Digital,
 * Applied Mathematics, Norway.
 *
 * Contact information: E-mail: tor.dokken@sintef.no                      
 * SINTEF Digital, Department of Mathematics and Cybernetics,                         
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
 * written agreement between you and SINTEF Digital. 
 */

#include "sisl.h"
#include <iostream>


using namespace std;


int main()
{
    int dim = 2;
    int kind = 1;
    double coefs[] = { -1, 1, 0, 0, 1, 0, 1, 1 };
    double knots[] = { 0, 0, 0, 0, 1, 1, 1, 1 };
    int num = 4;
    int order = 4;
    SISLCurve* sc = newCurve(num, order, knots, coefs, kind, dim, 1);

    double itpar;
    int stat;
    s1013(sc, 1.0, 0.01, 0.3, &itpar, &stat);
    double pt[4];
    int kleft;
    s1221(sc, 1, itpar, &kleft, pt, &stat);

    cout << itpar << ' ' << pt[2] << ' ' << pt[3] << endl;
}
