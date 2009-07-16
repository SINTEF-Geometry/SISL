//===========================================================================
//                                                                           
// File: s1013prob.C                                                         
//                                                                           
// Created: Wed Oct 11 11:20:30 2000                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: s1013prob.C,v 1.1 2001-09-05 14:30:11 afr Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================


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
