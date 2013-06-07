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

#include "sisl-copyright.h"

/*
 *
 * $Id: shchecktyp.c,v 1.2 2001-03-19 16:06:04 afr Exp $
 *
 */


#define SHCHECKTYPE

#include "sislP.h"                            


#if defined (SISLNEEDPROTOTYPES)
int
    shchecktype(SISLObject *pobj,double *parval)
#else
int shchecktype(pobj,parval)

     SISLObject *pobj;
     double *parval;
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : To classify (in the extremal sense) a onedimensional 
*              object's behaviour at a given parameter value.
*
*
*
* INPUT      : pobj     - Pointer to a surface or curve object.
*              parval   - Parameter value(s).
*
*
* OUTPUT     : shchecktype :
*                            0 - No extremal point
*                            1 - Maximum point
*                            2 - Minimum point
*                            3 - Saddle point (used in surface only)
*                            4 - Extremum, test inconclusive.
*                          < 0 - error
*
*
* METHOD     : Test are performed on the first and second derivatives.
*
*
* REFERENCES :
*
*-
* CALLS      : s1221      - Evaluate curve.
*              s1421      - Evaluate surface.
*
* WRITTEN BY :Ulf J. Krystad, January 1991.
*
*********************************************************************
*/                                     
{
   
   int kstat;              /* Local status variable.            */
   int kleft1 = 0;         /* Knot navigator                    */
   int kleft2 = 0;         /* Knot navigator                    */
   int kder   = 2;         /* Flag, compute 2. derivative.      */
   int kdim   = 1;         /* Dimension is one !                */
   double sval[9];         /* Storing uppdated parametervalues. */
   double sval1[9];        /* Storing uppdated parametervalues. */
   double *snorm = sval+6; /* Dummy normal pointer              */
   double tmax;        	   /* Size of derivatives               */
   double tdet;            /* Size of Hessian determinant       */
   int mult = 0;           /* Knot multiplicity                 */
   double ttol = 1000000.0*REL_COMP_RES;
   /* --------------------------------------------------------- */
   
   if (pobj == SISL_NULL ||
       (pobj->iobj != SISLCURVE && pobj->iobj != SISLSURFACE))
      return -1;
   
   if (pobj->iobj == SISLCURVE )
   {
      /* Curve case */
      if (pobj->c1->idim != kdim) return -1;
      
      /* Get function values */
      s1221(pobj->o1->c1,kder,parval[0],&kleft1,sval,&kstat);
      if (kstat < 0) return -2;
      
      mult = s6knotmult(pobj->o1->c1->et,pobj->o1->c1->ik,
		     pobj->o1->c1->in,&kleft1,parval[0],&kstat);
      if (kstat < 0) return -2;
      
      if (mult >= pobj->o1->c1->ik - 1)
      {	 
	 /* Get left side function values */
	 s1227(pobj->o1->c1,kder,parval[0],&kleft1,sval1,&kstat);
	 if (kstat < 0) return -2;
	 /* Test function values */
	 if (sval[1] < -ttol && sval1[1] >ttol)       return 1;
	 else if (sval[1] > ttol && sval1[1] < -ttol) return 2;
	 else                                       return 4;
      }
      else
      {
	 /* Test if first derivative is zero */
	 if (fabs(sval[1]) > ttol) return 0;
	 
	 /* Test if max, min or inconclusive point */
	 if      (sval[2] < -ttol) return 1;
	 else if (sval[2] >  ttol) return 2;
	 else                      return 4;
      }
   }
   else
   {
      /* Surface case */
      if (pobj->s1->idim != kdim) return -1;
      
      /* Get function values */
      
      s1421(pobj->o1->s1,kder,parval,&kleft1,&kleft2,
	    sval,snorm,&kstat);
      if (kstat < 0) return -2;
      
      /* Test function values */
      
      /* Test if first derivative is zero */
      tmax = sqrt(sval[1]*sval[1] + sval[2]*sval[2]);
      if (tmax > ttol) return 0;
      
      /* Test if max, min ,saddle or inconclusive point */
      tdet = (sval[3]*sval[5] - sval[4]*sval[4]);
      if      (tdet < -ttol) return 3;
      else if (tdet <  ttol) return 4;
      else if (sval[3] < DZERO) return 1;
      else return 2;
      
   }
   
}	 
   
