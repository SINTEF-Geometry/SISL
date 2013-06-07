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
 * $Id: s1902.c,v 1.7 2001-03-19 15:58:55 afr Exp $
 *
 */


#define S1902

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
s1902 (double epar[], int in, int ik, int cuopen, double **eknots, int *jstat)

#else
void
s1902 (epar, in, ik, cuopen, eknots, jstat)
     double epar[];
     int in;
     int ik;
     int cuopen;
     double **eknots;
     int *jstat;

#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    :	To produce the knot vector of a B-spline basis satisfying the
*		interpolation requirements reflected in the epar array.
*
* INPUT      :	epar	- Array containing a parametrization of the
*			  interpolation conditions. Each interpolation
*			  condition has got a distinct parameter value
*			  expect form the cases where several conditions
*			  are conflicting. In that case a multiple parameter
*                         value indicates the need of a multiple knot. The
*                         parameter values are sorted in increasing order.
*                         The dimension of the array is 'in' if the curve is
*			  open and 'in+1' if it is closed.
*               in     	- Number of interpolation conditions.
*               ik      - Order of B-spline basis.
*		cuopen	- Open/closed curve.
*
* OUTPUT     : eknots - The produced knot vector. The dimension of
*                       the array is in+ik if the curve is open.
*			If the curve is closed the dimension of the array
*			is in+2*ik-1 if the curve is even and in+2*ik if it is
*			odd.
*              jstat  - status messages
*
* METHOD     :
*
* REFERENCES :
*
* CALLS      :
*
* WRITTEN BY :	Vibeke Skytt, SI, 91-03
* REVISED BY :	Trond Vidar Stensby, SI, 91-06
*
*********************************************************************
*/
{
  int ki;			/* Counter used to traverse the knot vector.	*/
  int kpar;			/* Counter used to traverse the parametrization
	                      	   array.					*/
  int kn;			/* The number of conditions in epar. (closed)   */
  int kk2;			/* Half the order.				*/
  int kstop;			/* Control variable of loop.            	*/
  int smult;			/* Multiplisity of start parameter value        */
  int emult;			/* Multiplisity of end parameter value          */
  int kmult;
  double tprev;			/* Value of previous knot.                     	*/
  double th;			/* Value of current knot.                      	*/
  double tval1;			/* Start parameter value.			*/
  double tval2;			/* End parameter value.            		*/
  double tparint;		/* The parameter interval. (closed)		*/
  double tdum;			/* Help parameter used for parameter interval.	*/
  double dummy1, dummy2;
  double delta;
  int count1;

  /* Check if curve is closed or open. */

  if (cuopen == SISL_CRV_OPEN)
    {
      /* O P E N   C U R V E */

      *eknots = newarray (in +ik, DOUBLE);
      if (*eknots == SISL_NULL)
	goto err101;

      kk2 = ik / 2;
      tval1 = epar[0];
      tval2 = epar[in -1];

      /* Store a knot of multiplisity equal to the order at the start of the
	 curve The value of the knot is equal to the value of the start
	 parameter value. */

      for (ki = 0; ki < ik; ki++)
	(*eknots)[ki] = tval1;

      /* Find the multiplisity of the start parameter value. */

      for (smult = 1; epar[smult - 1] == epar[smult]; smult++) ;

      /* Find the multiplisity of the end parameter value. */

      for (emult = 1, count1 = in -2; epar[count1] == epar[count1 + 1];
	   count1--)
	emult++;


      if (ik % 2 == 0)
	{
	  /* The order is even.
	     Place the internal knots at the parameter values.  */

	  dummy1 = (double)0.5 * (epar[smult] + tval1);
	  dummy2 = (double)0.5 * (epar[in -emult - 1] +tval2);
	  if (dummy1 == dummy2)
	    {
	      dummy1 = (epar[smult] + tval1 + tval1)/(double)3;
	      dummy2 = (epar[in -emult - 1] + tval2 + tval2)/(double)3;
	    }
	  for (count1 = 0; count1 < smult - kk2; count1++, ki++)
	    (*eknots)[ki] = dummy1;

	  for (kpar = smult + MAX (0, kk2 - smult), kstop = in -emult -
	       MAX (0, kk2 - emult);
	       kpar < kstop; kpar++, ki++)
	    (*eknots)[ki] = epar[kpar];

	  for (count1 = 0; count1 < emult - kk2; count1++, ki++)
	    (*eknots)[ki] = dummy2;
	}
      else
	{
	  /* The order is odd.
	     Place the internal knots between the parameter values.  */

	  if (smult - kk2 > 0)
	    {
	      delta = (epar[smult] - tval1) / (smult - kk2 + 1);
	      dummy1 = tval1 + delta;
	      for (count1 = 0; count1 < smult - kk2;
		   count1++, dummy1 += delta, ki++)
		(*eknots)[ki] = dummy1;
	    }

	  for (kpar = smult + MAX (0, kk2 - smult), kstop = in -emult -
	       MAX (0, kk2 - emult);
	       kpar < kstop; kpar++, ki++)
	    (*eknots)[ki] = (double) 0.5 *(epar[kpar] + epar[kpar + 1]);

	  if (emult - kk2 > 0)
	    {
	      delta = (tval2 - epar[in -emult - 1]) / (emult - kk2 + 1);
	      dummy2 = epar[in -emult - 1] +delta;
	      for (count1 = 0; count1 < emult - kk2; count1++, dummy2 +=
		   delta, ki++)
		(*eknots)[ki] = dummy2;
	    }
	}

      /* Store a knot of multiplisity equal to the order at the end of
	 the curve. The value of the knot is equal to the value of the */


      for (ki = 0; ki < ik; ki++)
	(*eknots)[in +ki] = tval2;
    }
  else
    {
      /* C L O S E D   C U R V E */

      *eknots = newarray (in +2 * ik, DOUBLE);
      if (*eknots == SISL_NULL)
	goto err101;

      kn = in +1;
      kk2 = ik / 2;
      kstop = in +2 * ik - 1;
      tparint = epar[in] -epar[0];

      if (ik % 2 == 0)
	{
	  /* The order of the B-spline curve is even.
	     Make the ik/2 first knots as a shift of the last knots.  */

	  for (ki = 0, kpar = in -ik + ik / 2; ki < ik / 2; ki++, kpar++)
	    (*eknots)[ki] = epar[kpar] - tparint;
	  
	  /* Make the knots corresponding to the data points. */

	  for (kpar = 0; kpar < kn; ki++, kpar++)
	    (*eknots)[ki] = epar[kpar];


	  /* Make the ik+ik/2-2 last knots.  */

	  for (kpar = 1; ki < kstop; ki++, kpar++)
	    {
	      tdum = tparint;

	      /* We may risk that a double cyclic use of the parameter
	         values may result.          */

	      if (kpar > in)
		{
		  tdum += tparint;
		  kpar = 1;
		}
	      (*eknots)[ki] = epar[kpar] + tdum;
	    }
	}
      else
	{
	  /* The order of the B-spline curve is odd.
	     Make the ik/2+1 first knots.             */

	  for (ki = 0, kpar = in -ik + ik / 2; ki < ik / 2 + 1; ki++, kpar++)
	    (*eknots)[ki] = (double) 0.5 *(epar[kpar] + epar[kpar + 1]) -
	     tparint;

	  /* Make the in next knots.  */

	  for (kpar = 0; kpar < in; ki++, kpar++)
	    (*eknots)[ki] = (double) 0.5 *(epar[kpar] + epar[kpar + 1]);

	  /* Make the ik remaining knots.  */

	  for (kpar = 0; ki < kstop; ki++, kpar++)
	    {
	      tdum = tparint;

	      /* We may risk that a double cyclic use of the parameter
	         values may result.        */

	      if (kpar >= in)
		{
		  tdum += tparint;
		  kpar = 0;
		}
	      (*eknots)[ki] = (double) 0.5 *(epar[kpar] + epar[kpar + 1]) + tdum;
	    }
	}
    }


  /* Check that the produced knots are in increasing order and that
     the multiplicity is not greater than ik.                       */

  if (cuopen == SISL_CRV_OPEN)
    kstop = in +ik;

  for (ki = 1, tprev = (*eknots)[0], kmult = 0; ki < kstop; ki++, tprev = th)
    {
       th = (*eknots)[ki];
      kmult++;
      if (tprev > th)
	goto err112;		/* Decreasing parameter value. */
      if (tprev < th)
	kmult = 1;
      if (kmult > ik)
	goto err112;		/* Knot multiplisity greater than order. */
    }

  /* The knot vector is produced.  */

  *jstat = 0;
  goto out;

  /* Error in scratch allocation. */

err101:
  *jstat = -101;
  s6err ("s1902", *jstat, 0);
  goto out;

  /* Error in the knot vector.  */

err112:*jstat = -112;
  s6err ("s1902", *jstat, 0);
  goto out;

out:
  return;
}
