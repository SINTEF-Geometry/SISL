/*****************************************************************************/
/*                                                                           */
/*                                                                           */
/* (c) Copyright 1989,1990,1991,1992 by                                      */
/*     Senter for Industriforskning, Oslo, Norway                            */
/*     All rights reserved. See the sisl-copyright.h for more details.            */
/*                                                                           */
/*****************************************************************************/

#include "sisl-copyright.h"

/*
 *
 * $Id: s1910.c,v 1.2 2001-03-19 15:58:56 afr Exp $
 *
 */


#define S1910

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
s1910 (double econd[], int ntype[], int inpt, int idim, int iopen,
       double astpar, double *cendpar, double *epar1[],
       double *epar2[], int *jstat)
#else
void
s1910 (econd, ntype, inpt, idim, iopen, astpar, cendpar, epar1, epar2, jstat)
     double econd[];
     int ntype[];
     int inpt;
     int idim;
     int iopen;
     double astpar;
     double *cendpar;
     double *epar1[];
     double *epar2[];
     int *jstat;
#endif
/*
*********************************************************************
*
* PURPOSE    : Make uniform parametrization of interpolation conditions.
*	       Derivative conditions are given the same value
*	       as the corresponding point condition.
*
* INPUT      : econd  - Array of interpolation conditions. The last
*                       condition must be a point. Dimension is inpt*idim.
*              ntype  - Array containing kind of condition. Dimension
*                       is inpt.
*                       =  0 : A point is given.
*                       =  d : The d'th derivatative condition to the
*                              previous point is given.
*                       = -d : The d'th derivatative condition to the
*                              next point is given.
*              inpt   - Number of interpolation conditions.
*              idim   - Dimension of geometry space.
*	       iopen  - Open/closed curve.
*              astpar - Start parameter of parametrization.
*
* OUTPUT     : cendpar - End parameter of parametrization.
*              epar1   - Parametrization array. Derivative conditions has
*                        got the same parameter value as the corresponding
*                        positional condition.
*              epar2   - Parametrization array to use when making a knot
*                        vector. All entries are different.
*              jstat  - status messages
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*
* METHOD     :
*
* REFERENCES :
*
* CALLS      : s6dist,s6err.
*
* WRITTEN BY : Trond Vidar Stensby, SI, 1991-07
*
*********************************************************************
*/
{
  int count1, count2;	/* Loop control variables. */
  int prev;
  int knpt;			/* Number of parameter values to
				   be produced. */
  double sum, cord;
  int kpos = 0;

  *jstat = 0;

  if (iopen != SISL_CRV_OPEN)
    knpt = inpt + 1;
  else
    knpt = inpt;


  /* Allocate arrays. */

  *epar1 = newarray (knpt, DOUBLE);
  if (*epar1 == SISL_NULL)
    goto err101;

  *epar2 = newarray (knpt, DOUBLE);
  if (*epar2 == SISL_NULL)
    goto err101;


  /* Find average cord_length. */

  sum = (double) 0.0;
  prev = -1;
  for (count1 = 0, count2 = 0; count1 < inpt; count1++)
    {
      if (ntype[count1] == 0)
	{
	  if (prev >= 0)
	    sum += s6dist (&econd[count1 * idim], &econd[prev * idim],
			   idim);
	  prev = count1;
	  count2++;
	}
    }
  cord = sum / (double) (count2 -  1.0);


  for (count1 = 0; count1 < inpt; count1++)
    {
      if (ntype[count1] > 0)
	(*epar1)[count1] = astpar - cord;
      else if (ntype[count1] < 0)
	(*epar1)[count1] = astpar;
      else
	{
	  (*epar1)[count1] = astpar;
	  astpar += cord;
	}
    }

  if (iopen != SISL_CRV_OPEN)
    (*epar1)[inpt] = astpar;

  *cendpar = (*epar1)[knpt - 1];


  /* Find distinct values. */

  count2 = 1;
  (*epar2)[0] = (*epar1)[0];

  for (count1 = 1; count1 < knpt; count1++)
    {
      if ((*epar1)[count1 - 1] < (*epar1)[count1])
	{
	  (*epar2)[count2] = (*epar1)[count1];
	  count2++;
	}
    }

  *epar2 = increasearray (*epar2, count2, DOUBLE);
  if (*epar2 == SISL_NULL)
    goto err101;


  /* Parametrization computed.  */

  goto out;


  /* Error in scratch allocation. */

err101:
  *jstat = -101;
  s6err ("s1910", *jstat, kpos);
  goto out;

out:

  return;
}
