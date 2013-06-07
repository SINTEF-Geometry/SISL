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
 * $Id: s1700.c,v 1.1 1994-04-21 12:10:42 boh Exp $
 *
 */


#define S1700

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s1700(int imy,int ik,int in,int iv,
	   int *jpl,int *jfi,int *jla,double *et,double apar,double *galfa,int *jstat)
#else
void s1700(imy,ik,in,iv,jpl,jfi,jla,et,apar,galfa,jstat)
     int    imy;
     int    ik;
     int    in;
     int    iv;
     int    *jpl;
     int    *jfi;
     int    *jla;
     double *et;
     double apar;
     double *galfa;
     int    *jstat;
#endif
/*
********************************************************************
*
*********************************************************************
*
* PURPOSE  : To compute in a compact format a line in the discrete
*            B-spline matrix converting between an orginal basis
*            "et" and a new basis consisting of "et" and one p-tuppel
*            knot, where 0 < p <= ik.
*            To insert one p-tuppel knot the function must be called
*            ik - (multiplicity of the knot) + 2*(p-1) times.
*            The function can also be used for inserting knots at
*            more than one value if at least ik-1 knots exist in "et"
*            between the values in question.
*            In this case apar must be changed and the index of the
*            new vertices must be adjusted in the proper manner
*            by the calling function.
*
*
*
* INPUT    : imy   - If one p-tuppel knot are to be inserted, and
*                    the index of the new vertice are j then
*                    imy = j when et[j] < apar and then unchanged
*                    for the rest of j.
*                    If knots are to be inserted in more than one
*                    value, index the new vertices as it is only
*                    one value to be insert, and adjust the index
*                    after calling this function.
*            ik    - The order of the B-spline.
*            in    - The number of the orginal vertices.
*            iv    - The number of new knots to insert in the
*                    area between knot number j and knot number j+ik
*                    in the new knot vector.
*            et    - The knot vector.
*            apar  - The value where to insert new knots.
*
*
*
* OUTPUT   : jpl   - The negativ difference between the index in galfa
*                    and the real knot inserten matrix.
*            jfi   - The index of the first element in the line j in the
*                    the real knot inserten matrix whice is not zero.
*                    The element with the index (jfi+jpl) in galfa
*                    is the same as the element with index jfi in
*                    the real line j in the knot inserten matrix.
*            jla   - The index of the last element in the line j in the
*                    real knot inserten matrix whice is not zero.
*                    The element with the index (jla+jpl) in galfa
*                    is the same as the element with index jla in
*                    the real line j in the knot inserten matrix.
*            galfa - A compressed line in the knot inserten matrix.
*            jstat - status messages
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*
*
* METHOD     : Using the Oslo-algorithm
*
*
* REFERENCES : Making The Oslo algorithm more efficient.
*              by  T.Lyche and K.Moerken.
*              SIAM J.NUMER.ANAL  Vol. 23, No. 3, June 1986.
*
*-
* CALLS      :
*
* WRITTEN BY : Arne Laksaa, SI, 88-06.
*
**********************************************************************/
{
  int kpos=0;              /* Posisjon of error.           */
  int kj,kv;               /* Help variable                */
  int kp;                  /* Control variable in loop.    */
  double *salfa;           /* Help pointer to galfa.       */
  double tbeta,tbeta1;     /* Help variabels               */
  double td1,td2;          /* Help variabels               */
  double *t1,*t2;          /* Pointers to the knot vector. */
  
  
  /* Check that the number of knots we insert is not to large */
  
  if (iv >= ik) goto err152;
  
  
  /* Compute the negativ difference between the index in galfa and
     the real knot inserten matrix. */
  
  *jpl=ik-imy-1;
  
  
  /* Changing the galfa so we may use the index in the real matrix. */
  
  galfa += *jpl;
  
  
  /* Initialise the last element. */
  
  galfa[imy] = 1;
  
  
  /* Here we go one time for each new knot we insert. */
  
  for (kj=in+iv-2,in+=ik-1,kv=ik-iv,kp=0; kp<iv; kp++,kv++)
    {
      /* The initialising:  The two first are not changing.
	 kj = in+iv-2, minus the maximum of kp it
	 gives the index of the last
	 orginal vertices.
	 in = in+ik-1, the index of the last element in et.
	 kv = ik-iv ,  the nuber of old knots in the field.
	 This variabel is counting up to ik
	 (the order) during the loops. */
      
      
      /* Here we note the special case where we are at the
	 start of the matrix and we does not have a k-touple
	 knot at this end. */
      
      if (kp>=imy) tbeta1=(apar - *et)* *galfa/(et[kv] - *et);
      else         tbeta1=(double)0.0;
      
      
      *jfi=max(1,imy-kp); *jla=min(imy,kj-kp);
      
      
      /* For details about this loop look in the reference. */
      
      for (salfa=galfa+*jfi,t1=et+*jfi,t2=et+*jla; t1<=t2; t1++,salfa++)
	{
	  td1 = apar - *t1;
	  td2 = t1[kv] - apar;
	  tbeta = *salfa/(td1 + td2);
	  salfa[-1] = td2*tbeta + tbeta1;
	  tbeta1 = td1*tbeta;
	}
      
      
      /* Here we note the special case where we are at the
	 end of the matrix and we does not have a k-touple
	 knot at this end. */
      
      if (*jla<imy)
	{
	  t1 = et + in;
	  *(salfa-1) = tbeta1+(*t1-apar)* *salfa/(*t1 - *(t2+1));
	} else  *(salfa-1) = tbeta1;
    }
  
  
  /* Adjusting the index of first and last in galfa. */
  
  if (iv) (*jfi)--;
  else   *jfi = *jla = imy;
  
  
  /* Updating output. */
  
  *jstat = 0;
  goto out;
  
  
  /* Error, to many insertions knots. */
  
 err152:
  *jstat = -152;
  s6err("s1700",*jstat,kpos);
  goto out;
  
 out: 
  return;
}
