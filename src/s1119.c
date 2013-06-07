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
 * $Id: s1119.c,v 1.2 2001-03-19 15:58:41 afr Exp $
 *
 */


#define S1119

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s1119(double *ecoef,double *et1,double *et2,int ik1,int in1,int ik2,
	   int in2,int *jsimple,int *jind1,int *jind2,int *jstat)
#else
void s1119(ecoef,et1,et2,ik1,in1,ik2,in2,jsimple,jind1,jind2,jstat)
     double *ecoef;
     double *et1;
     double *et2;
     int    ik1;
     int    in1;
     int    ik2;
     int    in2;
     int    *jsimple;
     int    *jind1;
     int    *jind2;
     int    *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : Check if a one-dimensional surface can have only one
*              single maximal-point.
*
*
*
* INPUT      : ecoef  - Vertices of surface.
*              et1    - Knots, parameterdirection one.
*              et2    - Knots, parameterdirection two.
*              ik1    - Order of surface in first parameter direction.
*              in1    - Number of vertices in first parameter direction.
*              ik2    - Order of surface in second parameter direction.
*              in2    - Number of vertices in second parameter direction.
*
*
*
* OUTPUT     : jsimple - Indicates if single case
*                        = 0 : No interior max possible.
*                        = 1 : Only one maximal point or curve possible.
*                        = 2 : More than one maximal point possible. 
*              jind1   - Index to an interior knot with multiplicity ik1 in et1.
*                        = 0 : No knot with interior multiplicity.
*              jind1   - Index to an interior knot with multiplicity ik2 in et2.
*                        = 0 : No knot with interior multiplicity.
*              jstat   - status messages  
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*
*
* METHOD     : Count number of times the control-polygon changes
*              direction in each parameter direction. If maximum
*              times of changing is greater than one, this is not
*              a simple case (possibility of several maxima).
*
*
* REFERENCES :
*
*-
* CALLS      :
*
* WRITTEN BY : UJK, SI, 89-06.
*
*********************************************************************
*/
{
  int ki,kj;     /* Counters.                                          */
  int ksimple;   /* Indicates if simple case.                          */
  int ksimple1;  /* Indicates if simple case.                          */
  int ksimple2;  /* Indicates if simple case.                          */
  int ksign;     /* Number of direction changes in line/column.        */
  int kconvex1;  /* Flag, if true, we have no interior min 
		    in first direction.*/
  int kconcav1;  /* Flag, if true, we have no interior max 
		    in first direction.*/
  int kconvex2;  /* Flag, if true, we have no interior min 
		    in second direction.*/
  int kconcav2;  /* Flag, if true, we have no interior max 
		    in second direction.*/
  int kbez;      /* Flagging for bezier case.                          */
  double tfirst; /* First non-zero difference between two adjacent vertices. */
  double tprev;  /* Previous difference between two adjacent vertices. */
  double tdiff;  /* Current difference between two adjacent vertices.  */
  double *s1;    /* Pointer used to traverse array of vertices.        */
  
  /* First we search for interior knotmultiplicity in 
     both parameter directions*/
  *jind1    = 0;
  ksimple1 = 1;
  if (in1 > 1)
    for (ki=ik1+1; ki<in1 && ksimple1; ki++) 
      {  
	if (et1[ki] == et1[ki+ik1-1]) 
	  {
	    *jind1    = ki;
	    ksimple1  =  0;
	  }
      }
  
  *jind2   = 0;
  ksimple2 = 1;
  if (in2 > 1)
    for (ki=ik2+1; ki<in2 && ksimple2; ki++) 
      {  
	if (et2[ki] == et2[ki+ik2-1]) 
	  {
	    *jind2    = ki;
	    ksimple2  =  0;
	  }
      }
  
  
  ksimple = ksimple1 && ksimple2;
  kbez = (((ik1 == in1) && (ik2 == in2)) ? 1 : 0);
  
  /* Count number of direction changes in first parameter direction. */
  /* Notify that we cannot accept equal coeffisient neighbours when 
     we are in a none-bezier case.                                   */
  
  kconcav1 = kconvex1 = 1;
  
  if (in1 > 1)
    for (s1=ecoef,kj=0; kj<in2 && ksimple; kj++)
      {
	ksign = 0;
	tfirst = DZERO;
	
	for (ki=0; ki<in1-1 && ksimple; ki++,s1++)
	  {
	    tdiff = *(s1+1) - *s1;
	    if (DEQUAL(tdiff,DZERO) )
	      { 
		if (kbez == 0) ksimple = 0;
	      }
	    else if (DEQUAL(tfirst,DZERO) )
	      {
		/* First none-zero vector, save it. */
		tfirst = tdiff;
		tprev  = tdiff;
	      }
	    
	    else if (tprev*tdiff < DZERO)
	      {
		tprev = tdiff;
		ksign++;
		if (ksign > 1) ksimple = 0;
	      }
	  }
	
	
	if (kbez == 0)
	  {
	    /* We permit status simple case only in bezier case. 
	       However, if the surface is strictly concav in one 
	       parameter direction, we have found 
	       the max on the edges. */ 
	    kconvex1 = 0;
	    kconcav1 = ((tfirst < DZERO) && kconcav1); 
	  }
	else
	  {
	    
	    kconvex1 = (((ksign == 0) || 
			 (ksign == 1 && tfirst >= DZERO)) && kconvex1); 
	    kconcav1 = (((ksign == 0) || 
			 (ksign == 1 && tfirst <= DZERO)) && kconcav1); 
	  }
	
	ksimple  = ((kconvex1 || kconcav1) && ksimple);
	s1++;
	
      }
  
  /* Count number of direction changes in second parameter direction. */
  kconcav2 = kconvex2 = 1;    
  if (in2 > 1)
    for (kj=0; kj<in1 && ksimple; kj++)
      {
	ksign = 0;
	tfirst = DZERO;
	s1 = ecoef + kj;
	
	for (ki=0; ki<in2-1 && ksimple; ki++,s1+=in1)
	  {
	    tdiff = *(s1+in1) - *s1;
	    if (DEQUAL(tdiff,DZERO) )
	      { 
		if (kbez == 0) ksimple = 0;
	      }
	    else if (DEQUAL(tfirst,DZERO) )
	      {
		/* First none-zero vector, save it. */
		tfirst = tdiff;
		tprev  = tdiff;
	      }
	    
	    else if (tprev*tdiff < DZERO)
	      {
		tprev = tdiff;
		ksign++;
		if (ksign > 1) ksimple = 0;
	      }
	  }
	
	if (kbez == 0)
	  {
	    /* We permit status simple case only in bezier case. 
	       However, if the surface is strictly concav in one 
	       parameter direction, we have found 
	       the max on the edges. */ 
	    kconvex2 = 0;
	    kconcav2 = ((tfirst < DZERO) && kconcav2); 
	  }
	else
	  {
	    
	    kconvex2 = (((ksign == 0) || 
			 (ksign == 1 && tfirst >= DZERO)) && kconvex2); 
	    kconcav2 = (((ksign == 0) || 
			 (ksign == 1 && tfirst <= DZERO)) && kconcav2); 
	  }
	ksimple  = ((kconvex2 || kconcav2) && ksimple);
      }
  
  /* Simple case test performed.  */
  
  if (ksimple)
    {
      if (kconvex1 && kconvex2)
	*jsimple = 1;
      else
	*jsimple = 0;	
    }
  else
    *jsimple = 2;
  *jstat = 0;
  
  return;
}
                              
