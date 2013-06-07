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
 * $Id: tstcyclknt.c,v 1.1 1994-04-21 12:10:42 boh Exp $
 *
 */


#define TEST_CYCLIC_KNOTS

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
  test_cyclic_knots(double et[],int in,int ik,int *jstat)
#else
void test_cyclic_knots(et,in,ik,jstat)
     double et[];
     int    in;
     int    ik;
     int    *jstat;
#endif
     /*
      *********************************************************************
      *                                                                   
      * PURPOSE    : To test if a knot vector is cyclic.
      *
      *             
      *
      * INPUT      : et     - Knot vector
      *              in     - Number of knots
      *              ik     - Polynomial order
      *
      * OUTPUT     : 
      *              jstat  - status messages  
      *                                        = 2 : Cyclic full freedom.
      *                                        = 1 : Cyclic NOT full freedom.
      *                                        = 0 : NOT cyclic.
      *                                        < 0 : Error.
      *
      * METHOD     : 1. Check multiplicity of start and end parameter value
      *              2. Check if knots before start parameter value and after
      *                 end parameter value are repeated in a cyclic way.
      *              3. Check that number of basis function not repeated
      *                 are at least ik.
      *
      *
      * REFERENCES :
      *
      *-                                      
      * CALLS      : 
      *
      * WRITTEN BY : Tor Dokken, SI, Oslo, Norway, April 1992
      *
      *********************************************************************
      */
{
  int    kleft;          /* Pointer into knot interval             */
  int    kmult1;         /* Multiplicity of start parameter value  */
  int    kmult2;         /* Multiplicity of end  parameter value   */
  int    ki;             /* Control variable in loop               */
  int    kpos = 1;       /* Position of error                      */
  int    kant;           /* Number of knots before start parameter value */
  int    kcyclic;        /* Flag telling if cyclic basis           */
  int    kstat;          /* Local status variable                  */
  
  double tperiode;       /* Periode of basis                       */
  
  /* Find multiplicity of et[ik-1] and et[in] */
  
  kleft = ik-1;
  
  kmult1 = s6knotmult(et,ik,in,&kleft,et[ik-1],&kstat);
  if(kstat<0) goto error;
  
  kleft = in;
  
  kmult2 = s6knotmult(et,ik,in,&kleft,et[in],&kstat);
  if(kstat<0) goto error;
  
  if (kmult1 != kmult2 || kmult1 == ik) goto noncyclic;
  
  kant = ik - kmult1;
  tperiode = et[in] -et[ik-1];  
  
  /* Test that the first kant knots are repetitions of the knots in-kant,...,in-1 */
  
  for (ki=0, kcyclic=1; ki<kant ; ki++)
    if (DNEQUAL((et[ki]+tperiode),et[in-kant+ki])) kcyclic = 0;
  
  /* Test that the last kant knots are repetions of knots ik,..,ik+kant-1 */
  
  for (ki=0; ki<kant ; ki++)
    if (DNEQUAL((et[ik+ki]+tperiode),et[in+kmult1+ki])) kcyclic = 0;
  
  if (kcyclic == 0) goto noncyclic;
  
  /* The basis should have at least kant+ik degrees of freedom to allow for
     a proper cyclic curve with ik degrees of real freadom since kant vertices
     are repeated */
  
  if (in<ik+kant) goto missing_freedom; 
  
  /* Cyclic with enough degrees of freedom */
  
  *jstat = 2;
  goto out;
  
  /* Cyclic with less than ik+kant degrees of freedom */
 missing_freedom:
  
  *jstat = 1;
  goto out;
  
  /* Noncyclic basis */
 noncyclic:
  
  *jstat = 0;
  goto out;
  
 error:
  *jstat = kstat;
  s6err("test_cyclic_knots",*jstat,kpos);
  goto out;
  
 out:
  
  return;
  
}


