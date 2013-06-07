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
 * $Id: s6degnorm.c,v 1.1 1994-04-21 12:10:42 boh Exp $
 *
 */


#define S6DEGNORM

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
   s6degnorm(SISLSurf *ps1,int ider,double epar[],double eder[],
		double utang[],double vtang[],double enorm[],int *jstat)
#else
void s6degnorm(ps1,ider,epar,eder,utang,vtang,enorm,jstat)
     SISLSurf   *ps1;
     int    ider;
     double epar[];
     double eder[];
     double utang[];
     double vtang[];
     double enorm[];
     int    *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : Evaluate the (3d) surface tangents and normal
*              (at an edge or corner) when the surface is
*              degenerate (at an edge or corner).
*              All three vectors (if found) are normalised to length 1.
*         
*
*
* INPUT      : ps1    - Pointer to the surface to evaluate.
*              ider   - Number of derivatives calculated.
*                       Must be >=2.
*              epar   - Parameter-value at which to calculate. Dimension
*                       of epar is 2.
*              eder   - Array where the derivative of the curve in
*                       apar is placed. The sequence is position,
*                       first derivative in first parameter direction,
*                       first derivative in second parameter direction,
*                       (2,0) derivative, (1,1) derivative, (0,2) 
*                       derivative, etc. Dimension of eder is 
*                       idim*(1+2+...+(ider+1)).
*             
*
* OUTPUT     : utang -  Tangent to surface in u direction, length 1.
*              vtang -  Tangent to surface in v direction, length 1.
*              enorm  - Normal of surface, length 1.
*              jstat  - status messages  
*                             = 4      : Nothing found.
*                                        both tangents and normal are 0.
*                             = 3      : Only v tangent found.
*                                        u tangent and normal are 0.
*                             = 2      : Only u tangent found.
*                                        v tangent and normal are 0.
*                             = 1      : Only the tangents found.
*                                        Normal is 0.
*                             = 0      : OK. Both tangents and
*                                        normal are length 1.
*                             < 0      : error
*                      
*
* METHOD     :  Find the smallest non-zero derivatives and
*               use them in a Taylor expansion.
*
* REFERENCES :
*
*-
* CALLS      : 
*
* WRITTEN BY : Michael Floater, 17/1/92.
* Revised by : Christophe Rene Birkeland, SINTEF OSLO, June 1993.
*   
*********************************************************************
*/                                     
{
  int kstat=0;        /* Local status variable.                          */
  int kpos=0;         /* Position of error.                              */
  int kdim;           /* Dimension of the space in which the surface lies. */
  int ki;             /* Control variables in for loop                   */
  double *et1,*et2;   /* Local pointer to knot vectors. */
  int in1,in2;        /* Number of points in ps1 in the 2 direcs. */
  int ik1,ik2;        /* Degree of ps1 in the 2 direcs. */
  double upar,vpar;   /* Parameter values. */
  double *xu,*xv;         /* Pointers to first derivatives. */
  double *xuu,*xuv,*xvv;  /* Pointers to second derivatives. */
  double len;        /* Vector length. */
  int ius,ivs,is;    /* Flags. u=min => ius = 1, u=max => ius = -1. */
  int endu,endv;     /* Flags for whether u or v are extreme. */
  int iu,iv;         /* Which first derivs are zero? */
  int iuu,iuv,ivv;   /* Which second derivs are zero? */
  double vec[3];     /* Temporary vector. */
  double vec1[3],vec2[3];  /* Temporary vectors. */
  double normal[3];  /* Temporary normal. */
  int usuccess;      /* Flag if u tangent found. */
  int vsuccess;      /* Flag if v tangent found. */
  
  
  /* Set up local variables. */

  kdim = ps1 -> idim;
  et1 = ps1 -> et1;
  et2 = ps1 -> et2;
  in1 = ps1 -> in1;
  in2 = ps1 -> in2;
  ik1 = ps1 -> ik1;
  ik2 = ps1 -> ik2;

  /* Check input. */

  if(kdim != 3) goto err101;
  if(ider < 2) goto err101;
  

  upar = epar[0];
  vpar = epar[1];

  xu  = eder + kdim;
  xuu = xu   + kdim;
  xv  = xuu  + kdim;
  xuv = xv   + kdim;
  xvv = xuv  + kdim + kdim;

  /* Find out whether (u,v) is at a corner, edge or in the
     middle of the surface ps1. */

  ius = 0;
  ivs = 0;
  is = 0;

  if(upar == et1[ik1-1])
  {
      endu = TRUE;
      ius = 1;
  }
  else if(upar == et1[in1])
  {
      endu = TRUE;
      ius = -1;
  }
  else
  {
      endu = FALSE;
  }

  if(vpar == et2[ik2-1])
  {
      endv = TRUE;
      ivs = 1;
  }
  else if(vpar == et2[in2])
  {
      endv = TRUE;
      ivs = -1;
  }
  else
  {
      endv = FALSE;
  }

  if(endu && endv) is = ius * ivs;

  if(!endu && !endv) goto err101;

  /* For each derivative, set flag to 0 or 1 according to
     whether the length is 0 or non-zero. */

  len = s6length(xu,kdim,&iu);
  len = s6length(xv,kdim,&iv);
  len = s6length(xuu,kdim,&iuu);
  len = s6length(xuv,kdim,&iuv);
  len = s6length(xvv,kdim,&ivv);


  /* Calculate tangent in u using higher derivatives. */

  usuccess = FALSE;

  if(iu == 0)
  {
      if(endu && iuu == 1)
      {
          len = s6norm(xuu,kdim,vec,&kstat);
          for(ki=0; ki<kdim; ki++) vec[ki]*=ius;
          usuccess = TRUE;
      }
      else if(endv && iuv == 1)
      {
          len = s6norm(xuv,kdim,vec,&kstat);
          for(ki=0; ki<kdim; ki++) vec[ki]*=ivs;
          usuccess = TRUE;
      }
  }
  else
  {
      len = s6norm(xu,kdim,vec,&kstat);
      usuccess = TRUE;
  }


  if(usuccess)
  {
      /* u tangent found. Return result. */

      for(ki=0; ki<kdim; ki++) utang[ki] = vec[ki];
  }
  else
  {
      /* No u tangent found. Return zero and flag. */

      for(ki=0; ki<kdim; ki++) utang[ki] = (double)0.0;
  }

  /* Calculate tangent in v using higher derivatives. */

  vsuccess = FALSE;

  if(iv == 0)
  {
      if(endu && iuv == 1)
      {
          len = s6norm(xuv,kdim,vec,&kstat);
          for(ki=0; ki<kdim; ki++) vec[ki]*=ius;
          vsuccess = TRUE;
      }
      else if(endv && ivv == 1)
      {
          len = s6norm(xvv,kdim,vec,&kstat);
          for(ki=0; ki<kdim; ki++) vec[ki]*=ivs;
          vsuccess = TRUE;
      }
  }
  else
  {
      len = s6norm(xv,kdim,vec,&kstat);
      vsuccess = TRUE;
  }


  if(vsuccess)
  {
      /* v tangent found. Return result. */

      for(ki=0; ki<kdim; ki++) vtang[ki] = vec[ki];
  }
  else
  {
      /* No v tangent found. Return zero and flag. */

      for(ki=0; ki<kdim; ki++) vtang[ki] = (double)0.0;
  }


  /* Calculate normal using higher derivatives. */

  if(iu == 0)
  {
      if(iv == 0)
      {
	  if(endu && iuu == 1 && iuv == 1)
	  {
	      s6crss(xuu,xuv,vec);
	      len = s6norm(vec,kdim,normal,&kstat);
	      if(kstat == 1) goto normfound;
	  }
	  if(endv && iuv == 1 && ivv == 1)
	  {
	      s6crss(xuv,xvv,vec);
	      len = s6norm(vec,kdim,normal,&kstat);
	      if(kstat == 1) goto normfound;
	  }
	  if(endu && endv && iuu == 1 && ivv == 1)
	  {
	      s6crss(xuu,xvv,vec);
	      for(ki=0; ki<kdim; ki++) vec[ki]*=is;
	      len = s6norm(vec,kdim,normal,&kstat);
	      if(kstat == 1) goto normfound;
	  }
      }
      else
      {
	  if(endu && iuu == 1)
	  {
	      s6crss(xuu,xv,vec);
	      for(ki=0; ki<kdim; ki++) vec[ki]*=ius;
	      len = s6norm(vec,kdim,normal,&kstat);
	      if(kstat == 1) goto normfound;
	  }
	  if(endv && iuv == 1)
	  {
	      s6crss(xuv,xv,vec);
	      for(ki=0; ki<kdim; ki++) vec[ki]*=ivs;
	      len = s6norm(vec,kdim,normal,&kstat);
	      if(kstat == 1) goto normfound;
	  }
      }
  }
  else
  {
      if(iv == 0)
      {
	  if(endu && iuv == 1)
	  {
	      s6crss(xu,xuv,vec);
	      for(ki=0; ki<kdim; ki++) vec[ki]*=ius;
	      len = s6norm(vec,kdim,normal,&kstat);
	      if(kstat == 1) goto normfound;
	  }
	  if(endv && iuv == 1)
	  {
	      s6crss(xuv,xv,vec);
	      for(ki=0; ki<kdim; ki++) vec[ki]*=ivs;
	      len = s6norm(vec,kdim,normal,&kstat);
	      if(kstat == 1) goto normfound;
	  }
      }
      else
      {
	  if(endu && (iuu == 1 || iuv == 1))
	  {
	      s6crss(xuu,xv,vec1);
	      s6crss(xu,xuv,vec2);
	      for(ki=0; ki<kdim; ki++) vec[ki]=vec1[ki]+vec2[ki];
	      for(ki=0; ki<kdim; ki++) vec[ki]*=ius;
	      len = s6norm(vec,kdim,normal,&kstat);
	      if(kstat == 1) goto normfound;
	  }
	  if(endv && (iuv == 1 || ivv == 1))
	  {
	      s6crss(xuv,xv,vec1);
	      s6crss(xu,xvv,vec2);
	      for(ki=0; ki<kdim; ki++) vec[ki]=vec1[ki]+vec2[ki];
	      for(ki=0; ki<kdim; ki++) vec[ki]*=ivs;
	      len = s6norm(vec,kdim,normal,&kstat);
	      if(kstat == 1) goto normfound;
	  }
      }
  }

  /* No normal found. Return zero and flag. */

  for(ki=0; ki<kdim; ki++) enorm[ki] = (double)0.0;

  /* Set diagnostics flag. */
  if(usuccess)
  {
      if(vsuccess) *jstat = 1;
      else *jstat = 2;
  }
  else
  {
      if(vsuccess) *jstat = 3;
      else *jstat = 4;
  }
  goto out;

   /* Normal found and hence tangents found. Return result. */

normfound:

  for(ki=0; ki<kdim; ki++) enorm[ki] = normal[ki];
  *jstat = 0;
  goto out;


  /* Error in input. */

err101: *jstat = -101;
  s6err("s6degnorm",*jstat,kpos);
  goto out;
  
  
 out:
  
  return;
}

