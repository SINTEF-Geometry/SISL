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
 * $Id: s1611.c,v 1.2 2001-03-19 15:58:51 afr Exp $
 *
 */


#define S1611

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s1611 (double epoint[], int inbpnt, int idim, double eptyp[], int iopen,
       int ik, double astpar, double aepsge,
       double *cendpar, SISLCurve ** rc, int *jstat)
#else
void 
s1611 (epoint, inbpnt, idim, eptyp, iopen, ik, astpar, aepsge,
       cendpar, rc, jstat)
     double epoint[];
     int inbpnt;
     int idim;
     double eptyp[];
     int iopen;
     int ik;
     double astpar;
     double aepsge;
     double *cendpar;
     SISLCurve **rc;
     int *jstat;
#endif
/*
*************************************************************************
*
* PURPOSE: To approximate a conic arc with a b-spline curve in the two or
*          the three dimensional space. If two points are given a straight
*          line is produced, if three a circular arc and if four or five
*          a conic arc.
*
* INPUT:
*        Epoint - Array (length idim*inbpnt) containing the points/
*                 tangents to be interpolated.
*        Inbpnt - No. of points/tangents in the epoint array.
*        Idim   - The dimension of the space in which the points lie.
*        Eptyp  - Array (length inbpnt) containing type indicator for
*                 points/tangents :
*                  1 - Ordinary point.
*                  2 - Knuckle point. (Is treated as an ordinary point.)
*                  3 - Tangent to next point.
*                  4 - Tangent to prior point.
*        Iopen  - Flag telling if the curve should be open or closed:
*                  0 : Closed curve.
*                  1 : Open curve.
*                 NOT IN USE!
*        Ik     - The order of the B-spline curve to be produced.
*        Astpar - Parameter-value to be used at the start of the curve.
*        Aepsge - The geometry resolution.
*
* Output:
*        Cendpar - Parameter-value used at the end of the curve.
*        Rc      - Pointer to output-curve.
*        Jstat   - Status variable:
*                  < 0 : Error.
*                  = 0 : Ok.
*                  > 0 : Warning.
*
* Method: First we find the conic curve on which the points lie.
*         Then a b-spline curve describing the conic is produced.
*-
* Calls: s1614, s6crss, s6rotmat, s6inv4, s6mulvec, s1615, s1616, s1617,
*        s1385, s1602, s6err.
*
* Written by: A.M. Ytrehus, SI Oslo Oct.91.
* After FORTRAN (P19506), written by: T. Dokken  SI.
* REVISED BY : Christophe Birkeland, SI, July 1992 (eptyp is now DOUBLE)
* REVISED BY : Johannes Kaasa, SI, Aug. 1992 (transported status 105 from
*              s1616 upwards)
* REVISED BY : Johannes Kaasa, SI, Aug. 1992 (Fixed bug in production of
*              the rotational matrix, and made a seperate handling of
*              circular arcs)
*****************************************************************
*/
{
  double *spoint = SISL_NULL;
  int    *sptyp = SISL_NULL;
  double smatrix[16];
  double sinv[16];
  double *st = SISL_NULL;
  double sbeg[3], slutt[3], svect1[3], svect2[3], snorm[3];
  int knbpnt = 0;
  int ktyp;
  int kp, ki, kk, kn;
  double tdum, tpar;
  double tshape;
  double spnt[3], sdum[3];
  double start[3], stang[3], stop[3];
  double strans[3];
  double smul[3];
  double sconic[6];
  int kpos = 0;
  int krem = 0;
  int kstat = 0;
  int *ieptyp = SISL_NULL; /* The point type descriptions (eptyp) in integer form. */
  double centerpt[3], axis[3], sdir[3];
  double angle, sdot, tlength1, tlength2, tcos;
  int kstat1, kstat2;

  *jstat = 0;

  
  /* Copy eptyp to ieptyp. */

  ieptyp = newarray(inbpnt,INT);
  if (ieptyp == SISL_NULL)
     goto err101;
  
  for(ki=0; ki<inbpnt; ki++)
  {
      ieptyp[ki] = (int)eptyp[ki];
  }
  

  /* Test if legal input. */

  if (idim < 2 || idim > 3 || ik < 3)
    goto err104;


  /* Allocate new arrays. */

  spoint = newarray (inbpnt * idim, DOUBLE);
  if (spoint == SISL_NULL)
    goto err101;

  sptyp = newarray (inbpnt, INT);
  if (sptyp == SISL_NULL)
    goto err101;


  /* Check the types of the data/points. */

  s1614 (epoint, inbpnt, idim, ieptyp,spoint, &knbpnt, sptyp, &kstat);
  if (kstat < 0)
    goto error;


  /* Save the start- and end-points, in case of straight line. */

  kk = idim * (knbpnt - 1);

  for (kp = 0; kp < idim; kp++)
    {
      sbeg[kp] = spoint[kp];
      sdir[kp] = spoint[idim + kp];
      slutt[kp] = spoint[kk + kp];
    }


  /* Test if straight line is to be produced. */

  if (knbpnt == 2)
    goto straight;


  /* The calculations on the conic are performed two-dimensionally.
     Obtain a 2D description if the conic lies in a 3D space. */

  if (idim == 3)
    {

      /* Translate the points (and not the tangents) such that
         the first point lies in the origin. Remember translation. */

      for (kp = 0; kp < idim; kp++)
	strans[kp] = spoint[kp];

      for (ki = 0; ki < knbpnt; ki++)
	{
	  ktyp = sptyp[ki];
	  if (ktyp < 3)
	    {
	      kk = idim * ki;

	      for (kp = 0; kp < idim; kp++)
		spoint[kk + kp] = spoint[kk + kp] - strans[kp];
	    }
	}

      /* Produce rotational matrix. */

      kk = idim * (knbpnt - 1);

      for (kp = 0; kp < idim; kp++)
	spnt[kp] = spoint[kk + kp];

      for (ki = 1; ki < knbpnt - 1; ki++)
	{
	  kk = idim * ki;

	  for (kp = 0; kp < idim; kp++)
	    sdum[kp] = spoint[kk + kp];

	  for (kp = 0; kp < idim; kp++)
	    {
	      svect1[kp] = spnt[kp] - spoint[kp];
	      svect2[kp] = sdum[kp] - spoint[kp];
	    }

	  s6crss (svect1, svect2, snorm);

	  s6rotmat (spoint, spnt, snorm, smatrix, &kstat);
	  if (kstat == 0)
	    break;
	}

      /* If we are here with kstat< 0, we have not been able
         to produce rotational matrix. Make straight line. */

      if (kstat < 0)
	goto straight;


      /* Invert the matrix to prepare for rotation of points. */

      s6inv4 (smatrix, sinv, &kstat);
      if (kstat < 0)
	goto error;


      /* Use the inverted matrix on the points and tangents. */

      for (ki = 0; ki < knbpnt; ki++)
	{
	  kk = idim * ki;

	  for (kp = 0; kp < idim; kp++)
	    sdum[kp] = spoint[kk + kp];

	  s6mulvec (sinv, sdum, smul);

	  for (kp = 0; kp < idim; kp++)
	    spoint[kk + kp] = smul[kp];
	}
    }


  /* Test if the points lie on the same branch if the curve
     is a hyperbola. If kstat = 1, different branches, create straight line. */

  s1615 (spoint, knbpnt, idim, sptyp, &kstat);
  if (kstat < 0)
    goto error;
  if (kstat == 1)
    goto straight;


  /* Find the conic equation of the points. */

  s1616 (spoint, knbpnt, idim, sptyp, sconic, &kstat);
  if (kstat < 0)
    goto error;
  if (kstat > 0)
    krem = kstat;

  if (knbpnt == 3)
    {
       
       /* Circle segment */
       
       /* Establish centerpoint. */
       
       centerpt[0] = - sconic[3];
       centerpt[1] = - sconic[4];
       centerpt[2] = (double) 0.0;

      /* If the dimension is 3, rotate the centerpoint
         back into the 3D space and translate. */

      if (idim == 3)
        {
          s6mulvec (smatrix, centerpt, smul);

          for (kp = 0; kp < idim; kp++)
	    centerpt[kp] = smul[kp] + strans[kp];
        }


      /* Establish rotational angle and axis. */
	
      s6diff(sbeg, centerpt, idim, svect1);
      s6diff(slutt, centerpt, idim, svect2);

      sdot = s6scpr(svect1, svect2, idim);
  
      tlength1 = s6length(svect1, idim, &kstat1);
      tlength2 = s6length(svect2, idim, &kstat2);
  
      if (!kstat1 || !kstat2)
        angle = DZERO;
      else
      {
        tcos = sdot/(tlength1*tlength2);
        tcos = MIN((double)1.0,tcos);
        angle = acos(tcos);
      }

      if (idim == 3)
        s6crss(svect1, sdir, axis);
      else if ((svect1[0]*sdir[1] - svect1[1]*sdir[0]) < 0.)
	 angle = -angle;
      
      s6diff(slutt, sbeg, idim, svect1);
      sdot = 0.0;
      for (kp = 0; kp < idim; kp++)
	 sdot += sdir[kp]*svect1[kp];
      if (sdot < 0.0)
      {
	 if (idim == 3)
	    angle = TWOPI - angle;
	 else
	 {
	 if (angle < 0)
	   angle = - (angle + TWOPI);
	 else
	   angle = TWOPI - angle;
	 }
      }

      s1303 (sbeg, aepsge, angle, centerpt, axis, idim, rc, &kstat);
      if (kstat < 0)
        goto error;
    
    }
  else
    {
       
      /* Produce knots for interpolation. */

      s1617 (spoint, knbpnt, idim, sptyp, aepsge, sconic,
	     start, stang, stop, &tshape, &kstat);
      if (kstat < 0)
        goto error;
      if (kstat == 1)
	goto straight;


      /* If the dimension is 3, rotate the points
         back into the 3D space and translate. */

      if (idim == 3)
        {
          s6mulvec (smatrix, start, smul);

          for (kp = 0; kp < idim; kp++)
	    start[kp] = smul[kp] + strans[kp];

          s6mulvec (smatrix, stang, smul);

          for (kp = 0; kp < idim; kp++)
	    stang[kp] = smul[kp] + strans[kp];

          s6mulvec (smatrix, stop, smul);

          for (kp = 0; kp < idim; kp++)
	    stop[kp] = smul[kp] + strans[kp];
        }


      /* Make description of conic curve as
         B-spline within user defined tolerance */

      s1385 (start, stang, stop, tshape, idim, aepsge, rc, &kstat);
      if (kstat < 0)
        goto error;
    
    }


  /* Adjust knot vector to match input start parameter value */

  kk = (*rc)->ik;
  kn = (*rc)->in;

  st = (*rc)->et;

  tdum = astpar - st[kk - 1];

  for (ki = 0; ki < kn + kk; ki++)
    st[ki] += tdum;


  /* Find end parameter value of curve */

  *cendpar = st[kn];

  goto out;

straight:

  /* Create straight line, based on first and last real point. */

  s1602 (sbeg, slutt, ik, idim, astpar, &tpar, rc, &kstat);
  if (kstat < 0)
    goto error;


  /* Find end parameter value of curve */

  *cendpar = tpar;

  goto out;


  /* Error in allocation. */

err101:
  *jstat = -101;
  s6err ("s1611", *jstat, kpos);
  goto out;

  /* Dimension error. */

err104:
  *jstat = -104;
  s6err ("s1611", *jstat, kpos);
  goto out;

  /* Error in lower level routine. */

error:
  *jstat = kstat;
  s6err ("s1611", *jstat, kpos);
  goto out;

out:

  /* Free space occupied by local arrays. */

  if (spoint != SISL_NULL)
    freearray (spoint);
  if (sptyp != SISL_NULL)
    freearray (sptyp);
  if (ieptyp != SISL_NULL)
     freearray (ieptyp);

  if (*jstat >= 0)
     *jstat = krem;
  return;
}
