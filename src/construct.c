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
 * $Id: construct.c,v 1.2 2001-03-19 15:58:39 afr Exp $
 *
 */


#define CONSTRUCT

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
SISLIntpt *copyIntpt (SISLIntpt * ppt)
#else
SISLIntpt *
copyIntpt (ppt)
     SISLIntpt *ppt;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Make a copy of an instance of Intpt, but do not
*              copy the contents of pcurve.
*
*
*
* INPUT      : ppt    - Pointer to the intersection point that is
*                       to be copied.
*
*
*
* OUTPUT     : copyIntpt - New intersection point.
*
*
* METHOD     :
*
*
* REFERENCES :
*
*-
* CALLS      : newIntpt - Make new intersection point.
*
* WRITTEN BY : Vibeke Skytt, SI, 88-06.
*
*********************************************************************
*/
{
  SISLIntpt *qcopy;		/* Local pointer to copied intersection point. */

  /* Create copy.  */

  qcopy = newIntpt (ppt->ipar, ppt->epar, ppt->adist);
  if (qcopy == SISL_NULL)
    goto err101;

  /* Set remaining parameter.  */

  qcopy->iinter = ppt->iinter;

  /* Copy made.  */

  goto out;

  /* Error in space allocation. Return zero.  */

err101:goto out;

out:return (qcopy);
}

#if defined(SISLNEEDPROTOTYPES)
SISLIntpt *hp_copyIntpt (SISLIntpt * ppt)
#else
SISLIntpt *
hp_copyIntpt (ppt)
     SISLIntpt *ppt;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Make a copy of an instance of Intpt, but do not
*              copy the contents of pcurve.
*
*
*
* INPUT      : ppt    - Pointer to the intersection point that is
*                       to be copied.
*
*
*
* OUTPUT     : copyIntpt - New intersection point.
*
*
* METHOD     :
*
*
* REFERENCES :
*
*-
* CALLS      : hp_newIntpt - Make new intersection point.
*
* WRITTEN BY : Ulf J. Krystad, SI,  90-06.
*
*********************************************************************
*/
{
  SISLIntpt *qcopy;		/* Local pointer to copied intersection point. */

  /* Create copy.  */

  qcopy = hp_newIntpt (ppt->ipar, ppt->epar, ppt->adist, ppt->iinter,
		       ppt->left_obj_1[0], ppt->right_obj_1[0],
		       ppt->left_obj_2[0], ppt->right_obj_2[0],
		       ppt->size_1, ppt->size_2, ppt->geo_data_1,
		       ppt->geo_data_2);
  if (qcopy == SISL_NULL)
    goto err101;

  /* Copy made.  */

  goto out;

  /* Error in space allocation. Return zero.  */

err101:goto out;

out:return (qcopy);
}

#if defined(SISLNEEDPROTOTYPES)
SISLbox *
newbox (int idim)
#else
SISLbox *
newbox (idim)
     int idim;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Create and initialize a curve/surface surrounded
*	       SISLbox -instance.
*
*
*
* INPUT      : idim  -  Dimension of geometry space.
*
*
*
* OUTPUT     : newbox - Pointer to new SISLbox structure. If there is
*                       impossible to allocate space for the structure,
*                       newbox returns zero.
*
*
* METHOD     :
*
*
* REFERENCES :
*
*-
* CALLS      :
*
* WRITTEN BY : Arne Laksaa, SI, 89-03.
* REWISED BY : Vibeke Skytt, SI, 91-01.
*
*********************************************************************
*/
{
  SISLbox *qnew;		/* Local pointer to new direction structure.*/
  int ki;			/* Counter.                                 */
  int knum;			/* Number of corners in the box.	          */


  /* Initialise number of corners. */

  if (idim == 3)
    knum = 12;
  else if (idim == 2)
    knum = 4;
  else
    knum = idim;

  /* Allocate space for SISLbox structure.  */

  if ((qnew = newarray (1, SISLbox)) != SISL_NULL)
    {
      /* Initialise new direction structure. */

      qnew->imin = 0;
      qnew->imax = 0;

      /* Initialize arrays.  */

      for (ki = 0; ki < 3; ki++)
	{
	  qnew->e2max[ki] = SISL_NULL;
	  qnew->e2min[ki] = SISL_NULL;
	  qnew->etol[ki] = DZERO;
	}

      if ((qnew->emax = newarray (knum, double)) == SISL_NULL)
	{
	  freearray (qnew);
	  qnew = SISL_NULL;
	}
      else if ((qnew->emin = newarray (knum, double)) == SISL_NULL)
	{
	  freearray (qnew->emax);
	  freearray (qnew);
	  qnew = SISL_NULL;
	}
    }
  return (qnew);
}

#if defined(SISLNEEDPROTOTYPES)
SISLCurve *newCurve (int in, int ik, double *et, double *ecoef,
	  int ikind, int idim, int icopy)
#else
SISLCurve *
newCurve (in, ik, et, ecoef, ikind, idim, icopy)
     int in;
     int ik;
     double *et;
     double *ecoef;
     int ikind;
     int idim;
     int icopy;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Create and initialize a Curve-instance.
*
*
*
* INPUT      : in     - Number of vertices in new curve.
*              ik     - Order of curve.
*              et     - Knotvector of curve.
*              ecoef  - Vertices of curve.
*              ikind  - Kind of curve
*                        = 1 : Polynomial B-spline curve.
*                        = 2 : Rational B-spline curve.
*                        = 3 : Polynomial Bezier curve.
*                        = 4 : Rational Bezier curve.
*              idim   - Dimension of the space in which the curve lies.
*              icopy  - Flag
*                       = 0 : Set pointer to input arrays.
*                       = 1 : Copy input arrays.
*                       = 2 : Set pointer and remember to free arrays.
*
*
*
* OUTPUT     : newCurve - Pointer to new curve. If there is impossible
*                         to allocate space for the curve, newCurve
*                         returns zero.
*
*
* METHOD     :
*              If curve is rational,
*              points input in ecoef should have the form (w*P1,...,w*Pidim,w)
*              The new curve will then have two arrays:
*              qnew -> ecoef with points of the form (P1,...,Pidim)
*              qnew -> rcoef with points of the form (w*P1,...w*Pidim,w).
*
*
* REFERENCES :
*
*-
* CALLS      :
*
* WRITTEN BY  : Vibeke Skytt, SI, 88-05.
* MODIFIED BY : Mike Floater, SI, 91-01 for rational curves.
*
*             : UJK, 92.03.27 Default value set for cuopen SISL_CRV_OPEN
*             : VSK, 92.09.04 Remove superfluous knots and vertices
*                             in the start and end of the curve.
*
*********************************************************************
*/
{
  SISLCurve *qnew;		/* Local pointer to new curve.  */
  int i, j, J, jj, k;		/* loop variables               */
  int k1,k2;                    /* Superflous knots in the ends. */
  int kdim;			/* Dimension of space (also including potential
				   homogenous coordinate        */
  double *st = SISL_NULL;		/* Copy of knotvector.          */
  double *rcoef = SISL_NULL;		/* Copy of vertices in rational case.  */
  double *scoef = SISL_NULL;		/* Copy of vertices.            */


  /* Allocate space for curve.  */

  if ((qnew = newarray (1, SISLCurve)) == SISL_NULL)
    goto err101;

  if (ikind == 2 || ikind == 4)
    kdim = idim + 1;
  else
    kdim = idim;

  /* Count superflous knots in the start.  */
  
  for (k1=0; k1<in; k1++)
     if (et[ik-1] < et[ik+k1]) break;
  
  /* Count superflous knots in the end.  */
  
  for (k2=0; k2<in; k2++)
     if (et[in] > et[in-1-k2]) break;
  
  /* Reduce knots and vertices according to k1 and k2.  */
  
  if (k1 > 0)
  {
     memcopy(ecoef,ecoef+k1*kdim,(in-k1)*kdim,DOUBLE);
     memcopy(et,et+k1,in+ik-k1,DOUBLE);
  }
  in -= (k1+k2);

  /* Check if the curve is still valid. Otherwise return zero. */
  
  if (in < ik) goto err101;
     
  if (icopy == 1)
    {

      /* Copy input arrays. First allocate space for new arrays. */

      if ((st = newarray (in +ik, DOUBLE)) == SISL_NULL ||
	  (scoef = newarray (in *kdim, DOUBLE)) == SISL_NULL)
	goto err101;

      /* Copy contents of arrays.  */

      memcopy (st, et, in +ik, double);
      memcopy (scoef, ecoef, in *kdim, double);
    }
  else
    {
      st = et;
      scoef = ecoef;
    }

  /* Initialize new curve.  */

  qnew->in = in;
  qnew->ik = ik;
  qnew->ikind = ikind;
  qnew->idim = idim;
  qnew->icopy = icopy;
  qnew->et = st;
  qnew->pdir = SISL_NULL;
  qnew->pbox = SISL_NULL;

  if (ikind == 2 || ikind == 4)
    {
      /* Calculate the weighted control points if the object is rational  */
      rcoef = newarray (in *idim, DOUBLE);
      if (rcoef == SISL_NULL)
	goto err101;
      for (i = 0, j = 0, J = 0, k = idim; i < in; i++, k += kdim)
	{
	  for (jj = 0; jj < idim; jj++, j++, J++)
	    {
	      rcoef[J] = scoef[j] / scoef[k];
	    }
	  j++;
	}
      qnew->ecoef = rcoef;
      qnew->rcoef = scoef;
    }
  else
    {
      qnew->ecoef = scoef;
      qnew->rcoef = SISL_NULL;
    }


  /* UJK, 92.03.27 Default value must be set for cuopen */
  qnew->cuopen = SISL_CRV_OPEN;

  /* Task done. */
  goto out;

  /* Error in space allocation. Return zero. */

err101:if (qnew != SISL_NULL)
          { freearray (qnew);  qnew = SISL_NULL;}  
  if (st != SISL_NULL)
    freearray (st);
  if (rcoef != SISL_NULL)
    freearray (rcoef);
  if (scoef != SISL_NULL)
    freearray (scoef);
  goto out;

out:return (qnew);
}

#if defined(SISLNEEDPROTOTYPES)
SISLdir *
newdir (int idim)
#else
SISLdir *
newdir (idim)
     int idim;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Create and initialize a curve/surface direction-instance.
*
*
*
* INPUT      : idim   - Dimension of the space in which the object lies.
*
*
*
* OUTPUT     : newdir - Pointer to new direction structure. If there is
*                       impossible to allocate space for the structure,
*                       newdir returns zero.
*
*
* METHOD     :
*
*
* REFERENCES :
*
*-
* CALLS      :
*
* WRITTEN BY : Arne Laksaa, SI, 89-03.
*
*********************************************************************
*/
{
  SISLdir *qnew;		/* Local pointer to new direction structure.*/

  /* Allocate space for direction structure.  */

  if ((qnew = newarray (1, SISLdir)) != SISL_NULL)
    {
      /* Initialise new direction structure. */

      qnew->igtpi = 0;
      qnew->esmooth = SISL_NULL;
      if ((qnew->ecoef = newarray (idim, double)) == SISL_NULL)
	freearray (qnew);
    }
  return (qnew);
}

#if defined(SISLNEEDPROTOTYPES)
SISLEdge *newEdge (int iedge)
#else
SISLEdge *
newEdge (iedge)
     int iedge;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Create and initialize to zero an instance of the
*              structure Edge. The instance keeps track of intersection
*              points on the edges/endpoints of an object.
*
*
*
* INPUT      : iedge  - Number of edges/endpoints of the object to
*                       which it corresponds.
*                       A curve has iedge = 2.
*                       A tensor-product surface has iedge = 4.
*
*
*
* OUTPUT     : newEdge - A pointer to the instance created.
*
*
* METHOD     :
*
*
* REFERENCES :
*
*-
* CALLS      :
*
* WRITTEN BY : Vibeke Skytt, SI, 88-05.
*
*********************************************************************
*/
{
  int ki;			/* Counter.                          */
  SISLEdge *pnew;		/* Local pointer to the instance.    */
  SISLPtedge *(*ppt);		/* Pointer to traverst array prpt of
			       pointers to Ptedge-instances.     */

  /* Allocate space for instance.  */

  pnew = newarray (1, SISLEdge);
  if (pnew == SISL_NULL)
    goto err101;

  /* Allocate space for array of pointers to Ptedge. */

  pnew->prpt = newarray (iedge, SISLPtedge *);
  if (pnew->prpt == SISL_NULL)
    goto err101;

  /* Initiate the variables of the instance. */

  pnew->iedge = iedge;
  pnew->ipoint = 0;

  ppt = pnew->prpt;
  for (ki = 0; ki < iedge; ki++)
    *(ppt++) = SISL_NULL;

  /* Task done. */

  goto out;

  /* Error in space allocation. Retrun zero. */

err101:pnew = SISL_NULL;
  goto out;

out:return (pnew);
}

#if defined(SISLNEEDPROTOTYPES)
SISLIntcurve *newIntcurve (int ipoint, int ipar1, int ipar2,
			   double *epar1, double *epar2, int itype)
#else
SISLIntcurve *
newIntcurve (ipoint, ipar1, ipar2, epar1, epar2, itype)
     int ipoint;
     int ipar1;
     int ipar2;
     double *epar1;
     double *epar2;
     int itype;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Create and initialize a new instance of the structure
*              Intlist.
*
*
*
* INPUT      : ipoint - Number of points that define the curve.
*              ipar1  - Number of parameter directions of first object
*                       involved in the intersection.
*              ipar2  - Number of parameter direction of second object.
*              epar1  - Parameter values of points in the parameter
*                       area of the first object.
*                       NB! The epar1-pointer is set to point to this
*                           array. The values are not copied.
*              epar2  - Parameter values of points in the parameter
*                       area of the second object.
*                       NB! The epar2-pointer is set to point to this
*                           array. The values are not copied.
*              itype  - Kind of curve. (See intcurve.dcl).
*
*
*
* OUTPUT     : newIntcurve - An instance of the structure Intlist.
*                            If it is not possible to allocate the
*                            necessary space, zero is returned.
*              istat  - status messages
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*
*
* METHOD     :
*
*
* REFERENCES :
*
*-
* CALLS      :
*
* WRITTEN BY : Vibeke Skytt, SI, 88-05.
*
*********************************************************************
*/
{
  SISLIntcurve *qnew;

  /* Allocate space for the new Intcurve.  */

  qnew = newarray (1, SISLIntcurve);
  if (qnew == SISL_NULL)
    goto err101;

  /* Set variables of the intersection curve.  */

  qnew->ipoint = ipoint;
  qnew->ipar1 = ipar1;
  qnew->ipar2 = ipar2;
  qnew->epar1 = epar1;
  qnew->epar2 = epar2;
  qnew->pgeom = SISL_NULL;
  qnew->ppar1 = SISL_NULL;
  qnew->ppar2 = SISL_NULL;
  qnew->itype = itype;

  /* Task done.  */

  goto out;

  /* Error in space allocation. Return zero.  */

err101:qnew = SISL_NULL;
  goto out;

out:return (qnew);
}

#if defined(SISLNEEDPROTOTYPES)
SISLIntdat *newIntdat (void)
#else
SISLIntdat *
newIntdat ()
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Create and initialize to zero an instance of the
*              structure Intdat.
*
*
*
* INPUT
*
*
*
* OUTPUT     : newIntdat - A pointer to the instance created.
*
*
* METHOD     :
*
*
* REFERENCES :
*
*-
* CALLS      :
*
* WRITTEN BY : UJK         , SI, 89-04.
*
*********************************************************************
*/
{
  SISLIntdat *pnew = SISL_NULL;	/* Local pointer to the instance.       */

  /* Allocate space for instance.                                      */

  if ((pnew = newarray (1, SISLIntdat)) != SISL_NULL)
    {
      /* Initiate the variables of the instance.                   */
      pnew->ipmax = 20;
      pnew->ilmax = 10;
      pnew->ipoint = 0;
      pnew->ilist = 0;

      /* Allocate space for array of pointers to Intlist.          */

      if ((pnew->vlist = new0array (pnew->ilmax, SISLIntlist *)) != SISL_NULL)
	{
	  /* Allocate space for array of pointers to SISLIntpt     */
	  if ((pnew->vpoint = new0array (pnew->ipmax, SISLIntpt *))
	      != SISL_NULL) ;

	  /* Task done.                                        */

	  else
	    {
	      /* Error in space allocation of pnew->vpoint.*/
	      freearray (pnew->vlist);
	      freearray (pnew);
	    }
	}
      else
	/* Error in space allocation of pnew->vlist.	     */
	freearray (pnew);
    }
  return pnew;
}

#if defined(SISLNEEDPROTOTYPES)
SISLIntlist *newIntlist (SISLIntpt * pfirst, SISLIntpt * plast, int itype)
#else
SISLIntlist *
newIntlist (pfirst, plast, itype)
     SISLIntpt *pfirst;
     SISLIntpt *plast;
     int itype;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Create and initialize an instance of the structure
*              SISLIntlist that contains a list of intersection points
*              defining an intersection curve.
*
*
*
* INPUT      : pfirst - Pointer to first intersection point in the list.
*              plast  - Pointer to last intersection point in the list.
*              itype  - Kind of curve. See definition of the structure
*                       Intlist.
*
*
*
* OUTPUT     : newIntlist - Pointer to the new instance of Intlist.
*
*
* METHOD     :
*
*
* REFERENCES :
*
*-
* CALLS      :
*
* WRITTEN BY : Vibeke Skytt, SI, 88-05.
*
*********************************************************************
*/
{
  SISLIntlist *pnew;		/* Local pointer to this instance. */

  /* Allocate space for the instance.  */

  pnew = newarray (1, SISLIntlist);
  if (pnew == SISL_NULL)
    goto err101;

  /* Initialize.  */

  pnew->pfirst = pfirst;
  pnew->plast = plast;
  pnew->itype = itype;
  pnew->inumb = 2;

  /* Tast done. */

  goto out;

  /* Error in space allocation. Return zero.  */

err101:pnew = SISL_NULL;
  goto out;

out:return (pnew);
}

#if defined(SISLNEEDPROTOTYPES)
SISLIntpt *newIntpt (int ipar, double *epar, double adist)
#else
SISLIntpt *
newIntpt (ipar, epar, adist)
     int ipar;
     double *epar;
     double adist;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Create and initializes a new instance of the structure
*              Intpt. A pointer to this instance is returned.
*
*
*
* INPUT      : ipar   - Number of parameter directions in the intersection
*                       problem, i.e. dimension of the array epar.
*              epar   - Parameter value of intersection point, possibly
*                       in two objects.
*              adist  - Distance between the objects in the point if
*                       there is two objects in the problem.
*
*
*
* OUTPUT     : newIntpt - Pointer to the new instance of Intpt.
*
*
* METHOD     :
*
*
* REFERENCES :
*
*-
* CALLS      :
*
* WRITTEN BY : Vibeke Skytt, SI, 88-05.
*
*********************************************************************
*/
{
  SISLIntpt *pnew;		/* Local pointer to instance to create. */
  int ki;			/* Counter.                             */

  /* Allocate space for instance of Intpt. */

  pnew = newarray (1, SISLIntpt);
  if (pnew == SISL_NULL)
    goto err101;

  /* Initialize instance. First allocate space for parameter array. */

  if (ipar > 0)
    {
      pnew->epar = newarray (ipar, DOUBLE);
      if (pnew->epar == SISL_NULL)
	goto err101;
    }


  /* Initialize the variables of the instance. */

  pnew->ipar = ipar;
  for (ki = 0; ki < ipar; ki++)
    pnew->epar[ki] = epar[ki];
  pnew->adist = adist;
  pnew->pcurve = SISL_NULL;
  pnew->iinter = 0;

  /* Set intersection atributes to SISL_NULL */
  pnew->no_of_curves_alloc = 0;
  pnew->no_of_curves = 0;

  pnew->pnext = SISL_NULL;
  pnew->curve_dir = SISL_NULL;
  pnew->left_obj_1 = SISL_NULL;
  pnew->left_obj_2 = SISL_NULL;
  pnew->right_obj_1 = SISL_NULL;
  pnew->right_obj_2 = SISL_NULL;
  pnew->geo_data_1 = SISL_NULL;
  pnew->size_1 = 0;
  pnew->geo_data_2 = SISL_NULL;
  pnew->size_2 = 0;

  pnew->trim[0] = SISL_NULL;
  pnew->trim[1] = SISL_NULL;

  /* Task done.  */


  goto out;

  /* Error in space allocation. Return zero. */

err101:pnew = SISL_NULL;
  goto out;

out:return (pnew);
}

#if defined(SISLNEEDPROTOTYPES)
SISLIntsurf *newIntsurf (SISLIntlist * pintlist)
#else
SISLIntsurf *
newIntsurf (pintlist)
     SISLIntlist *pintlist;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Create a new instance of the structure
*              Intsurf.
*
*
*
* INPUT      : pintlist - Pointer to the surround curve.
*
*
* OUTPUT     : newIntsurf - Pointer to the new instance,
*
*
* METHOD     :
*
*
* REFERENCES :
*
*-
* CALLS      :
*
* WRITTEN BY : Ulf J. Krystad, SI, 91-10.
*
*********************************************************************
*/
{
  SISLIntsurf *pnew = SISL_NULL;	/* Local pointer to instance to create. */
  SISLIntpt *qpt = SISL_NULL;	/* Local help pointer                   */
  SISLIntpt *qpfirst = SISL_NULL;	/* Local help pointer                   */
  SISLIntpt *qplast = SISL_NULL;	/* Local help pointer                   */
  SISLIntpt *qprev = SISL_NULL;	/* Local help pointer                   */
  SISLIntpt *qnext = SISL_NULL;	/* Local help pointer                   */
  int index, ipar, ipoint;
  int ki, kk, kdir;		/* Counter.                             */
  int dummy,kstat;
  double *stpar1, *stpar2;
  /* ------------------------------------------------------------------ */

  if (pintlist == SISL_NULL)
    goto out;

  qpfirst = pintlist->pfirst;
  qplast  = pintlist->plast;
  ipoint  = pintlist->inumb - 1;
  ipar    = qpfirst->ipar;
  index   = pintlist->ind_first;


  if (ipar <= 0)
    goto out;
  if (ipoint <= 1)
    goto out;

  /* Allocate space for instance of Intsurf. */
  pnew = newarray (1, SISLIntsurf);
  if (pnew == SISL_NULL)
    goto err101;

  pnew->ipar = ipar;
  pnew->ipoint = ipoint;

  /* First allocate space for parameter array. */
  pnew->epar = stpar1 = newarray (ipar * ipoint, DOUBLE);
  if (pnew->epar == SISL_NULL)
    goto err101;

  /* Allocate space for constant direction array. */
  /* UJK, sept 92 */
  /* pnew->const_par = newarray (ipar, int); */
  pnew->const_par = newarray (ipoint, int);
  if (pnew->const_par == SISL_NULL)
    goto err101;

  /* Fill in arrays */

  qpt = qprev = qpfirst;
  qnext = qpt->pnext[index];

    for (ki = 0; ki < ipoint; ki++)
    {
      qpt->marker = -99;
      stpar2 = qpt->epar;
      for (kk = 0; kk < ipar; kk++)
	*(stpar1++) = *(stpar2++);

      for (kdir = 0; kdir < ipar; kdir++)
	if (qpt->curve_dir[index] &
	    (1 << (kdir + 1)))
	  break;

      pnew->const_par[ki] = kdir;

      /* Next point */
      qprev = qpt;
      qpt = qnext;
      sh6getother (qpt, qprev, &qnext, &kstat);

      sh6getlist (qpt, qnext, &index, &dummy, &kstat);
    }

  /* Task done.  */


  goto out;

  /* Error in space allocation. Return zero. */

err101:pnew = SISL_NULL;
  goto out;

out:return (pnew);
}





#if defined(SISLNEEDPROTOTYPES)
SISLTrimpar *newTrimpar (int pt, int par)
#else
SISLTrimpar *
newTrimpar (pt, par)
     int pt;
     int par;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Create and initializes a new instance of the structure
*              SISLTrimpar. A pointer to this instance is returned.
*
*
*
* INPUT      : pt     - Index in pnext to a pointer to an adjacent
*                       point in the trim curve.
*              par    - An integer between 0 and 3 indicating
*                          which parameter is constant on the
*                          trim curve joining ptindex.
*
*
* OUTPUT     : newTrimpar - Pointer to the new instance of Trimpar.
*
*
* METHOD     :
*
*
* REFERENCES :
*
*-
* CALLS      :
*
* WRITTEN BY : Michael Floater, SI, 91-10.
*
*********************************************************************
*/
{
  SISLTrimpar *newtrim;		/* Local pointer to instance to create. */

  /* Allocate space for instance of Trimpar. */

  newtrim = newarray (1, SISLTrimpar);
  if (newtrim == SISL_NULL)
    goto err101;

  /* Initialize the variables of the instance. */

  newtrim->ptindex = pt;
  newtrim->parindex = par;

  /* Task done.  */


  goto out;

  /* Error in space allocation. Return zero. */

err101:newtrim = SISL_NULL;
  goto out;

out:return (newtrim);
}

#if defined(SISLNEEDPROTOTYPES)
SISLIntpt *
hp_newIntpt (int ipar, double *epar, double adist, int itype,
	     int ileft1, int iright1, int ileft2, int iright2,
	     int size_1, int size_2, double egeom1[], double egeom2[])
#else
SISLIntpt *
hp_newIntpt (ipar, epar, adist, itype, ileft1, iright1, ileft2, iright2,
	     size_1, size_2, egeom1, egeom2)
     int ipar;
     double *epar;
     double adist;
     int itype;
     int ileft1;
     int iright1;
     int ileft2;
     int iright2;
     int size_1;
     int size_2;
     double *egeom1;
     double *egeom2;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Create and initializes a new instance of the structure
*              Intpt. A pointer to this instance is returned.
*
*
*
* INPUT      : ipar   - Number of parameter directions in the intersection
*                       problem, i.e. dimension of the array epar.
*              epar   - Parameter value of intersection point, possibly
*                       in two objects.
*              adist  - Distance between the objects in the point if
*                       there is two objects in the problem.
*              itype  - Type of intersection point.
*              ileft1, iright1, ileft2, iright2 - Pre-topology information
*                       conserning interpolation point.
*              size_1 - length of egeom1.
*              size_2 - length of egeom2.
*              egeom1 - Geometry information of intersection point in
*                       first object.
*              egeom2 - Geometry information of intersection point in
*                       second object.
*
*
*
* OUTPUT     : hp_newIntpt - Pointer to the new instance of Intpt.
*
*
* METHOD     :
*
*
* REFERENCES :
*
*-
* CALLS      :
*
* WRITTEN BY : Ulf J. Krystad, SI, 91-06.
*
*********************************************************************
*/
{
  SISLIntpt *pnew;		/* Local pointer to instance to create. */
  int ki;			/* Counter.                             */


  /* Allocate space for instance of Intpt. */

  pnew = new0array (1, SISLIntpt);
  if (pnew == SISL_NULL)
    goto err101;

  /* Initialize parameters concerning the size of arrays describing
     the intersection curves in which this points lie, and allocate
     scratch for the arrays.      */

  pnew->no_of_curves_alloc = 4;
  pnew->no_of_curves = 0;

  if ((pnew->pnext = newarray (pnew->no_of_curves_alloc, SISLIntpt *)) == SISL_NULL)
    goto err101;
  if ((pnew->curve_dir = newarray (pnew->no_of_curves_alloc, INT)) == SISL_NULL)
    goto err101;
  if ((pnew->left_obj_1 = newarray (pnew->no_of_curves_alloc, INT)) == SISL_NULL)
    goto err101;
  if ((pnew->left_obj_2 = newarray (pnew->no_of_curves_alloc, INT)) == SISL_NULL)
    goto err101;
  if ((pnew->right_obj_1 = newarray (pnew->no_of_curves_alloc, INT)) == SISL_NULL)
    goto err101;
  if ((pnew->right_obj_2 = newarray (pnew->no_of_curves_alloc, INT)) == SISL_NULL)
    goto err101;

  /* Initialize instance. First allocate space for parameter array. */

  pnew->epar = SISL_NULL;
  if (ipar > 0)
    {
      pnew->epar = newarray (ipar, DOUBLE);
      if (pnew->epar == SISL_NULL)
	goto err101;
    }


  /* Initialize the variables of the instance. */

  pnew->ipar = ipar;
  for (ki = 0; ki < ipar; ki++)
    pnew->epar[ki] = epar[ki];
  pnew->adist = adist;
  pnew->pcurve = SISL_NULL;
  pnew->edge_1 = 0;
  pnew->edge_2 = 0;
  pnew->iinter = itype;
  pnew->marker = 0;
  pnew->evaluated = 0;

  if (size_1 > 0)
    {
      pnew->geo_data_1 = newarray (size_1, DOUBLE);
      pnew->size_1 = size_1;
      memcopy (pnew->geo_data_1, egeom1, size_1, double);
    }
  else
    {
      pnew->geo_data_1 = SISL_NULL;
      pnew->size_1 = 0;
    }

  if (size_2 > 0)
    {
      pnew->geo_data_2 = newarray (size_2, DOUBLE);
      pnew->size_2 = size_2;
      memcopy (pnew->geo_data_2, egeom2, size_2, double);
    }
  else
    {
      pnew->geo_data_2 = SISL_NULL;
      pnew->size_2 = 0;
    }


  *(pnew->left_obj_1) = ileft1;
  *(pnew->left_obj_2) = ileft2;
  *(pnew->right_obj_1) = iright1;
  *(pnew->right_obj_2) = iright2;

  for (ki = 0; ki < pnew->no_of_curves_alloc; ki++)
    pnew->pnext[ki] = SISL_NULL;

  /*  for (ki=0; ki<6; ki++) pnew -> geo_aux[ki] = DZERO; */

  pnew->trim[0] = SISL_NULL;
  pnew->trim[1] = SISL_NULL;

  /* Init the left/right evaluator to default value zero. */
  pnew->iside_1 = 0;
  pnew->iside_2 = 0;
  
  /* Task done.  */

  goto out;

  /* Error in space allocation. Return zero. */

err101:pnew = SISL_NULL;
  goto out;

out:return (pnew);
}

#if defined(SISLNEEDPROTOTYPES)
SISLTrack *
newTrack (SISLSurf * psurf_1, SISLSurf * psurf_2, SISLCurve * pcurve_3d,
	  SISLCurve * pcurve_2d_1, SISLCurve * pcurve_2d_2,
	int ideg, double eimpli[], int sing_start, int sing_end, int turned)
#else
SISLTrack *
newTrack (psurf_1, psurf_2, pcurve_3d,
	  pcurve_2d_1, pcurve_2d_2,
	  ideg, eimpli, sing_start, sing_end, turned)
     SISLSurf *psurf_1;
     SISLSurf *psurf_2;
     SISLCurve *pcurve_3d;

     SISLCurve *pcurve_2d_1;
     SISLCurve *pcurve_2d_2;

     int ideg;
     double eimpli[];
     int sing_start;
     int sing_end;
     int turned;

#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Create and initializes a new instance of the structure
*              SISLTrack. A pointer to this instance is returned.
*
*
*
* INPUT      : Compatible with the attributes in the structure track
*
*
*
* OUTPUT     : newTrack - Pointer to the new instance of SISLTrack
*
*
* REFERENCES :
*
*-
* CALLS      :
*
* WRITTEN BY : Ulf J. Krystad, SI, 91-06.
*
*********************************************************************
*/
{
  SISLTrack *pnew;		/* Local pointer to instance to create. */


  /* Allocate space for instance of Intpt. */

  pnew = new0array (1, SISLTrack);
  if (pnew == SISL_NULL)
    goto err101;

  /* Set atributes */


  pnew->psurf_1 = psurf_1;
  pnew->psurf_2 = psurf_2;
  pnew->pcurve_3d = pcurve_3d;
  pnew->pcurve_2d_1 = pcurve_2d_1;
  pnew->pcurve_2d_2 = pcurve_2d_2;
  pnew->ideg = ideg;
  if (pnew->ideg)
    memcopy (pnew->eimpli, eimpli, 16, DOUBLE);
  pnew->sing_start = sing_start;
  pnew->sing_end = sing_end;
  pnew->turned = turned;

  /* Task done.  */

  goto out;

  /* Error in space allocation. Return zero. */

err101:pnew = SISL_NULL;
  goto out;

out:return (pnew);
}

#if defined(SISLNEEDPROTOTYPES)
SISLObject *newObject (int iobj)
#else
SISLObject *
newObject (iobj)
     int iobj;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Create and initialize partly an Object-instance.
*
*
*
* INPUT      : iobj   - Kind of object.
*
*
*
* OUTPUT     : newObject - Pointer to new object. If there is impossible
*                          to allocate space for the object, newObject
*                          returns zero.
*              istat     - status messages
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*
*
* METHOD     :
*
*
* REFERENCES :
*
*-
* CALLS      :
*
* WRITTEN BY : Vibeke Skytt, SI, 88-05.
* MODIFIED BY: UJK               89-04.
*              The structure SISLPoint is included.
*              The pointer to structure object itself is included.
*
*********************************************************************
*/
{
  SISLObject *qnew;		/* Local pointer to new object.  */

  qnew = newarray (1, SISLObject);
  if (qnew == SISL_NULL)
    goto out;

  qnew->iobj = iobj;

  qnew->p1 = SISL_NULL;
  qnew->c1 = SISL_NULL;
  qnew->s1 = SISL_NULL;
  qnew->o1 = SISL_NULL;
  qnew->edg[0] = SISL_NULL;
  qnew->edg[1] = SISL_NULL;
  qnew->edg[2] = SISL_NULL;
  qnew->edg[3] = SISL_NULL;
  qnew->psimple = SISL_NULL;

  /* Task done.  */

out:return qnew;
}

#if defined(SISLNEEDPROTOTYPES)
SISLPoint *newPoint (double *ecoef, int idim, int icopy)
#else
SISLPoint *
newPoint (ecoef, idim, icopy)
     double *ecoef;
     int idim;
     int icopy;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Create and initialize a Point-instance.
*
*
*
* INPUT      : ecoef  - Coordinates of point.
*              idim   - Dimension of the space in which the point lies.
*              icopy  - Flag
*                       = 0 : Set pointer to input array.
*                       = 1 : Copy input arrays.
*                       = 2 : Set pointer and remember to free array.
*
*
*
* OUTPUT     : newPoint - Pointer to new point. If there is impossible
*                         to allocate space for the point, newPoint
*                         returns zero.
*
*
* METHOD     :
*
*
* REFERENCES :
*
*-
* CALLS      :
*
* WRITTEN BY : UJK, SI. 89-04.
*
*********************************************************************
*/
{
  SISLPoint *qnew;		/* Local pointer to new point.  */
  double *scoef = SISL_NULL;		/* Copy of coordinates          */


  /* Allocate space for point.  */

  if ((qnew = newarray (1, SISLPoint)) == SISL_NULL)
    goto err101;


  if (icopy == 1)
    {

      /* Copy input array. First allocate space for new array. */

      if (idim < 4) 	scoef = qnew->ec;
      else if ((scoef = newarray (idim, double)) == SISL_NULL) goto err101;

      /* Copy contents of arrays.  */

      memcopy (scoef, ecoef, idim, double);
    }
  else
    {
      scoef = ecoef;
    }

  /* Initialize new point.  */

  qnew->idim = idim;
  qnew->icopy = icopy;
  qnew->ecoef = scoef;
  qnew->pbox = SISL_NULL;

  /* Task done. */

  goto out;

  /* Error in space allocation. Return zero. */

err101:if (qnew != SISL_NULL)
    freearray (qnew);
  qnew = SISL_NULL;
  goto out;
out:return qnew;
}

#if defined(SISLNEEDPROTOTYPES)
SISLPtedge *newPtedge (SISLIntpt * ppt)
#else
SISLPtedge *
newPtedge (ppt)
     SISLIntpt *ppt;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Create and initialize an instance of the structure
*              Ptedge, that contains elements in a list where each
*              element has a pointer to an intersection point.
*
*
*
* INPUT      : ppt    - The intersection point at which this instance
*                       point.
*
*
*
* OUTPUT     : newPtedge - Pointer to the crated instance.
*
*
* METHOD     :
*
*
* REFERENCES :
*
*-
* CALLS      :
*
* WRITTEN BY : Vibeke Skytt, SI, 88-05.
*
*********************************************************************
*/
{
  SISLPtedge *pnew;		/* Local pointer to the instance to create. */

  /* Allocate space for the instance. */

  pnew = newarray (1, SISLPtedge);
  if (pnew == SISL_NULL)
    goto err101;

  /* Initialize instance.  */

  pnew->ppt = ppt;
  pnew->pnext = SISL_NULL;

  /* Task done.  */

  goto out;

  /* Error in space allocation. Return zero. */

err101:pnew = SISL_NULL;
  goto out;

out:return (pnew);
}

#if defined(SISLNEEDPROTOTYPES)
SISLSurf *
newSurf (int in1, int in2, int ik1, int ik2, double *et1, double *et2,
	 double *ecoef, int ikind, int idim, int icopy)
#else
SISLSurf *
newSurf (in1, in2, ik1, ik2, et1, et2, ecoef, ikind, idim, icopy)
     int in1;
     int in2;
     int ik1;
     int ik2;
     double *et1;
     double *et2;
     double *ecoef;
     int ikind;
     int idim;
     int icopy;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Create and initialize a surface (instance of SISLSurf).
*
*
*
* INPUT      : in1    - Number of vertices in first parameter direction
*                       of new surface.
*            : in2    - Number of vertices in second parameter direction
*                       of new surface.
*              ik1    - Order of surface in first parameter direction.
*              ik2    - Order of surface in second parameter direction.
*              et1    - Knotvector of surface in first parameter direction.
*              et2    - Knotvector of surface in second parameter direction.
*              ecoef  - Vertices of surface.
*              ikind  - Kind of surface
*                        = 1 : Polynomial B-spline surface.
*                        = 2 : Rational B-spline surface.
*                        = 3 : Polynomial Bezier surface.
*                        = 4 : Rational Bezier surface.
*              idim   - Dimension of the space in which the surface lies.
*              icopy  - Flag
*                       = 0 : Set pointer to input arrays.
*                       = 1 : Copy input arrays.
*                       = 2 : Set pointer and mark to free.
*
*
*
* OUTPUT     : newSurf  - Pointer to new surface. If there is impossible
*                         to allocate space for the surface, newSurf
*                         returns zero.
*
*
* METHOD     :
*              If surface is rational,
*              points input in ecoef should have the form (w*P1,...,w*Pidim,w)
*              The new surface will then have two arrays:
*              qnew -> ecoef with points of the form (P1,...,Pidim)
*              qnew -> rcoef with points of the form (w*P1,...w*Pidim,w).
*
*
* REFERENCES :
*
*-
* CALLS      :
*
* WRITTEN BY : Vibeke Skytt, SI, 88-05.
* MODIFIED BY : Mike Floater, SI, 91-01 for rational surfaces.
*
*********************************************************************
*/
{
  SISLSurf *qnew;		/* Local pointer to new surface.       */
  int i, j, J, jj, k;		/* loop variables                      */
  int k1, k2;                   /* Superfluous knots at the ends.      */
  int kdim;			/* Dimension indicator.                */
  double *st1 = SISL_NULL, *st2 = SISL_NULL;	/* Copy of knotvectors.        */
  double *rcoef = SISL_NULL;		/* Copy of vertices in rational case.  */
  double *scoef = SISL_NULL;		/* Copy of vertices.                   */
  double *ucoef = SISL_NULL;         /* Utility coefficient array.          */

  /* Allocate space for surface.  */

  if ((qnew = newarray (1, SISLSurf)) == SISL_NULL)
    goto err101;

  if (ikind == 2 || ikind == 4)
    kdim = idim + 1;
  else
    kdim = idim;

  /* Count superfluous knots at ends, first in u parameter directions. */

  if (ik1 == 0)
  {
    k1 = k2 = 0;
  }
  else
  {
    /* Count superfluous knots in the start. */

    for (k1 = 0; k1 < in1; k1++)
      if (et1[ik1 - 1] < et1[ik1 + k1]) break;

    /* Count superfluous knots in the end. */

    for (k2 = 0; k2 < in1; k2++)
      if (et1[in1] > et1[in1 - 1 - k2]) break;
  }

  /* Reduce knots and vertices according to k1 and k2. */

  if (k1 > 0 || k2 > 0)
  {
     ucoef = newarray(in1*in2*kdim, DOUBLE);
     s6chpar(ecoef, in1, in2, kdim, ucoef);
  }
  if (k1 > 0)
  {
     memcopy(ucoef, ucoef + k1*in2*kdim, (in1 - k1)*in2*kdim, DOUBLE);
     memcopy(et1, et1 + k1, in1 + ik1 - k1, DOUBLE);
  }
  in1 -= (k1 + k2);
  if (k1 > 0 || k2 > 0)
  {
     s6chpar(ucoef, in2, in1, kdim, ecoef);
     if (ucoef != SISL_NULL) freearray(ucoef);
  }

  /* Count superfluous knots at ends in v parameter directions. */

  if (ik2 == 0)
  {
    k1 = k2 = 0;
  }
  else
  {
    /* Count superfluous knots in the start. */

    for (k1 = 0; k1 < in2; k1++)
     if (et2[ik2 - 1] < et2[ik2 + k1]) break;

    /* Count superfluous knots in the end. */

    for (k2 = 0; k2 < in2; k2++)
      if (et2[in2] > et2[in2 - 1 - k2]) break;
  }
  /* Reduce knots and vertices according to k1 and k2. */

  if (k1 > 0)
  {
     memcopy(ecoef, ecoef + k1*in1*kdim, (in2 - k1)*in1*kdim, DOUBLE);
     memcopy(et2, et2 + k1, in2 + ik2 - k1, DOUBLE);
  }
  in2 -= (k1 + k2);

  if (icopy == 1)
    {

      /* Copy input arrays. First allocate space for new arrays. */

      st1 = newarray (in1 + ik1, DOUBLE);
      st2 = newarray (in2 + ik2, DOUBLE);
      scoef = newarray (in1 * in2 * kdim, DOUBLE);
      if (st1 == SISL_NULL || st2 == SISL_NULL || scoef == SISL_NULL)
	goto err101;

      /* Copy contents of arrays.  */
      memcopy (st1, et1, in1 + ik1, double);
      memcopy (st2, et2, in2 + ik2, double);
      memcopy (scoef, ecoef, in1 * in2 * kdim, double);
    }
  else
    {
      st1 = et1;
      st2 = et2;
      scoef = ecoef;
    }

  /* Initialize new surface. */

  qnew->in1 = in1;
  qnew->in2 = in2;
  qnew->ik1 = ik1;
  qnew->ik2 = ik2;
  qnew->ikind = ikind;
  qnew->idim = idim;
  qnew->icopy = icopy;
  qnew->et1 = st1;
  qnew->et2 = st2;
  qnew->pdir = SISL_NULL;
  qnew->pbox = SISL_NULL;

  if (ikind == 2 || ikind == 4)
    {
      /* Calculate the weighted control points if the object is rational  */
      rcoef = newarray (in1 * in2 * idim, DOUBLE);
      if (rcoef == SISL_NULL)
	goto err101;
      for (i = 0, j = 0, J = 0, k = idim; i < in1 * in2; i++, k += kdim)
	{
	  for (jj = 0; jj < idim; jj++, j++, J++)
	    {
	      rcoef[J] = scoef[j] / scoef[k];
	    }
	  j++;
	}
      qnew->ecoef = rcoef;
      qnew->rcoef = scoef;
    }
  else
    {
      qnew->ecoef = scoef;
      qnew->rcoef = SISL_NULL;
    }
  
  /* UJK, 92.05.05 Default value must be set for cuopen */
  qnew->cuopen_1 = SISL_SURF_OPEN;
  qnew->cuopen_2 = SISL_SURF_OPEN;

  /* Task done. */

  goto out;

  /* Error in space allocation. Return zero. */

err101:if (qnew != SISL_NULL)
    freearray (qnew);
  if (st1 != SISL_NULL)
    freearray (st1);
  if (st2 != SISL_NULL)
    freearray (st2);
  if (rcoef != SISL_NULL)
    freearray (rcoef);
  if (scoef != SISL_NULL)
    freearray (scoef);
  goto out;

out:return (qnew);
}

#if defined(SISLNEEDPROTOTYPES)
SISLCurve *copyCurve (SISLCurve * pcurve)
#else
SISLCurve *
copyCurve (pcurve)
     SISLCurve *pcurve;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Make a copy of a SISLCurve.
*
*
*
* INPUT      : pcurve  -  Curve to be copied.
*
*
*
* OUTPUT     : copyCurve - The new curve.
*
*
* METHOD     :
*
*
* REFERENCES :
*
*-
* CALLS      : newCurve - Make new SISLCurve.
*
* WRITTEN BY : Vibeke Skytt, SI, 93-11.
*
*********************************************************************
*/
{
  int kstat = 0;
  SISLCurve *qc = SISL_NULL;
  int knum;
  int ki;
  
  /* Make new curve. */

  if (pcurve->ikind == 2 || pcurve->ikind == 4)
  {
     if ((qc = newCurve(pcurve->in, pcurve->ik, pcurve->et, pcurve->rcoef,
			pcurve->ikind, pcurve->idim, 1)) == SISL_NULL) goto err101;
  }
  else
  {
     if ((qc = newCurve(pcurve->in, pcurve->ik, pcurve->et, pcurve->ecoef,
			pcurve->ikind, pcurve->idim, 1)) == SISL_NULL) goto err101;
  }
  
  /* Set the open/closed flag. */
  
  qc->cuopen = pcurve->cuopen;
  
  /* Copy box and  cone information. */
  
  if (pcurve->pbox != SISL_NULL)
  {
     /* Copy box information. */
     
     if ((qc->pbox = newbox(pcurve->idim)) == SISL_NULL) 
     {
	freeCurve(qc);
	qc = SISL_NULL;
	goto err101;
     }
     if (pcurve->idim == 3) knum = 9;
     else if (pcurve->idim == 2) knum = 4;
     else knum = pcurve->idim;
     
     memcopy(qc->pbox->emin, pcurve->pbox->emin, knum, DOUBLE);
     memcopy(qc->pbox->emax, pcurve->pbox->emax, knum, DOUBLE);
     memcopy(qc->pbox->etol, pcurve->pbox->etol, 3, DOUBLE);

     for (ki=0; ki<3; ki++)
     {
	if (s6existbox(pcurve->pbox, ki, pcurve->pbox->etol[ki]))
	{
	   s6newbox(qc->pbox, knum, ki, pcurve->pbox->etol[ki], &kstat);
	   if (kstat < 0)
	   {
	      freeCurve(qc);
	      qc = SISL_NULL;
	      goto err101;
	   }
	   memcopy(qc->pbox->e2min[ki], pcurve->pbox->e2min[ki], knum,
		   DOUBLE);
	   memcopy(qc->pbox->e2max[ki], pcurve->pbox->e2max[ki], knum,
		   DOUBLE);
	}
     }
  }
  
  if (pcurve->pdir != SISL_NULL)
  {
     /* Copy cone information. */
     
     if ((qc->pdir = newdir(pcurve->idim)) == SISL_NULL)
     {
	freeCurve(qc);
	qc = SISL_NULL;
	goto err101;
     }
     qc->pdir->igtpi = pcurve->pdir->igtpi;
     qc->pdir->aang = pcurve->pdir->aang;
     memcopy(qc->pdir->ecoef, pcurve->pdir->ecoef, pcurve->idim, DOUBLE);
     
     if (pcurve->pdir->esmooth != SISL_NULL)
     {
	if ((qc->pdir->esmooth = newarray(qc->in*qc->idim,DOUBLE)) == SISL_NULL)
	{
	   freeCurve(qc);
	   qc = SISL_NULL;
	   goto err101;
	}
	memcopy(qc->pdir->esmooth, pcurve->pdir->esmooth, qc->in*qc->idim,
		DOUBLE);
     }
  }
  
  goto out;
  
  /* Error in allocation. */
  
  err101 : 
     goto out;
  
  out :
     return qc;
}

#if defined(SISLNEEDPROTOTYPES)
SISLSurf *copySurface (SISLSurf * psurf)
#else
SISLSurf *
copySurface (psurf)
     SISLSurf *psurf;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Make a copy of a SISLSurface.
*
*
*
* INPUT      : psurf  -  Surface to be copied.
*
*
*
* OUTPUT     : copySurface - The new curve.
*
*
* METHOD     :
*
*
* REFERENCES :
*
*-
* CALLS      : newSurf - Make new SISLSurface.
*
* WRITTEN BY : Vibeke Skytt, SI, 93-11.
*
*********************************************************************
*/
{
  int kstat = 0;
  SISLSurf *qs = SISL_NULL;
  int knum;
  int ki;
  
  /* Make new curve. */

  if (psurf->ikind == 2 || psurf->ikind == 4)
  {
     if ((qs = newSurf(psurf->in1, psurf->in2, psurf->ik1, psurf->ik2,
		       psurf->et1, psurf->et2, psurf->rcoef,
		       psurf->ikind, psurf->idim, 1)) == SISL_NULL) goto err101;
  }
  else
  {
     if ((qs = newSurf(psurf->in1, psurf->in2, psurf->ik1, psurf->ik2,
		       psurf->et1, psurf->et2, psurf->ecoef,
		       psurf->ikind, psurf->idim, 1)) == SISL_NULL) goto err101;
  }
  
  /* Set the open/closed flag. */
  
  qs->cuopen_1 = psurf->cuopen_1;
  qs->cuopen_2 = psurf->cuopen_2;
  
  /* Copy box and  cone information. */
  
  if (psurf->pbox != SISL_NULL)
  {
     /* Copy box information. */
     
     if ((qs->pbox = newbox(psurf->idim)) == SISL_NULL) 
     {
	freeSurf(qs);
	qs = SISL_NULL;
	goto err101;
     }
     if (psurf->idim == 3) knum = 9;
     else if (psurf->idim == 2) knum = 4;
     else knum = psurf->idim;
     
     memcopy(qs->pbox->emin, psurf->pbox->emin, knum, DOUBLE);
     memcopy(qs->pbox->emax, psurf->pbox->emax, knum, DOUBLE);
     memcopy(qs->pbox->etol, psurf->pbox->etol, 3, DOUBLE);

     for (ki=0; ki<3; ki++)
     {
	if (s6existbox(psurf->pbox, ki, psurf->pbox->etol[ki]))
	{
	   s6newbox(qs->pbox, knum, ki, psurf->pbox->etol[ki], &kstat);
	   if (kstat < 0)
	   {
	      freeSurf(qs);
	      qs = SISL_NULL;
	      goto err101;
	   }
	   memcopy(qs->pbox->e2min[ki], psurf->pbox->e2min[ki], knum,
		   DOUBLE);
	   memcopy(qs->pbox->e2max[ki], psurf->pbox->e2max[ki], knum,
		   DOUBLE);
	}
     }
  }
  
  if (psurf->pdir != SISL_NULL)
  {
     /* Copy cone information. */
     
     if ((qs->pdir = newdir(psurf->idim)) == SISL_NULL)
     {
	freeSurf(qs);
	qs = SISL_NULL;
	goto err101;
     }
     qs->pdir->igtpi = psurf->pdir->igtpi;
     qs->pdir->aang = psurf->pdir->aang;
     memcopy(qs->pdir->ecoef, psurf->pdir->ecoef, psurf->idim, DOUBLE);
     
     if (psurf->pdir->esmooth != SISL_NULL)
     {
	if ((qs->pdir->esmooth = newarray(qs->in1*qs->in2*qs->idim,DOUBLE)) 
	    == SISL_NULL)
	{
	   freeSurf(qs);
	   qs = SISL_NULL;
	   goto err101;
	}
	memcopy(qs->pdir->esmooth, psurf->pdir->esmooth, 
		qs->in1*qs->in2*qs->idim, DOUBLE);
     }
  }
  
  goto out;
  
  /* Error in allocation. */
  
  err101 : 
     goto out;
  
  out :
     return qs;
}   
