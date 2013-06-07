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
 * $Id: destruct.c,v 1.3 2001-03-19 15:58:40 afr Exp $
 *
 */


#define DESTRUCT

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
     freeCurve(SISLCurve *pcurve)
#else
void freeCurve(pcurve)
     SISLCurve *pcurve;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Free the space occupied by the curve pointed at by
*              pcurve.
*
*
*
* INPUT      : pcurve - Pointer to curve.
*
*
*
* OUTPUT     :
*
*
* METHOD     :
*
*
* REFERENCES :
*
*-
* CALLS      : freearray - Free space occupied by a given array.
*
* WRITTEN BY : Vibeke Skytt, SI, 88-05. Arne Laksaa, SI, 89-03.
*              Vibeke Skytt, SI, 91-01.
* Revised by : Paal Fugelli, SINTEF, Oslo, Norway, 94-09.  Added free'ing
*              of ecoef for rationals when icopy==0.  ecoef is ALWAYS
*              allocated for rationals, regardless of the icopy flag value.
*
*********************************************************************
*/
{
  int ki;         /* Counter.  */


  if ( pcurve->icopy != 0 )
  {

    /* Free arrays.  */

    freearray(pcurve->et);
    freearray(pcurve->ecoef);
    if ( pcurve->rcoef != SISL_NULL )  freearray(pcurve->rcoef);
  }
  else if ( pcurve->ikind == 2  ||  pcurve->ikind == 4 )
  {
    /* VSK, 940902. The array pcurve->ecoef is ALWAYS allocated in the
       constructor for rationals and must be free'ed.                  */

    freearray(pcurve->ecoef);
  }

  if ( pcurve->pdir )
  {

    /* Free direction structure. */

    if ( pcurve->pdir->ecoef )  freearray(pcurve->pdir->ecoef);
    if ( pcurve->pdir->esmooth != SISL_NULL )  freearray(pcurve->pdir->esmooth);
    freearray(pcurve->pdir);
  }

  if ( pcurve->pbox )
  {

    /* Free surrounded SISLbox structure. */

    if ( pcurve->pbox->emax )  freearray(pcurve->pbox->emax);
    if ( pcurve->pbox->emin )  freearray(pcurve->pbox->emin);

    for ( ki=0;  ki < 3;  ki++ )
    {
      if ( pcurve->pbox->e2max[ki] )  freearray(pcurve->pbox->e2max[ki]);
      if ( pcurve->pbox->e2min[ki] )  freearray(pcurve->pbox->e2min[ki]);
    }
    freearray(pcurve->pbox);
  }

  /* Free instance of curve. */

  freearray(pcurve);

  return;
}

#if defined(SISLNEEDPROTOTYPES)
void
     freeEdge(SISLEdge *pedge)
#else
void freeEdge(pedge)
     SISLEdge *pedge;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Free the space occupied by an instance of the structure
*              Edge. Also free space occupied by instances of Ptedge
*              that belongs to pedge.
*
*
*
* INPUT      : pedge  - Pointer to the instance.
*
*
*
* OUTPUT     :
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

  SISLPtedge *p1,*p2;  /* Pointers to traverse lists of Ptedge-elements.*/
  SISLPtedge *(*pel);  /* Pointer to an array element.      */
  int ki;                 /* Counter.                          */

  /* First free the space occupied by the lists pointed at by
     the array prpt.                                           */

  pel = pedge -> prpt;
  for (ki=0; ki<pedge->iedge; ki++)
    {

      /* Traverse the list connected to edge nr ki and free the elements. */

      p1 = *pel;
      while (p1 != SISL_NULL)
	{
	  p2 = p1 -> pnext;
	  freePtedge(p1);
	  p1 = p2;
	}
      pel++;
    }

  /* Free the space occupied by the prpt array. */

  freearray(pedge -> prpt);

  /* Free the space occupied by the instance. */

  freearray(pedge);

  return;
}

#if defined(SISLNEEDPROTOTYPES)
void
     freeIntcrvlist(SISLIntcurve **viclist, int icrv)
#else
void freeIntcrvlist(viclist,icrv)
     SISLIntcurve **viclist;
     int      icrv;
#endif
/*
************************************************************************
*
* Purpose	: Free a list of Intcurve.
*
* Input		: viclist  - Array of pointers to pointers to instance
*                            of Intcurve.
*                 icrv     - # of SISLIntcurve in the list.
*
* Output	: None.
*
* Method	:
*
* References	:
*-
* Calls		: freeIntcurve - free instance of Intcurve.
*
*
* Written by	: B. O. Hoset, SI, Oslo, Norway, aug. 1989
*
************************************************************************
*/
{
  /*
   * Local declarations
   * ------------------
   */

  int         ki;

  /*
   * Free each SISLIntcurve from bottom of the list and u to the top.
   * ------------------------------------------------------------
   */

  if (viclist)
    {
      for ( ki = icrv - 1; ki >= 0; ki --)
	{
	  if (viclist[ki])
	    {
	      freeIntcurve(viclist[ki]);
	      viclist[ki] = SISL_NULL;
	    }
	}
      freearray(viclist);
      viclist = SISL_NULL;
    }
}

#if defined(SISLNEEDPROTOTYPES)
void
     freeIntcurve(SISLIntcurve *pintc)
#else
void freeIntcurve(pintc)
     SISLIntcurve *pintc;
#endif

/*
****************************************************************************
*
* Purpose    : Free the space occupied by the intersection curve structure.
*
* Input      : pintc   - Pointer to Intcurve.
*
*
*
* Output     :
*
*
* Method     :
*
*
* References :
*
*-
* Calls      : freeCurve - Free space occupied by a curve.
*
* Written by : Bjoern Olav Hoset, SI, 06-89.
*
****************************************************************************
*/
{
  /* Free whatever is allocated */

  if (pintc)
    {
      if (pintc->epar1) freearray(pintc->epar1);
      if (pintc->epar2) freearray(pintc->epar2);
      if (pintc->pgeom) freeCurve(pintc->pgeom);
      if (pintc->ppar1) freeCurve(pintc->ppar1);
      if (pintc->ppar2) freeCurve(pintc->ppar2);
      freearray(pintc);
    }
  return;
}

#if defined(SISLNEEDPROTOTYPES)
void
     freeIntdat(SISLIntdat *pintdat)
#else
void freeIntdat(pintdat)
     SISLIntdat *pintdat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Free the space occupied by an instance of the structure
*              Intdat. Also free space occupied by instances of Intlist
*              and SISLIntpt that belongs to intdat.
*
*
*
* INPUT      : pintdat - Pointer to the instance.
*
*
*
* OUTPUT     :
*
*
* METHOD     :
*
*
* REFERENCES :
*
*-
* CALLS      : freeIntpt
*              freeIntlist
*
* WRITTEN BY : UJK         , SI, 89-04.
*
*********************************************************************
*/
{
  int ki;                    /* Counter.                                     */

  if (pintdat == SISL_NULL) goto out;

  /* First free the space occupied by the Intpt's pointed at by
     the array vpoint.                                                       */

  for (ki=0; ki<pintdat->ipoint; ki++)
    if (pintdat -> vpoint[ki])   freeIntpt(pintdat -> vpoint[ki]);

  /* Free the space occupied by the vpoint array. */

  freearray(pintdat -> vpoint);

  /* Next free the space occupied by the Intlists pointed at by
     the array vlist.                                                       */

  for (ki=0; ki<pintdat->ilist; ki++)
    if    (pintdat -> vlist[ki])   freeIntlist(pintdat -> vlist[ki]);

  /* Free the space occupied by the vpoint array. */

  freearray(pintdat -> vlist);

  /* Free the space occupied by the instance. */

  freearray(pintdat);

 out:
  return;
}

#if defined(SISLNEEDPROTOTYPES)
void
     freeIntlist(SISLIntlist *plist)
#else
void freeIntlist(plist)
     SISLIntlist *plist;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Free the space occupied by an instance of the structure
*              Intlist. The instance contains a list of intersection
*              points defining an intersection curve.
*
*
*
* INPUT      : plist  - The space occupied by this list is to be free.
*
*
*
* OUTPUT     :
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
  /* Free space. */

  freearray(plist);

  return;
}

#if defined(SISLNEEDPROTOTYPES)
void
     freeIntpt(SISLIntpt *ppt)
#else
void freeIntpt(ppt)
     SISLIntpt *ppt;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Free the space allocated to keep the instance of Intpt
*              pointed at by ppt.
*
*
*
* INPUT      : ppt    - Pointer to instance.
*
*
*
* OUTPUT     :
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
  /* Free the arrays contained in the instance. */

  if (ppt->ipar)
    freearray(ppt -> epar);
  if (ppt->pnext)       freearray(ppt->pnext);
  if (ppt->curve_dir)   freearray(ppt->curve_dir);
  if (ppt->left_obj_1)  freearray(ppt->left_obj_1);
  if (ppt->left_obj_2)  freearray(ppt->left_obj_2);
  if (ppt->right_obj_1) freearray(ppt->right_obj_1);
  if (ppt->right_obj_2) freearray(ppt->right_obj_2);
  if (ppt->geo_data_1)  freearray(ppt->geo_data_1);
  if (ppt->geo_data_2)  freearray(ppt->geo_data_2);

  if(ppt->trim[0] != SISL_NULL) freeTrimpar(ppt->trim[0]);
  if(ppt->trim[1] != SISL_NULL) freeTrimpar(ppt->trim[1]);

  /* Free the instance pointed at by ppt. */

  freearray(ppt);
}

#if defined(SISLNEEDPROTOTYPES)
void
     freeIntsurf(SISLIntsurf *intsurf)
#else
void freeIntsurf(intsurf)
     SISLIntsurf *intsurf;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Free the space allocated to keep the instance of Intsurf
*              pointed at by intsurf.
*
*
*
* INPUT      : intsurf    - Pointer to instance.
*
*
*
* OUTPUT     :
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
* WRITTEN BY : Michael Floater, SI, 91-11.
*
*********************************************************************
*/
{

  /* Free the arrays if not SISL_NULL. */

  if(intsurf->epar != SISL_NULL) freearray(intsurf->epar);
  if(intsurf->const_par != SISL_NULL) freearray(intsurf->const_par);

  /* Free the instance pointed at by intsurf. */

  freearray(intsurf);

  return;
}

#if defined(SISLNEEDPROTOTYPES)
void
     freeTrimpar(SISLTrimpar *trimpar)
#else
void freeTrimpar(trimpar)
     SISLTrimpar *trimpar;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Free the space allocated to keep the instance of Trimpar
*              pointed at by trimpar.
*
*
*
* INPUT      : trimpar    - Pointer to instance.
*
*
*
* OUTPUT     :
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


  /* Free the instance pointed at by trimpar. */

  freearray(trimpar);

  return;
}

#if defined(SISLNEEDPROTOTYPES)
void
     freeTrack(SISLTrack *ppt)
#else
void freeTrack(ppt)
     SISLTrack *ppt;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Free the space allocated to keep the instance of Track
*              pointed at by ppt.
*              Curves and surfaces are also removed.
*
*
* INPUT      : ppt    - Pointer to instance.
*
*
*
* OUTPUT     :
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
  /* Free the arrays contained in the instance. */

   /*  if (ppt->psurf_1) freeSurf(ppt->psurf_1);
      if (ppt->ideg == 0 && ppt->psurf_2) freeSurf(ppt->psurf_2); */
  if (ppt->pcurve_3d) freeCurve(ppt->pcurve_3d);
  if (ppt->pcurve_2d_1) freeCurve(ppt->pcurve_2d_1);
  if (ppt->ideg==0 && ppt->pcurve_2d_2) freeCurve(ppt->pcurve_2d_2);

  /* Free the instance pointed at by ppt. */

  freearray(ppt);
}

#if defined(SISLNEEDPROTOTYPES)
void
     freeObject(SISLObject *pobj)
#else
void freeObject(pobj)
     SISLObject *pobj;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Free the space occupied by the object pointed at by
*              pobj. Also free the point, curve or surface that the
*              object represents.
*
*
*
* INPUT      : pobj   - Pointer to object.
*
*
*
* OUTPUT     :
*
*
* METHOD     :
*
*
* REFERENCES :
*
*-
* CALLS      : freePoint - Free space occupied by a point.
*              freeCurve - Free space occupied by a curve.
*              freeSurf  - Free space occupied by a surface.
*
* WRITTEN BY : Vibeke Skytt, SI, 88-05.
*
* MODIFIED   : UJK               89-04.
*              The structure SISLPoint is included.
*
*********************************************************************
*/
{
  int ki;

  /* Free point, curve or surface represented by pobj.                */

  if (pobj -> iobj == SISLPOINT)
    { if (pobj -> p1 != SISL_NULL) freePoint(pobj -> p1); }
  else if (pobj -> iobj == SISLCURVE)
    { if (pobj -> c1 != SISL_NULL) freeCurve(pobj -> c1); }
  else if (pobj -> iobj == SISLSURFACE)
    { if (pobj -> s1 != SISL_NULL) freeSurf(pobj -> s1);  }

  for (ki=0; ki<4; ki++)
    if (pobj->edg[ki] != SISL_NULL) freeObject(pobj->edg[ki]);

  /* Free instance of object. */

  freearray(pobj);

  return;
}

#if defined(SISLNEEDPROTOTYPES)
void
     freePoint(SISLPoint *ppoint)
#else
void freePoint(ppoint)
     SISLPoint *ppoint;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Free the space occupied by the point pointed at by
*              ppoint.
*
*
*
* INPUT      : ppoint - Pointer to point.
*
*
*
* OUTPUT     :
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
* WRITTEN BY : UJK, SI, 89-04.   Vibeke Skytt, SI, 91-01.
*
*********************************************************************
*/
{
   int ki;   /* Counter.  */

  if (ppoint != SISL_NULL)
    {
      if (ppoint -> pbox != SISL_NULL)
      {
	 if (ppoint->pbox->emax) freearray(ppoint->pbox->emax);
	 if (ppoint->pbox->emin) freearray(ppoint->pbox->emin);

	 for (ki=0; ki<3; ki++)
	 {
	    if (ppoint->pbox->e2max[ki]) freearray(ppoint->pbox->e2max[ki]);
	    if (ppoint->pbox->e2min[ki]) freearray(ppoint->pbox->e2min[ki]);
	 }
	 freearray(ppoint->pbox);
	}

      if ((ppoint -> idim > 3) &&
	  (ppoint -> icopy != 0) &&
	  (ppoint -> ecoef != SISL_NULL))
	freearray(ppoint -> ecoef);	/* Free array.  */

      freearray(ppoint);		/* Free instance of point. */
    }
  return;
}

#if defined(SISLNEEDPROTOTYPES)
void
     freePtedge(SISLPtedge *p1)
#else
void freePtedge(p1)
     SISLPtedge *p1;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Free the space allocated by an instance of the structure
*              Ptedge.
*
*
*
* INPUT      : p1     - A pointer to the instance.
*
*
*
* OUTPUT     :
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
  /* Free the space that p1 occupies. */

  freearray(p1);

  return;
}

#if defined(SISLNEEDPROTOTYPES)
void
     freeSurf(SISLSurf *psurf)
#else
void freeSurf(psurf)
     SISLSurf *psurf;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Free the space occupied by the surface pointed at by
*              psurf.
*
*
*
* INPUT      : psurf - Pointer to surface.
*
*
*
* OUTPUT     :
*
*
* METHOD     :
*
*
* REFERENCES :
*
*-
* CALLS      : freearray - Free space occupied by a given array.
*
* WRITTEN BY : Vibeke Skytt, SI, 88-05. Arne Laksaa, SI, 89-03.
*              Vibeke Skytt, SI, 91-01.
* Revised by : Paal Fugelli, SINTEF, Oslo, Norway, 94-09.  Added free'ing
*              of ecoef for rationals when icopy==0.  ecoef is ALWAYS
*              allocated for rationals, regardless of the icopy flag value.
*
*********************************************************************
*/
{
  int ki;     /* Counter.  */


  if ( psurf->icopy != 0 )
  {

    /* Free arrays.  */

    freearray(psurf->et1);
    freearray(psurf->et2);
    freearray(psurf->ecoef);
    if ( psurf->rcoef != SISL_NULL )  freearray(psurf->rcoef);
  }
  else if ( psurf->ikind == 2  ||  psurf->ikind == 4 )
  {
    /* VSK, 940902. The array pcurve->ecoef is ALWAYS allocated in the
       constructor for rationals and must be free'ed.                  */

    freearray(psurf->ecoef);
  }


  if ( psurf->pdir )
  {

    /* Free direction structure. */

    if ( psurf->pdir->ecoef )  freearray(psurf->pdir->ecoef);
    if ( psurf->pdir->esmooth != SISL_NULL )  freearray(psurf->pdir->esmooth);
    freearray(psurf->pdir);
  }

  if ( psurf->pbox )
  {

    /* Free surrounded SISLbox structure. */

    if ( psurf->pbox->emax )  freearray(psurf->pbox->emax);
    if ( psurf->pbox->emin )  freearray(psurf->pbox->emin);

    for ( ki=0;  ki < 3;  ki++ )
    {
      if ( psurf->pbox->e2max[ki] )  freearray(psurf->pbox->e2max[ki]);
      if ( psurf->pbox->e2min[ki] )  freearray(psurf->pbox->e2min[ki]);
    }
    freearray(psurf->pbox);
  }

  /* Free instance of surface. */

  freearray(psurf);

  return;
}
