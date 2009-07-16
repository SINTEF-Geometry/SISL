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
 * $Id: s1440.c,v 1.3 2001-03-19 15:58:49 afr Exp $
 *
 */


#define S1440

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
s1440(SISLSurf *ps1,SISLSurf **rs2,int *jstat)
#else
void s1440(ps1,rs2,jstat)
     SISLSurf *ps1;
     SISLSurf **rs2;
     int  *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Interchange the two parameter directions used in the
*              mathematical description of a surface and thereby
*              change the direction of the normal vector of the surface.
*
*
*
* INPUT      : ps1    - Pointer to the original surface.
*
*
*
* OUTPUT     : rs2    - Pointer to the surface with interchanged
*                       parameter directions.
*              jstat  - status messages
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
* CALLS      : s6chpar  - Change parameter directions of the vertices
*                         of the surface.
*              newSurf  - Create new surface.
*
* WRITTEN BY : Vibeke Skytt, SI, 88-11.
* REVISED BY : Johannes Kaasa, SI, 91-09 (Introduced NURBS).
* Revised by : Paal Fugelli, SINTEF, Oslo, Norway, Nov. 1994.  Added
*              handling of 'cuopen' flags.
*
*********************************************************************
*/
{
  int kpos = 0;          /* Position of error.                  */
  double *ssurf = SISL_NULL;  /* Pointer to vertices of new surface. */
  int kdim;              /* Local (rational) dimension.         */
  double *vert;          /* Pointer to vertices.                */

  /* Check for rational surface. */

  if (ps1->ikind == 2 || ps1->ikind == 4)
    {
      kdim = ps1->idim + 1;
      vert = ps1->rcoef;
    }
  else
    {
      kdim = ps1->idim;
      vert = ps1->ecoef;
    }

  /* Allocate scratch for vertices of new surface.  */

  ssurf = newarray(ps1->in1*ps1->in2*kdim,double);
  if (ssurf == SISL_NULL) goto err101;

  /* Change parameter directions of vertices.  */

  s6chpar(vert,ps1->in1,ps1->in2,kdim,ssurf);

  /* Create output surface.  */

  *rs2 = SISL_NULL;
  if ((*rs2 = newSurf(ps1->in2,ps1->in1,ps1->ik2,ps1->ik1,ps1->et2,
		      ps1->et1,ssurf,ps1->ikind,ps1->idim,1)) == SISL_NULL) goto err101;

  /* Set periodicity flag */

  (*rs2)->cuopen_1 = ps1->cuopen_2;
  (*rs2)->cuopen_2 = ps1->cuopen_1;

  /* Parameter directions changed.  */

  *jstat = 0;
  goto out;

  /* Error in space allocation.  */

 err101: *jstat = -101;
  s6err("s1440",*jstat,kpos);
  goto out;

 out:

  /* Free space occupied by local array.  */

  if (ssurf != SISL_NULL) freearray(ssurf);

  return;
}
