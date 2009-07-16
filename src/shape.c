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
 * $Id: shape.c,v 1.1 1994-04-21 12:10:42 boh Exp $
 *
 */


#define SHAPE

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
      shape(double emid[],double etang[],int idim,int iedge,int *jstat)
#else	 
void shape(emid,etang,idim,iedge,jstat)
     int idim,iedge,*jstat;
     double emid[],etang[];
#endif     
/*
*********************************************************************
*                                                                   
* PURPOSE    : This routine gives a possibility for the application
*              to adjust the value and the derivatives in the midpoint
*              of the vertex region, i.e. the point in which the region
*              is divided.
*
*
* INPUT      : idim    - Dimension of geometry space.
*              iedge   - Number of edges of the vertex region.
*
*
* INPUT/OUTPUT  : emid    - The value in the midpoint of the region.
*                           Dimension is idim.
*                 etang   - The tangents of the blending surfaces in
*                           the midpoint of the region, along the 
*                           curves which divides the region into 4-sided
*                           blending surfaces. Dimension is iedge*idim.
*                       
*
* OUTPUT     : jstat   - status messages  
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*
*
* USE        : The input of the arrays emid and etang may be changed.
*              This will effect the geometry of the blend strongly. Make
*              sure to return sensible midpoint and tangents. The midpoint
*              should lie close to the real midpoint of the region, pulling
*              it close to an edge, may result in blending surfaces with
*              cusps. Also if the tangents is very long, cusps may occur.
*              The numbering of the tangents must be the same as the numbering
*              of the edges, otherwise the G1-continuity will be lost.
*
*-
* CALLS      : 
*
* WRITTEN BY : Vibeke Skytt, SI, 05.90.
*
*********************************************************************
*/
{
  int ki;
  double tfac = 1.0;
  
  /*ü printf("Give factor with which to multiply tangent : ");*/
  /*ü scanf("%lf",&tfac);*/
  
  for (ki=0; ki<iedge*idim; ki++) etang[ki] *= tfac;

  *jstat = 0;

  return;
}
