/* SISL - SINTEF Spline Library version 4.4.                              */
/* Definition and interrogation of NURBS curves and surface.              */
/* Copyright (C) 1978-2005, SINTEF ICT, Applied Mathematics, Norway.      */

/* This program is free software; you can redistribute it and/or          */
/* modify it under the terms of the GNU General Public License            */
/* as published by the Free Software Foundation version 2 of the License. */

/* This program is distributed in the hope that it will be useful,        */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of         */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          */
/* GNU General Public License for more details.                           */

/* You should have received a copy of the GNU General Public License      */
/* along with this program; if not, write to the Free Software            */
/* Foundation, Inc.,                                                      */
/* 59 Temple Place - Suite 330,                                           */
/* Boston, MA  02111-1307, USA.                                           */

/* Contact information: e-mail: tor.dokken@sintef.no                      */
/* SINTEF ICT, Department of Applied Mathematics,                         */
/* P.O. Box 124 Blindern,                                                 */
/* 0314 Oslo, Norway.                                                     */

/* SISL commercial licenses are available for:                            */
/* - Building commercial software.                                        */
/* - Building software whose source code you wish to keep private.        */

#include <stdio.h>
#include <math.h>

// #ifdef _MSC_VER
// #  include <vector>
// #else
// #  include <vector.h>
// #endif
#include <vector>

#ifdef _MSC_VER
#  define max(a, b) ((a)>(b)? (a) : (b))
#endif

#if defined(_MSC_VER) || defined(LINUX)
#  include <sys/timeb.h>
#endif

#ifndef _MSC_VER
#  include <strings.h>
#  include <unistd.h>
#else
#  include <stdlib.h>
#endif

#include <time.h>




double jon_timer;


static double jon_sec(void)
{
#ifndef _MSC_VER
#ifdef LINUX
  struct timeb jon_ts;
#else
  struct timespec jon_ts;
#endif
#else
  struct _timeb jon_ts;
#endif
  
#ifndef _MSC_VER
#ifdef LINUX
  ftime(&jon_ts);
  return (jon_ts.time + 0.001*jon_ts.millitm);
#else
  clock_gettime(CLOCK_REALTIME, &jon_ts);
  return (int)(jon_ts.tv_sec) + ((double)(jon_ts.tv_nsec)/1e9);
#endif
#else
  _ftime(&jon_ts);
  return (jon_ts.time + 0.001*jon_ts.millitm);
#endif
}


void tic(void)
{
  jon_timer=jon_sec();
}


void toc(void)
{
  printf("Time elapsed: %f\n", jon_sec()-jon_timer);
}
