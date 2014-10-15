//
//  GeomTools.h
//  
//
//  Created by Xujun ZHao on 10/14/14.
//
//

#ifndef _GeomTools_h
#define _GeomTools_h

// C++ Includes
#include <iostream>
#include <algorithm>
#include <vector>
#include <math.h>


// LibMesh library includes
#include "libmesh/point.h"

// Bring in everything from the libMesh namespace
using namespace libMesh;



/*
 this class defines basic tools used in our codes
 */
class GeomTools
{
public:
  
  // return the norm of a point |x| = sqrt(x^2 + y^2 + z^2)
  static Real point_norm(const Point& pt);
  
  // quadratic function used for applying BC to avoid singularities at corners
  static Real quadratic_function(const Real& yb);
};



// =============================================================================================
Real GeomTools::point_norm(const Point& pt)
{
  Real val = 0.;
  for(unsigned int i=0; i<LIBMESH_DIM; ++i)
    val += pt(i)*pt(i);
  
  return std::sqrt(val);
}


// =============================================================================================
Real GeomTools::quadratic_function(const Real& yb)
{
  // a controls the magnitude of velocity, b is the mag of geometry(height of channal)
  Real a = 1.0, b = 0.5, y0 = 0.0;
  return b*b - a*(yb - y0)*(yb - y0);
}




#endif
