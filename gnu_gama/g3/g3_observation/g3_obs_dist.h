/*  
    GNU Gama -- adjustment of geodetic networks
    Copyright (C) 2003  Ales Cepek <cepek@fsv.cvut.cz>

    This file is part of the GNU Gama C++ library.
    
    This library is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this library; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

/*
 *  $Id: g3_obs_dist.h,v 1.1 2003/03/28 22:07:31 cepek Exp $
 */

#include <gnu_gama/g3/g3_observation/g3_obs_base.h>


#ifndef GNU_gama__g3_obs_distance_h_gnugamag3obs_distanceh___gnu_gama_g3obs
#define GNU_gama__g3_obs_distance_h_gnugamag3obs_distanceh___gnu_gama_g3obs


namespace GNU_gama {  namespace g3 {


  class Distance : public Observation {
  public:  

    Point::Name name[2];

    Distance() : Observation(6) {}

    double parlist_value() const;
    void   parlist_init (Model*);
    double derivative   (Parameter*);
  };


  class DistanceAnalyticalDerivative {
  public:
    
    virtual double analytical_derivative(Distance*) = 0;
  };

}}


#endif
