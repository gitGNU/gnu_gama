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
 *  $Id: g3_model.h,v 1.5 2003/04/08 16:41:51 cepek Exp $
 */

#include <gnu_gama/pointbase.h>
#include <gnu_gama/obsdata.h>
#include <gnu_gama/list.h>
#include <gnu_gama/ellipsoids.h>
#include <gnu_gama/g3/g3_point.h>
#include <gnu_gama/g3/g3_observation.h>


#ifndef GNU_gama__g3_model_h_gnugamag3modelh___gnu_gama_g3model
#define GNU_gama__g3_model_h_gnugamag3modelh___gnu_gama_g3model


namespace GNU_gama {  namespace g3 {

  
  class Model {
  public:
    
    typedef GNU_gama::PointBase<g3::Point>              PointBase;
    typedef GNU_gama::ObservationData<g3::Observation>  ObservationData;
    
    PointBase           points;
    ObservationData     obs;
    
    GNU_gama::Ellipsoid ellipsoid;


    Model();
    ~Model();

    Point* get_point(const Point::Name&);
  };


}}

#endif
