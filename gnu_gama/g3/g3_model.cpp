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
 *  $Id: g3_model.cpp,v 1.6 2003/05/17 17:07:08 cepek Exp $
 */

#include <gnu_gama/g3/g3_model.h>


using namespace GNU_gama::g3;


Model::Model() 
{ 
  using namespace GNU_gama;

  points = new PointBase;
  obs    = new ObservationData;

  points->set_common_data(this); 
  set(&ellipsoid, ellipsoid_wgs84);
}


Model::~Model()
{
  delete points;
  delete obs;

}


Point* Model::get_point(const Point::Name& name)
{
  Point* p = points->find(name);
  if (p == 0)
    {
      p = new Point;
      p->name = name;
      points->put(p);
    }

  return p;
}
