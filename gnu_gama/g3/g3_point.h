/*  
    GNU Gama -- adjustment of geodetic networks
    Copyright (C) 2003  Ales Cepek <cepek@fsv.cvut.cz>

    This file is part of the GNU Gama library.
    
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
 *  $Id: g3_point.h,v 1.3 2003/03/19 10:58:57 cepek Exp $
 */

#include <gnu_gama/g3/g3_parameter.h>
#include <gamalib/pointid.h>


#ifndef GNU_gama__g3_point_h_gnugamag3pointh___gnu_gama_g3point
#define GNU_gama__g3_point_h_gnugamag3pointh___gnu_gama_g3point


namespace GNU_gama {  namespace g3 {

  class Model;


  class Parameter_position : public Parameter {
  public:
  };


  class Parameter_height   : public Parameter {
  public:
  };



  class Point {
  public:
  
    typedef GaMaLib::PointID Name;
    typedef Model            Common;

    Name                 name;
    Common*              common; 

    Parameter_position*  B;
    Parameter_position*  L;
    Parameter_height*    H;

    ParameterList        parlist;


    Point();
    ~Point();



  private:

    Point(const Point&);
    Point& operator=(const Point&);

  };


}}

#endif
