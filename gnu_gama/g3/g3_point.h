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
 *  $Id: g3_point.h,v 1.9 2003/03/25 12:38:33 cepek Exp $
 */

#include <gnu_gama/g3/g3_parameter.h>
#include <gamalib/pointid.h>


#ifndef GNU_gama__g3_point_h_gnugamag3pointh___gnu_gama_g3point
#define GNU_gama__g3_point_h_gnugamag3pointh___gnu_gama_g3point


namespace GNU_gama {  namespace g3 {

  class Model;


  class Parameter_position : public Parameter {
  public:

    Parameter_position* clone() { return new Parameter_position(*this); }
  };


  class Parameter_height : public Parameter {
  public:

    Parameter_height* clone() { return new Parameter_height(*this); }
  };



  class Point {
  public:
  
    typedef GaMaLib::PointID Name;
    typedef Model            Common;

    Name                 name;
    Common*              common; 

    Point();
    Point(const Point&);
    Point& operator=(const Point&);
    ~Point();

    // -----------------------------

    Parameter_position*  B;
    Parameter_position*  L;
    Parameter_height*    H;

    ParameterList        parlist;


    enum {
      unused               = 0,
      fixed_position       = 1,
      fixed_height         = 2,
      fixed_point          = fixed_position + fixed_height,
      free_position        = 4,
      free_height          = 8,
      free_point           = free_position  + free_height,
      constrained_position = 16 + free_position,
      constrained_height   = 32 + free_height,
      constrained_point    = constrained_position + constrained_height,
      position             = fixed_position + free_position, 
      height               = fixed_height   + free_height,
      point                = position + height  
    };

    int  state(int s) const { return s & state_; }
    void set_state(int s)   { state_  = s;       }

  private:

    int state_;

  };


}}

#endif
