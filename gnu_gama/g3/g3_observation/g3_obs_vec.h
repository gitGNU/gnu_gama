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
 *  $Id: g3_obs_vec.h,v 1.10 2003/12/27 21:00:58 uid66336 Exp $
 */

#include <gnu_gama/g3/g3_observation/g3_obs_base.h>

#ifndef GNU_gama__g3_obs_vector_h_gnugamag3obs_vectorh___gnu_gama_g3obs__vec
#define GNU_gama__g3_obs_vector_h_gnugamag3obs_vectorh___gnu_gama_g3obs__vec



namespace GNU_gama {  namespace g3 {


  class Vector : public Observation {
  public:  
    
    Point::Name from;
    Point::Name to;
    
    Vector()
    {
    }
    Vector(double x, double y, double z)
    {
      dx_ = x; dy_ = y; dz_ = z;
    }
    Vector(const Vector& v)
    {
      from = v.from;
      to   = v.to;
      
      dx_ = v.dx_; dy_ = v.dy_; dz_ = v.dz_;
    }
    
    int  dimension() const { return 3; }
    void set_dxyz(double x, double y, double z)
    {
      dx_ = x; dy_ = y; dz_ = z;
    }
    
    double dx() const { return dx_; }
    double dy() const { return dy_; }
    double dz() const { return dz_; }
    
    bool revision_accept(ObservationVisitor* visitor)
    {
      if (Revision<Vector>* rv = dynamic_cast<Revision<Vector>*>(visitor))
        {
          return rv->revision_visit(this);
        }
      else
        return false;
    }
  
    void linearization_accept(ObservationVisitor* visitor)
    {
      if (Linearization<Vector>* 
          lv = dynamic_cast<Linearization<Vector>*>(visitor))
        {
          lv->linearization_visit(this);
        }
    }
  
    
  private:
    
    double dx_, dy_, dz_;
  };
  
}}


#endif
