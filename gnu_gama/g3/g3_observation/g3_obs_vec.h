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
 *  $Id: g3_obs_vec.h,v 1.8 2003/12/24 11:34:11 uid66336 Exp $
 */

#include <gnu_gama/g3/g3_observation/g3_obs_base.h>

#ifndef GNU_gama__g3_obs_vector_h_gnugamag3obs_vectorh___gnu_gama_g3obs__vec
#define GNU_gama__g3_obs_vector_h_gnugamag3obs_vectorh___gnu_gama_g3obs__vec



namespace GNU_gama {  namespace g3 {


  class Vector : private Observation {
  public:  
    
    Point::Name name[2];
    
    Vector() : Observation(6), select(0) 
    {
    }
    Vector(double x, double y, double z) : Observation(6), select(0)
    {
      dxyz_[0] = x; dxyz_[1] = y; dxyz_[2] = z;
    }
    Vector(const Vector& v) : Observation(6), select(0)
    {
      name[0] = v.name[0];
      name[1] = v.name[1];

      dxyz_[0] = v.dxyz_[0]; dxyz_[1] = v.dxyz_[1]; dxyz_[2] = v.dxyz_[2];
    }

    int  dimension() const { return 3; }
    void set_dxyz(double x, double y, double z)
    {
      dxyz_[0] = x; dxyz_[1] = y; dxyz_[2] = z;
    }

    double dx() const { return dxyz_[0]; }
    double dy() const { return dxyz_[1]; }
    double dz() const { return dxyz_[2]; }

    void   parlist_init (Model*);

  private:

    int select;

    double dxyz_[3];

    bool revision_accept(ObservationVisitor* visitor)
    {
      if (Revision<Vector>* rv = dynamic_cast<Revision<Vector>*>(visitor))
        {
          return rv->revision_visit(this);
        }
      else
        return false;
    }

    friend class Diff;
  };



  class Diff : public Observation {
  public: 

    Point::Name name[2];
    
    Diff(Vector* v, int r) : Observation(6), vec(v), sel(r)
    {
      name[0] = vec->name[0];
      name[1] = vec->name[1];
    }

    double obs() const 
    { 
      return vec->dxyz_[sel]; 
    }
    bool   revision_accept(ObservationVisitor* visitor) 
    { 
      return vec->revision_accept(visitor); 
    }

    Vector* vector() const { return vec; }

  private:
    
    Vector*   vec;
    const int sel;

  };



  class DiffX : public Diff { public: DiffX(Vector* v) : Diff(v, 0) {} };
  class DiffY : public Diff { public: DiffY(Vector* v) : Diff(v, 1) {} };
  class DiffZ : public Diff { public: DiffZ(Vector* v) : Diff(v, 2) {} };

}}


#endif
