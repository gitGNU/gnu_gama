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
 *  $Id: g3_obs_vec.h,v 1.2 2003/05/06 18:16:34 cepek Exp $
 */

#include <gnu_gama/g3/g3_observation/g3_obs_base.h>


#ifndef GNU_gama__g3_obs_vector_h_gnugamag3obs_vectorh___gnu_gama_g3obs__vec
#define GNU_gama__g3_obs_vector_h_gnugamag3obs_vectorh___gnu_gama_g3obs__vec




namespace GNU_gama {  namespace g3 {


  class Vector : private Observation {
  public:  
    
    Point::Name name[2];
    
    Vector() : Observation(6), select(0) {}

  private:

    int select;

    double parlist_value() const;    
    void   parlist_init (Model*);
    double derivative   (Parameter*);
    void   prepare_to_linearization();

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

    void   parlist_init (Model*);
    double parlist_value() const 
    { 
      vec->select = sel;   return vec->parlist_value();
    }
    double derivative(Parameter* p)
    {
      vec->select = sel;   return vec->derivative(p);
    }
    void linearization(GNU_gama::SparseVector<>& row)
    {
      vec->select = sel;   vec->linearization(row);
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
