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
 *  $Id: g3_obs_base.h,v 1.8 2003/04/11 09:38:26 cepek Exp $
 */

#include <gnu_gama/g3/g3_parameter.h>
#include <gnu_gama/sparse/svector.h>
#include <gnu_gama/ellipsoid.h>

#ifndef GNU_gama__g3_obs_base_h_gnugamag3obs_baseh___gnu_gama_g3obs
#define GNU_gama__g3_obs_base_h_gnugamag3obs_baseh___gnu_gama_g3obs

#include <gnu_gama/g3/g3_point.h>


namespace GNU_gama {  namespace g3 {


  class Observation {
  public:

    Observation(int n) : parlist(n), model(0), time(0) {}
    virtual ~Observation() {}

    double  obs() const { return 0; }

    virtual double parlist_value()  const    = 0;      
    virtual void   parlist_init (Model*)     = 0;
    virtual double derivative   (Parameter*) = 0;

    bool   active() const     { return active_; }
    void   set_active(bool b) { active_ = b;    }

    double numerical_derivative (Parameter*);
    void   linearization        (GNU_gama::SparseVector<>&);

    ParameterList  parlist;

  protected:  

    GNU_gama::g3::Model* model;
    double  time;

    virtual void prepare_to_linearization() {}

  private:

    bool active_;
  };

}}


#endif
