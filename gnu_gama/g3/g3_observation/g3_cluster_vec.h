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
 *  $Id: g3_cluster_vec.h,v 1.1 2003/05/21 17:14:02 cepek Exp $
 */

#include <gnu_gama/g3/g3_observation/g3_obs_vec.h>
#include <gnu_gama/g3/g3_model.h>

#ifndef GNU_gama__g3_cluster_vector_gnugamag3clustervectorhgnugamag3clstrvec
#define GNU_gama__g3_cluster_vector_gnugamag3clustervectorhgnugamag3clstrvec


namespace GNU_gama {  namespace g3 {

  class Vectors : public GNU_gama::Cluster<Observation> {
  public:

    GNU_gama::List<Vector*> vectors;

    Vectors(const Model::ObservationData* obs) : Cluster<Observation>(obs) 
      {
      }
    Vectors* clone(const Model::ObservationData*) const 
      { 
        throw GNU_gama::Exception::string("Vector::clone() not implemented");
        return 0; 
      }
    void add(Vector*);
  };

}}


#endif
