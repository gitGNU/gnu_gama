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
 *  $Id: g3_cluster.h,v 1.5 2003/12/23 19:52:49 uid66336 Exp $
 */


// #include <gnu_gama/g3/g3_observation/g3_cluster_vec.h>

#ifndef GNU_gama__g3_cluster___g3clustergnugamag3cluster__gnu_gama_g3_cluster
#define GNU_gama__g3_cluster___g3clustergnugamag3cluster__gnu_gama_g3_cluster

#include <gnu_gama/model.h>
#include <gnu_gama/g3/g3_observation/g3_obs_vec.h>
#include <gnu_gama/g3/g3_model.h>
#include <iostream>

namespace GNU_gama { namespace g3 {


  class g3Cluster :  public GNU_gama::Cluster<Observation> {
  public:
  
      g3Cluster(const Model::ObservationData* obs) : Cluster<Observation>(obs) 
      {
      }

      g3Cluster* clone(const Model::ObservationData*) const 
        { 
          throw 
            GNU_gama::Exception::string("g3Cluster::clone() not implemented");
          return 0; 
        }

      virtual void write_xml(std::ostream&) const = 0;
      virtual void parlist_init(Model*) {}
 };

}}

#endif
