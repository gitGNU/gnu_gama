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
 *  $Id: g3_cluster_vec.h,v 1.4 2003/06/08 08:11:13 cepek Exp $
 */

#include <gnu_gama/g3/g3_cluster.h>

#ifndef GNU_gama__g3_cluster_vector_gnugamag3clustervectorhgnugamag3clstrvec
#define GNU_gama__g3_cluster_vector_gnugamag3clustervectorhgnugamag3clstrvec


namespace GNU_gama {  namespace g3 {

  class Vectors :  public g3Cluster {
  public:

    GNU_gama::List<Vector*> vectors;

    Vectors(const Model::ObservationData* obs) : g3Cluster(obs) 
      {
      }
    ~Vectors();

    void add(Vector*);
    void write_xml(std::ostream&) const;
    void parlist_init(Model*);
  };

}}


#endif
