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
 *  $Id: g3_model.h,v 1.12 2003/11/23 14:40:27 cepek Exp $
 */

#include <gnu_gama/model.h>
#include <gnu_gama/pointbase.h>
#include <gnu_gama/obsdata.h>
#include <gnu_gama/list.h>
#include <gnu_gama/ellipsoids.h>
#include <gnu_gama/g3/g3_point.h>
#include <gnu_gama/g3/g3_observation.h>


#ifndef GNU_gama__g3_model_h_gnugamag3modelh___gnu_gama_g3model
#define GNU_gama__g3_model_h_gnugamag3modelh___gnu_gama_g3model


namespace GNU_gama {  namespace g3 {

  
  class Model 
    : public Revision<Distance>
    {
    public:
    
      typedef GNU_gama::PointBase<g3::Point>              PointBase;
      typedef GNU_gama::ObservationData<g3::Observation>  ObservationData;
      
      PointBase           *points;
      ObservationData     *obs;
      
      GNU_gama::Ellipsoid  ellipsoid;
      
      
      Model();
      ~Model();
      

      Point* get_point(const Point::Name&);
      void   write_xml(std::ostream& out) const;
      void   pre_linearization();
      
      void reset()               { state_ = init_; }
      void reset_points()        { if (points_ < state_) state_ = points_; }
      void reset_observations()  { if (obsrvs_ < state_) state_ = obsrvs_; }
      void reset_linearization() { if (linear_ < state_) state_ = linear_; }
      void reset_adjustment()    { if (adjust_ < state_) state_ = adjust_; }
      
      bool check_points()        const { return state_ > points_; }
      bool check_observations()  const { return state_ > obsrvs_; }
      bool check_linearization() const { return state_ > linear_; }
      bool check_adjustment()    const { return state_ > adjust_; }
      
      void update_points();
      void update_observations();
      void update_linearization();
      void update_adjustment();
      
      
      // virtual functions derived from template class Revision<>
      bool revision(Distance*);


    private:   /*-----------------------------------------------------------*/
      
      Model(const Model&);
      Model& operator=(const Model&);
      
      // active observations' list (observations used in the adjustment)
      GNU_gama::List<Observation*>  active_obs;
      
      // basic revision steps 
      enum State_{init_, points_, obsrvs_, linear_, adjust_, ready_} state_;
      
      void next_state_(int s) { state_ = State_(++s); }
      bool check_init() const { return state_ > init_; }
      void update_init();
      
    };
  
}}

#endif
