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
 *  $Id: g3_model.h,v 1.20 2004/02/19 17:21:23 cepek Exp $
 */

#include <gnu_gama/model.h>
#include <gnu_gama/pointbase.h>
#include <gnu_gama/obsdata.h>
#include <gnu_gama/list.h>
#include <gnu_gama/ellipsoids.h>
#include <gnu_gama/g3/g3_point.h>
#include <gnu_gama/g3/g3_observation.h>
#include <gnu_gama/sparse/smatrix.h>
#include <gnu_gama/sparse/sbdiagonal.h>
#include <gnu_gama/adj/adj.h>

#ifndef GNU_gama__g3_model_h_gnugamag3modelh___gnu_gama_g3model
#define GNU_gama__g3_model_h_gnugamag3modelh___gnu_gama_g3model


namespace GNU_gama {  namespace g3 {

  
  class Model : 
    public GNU_gama::Model<g3::Observation>,
    public Revision     <Distance>,
    public Linearization<Distance>,
    public Revision     <Vector>,
    public Linearization<Vector>
  {
  public:
    
    typedef GNU_gama::ObservationData<g3::Observation>  ObservationData;
    typedef GNU_gama::List<Observation*>                ObservationList; 
    typedef GNU_gama::PointBase<g3::Point>              PointBase;
    typedef GNU_gama::List<Parameter*>                  ParameterList;
    
    PointBase           *points;    
    GNU_gama::Ellipsoid  ellipsoid;
    
    
    Model();
    virtual ~Model();
    
    
    Point* get_point(const Point::Name&);
    void   write_xml(std::ostream& out) const;
    
    void reset()               { state_ = init_; }
    void reset_parameters()    { if (params_ < state_) state_ = params_; }
    void reset_observations()  { if (obsrvs_ < state_) state_ = obsrvs_; }
    void reset_linearization() { if (linear_ < state_) state_ = linear_; }
    void reset_adjustment()    { if (adjust_ < state_) state_ = adjust_; }
    
    bool check_parameters()    const { return state_ > params_; }
    bool check_observations()  const { return state_ > obsrvs_; }
    bool check_linearization() const { return state_ > linear_; }
    bool check_adjustment()    const { return state_ > adjust_; }
    
    void update_parameters();
    void update_observations();
    void update_linearization();
    void update_adjustment();
        
    bool revision_visit     (Distance*);
    void linearization_visit(Distance*);
    bool revision_visit     (Vector*  );
    void linearization_visit(Vector*  );

    void write_xml_adjustment_input_data(std::ostream&);

  private:   /*-----------------------------------------------------------*/
      
    Model(const Model&);
    Model& operator=(const Model&);
    
    // active observations list (active observations used in the adjustment)
    ObservationList*  active_obs;

    // parameter list
    ParameterList*  par_list;
    void update_index(Parameter&);
    
    // basic revision steps 
    enum State_{init_, params_, obsrvs_, linear_, adjust_, ready_} state_;
    
    void next_state_(int s) { state_ = State_(++s); }
    bool check_init() const { return state_ > init_; }
    void update_init();
        

    // design matrix
    int dm_rows, dm_cols, dm_floats;
    SparseMatrix <>*  A;
    Vec               rhs;
    int               rhs_ind;
    BlockDiagonal<>*  B;
    GNU_gama::AdjInputData  adj_input_data;

  };
  
}}

#endif
