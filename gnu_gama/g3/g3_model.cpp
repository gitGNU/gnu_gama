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
 *  $Id: g3_model.cpp,v 1.23 2004/01/25 11:07:13 cepek Exp $
 */

#include <gnu_gama/g3/g3_model.h>
#include <gnu_gama/g3/g3_cluster.h>
#include <gnu_gama/outstream.h>


using namespace std;
using namespace GNU_gama::g3;


Model::Model() 
{ 
  using namespace GNU_gama;

  points     = new PointBase;
  active_obs = new ObservationList;
  par_list    = new ParameterList;

  A = 0;
  B = 0;

  points->set_common_data(this); 
  set(&ellipsoid, ellipsoid_wgs84);

  reset();
}


Model::~Model()
{
  delete points;
  delete active_obs;
  delete par_list;
  delete A;
  delete B;
}


Point* Model::get_point(const Point::Name& name)
{
  Point* p = points->find(name);
  if (p == 0)
    {
      p = new Point;
      p->name = name;
      points->put(p);
    }

  return p;
}


void Model::write_xml(std::ostream& out) const
{
  using namespace std;
 
  GNU_gama::SaveFlags sf(out);
  out.setf(ios_base::fixed, ios_base::floatfield);
  out.precision(5);

  out << "<g3-model>\n";

  {
    out << "\n";
    for (Model::PointBase::const_iterator
           b = points->begin(), e = points->end(); b != e; ++b)
      {
        const Point *p = *b;
        out << "<point>\t<id>" << p->name.c_str() << "</id>";
        
        if (p->has_xyz())
          out << "\n\t"
              << "<x>" << p->X.value() << "</x> "
              << "<y>" << p->Y.value() << "</y> "
              << "<z>" << p->Z.value() << "</z>";
        
        if (p->has_height())
          out << "\n\t"
              << "<height>" << p->height.value() << "</height>";
        
        if (p->unused())
          out << "\n\t"; // <unused/>";
        else if (p->fixed_position())   
          out << "\n\t<fixed/>";
        else if (p->free_position() && !p->constr_position())    
          out << "\n\t<free/>";
        else if (p->constr_position())  
          out << "\n\t<constr/>";
        else {
          if (p->fixed_horizontal_position())  
            out << "\n\t<fixed-position/>";
          if (p->free_horizontal_position()&&!p->constr_horizontal_position())
            out << "\n\t<free-position/>";
          if (p->constr_horizontal_position()) 
            out << "\n\t<constr-position/>";
          if (p->fixed_height())               
            out << "\n\t<fixed-height/>";
          if (p->free_height() && !p->constr_height())                
            out << "\n\t<free-height/>";
          if (p->constr_height())              
            out << "\n\t<constr-height/>";
        }

        out << "</point>\n";
      }
  }
  
  {
    ClusterList::const_iterator i = obsdata.CL.begin();
    ClusterList::const_iterator e = obsdata.CL.end();
    while (i != e)
      {
        if (const g3Cluster* c = dynamic_cast<const g3Cluster*>(*i))
          {
            c->write_xml(out);
          }
        ++i;
      }
  }

  out << "\n</g3-model>\n";
}


void Model::update_init()
{
  return next_state_(init_);
}


void Model::update_parameters()
{
  if (!check_init()) update_init();

  for (Model::PointBase::iterator 
         i=points->begin(), e=points->end(); i!=e; ++i)
    {
      Point* point = (*i);
      point->N_.set_index(0);
      point->E_.set_index(0);
      point->U_.set_index(0);
    }

  return next_state_(params_);
}


void Model::update_observations()
{
  if (!check_parameters()) update_parameters();

  par_list->clear();
  active_obs->clear();
  dm_rows = dm_cols = dm_floats = 0;   // dimension and size of design matrix

  for (Model::ObservationData::iterator 
         i=obsdata.begin(), e=obsdata.end(); i!=e; ++i)
    {
      (*i)->revision_accept(this);
    }

  return next_state_(obsrvs_);
}


void Model::update_linearization()
{
  if (!check_observations()) update_observations();

  cerr << "\n###  update_linearization() : dm_floats, dm_rows, dm_cols " 
       << dm_floats << " " << dm_rows << " " << dm_cols<< "\n";

  // delete A;   A = new SparseMatrix<>(dm_floats, dm_rows, dm_cols);
  // delete B;   B = new BlockDiagonal<>;

  for (ObservationList::iterator 
         i=active_obs->begin(), e=active_obs->end(); i!=e; ++i)
    {
      (*i)->linearization_accept(this);
    }

  return next_state_(linear_);
}


void Model::update_adjustment()
{
  if (!check_linearization()) update_linearization();

  return next_state_(adjust_);
}
