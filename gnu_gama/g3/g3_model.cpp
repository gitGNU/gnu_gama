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
 *  $Id: g3_model.cpp,v 1.28 2004/04/23 22:01:31 cepek Exp $
 */

#include <gnu_gama/g3/g3_model.h>
#include <gnu_gama/g3/g3_cluster.h>
#include <gnu_gama/outstream.h>
#include <gnu_gama/adj/adj.h>


using namespace std;
using namespace GNU_gama::g3;


Model::Model() 
{ 
  using namespace GNU_gama;

  points     = new PointBase;
  active_obs = new ObservationList;
  par_list   = new ParameterList;
  adj        = new Adj;

  points->set_common_data(this); 
  set(&ellipsoid, ellipsoid_wgs84);

  reset();
}


Model::~Model()
{
  delete  points;
  delete  active_obs;
  delete  par_list;
  delete  adj;
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
              << "<x>" << p->X() << "</x> "
              << "<y>" << p->Y() << "</y> "
              << "<z>" << p->Z() << "</z>";
        
        if (p->has_height())
          out << "\n\t"
              << "<height>" << p->height() << "</height>";
        
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


void Model::update_index(Parameter& p)
{
  if (!p.has_index()) 
    {
      if (p.free())
        p.set_index(++dm_cols);
      else
        p.set_index(1);

      par_list->push_back(&p);
    }
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
 
      point->N.set_index(0);
      point->E.set_index(0);
      point->U.set_index(0);


      if (point->fixed_horizontal_position())
        {
          point->N.set_fixed();
          point->E.set_fixed();
        }
      else if (point->constr_horizontal_position())
        {
          point->N.set_constr();
          point->E.set_constr();
        }
      else if (point->free_horizontal_position())
        {
          point->N.set_free();
          point->E.set_free();
        }
      else
        {
          point->N.set_unused();
          point->E.set_unused();
        }
      

      if (point->fixed_height())
        {
          point->U.set_fixed();
        }
      else if (point->constr_height())
        {
          point->U.set_constr();
        }
      else if (point->free_height())
        {
          point->U.set_free();
        }
      else
        {
          point->U.set_unused();
        }
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

  for (Model::ClusterList::iterator
         ci = obsdata.CL.begin(), ce = obsdata.CL.end(); ci!=ce; ++ci)
    {
      (*ci)->update();
    }
  
  return next_state_(obsrvs_);
}


void Model::update_linearization()
{
  if (!check_observations()) update_observations();

  cerr << "- - - - - - - - - - - - - - - - - - - - - - - .\n";
  adj_input_data = new AdjInputData;

  A = new SparseMatrix<>(dm_floats, dm_rows, dm_cols);
  rhs.reset(dm_rows);
  rhs_ind = 0;

  for (ObservationList::iterator 
         i=active_obs->begin(), e=active_obs->end(); i!=e; ++i)
    {
      (*i)->linearization_accept(this);
    }

  adj_input_data->set_mat(A);
  adj_input_data->set_rhs(rhs);

  {
    int minx=0;
    for (ParameterList::const_iterator 
           i=par_list->begin(), e=par_list->end(); i!=e; ++i)
      {
        if ((*i)->constr()) minx++;
      }
    
    if (minx)
      {
        GNU_gama::IntegerList<>* minlist = new GNU_gama::IntegerList<>(minx);
        GNU_gama::IntegerList<>::iterator m = minlist->begin();
        for (ParameterList::const_iterator 
               i=par_list->begin(), e=par_list->end(); i!=e; ++i)
          {
            const Parameter* p = *i;
            if (p->constr())
              {
                *m = p->index();
                ++m;
              }
          }
        
        adj_input_data->set_minx(minlist);
      }
  }
  
  {
    int nonzeroes=0, blocks=0;
    for (ClusterList::iterator
           ci = obsdata.CL.begin(), ce = obsdata.CL.end(); ci!=ce; ++ci)
      {
        if (int n = (*ci)->activeNonz())
          {
            nonzeroes += n;
            blocks++;
          }
      }

    GNU_gama::BlockDiagonal<>* 
      bd = new GNU_gama::BlockDiagonal<>(blocks, nonzeroes);

    for (ClusterList::const_iterator
           ci = obsdata.CL.begin(), ce = obsdata.CL.end(); ci!=ce; ++ci)
      {
        Cov C = (*ci)->activeCov();
        if (C.dim())
          {
            bd->add_block(C.dim(), C.bandWidth(), C.begin());
          }
      }

    adj_input_data->set_cov(bd);
  }
    
  adj->set(adj_input_data);

  return next_state_(linear_);
}


void Model::update_adjustment()
{
  cerr << "\n************ void Model::update_adjustment()\n\n";
  if (!check_linearization()) update_linearization();

  cerr << "x() = " << adj->x();
  for (ParameterList::iterator 
         i=par_list->begin(), e=par_list->end(); i!=e; ++i)
    {
      Parameter* p = *i;
      if (int k = p->index()) p->set_correction(adj->x()(k));
    }
  

  cerr << "\n------------ void Model::update_adjustment()\n\n";
  return next_state_(adjust_);
}


void Model::write_xml_adjustment_input_data(std::ostream& out)
{
  if (!check_linearization()) update_linearization();

  /* see also : GNU_gama::DataParser::xml_start  *
   *            GNU_gama::DataParser::xml_end    */
  adj_input_data->write_xml(out);
}


void Model::write_xml_adjustment_results(std::ostream& out)
{
  if (!check_adjustment()) update_adjustment();

  out << 
    "\n<!-- adjustment results  :"
    " dn / de / du  are in millimeters -->\n\n";
  
  out << "<adjustment-results>\n";

  for (ParameterList::iterator 
         i=par_list->begin(), e=par_list->end(); i!=e; ++i)
    {
      (*i)->write_xml_init();
    }
  for (ParameterList::iterator 
         i=par_list->begin(), e=par_list->end(); i!=e; ++i)
    {
      (*i)->write_xml(out);
    }
  
  out << "\n</adjustment-results>\n";
}


