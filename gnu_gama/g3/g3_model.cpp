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
 *  $Id: g3_model.cpp,v 1.11 2003/06/08 08:11:13 cepek Exp $
 */

#include <gnu_gama/g3/g3_model.h>
#include <gnu_gama/g3/g3_cluster.h>
#include <gnu_gama/xml/dataparser.h>
#include <gnu_gama/outstream.h>


using namespace GNU_gama::g3;


Model::Model() 
{ 
  using namespace GNU_gama;

  points = new PointBase;
  obs    = new ObservationData;

  points->set_common_data(this); 
  set(&ellipsoid, ellipsoid_wgs84);
}


Model::~Model()
{
  delete points;
  delete obs;

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


void Model::write_xml(std::ostream& out) 
{
  using namespace std;
 
  GNU_gama::SaveFlags sf(out);
  out.setf(ios::fixed, ios::floatfield);
  out.precision(5);

  out 
    << DataParser::xml_start
    << "<g3-model>\n";

  {
    out << "\n";
    for (PointBase::const_iterator 
           b = points->begin(), e = points->end(); b != e; ++b)
      {
        const Point *p = *b;
        out << "<point>\t<id>" << p->name << "</id>";

        if (p->has_xyz())
          out << "\n\t"
              << "<x>" << p->X.value(0) << "</x> "
              << "<y>" << p->Y.value(0) << "</y> "
              << "<z>" << p->Z.value(0) << "</z>";

        if (p->has_height())
          out << "\n\t"
              << "<height>" << p->height.value(0) << "</height>";

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
    for (ObservationData::ClusterList::const_iterator
           b = obs->CL.begin(), e=obs->CL.end();  b != e;  ++b)
      if (const g3Cluster* c = dynamic_cast<const g3Cluster*>(*b))
      {
        
         c->write_xml(out);
      }

  }

  out 
    << "\n</g3-model>\n"
    << DataParser::xml_end;
;
}


void Model::pre_linearization()
{
  for (ObservationData::ClusterList::iterator
         i=obs->CL.begin(), e=obs->CL.end();  i != e;  ++i)
    if (g3Cluster* cluster = dynamic_cast<g3Cluster*>(*i))
    {
      cluster->parlist_init(this);

      for (GNU_gama::List<Observation*>::iterator
             obs = cluster->observation_list.begin(),
             end = cluster->observation_list.end();  obs != end;  ++obs)
        {
          (*obs)->parlist_init(this);
        }
    }
}
