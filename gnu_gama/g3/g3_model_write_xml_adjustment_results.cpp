/*  
    GNU Gama -- adjustment of geodetic networks
    Copyright (C) 2003  Ales Cepek <cepek@gnu.org>

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
    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/

/*
 *  $Id: g3_model_write_xml_adjustment_results.cpp,v 1.1 2005/10/19 16:12:02 cepek Exp $
 */

#include <gnu_gama/g3/g3_model.h>
#include <gnu_gama/g3/g3_cluster.h>
#include <gnu_gama/outstream.h>
#include <gnu_gama/adj/adj.h>
#include <iomanip>


using namespace std;
using namespace GNU_gama::g3;


void Model::write_xml_adjustment_results(std::ostream& out)
{
  if (!check_adjustment()) update_adjustment();

  const Vec<>& r = adj->r();

  out << "\n<adjustment-statistics>\n\n";

  Adj::algorithm alg = adj->get_algorithm();
  out << "<algorithm> ";
  if      (alg == Adj::gso)      cout << "gso";
  else if (alg == Adj::svd)      cout << "svd";
  else if (alg == Adj::cholesky) cout << "cholesky";
  else                           cout << "unknown";
  out << " </algorithm>\n\n";

  {
    gama_ellipsoid id = gama_ellipsoid(ellipsoid.id);
    
    out.setf(ios_base::fixed, ios_base::floatfield);
    out << "<ellipsoid> "
        << "<caption> " << gama_ellipsoid_caption[id] << " </caption>\n";
    out << "            <id>      " << gama_ellipsoid_id[id] 
        << "         </id>\n";
    out.precision(5);
    out << "            <a>       " << ellipsoid.a() << " </a>\n"; 
    out.precision(5);
    out << "            <b>       " << ellipsoid.b() << " </b>\n"; 
    out << "            </ellipsoid>\n\n"; 
  }

  out << "<parameters>" << setw(5) << dm_cols       << " </parameters>\n";
  out << "<equations> " << setw(5) << dm_rows       << " </equations>\n";
  out << "<defect>    " << setw(5) << adj->defect() << " </defect>\n";

  out << "<redundancy>" << setw(5) << redundancy    << " </redundancy\n\n";

  out.setf(ios_base::scientific, ios_base::floatfield);
  out.precision(5);
  double rtr = trans(r)*r;
  out << "<sum-of-squares>    " << rtr << " </sum-of-squares>\n";

  double sigma_apriori = apriori_sd*apriori_sd;
  out << "<sigma-apriori>     " << sigma_apriori << " </sigma-apriori>\n";

  double sigma_aposteriori = aposteriori_sd*aposteriori_sd;
  out << "<sigma-aposteriori> "<<sigma_aposteriori<<" </sigma-aposteriori>\n";
  
  out << "\n</adjustment-statistics>\n\n";

  // -----------------------------------------------------------------------

  out << "\n"
    "<!-- adjustment results    : dn / de / du  are in millimeters -->\n"
    "<!-- deflection of vertical: db / dl       are in arc seconds -->\n\n";
  
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

  out << "\n<!-- adjusted observations -->\n";

  Index index = 1;
  for (ObservationList::iterator 
         i=active_obs->begin(), e=active_obs->end(); i!=e; ++i)
    {
      (*i)->write_xml_adjusted(out, this, index);
      index += (*i)->dimension();
    }

  
  out << "\n</adjustment-results>\n";
}


void Model::write_xml_adjustment_results_points      (std::ostream& out)
{
}

void Model::write_xml_adjustment_results_observations(std::ostream& out)
{
}

