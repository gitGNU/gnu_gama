/*  
    Geodesy and Mapping C++ Library (GNU GaMa / GaMaLib)
    Copyright (C) 2002  Ales Cepek <cepek@fsv.cvut.cz>

    This file is part of the GNU GaMa / GaMaLib C++ Library.
    
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
 *  $Id: adj.cpp,v 1.2 2002/11/22 17:46:22 cepek Exp $
 */

#include <gamalib/adj/adj.h>
#include <vector>

using namespace GaMaLib;

void AdjInputData::read_gama_local_old_format(std::istream& inp)
{
  using namespace std;
  vector<long>   ind;
  vector<double> flt;

  double cols, rows;
  inp >> cols >> rows;

  minx.reset(rows);
  rhs .reset(rows);

}


// -----------------------------------------------------------------


void Adj::init(const AdjInputData* inp)
{
  delete data; 
  data = inp; 
  solved = false; 
  n_obs_ = n_par_ = 0;

  if (data)
    {
      n_obs_ = data->A.rows();
      n_par_ = data->A.columns();
    }
}
