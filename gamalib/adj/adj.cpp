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
 *  $Id: adj.cpp,v 1.3 2002/11/22 21:06:30 cepek Exp $
 */

#include <gamalib/adj/adj.h>
#include <vector>

using namespace GaMaLib;

void AdjInputData::read_gama_local_old_format(std::istream& inp)
{
  // see void LocalNetwork::project_equations(std::ostream& out)

  using namespace std;
  vector<long>   ind;
  vector<double> flt;

  long cols, rows;
  inp >> cols >> rows;                    // dimensions 

  minx.reset(rows);
  rhs .reset(rows);
  gMatVec::Vec<> c(rows);

  IntegerList<>::iterator m = minx.begin();

  long floats=0;
  for (long nonz, n, k, i=1; i<=rows; i++)
    {
      inp >> nonz;                        // number of nonzero elements
      *m++ = nonz;
      floats += nonz;
      for (k=1; k<=nonz; k++)
        {
          inp >> n;                       // indexes of nonzero elements
          ind.push_back(n);
        }

      double d;
      inp >> d;                           // i-th weight
      c(i) = 1/d;
      inp >> d;                           // i-th right-hand site element
      rhs(i) = d;
      for (k=1; k<=nonz; k++)
        {
          inp >> d;                       // nonzeroes elements
          flt.push_back(d);
        }
    }


  A.reset(floats, rows, cols);
  m = minx.begin();
  for (long r=1; r<=rows; r++)
    {
      A.new_row();
      long nonz = *m++;
      for (long i=0; i<nonz; i++)  A.add_element(flt[i], ind[i]);       
    }

  minx.reset();  // no regularization is defined for singular systems

  cov.reset(1, rows);
  cov.add_block(rows, 0, c.begin());
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
