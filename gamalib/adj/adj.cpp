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
 *  $Id: adj.cpp,v 1.7 2002/12/15 16:12:44 cepek Exp $
 */

#include <gamalib/adj/adj.h>
#include <vector>

using namespace GaMaLib;


void AdjInputData::write_xml(std::ostream& out) const
{
  const char* indent = "  ";

  out << "\n" << indent << "<adj-input-data>\n";

  out << "\n" << indent << "  <sparse-mat>\n";
  out << indent << "    "
      << "<rows>" << A.rows() << "<rows> "
      << "<cols>" << A.columns() << "</cols> "
      << "<nonz>" << A.nonzeroes() << "</nonz>\n";

  for (Index m, k=1; k<=A.rows(); k++)
    {
      out << indent << "      <row>";
      out << " <nonz>" << (A.end(k)-A.begin(k)) << "</nonz>";
 
      double* n = A.begin(k);
      double* e = A.end  (k);
      for(Index* i=A.ibegin(k) ; n!=e; n++, i++, m++)
        {
          out << "\n        "
              << "<int>" << *i << "</int>"
              << "<flt>" << *n << "</flt>";
        }
      out << "\n        </row>\n";
    }

  out << indent << "  </sparse-mat>\n";

  out << "\n" << indent << "</adj-input-data>\n";
}



void AdjInputData::read_xml(std::istream& inp)
{
}



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
  for (long k=0, r=1; r<=rows; r++)
    {
      A.new_row();
      long nonz = *m++;
      for (long i=1; i<=nonz; i++, k++)  A.add_element(flt[k], ind[k]);       
    }

  minx.reset();  // no regularization is defined for singular systems

  cov.reset(1, rows);
  cov.add_block(rows, 0, c.begin());
}


// -----------------------------------------------------------------

Adj::~Adj()
{ 
  delete data; 
  delete least_squares;
}



void Adj::init(const AdjInputData* inp)
{
  delete data; 
  data = inp; 
  least_squares = 0;
  solved = false; 
  n_obs_ = n_par_ = 0;

  if (data)
    {
      n_obs_ = data->A.rows();
      n_par_ = data->A.columns();
    }
}



void Adj::init_least_squares()
{
  delete least_squares;

  switch (algorithm_) 
    {
    case svd: 
      least_squares = new OLSsvd<Double, GaMaLib::MatVecException>;
      break;
    case gso: 
      least_squares = new OLSgso<Double, GaMaLib::MatVecException>;
      break;
    default:
      throw AdjException("### unknown algorithm");
    }

  A_dot.reset(data->A.rows(), data->A.columns());
  A_dot.set_zero();
  b_dot.reset(data->A.rows());

  for (Index k=1; k<=data->A.rows(); k++)
    {
      double* n = data->A.begin(k);
      double* e = data->A.end  (k);
      for(size_t* i=data->A.ibegin(k) ; n!=e; n++, i++)  A_dot(k,*i) = *n; 
     }

  for (Index i, j, dim, width, r=0, b=1; b<=data->cov.blocks(); b++, r += dim)
    {
      dim   = data->cov.dim(b);
      width = data->cov.width(b);
      Cov C(dim, width);

      const Double *p = data->cov.begin(b), *e = data->cov.end(b);
      Cov::iterator c = C.begin();
      while (p != e) *c++ = *p++;
      cholesky(C);

      Vec t(dim);
      for (i=1; i<=dim; i++) t(i) = data->rhs(r+i);
      forwardSubstitution(C, t);
      for (i=1; i<=dim; i++) b_dot(r+i) = t(i);

      for (j=1; j<=data->A.columns(); j++)
        {
          for (i=1; i<=dim; i++) t(i) = A_dot(r+i,j);
          forwardSubstitution(C, t);
          for (i=1; i<=dim; i++) A_dot(r+i,j) = t(i);
        }
    }

  least_squares->reset(A_dot, b_dot);
  solved = true;
}



void Adj::preferred_algorithm(Adj::algorithm alg)
{
  switch (alg)
    {
    case svd:
    case gso:
      solved = false;
      algorithm_ = alg;
      break;

    default:
      throw AdjException("### unknown algorithm");
    }
}



Vec Adj::x()
{
  if (!solved) 
    {
      init_least_squares();
      x_ = least_squares->solve();
    }

  return x_;
}


/* ######################################################################
 * functions cholesky() and forwardSubstitution are identical in Adj and
 * in LocalNetwork and shall be moved to a single class
 * ###################################################################### */

void Adj::cholesky(Cov& chol)
{
  chol.cholDec();

  using namespace std;
  const Index N = chol.rows();
  const Index b = chol.bandWidth();

  for (Index m, j, i=1; i<=N; i++)
    {
      double d = sqrt(chol(i,i));
      chol(i,i) = d;

      m = i+b;  if(N < m) m = N;    // m = min(N, i+b);
      
      for (j=i+1; j<=m; j++) chol(i,j) *= d;
    }
}

void Adj::forwardSubstitution(const Cov& chol, Vec& v)
{
  using namespace std;
  const Index N = chol.rows();
  const Index b = chol.bandWidth();

  for (Index m, i=1; i<=N; i++)
    {
      if (i > b+1) m = i - b;
      else         m = 1;
      for (Index j=m; j<i; j++) v(i) -= chol(i,j)*v(j);

      v(i) /= chol(i,i);
    }
}

