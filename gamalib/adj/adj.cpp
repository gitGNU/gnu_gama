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
 *  $Id: adj.cpp,v 1.10 2003/01/04 15:51:51 cepek Exp $
 */

#include <gamalib/adj/adj.h>
#include <gamalib/xml/dataparser.h>
#include <vector>
#include <cstddef>
#include <algorithm>

using namespace GaMaLib;
using namespace std;

AdjInputData::AdjInputData()
{
  A     = 0;    // SparseMatrix<>  *
  pcov  = 0;    // BlockDiagonal<> *
  pminx = 0;    // IntegerList<>   * 
}



AdjInputData::~AdjInputData()
{
  delete A;
  delete pcov;
  delete pminx;
}



void AdjInputData::swap(AdjInputData *data)
{
  std::swap(A     , data->A    );     
  std::swap(pcov  , data->pcov );  
  std::swap(prhs  , data->prhs );
  std::swap(pminx , data->pminx); 
}



void AdjInputData::write_xml(std::ostream& out) const
{
  const char* indent = "  ";

  out << "\n" << indent << "<adj-input-data>\n";

  out << "\n" << indent << "  <sparse-mat>\n";
  out << indent << "    "
      << "<rows>" << A->rows() << "<rows> "
      << "<cols>" << A->columns() << "</cols> "
      << "<nonz>" << A->nonzeroes() << "</nonz>\n";

  for (Index m, k=1; k<=A->rows(); k++)
    {
      out << indent << "      <row>";
      out << " <nonz>" << (A->end(k)-A->begin(k)) << "</nonz>";
 
      double* n = A->begin(k);
      double* e = A->end  (k);
      for(Index* i=A->ibegin(k) ; n!=e; n++, i++, m++)
        {
          out << "\n        "
              << "<int>" << *i << "</int>"
              << "<flt>" << *n << "</flt>";
        }
      out << "\n        </row>\n";
    }

  out << indent << "  </sparse-mat>\n";

  // =================================================================

  const long blocks = pcov->blocks();

  out << "\n" << indent << " <block-diagonal>\n"
      << indent << "    <blocks>" << blocks << "</blocks>\n";

  for (long b=1; b<=blocks; b++)
    {
      long dim   = pcov->dim(b);
      long width = pcov->width(b);

      out << indent << "      <block> <dim>" 
          << dim    << "</dim> <width>" 
          << width  << "</width>\n";

      const double *m = pcov->begin(b);
      const double *e = pcov->end(b);
      while (m != e)
        cout << indent << "      <flt>" << *m++ << "</flt>\n";

      out << indent << "      </block>\n";
    }

  out << indent << "  </block-diagonal>\n";

  // =================================================================

  out << "\n" <<  indent << "  <vector>\n"
      << indent 
      << "    <dim>" << prhs.dim() << "</dim>\n";

  for (long i=1; i<=prhs.dim(); i++)
        cout << indent << "      <flt>" << prhs(i) << "</flt>\n";

  out << indent << "  </vector>\n";

  // =================================================================
 
  if (pminx)
    {
      out << "\n" <<  indent << "  <array>\n"
          << indent 
          << "    <dim>" << pminx->dim() << "</dim>\n";
      
      const std::size_t *indx = pminx->begin();
      for (long i=1; i<=pminx->dim(); i++)
        cout << indent << "      <int>" << *indx++ << "</int>\n";
       
      out << indent << "  </array>\n";
    }

  // =================================================================

  out << "\n" << indent << "</adj-input-data>\n";
}



void AdjInputData::read_xml(std::istream& inp)
{
  string            line;
  list<DataObject*> objects;
  DataParser dp(objects);

  while (getline(inp, line))
    {
      line += '\n';
      dp.xml_parse(line.c_str(), line.length(), 0);
    }
  dp.xml_parse("", 0, 1);

  for (list<DataObject*>::const_iterator i=objects.begin(); 
       i!=objects.end(); ++i)
    {
      if (AdjInputDataObject *adj = dynamic_cast<AdjInputDataObject*>(*i))
        {
          // take over the data from DataObject
          swap(adj->data);
        }

      delete *i;
    }
}



void AdjInputData::read_gama_local_old_format(std::istream& inp)
{
  // see void LocalNetwork::project_equations(std::ostream& out)

  using namespace std;
  vector<long>   ind;
  vector<double> flt;

  long cols, rows;
  inp >> cols >> rows;                    // dimensions 

  delete pminx; 
  pminx = 0;        // no regularization is defined for singular systems
  prhs.reset(rows);

  gMatVec::Vec<> c(rows);

  IntegerList<> tmplist(rows);        
  IntegerList<>::iterator m = tmplist.begin();

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
      prhs(i) = d;
      for (k=1; k<=nonz; k++)
        {
          inp >> d;                       // nonzeroes elements
          flt.push_back(d);
        }
    }

  delete A;
  A = new SparseMatrix<>(floats, rows, cols);

  m = tmplist.begin();
  for (long k=0, r=1; r<=rows; r++)
    {
      A->new_row();
      long nonz = *m++;
      for (long i=1; i<=nonz; i++, k++)  A->add_element(flt[k], ind[k]);       
    }

  delete pcov;
  pcov = new BlockDiagonal<>(1, rows);
  pcov->add_block(rows, 0, c.begin());
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
      n_obs_ = data->A->rows();
      n_par_ = data->A->columns();
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

  A_dot.reset(data->A->rows(), data->A->columns());
  A_dot.set_zero();
  b_dot.reset(data->A->rows());

  for (Index k=1; k<=data->A->rows(); k++)
    {
      double* n = data->A->begin(k);
      double* e = data->A->end  (k);
      for(size_t* i=data->A->ibegin(k) ; n!=e; n++, i++)  A_dot(k,*i) = *n; 
     }

  for (Index i, j, dim, width, r=0, b=1; 
       b<=data->pcov->blocks(); b++, r += dim)
    {
      dim   = data->pcov->dim(b);
      width = data->pcov->width(b);
      Cov C(dim, width);

      const Double *p = data->pcov->begin(b), *e = data->pcov->end(b);
      Cov::iterator c = C.begin();
      while (p != e) *c++ = *p++;
      cholesky(C);

      Vec t(dim);
      for (i=1; i<=dim; i++) t(i) = data->prhs(r+i);
      forwardSubstitution(C, t);
      for (i=1; i<=dim; i++) b_dot(r+i) = t(i);

      for (j=1; j<=data->A->columns(); j++)
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

