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
 *  $Id: adj.h,v 1.7 2003/01/04 22:11:48 cepek Exp $
 */

#include <gamalib/exception.h>
#include <gamalib/sparse/smatrix.h>
#include <gamalib/sparse/sbdiagonal.h>
#include <gamalib/sparse/intlist.h>
#include <gamalib/ls/olssvd.h>
#include <gamalib/ls/olsgso.h>

#include <iostream>

#ifndef GaMaLib_Adj__adjustment_class__h
#define GaMaLib_Adj__adjustment_class__h


namespace GaMaLib {

  class AdjException : public Exception {
  public:

    AdjException(std::string s) : Exception(s) {} 
  };

  class AdjInputData;

  class Adj 
  {
  public:
    
    enum algorithm {gso, svd};
    
    Adj () : data(0), algorithm_(gso) { init(0); }
    virtual ~Adj();
    
    int n_obs() const { return n_obs_; }
    int n_par() const { return n_par_; }
    
    void set(const AdjInputData* inp) { init(inp); }
    void preferred_algorithm(Adj::algorithm);
  
    Vec x();
  
  private:
    
    const AdjInputData *data;
    BaseOLS<Double, GaMaLib::MatVecException> *least_squares;

    bool      solved;
    algorithm algorithm_;
    int       n_obs_, n_par_;
    Mat       A_dot;
    Vec       b_dot;
    Vec       x_;
  
    void init(const AdjInputData*);
    void init_least_squares();
    void cholesky(Cov& chol);                          // move it away!   
    void forwardSubstitution(const Cov& chol, Vec& v); // move it away!

  };
  

  class AdjInputData {
  public:

    AdjInputData();
    ~AdjInputData();
    
    void write_xml(std::ostream&) const;
    void read_xml(std::istream&);

    const SparseMatrix <> * mat () const { return A;     }
    const BlockDiagonal<> * cov () const { return pcov;  }
    const Vec               rhs () const { return prhs;  }
    const IntegerList  <> * minx() const { return pminx; } 

    void set_mat (SparseMatrix <> * p) { delete A;     A     = p; }
    void set_cov (BlockDiagonal<> * p) { delete pcov;  pcov  = p; }
    void set_rhs (Vec               p) {               prhs  = p; }
    void set_minx(IntegerList  <> * p) { delete pminx; pminx = p; } 

    void swap(AdjInputData *);

    /* Sparse project equations for uncorrelated observations. *
     * Defined here only for backward data compatibility       */
    void read_gama_local_old_format(std::istream&);


  private:

    friend class Adj;

    SparseMatrix <> * A;
    BlockDiagonal<> * pcov;
    Vec               prhs;
    IntegerList  <> * pminx;

  };
  
}  // namespace GaMaLib

#endif









