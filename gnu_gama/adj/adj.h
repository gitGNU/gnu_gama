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
 *  $Id: adj.h,v 1.2 2005/03/27 17:43:26 cepek Exp $
 */

#include <gnu_gama/matvec.h>
#include <gnu_gama/sparse/smatrix.h>
#include <gnu_gama/sparse/sbdiagonal.h>
#include <gnu_gama/sparse/intlist.h>
#include <gamalib/ls/olssvd.h>
#include <gamalib/ls/olsgso.h>

#include <iostream>

#ifndef GaMaLib_Adj__adjustment_class__h
#define GaMaLib_Adj__adjustment_class__h


namespace GNU_gama {


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
  
    Vec<> x();
  
  private:
    
    const AdjInputData *data;
    GaMaLib::BaseOLS<double, Exception::matvec> *least_squares;

    bool      solved;
    algorithm algorithm_;
    int       n_obs_, n_par_;
    Mat <>    A_dot;
    Vec <>    b_dot;
    Vec <>    x_;
  
    void init(const AdjInputData*);
    void init_least_squares();
    void cholesky(Cov<>& chol);                            // move it away!   
    void forwardSubstitution(const Cov<>& chol, Vec<>& v); // move it away!

  };
  

  class AdjInputData {
  public:

    AdjInputData();
    ~AdjInputData();
    
    void write_xml(std::ostream&) const;
    void read_xml(std::istream&);

    const GNU_gama::SparseMatrix <> * mat () const { return A;     }
    const GNU_gama::BlockDiagonal<> * cov () const { return pcov;  }
    const           Vec          <>   rhs () const { return prhs;  }
    const GNU_gama::IntegerList  <> * minx() const { return pminx; } 

    void set_mat (SparseMatrix <> * p) { delete A;     A     = p; }
    void set_cov (BlockDiagonal<> * p) { delete pcov;  pcov  = p; }
    void set_rhs (Vec          <>   p) {               prhs  = p; }
    void set_minx(IntegerList  <> * p) { delete pminx; pminx = p; } 

    void swap(AdjInputData *);

    /* Sparse project equations for uncorrelated observations. *
     * Defined here only for backward data compatibility       */
    void read_gama_local_old_format(std::istream&);


  private:

    friend class Adj;

    SparseMatrix <> * A;
    BlockDiagonal<> * pcov;
    Vec          <>   prhs;
    IntegerList  <> * pminx;

  };
  
}  // namespace GaMaLib

#endif









