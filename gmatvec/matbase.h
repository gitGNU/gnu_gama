/*  
    C++ Matrix/Vector templates (GNU GaMa / gMatVec 0.9.21)
    Copyright (C) 1999  Ales Cepek <cepek@fsv.cvut.cz>

    This file is part of the gMatVec C++ Matrix/Vector template library.
    
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
 *  $Id: matbase.h,v 1.8 2002/11/14 14:58:52 cepek Exp $
 *  http://www.gnu.org/software/gama/
 */

#ifndef gMatVec_MatBase__h_
#define gMatVec_MatBase__h_

#include <iostream>
#include <cstdarg>
#include <gmatvec/matvecbase.h>


namespace gMatVec {


template <class Float=double, class Exc=Exception>
class MatBase : public MatVecBase<Float, Exc> {

protected:

  Index row_;
  Index col_;

  MatBase() : row_(0), col_(0) {}
  MatBase(Index r, Index c, Index nsz) 
    : MatVecBase<Float, Exc>(nsz), row_(r), col_(c) {}
  MatBase(Index r, Index c, const MatBase& m) 
    : MatVecBase<Float, Exc>(m), row_(r), col_(c) {}
  virtual ~MatBase() {}

public:

  typedef MatVecBase<Float, Exc>::iterator       iterator;
  typedef MatVecBase<Float, Exc>::const_iterator const_iterator;

  Index rows() const { return row_; }
  Index cols() const { return col_; }

  Index min_rc() const { return row_ > col_ ? col_ : row_; }
  Index max_rc() const { return row_ < col_ ? col_ : row_; }

  virtual Float& operator()(Index r, Index c) = 0;
  virtual Float  operator()(Index r, Index c) const = 0;

  void reset() { row_ = col_ = 0; resize(0); }
  virtual void reset(Index r, Index c) {
    if (r != row_ || c != col_) {
      row_ = r; col_ = c; resize(r*c);
    }
  }

  void set_diagonal(Float d)
    {
      set_zero();
      for (Index i=1; i<=min_rc(); i++)
        this->operator()(i,i) = d;
    }
  void set_identity() { set_diagonal(1.0); }

  virtual void transpose() 
    {
      throw Exc(NotImplemented, "MatBase::transpose()");
    }

  virtual void invert() 
    {
      throw Exc(NotImplemented, "MatBase::invert()");
    }

  virtual std::istream& read(std::istream& inp) 
    {
      Index r, c;
      if (inp >> r >> c)
         {
            reset(r, c);
            for (Index i=1; i<=r; i++)
               for (Index j=1; j<=c; j++)
                  inp >> operator()(i,j);
         }
      return inp;
    }
  virtual std::ostream& write(std::ostream& out) const 
    {

      const int fw = out.width();
      out.width(fw);  
      out << rows() << " ";
      out.width(fw);  
      out << cols() << "\n\n";
      for (Index i=1; i<=rows(); i++)
        {
          for (Index j=1; j<=cols(); j++) {
            out.width(fw);  
            out << operator()(i,j) << " ";
          }
          out << '\n';
        }
      return out;
    }

};


template <class Float, class Exc> 
std::istream& operator>>(std::istream& inp, MatBase<Float, Exc>& M) 
  {
    return M.read(inp);
  }


template <class Float, class Exc> 
std::ostream& operator<<(std::ostream& out, const MatBase<Float, Exc>& M) 
  {
    return M.write(out);
  }


}   // namespace gMatVec

#endif
