/*  
    C++ Matrix/Vector templates (gMatVec 0.9.14)
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
 *  $Id: matvecbase.h,v 1.1 2001/12/07 11:59:46 cepek Exp $
 */

#ifndef gMatVec_MatVecBase__h_
#define gMatVec_MatVecBase__h_

#include <gmatvec/memrep.h>

 
namespace gMatVec {
  
template <class Float=double, class Exc=Exception>
class MatVecBase : public MemRep<Float, Exc> {

protected:

  MatVecBase() {}
  MatVecBase(Index nsz) : MemRep<Float, Exc>(nsz) {}
  MatVecBase(const MemRep<Float, Exc>& mem) : MemRep<Float, Exc>(mem) {}

  void mul(Float f, MatVecBase& X) const
    {
      if (size() != X.size())
        throw Exc(BadRank, "MatVecBase::mul(Float f, MatVecBase& X)");

      MemRep<Float, Exc>::const_iterator a = begin();
      MemRep<Float, Exc>::iterator x = X.begin();
      MemRep<Float, Exc>::iterator e = X.end();
      while (x != e)
        *x++ = *a++ * f;
    }
 
  void add(const MatVecBase& B, MatVecBase& X) const
    {
      if (size() != B.size() || size() != X.size())
        throw Exc(BadRank, "MatVecBase::add(const MatVecBase&, MatVecBase&)");

      MemRep<Float, Exc>::const_iterator a = begin();
      MemRep<Float, Exc>::const_iterator b = B.begin();
      MemRep<Float, Exc>::iterator x = X.begin();
      MemRep<Float, Exc>::iterator e = X.end();
      while (x != e)
        *x++ = *a++ + *b++;
    }
 
  void sub(const MatVecBase& B, MatVecBase& X) const
    {
      if (size() != B.size() || size() != X.size())
        throw Exc(BadRank, "MatVecBase::sub(const MatVecBase&, MatVecBase&)");

      MemRep<Float, Exc>::const_iterator a = begin();
      MemRep<Float, Exc>::const_iterator b = B.begin();
      MemRep<Float, Exc>::iterator x = X.begin();
      MemRep<Float, Exc>::iterator e = X.end();
      while (x != e)
        *x++ = *a++ - *b++;
    }

  // List initialiser helper class

  class ListInitialiser {

    MemRep<Float, Exc>::iterator x, e, first;

  public:

  ListInitialiser(MemRep<Float, Exc>::iterator begin, 
                  MemRep<Float, Exc>::iterator end) : x(begin), e(end) 
    {
      first = x;
      if (x != e) ++first;
    }
  ~ListInitialiser()
    {
      if (x != first && x != e)  
        throw Exc(BadRank, "ListInitialiser : "
                           "not enough elements in the initialisation list");
    }

  void add(Float p)
    {
      if (x == e)  throw Exc(BadRank, "ListInitialiser : "
                             "too many elements in the initialisation list");
      *x = p;
      ++x;
    }

  ListInitialiser& operator,(Float x)
    {
      add(x);
      return *this;
    }

  };

  ListInitialiser list_init(Float p)
    {
      ListInitialiser linit( begin(), end() );
      linit.add(p);
      return linit;
    }
 
public:

  void operator*=(Float f)
  {
    MemRep<Float, Exc>::iterator b = begin();
    MemRep<Float, Exc>::iterator e = end();
    while (b != e)
    *b++ *= f;
  }
  void operator/=(Float f) { operator*=(1/f); }

  void set_all(Float f)
  {
    MemRep<Float, Exc>::iterator b = begin();
    MemRep<Float, Exc>::iterator e = end();
    while (b != e)
      *b++ = f;
  }

  void set_zero() { set_all(0.0); }

};

}   // namespace gMatVec

#endif








