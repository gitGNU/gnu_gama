/*  
    C++ Matrix/Vector templates (GNU GaMa / gMatVec 0.9.18)
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
 *  $Id: vecbase.h,v 1.5 2002/04/10 16:19:59 cepek Exp $
 *  http://www.gnu.org/software/gama/
 */

#ifndef gMatVec_VecBase__h_
#define gMatVec_VecBase__h_

#include <iostream>
#include <cstdarg>
#include <cmath>
#include <gmatvec/matvecbase.h>

namespace gMatVec {

  
template <class Float=double, class Exc=Exception>
class VecBase : public MatVecBase<Float, Exc> {

protected:

  VecBase() {}
  VecBase(Index nsz) :  MatVecBase<Float, Exc>(nsz) {}

public:

  typedef MatVecBase<Float, Exc>::iterator       iterator;
  typedef MatVecBase<Float, Exc>::const_iterator const_iterator;

  Index dim() const { return size(); }

  Float& operator()(Index n) { 
    Float* m = begin(); return m[--n]; 
  }
  Float  operator()(Index n) const { 
    const Float* m = begin(); return m[--n]; 
  } 

  void reset(Index n=0) { resize(n); }

  Float dot(const VecBase<Float, Exc> &B) const;

  Float length_sq() const;
  Float length()    const { return length_sq(); }
  Float norm_L1()   const;
  Float norm_L2()   const { return length_sq(); }
  Float norm_Linf() const;

};


template <class Float, class Exc>
Float VecBase<Float, Exc>::length_sq() const
  {
    const_iterator a = begin();
    const_iterator e = end();
    
    Float sum = 0, x;
    while (a != e) { x = *a++; sum += x*x; }
    
    return sum;
  }


template <class Float, class Exc>
Float VecBase<Float, Exc>::dot(const VecBase<Float, Exc> &B) const
  {
    if (dim() != B.dim())
      throw Exc(BadRank, "Float VecBase::dot(const VecBase&) const");
    
    const_iterator a = begin();
    const_iterator e = end();
    const_iterator b = B.begin();
    
    Float sum = 0;
    while (a != e) sum += *a++ * *b++;
    
    return sum;
  }


template <class Float, class Exc>
Float VecBase<Float, Exc>::norm_L1() const
  {
    const_iterator a = begin();
    const_iterator e = end();
    
    Float sum = 0;
    while (a != e) { sum += *a >= 0 ? *a : -(*a); ++a; }
    
    return sum;
  }


template <class Float, class Exc>
Float VecBase<Float, Exc>::norm_Linf() const
  {
    const_iterator a = begin();
    const_iterator e = end();
    
    Float norm = 0, x;
    while (a != e) 
      { 
        x = *a >= 0 ?  *a : -(*a); 
        ++a; 
        if (x > norm) norm = x;
      }
    
    return norm;
  }


template <class Float, class Exc>
std::istream& operator>>(std::istream& inp, VecBase<Float, Exc>& v)
  {
    Index size;
    inp >> size;
    
    if (size != v.dim())
      v.reset(size);
    
    MatVecBase<Float, Exc>::iterator b = v.begin();
    MatVecBase<Float, Exc>::iterator e = v.end();
    while (b != e)
      {
        inp >> *b;
        ++b;
      }
    
    return inp;
  }


}   // namespace gMatVec

#endif
