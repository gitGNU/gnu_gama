/*  
    Geodesy and Mapping C++ Library (GNU GaMa / GaMaLib)
    Copyright (C) 1999  Ales Cepek <cepek@fsv.cvut.cz>

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
 *  $Id: matvec.h,v 1.2 2002/09/16 10:42:04 cepek Exp $
 */

#ifndef GaMaLib_Bod_Mer_MatVec_H
#define GaMaLib_Bod_Mer_MatVec_H

#include <gamalib/float.h>
#include <gmatvec/svd.h>
#include <gmatvec/bandmat2.h>

namespace GaMaLib {

  class MatVecException : public GaMaLib::Exception {
  public:
    const int error;
    MatVecException(int e, const char *s) : GaMaLib::Exception(s), error(e) {}
  };
  
  typedef gMatVec::Index Index;
  
  typedef gMatVec::Vec<double, MatVecException> Vec;
  typedef gMatVec::Mat<double, MatVecException> Mat;
  typedef gMatVec::SVD<double, MatVecException> SVD;
  
  // GaMa covariance matrix
  typedef gMatVec::BandMat2<double, MatVecException> Cov;
  
}      // GaMaLib

#endif




