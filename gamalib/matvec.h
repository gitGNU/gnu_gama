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
 *  $Id: matvec.h,v 1.4 2005/03/27 17:43:26 cepek Exp $
 */

#ifndef GaMaLib_Bod_Mer_MatVec_H
#define GaMaLib_Bod_Mer_MatVec_H

#include <gamalib/exception.h>
#include <gamalib/float.h>
#include <gmatvec/svd.h>
#include <gmatvec/bandmat2.h>

namespace GaMaLib {

  class MatVecException : public GaMaLib::Exception {
  public:
    const int error;
    MatVecException(int e, const char *s) : GaMaLib::Exception(s), error(e) {}
  };
  
  typedef GNU_gama::Index Index;
  
  typedef GNU_gama::Vec     <double, MatVecException>   Vec;
  typedef GNU_gama::Mat     <double, MatVecException>   Mat;
  typedef GNU_gama::SVD     <double, MatVecException>   SVD;
  typedef GNU_gama::BandMat2<double, MatVecException>   CovMat;  // covariances
  
}      // GaMaLib

#endif




