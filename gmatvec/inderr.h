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
 *  $Id: inderr.h,v 1.8 2002/11/14 14:58:52 cepek Exp $
 *  http://www.gnu.org/software/gama/
 */

#ifndef gMatVec__IndexErr__h_
#define gMatVec__IndexErr__h_

#include <cstddef>

namespace gMatVec {

  typedef size_t Index;

  enum Error 
  { 
    BadRank,  
    BadIndex, 
    Singular,
    BadRegularization,
    NoConvergence,  
    ZeroDivision,  
    NonPositiveDefinite,
    NotImplemented,
    StreamError 
  };
  
  struct Exception 
  {
    Error       error;
    const char* description;
    Exception(Error e, const char* t) : error(e), description(t) {}
  };
  

}   // namespace gMatVec

#endif




