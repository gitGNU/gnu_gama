/*  
    GNU Gama -- adjustment of geodetic networks
    Copyright (C) 2003  Ales Cepek <cepek@fsv.cvut.cz>

    This file is part of the GNU Gama C++ library.
    
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
 *  $Id: model.h,v 1.1 2003/11/23 19:22:31 cepek Exp $
 */


#ifndef GNU_gama__modelmodel_h_gnugamamodeltraits___gnu_gama_gmodel___h
#define GNU_gama__modelmodel_h_gnugamamodeltraits___gnu_gama_gmodel___h


namespace GNU_gama {


  template <class Observation> 
    class Revision
    {
    public: 
      
      virtual bool revision(Observation*) = 0;
    };


  template <class Observation> 
    class Derivative 
    {
    public:
      
      virtual double derivative(Observation*) = 0;
    };
    
}

#endif
