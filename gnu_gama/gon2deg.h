/*  
    Geodesy and Mapping C++ Library (GNU Gama)
    Copyright (C) 2004  Ales Cepek <cepek@fsv.cvut.cz>

    This file is part of the GNU Gama C++ Library.
    
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
 *  $Id: gon2deg.h,v 1.3 2004/04/03 11:06:37 cepek Exp $
 */

#include <string>

#ifndef GNU_gama_gons_to_degrees_h___GNU_Gama_gon2deg__gon2deg
#define GNU_gama_gons_to_degrees_h___GNU_Gama_gon2deg__gon2deg


namespace GNU_gama {

  // sign 0  conversion without sign
  //      1  sign left-padded
  //      2  sign right-padded 
  //      3  signed with leading spaces trimmed
  std::string gon2deg(double gon,  int sign, int prec);

  bool        deg2gon(std::string, double &);

}


#endif








