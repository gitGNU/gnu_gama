/*  
    Geodesy and Mapping C++ Library (GNU GaMa / GaMaLib)
    Copyright (C) 2001  Ales Cepek <cepek@fsv.cvut.cz>

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
 *  $Id: lcoords.h,v 1.1 2001/12/07 12:38:37 cepek Exp $
 */

#ifndef GaMaLib____Local_Coordinate_System___h____
#define GaMaLib____Local_Coordinate_System___h____


namespace GaMaLib {

  class LocalCoordinateSystem {
  public:

    enum CS
    {
                                            // orientation of axes x and y :
      EN= 1, NW= 2, SE= 4, WS=  8,          //   plane left-handed  systems
      NE=16, SW=32, ES=64, WN=128,          //   plane right-handed systems

      left_handed  = (EN | NW | SE | WS),
      right_handed = (NE | SW | ES | WN)

    } local_coordinate_system;

    LocalCoordinateSystem(CS cs=NE) 
      : local_coordinate_system(cs)
      {
      }
    bool left_handed_coordinates () const 
      { 
        return local_coordinate_system & left_handed; 
      }
    bool right_handed_coordinates() const 
      { 
        return local_coordinate_system & right_handed; 
      }
    
  };

}   // namespace GaMaLib



#endif

















