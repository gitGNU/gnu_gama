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
 *  $Id: orientation.h,v 1.1 2001/12/07 12:38:37 cepek Exp $
 */

#ifndef GaMaLib_Bod_Mer_VYPORPOS_H
#define GaMaLib_Bod_Mer_VYPORPOS_H

#include <functional>
#include <vector>
#include <gamalib/local/pobs/bearing.h>
#include <gamalib/local/gamadata.h>

namespace GaMaLib {

class Orientation {

   PointData&       PL;
   ObservationList& OL;

public:

   Orientation(PointData& p, ObservationList& o) : PL(p), OL(o) {}
   
   // L1 estimate of the standpoint orientation (iter. to the first direction)
   void orientation(ObservationList::const_iterator&, Double&, int&);
   
   // add all possible orientations for the observation list
   void add_all();
};

}
#endif
