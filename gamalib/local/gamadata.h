/*  
    Geodesy and Mapping C++ Library (GNU GaMa / GaMaLib)
    Copyright (C) 2000  Ales Cepek <cepek@fsv.cvut.cz>

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
 *  $Id: gamadata.h,v 1.4 2004/04/04 11:19:29 cepek Exp $
 */ 

#ifndef GaMaLib____ObservationData____and_other_data_objects___h_____
#define GaMaLib____ObservationData____and_other_data_objects___h_____

#include <gnu_gama/obsdata.h>

#include <gamalib/pointid.h>
#include <gamalib/local/lpoint.h>
#include <gamalib/observation.h>
#include <gamalib/cluster.h>
#include <gamalib/angobs.h>
#include <gamalib/local/lcoords.h>

#include <map>
#include <list>
#include <algorithm>

namespace GaMaLib {

  typedef std::list<PointID>    PointIDList;

  class PointData : public std::map <PointID, LocalPoint>,
                    public LocalCoordinateSystem,
                    public AngularObservations
    {
    };


  std::ostream& operator << (std::ostream&,     PointData&);
  std::ostream& operator << (std::ostream& str, ObservationData&);

  inline bool Consistent(const PointData& lcs)
    {
      return 
        (lcs.right_handed_coordinates() &&  lcs.right_handed_angles) ||
        (lcs. left_handed_coordinates() && !lcs.right_handed_angles) ; 
    }


}   // namepsace GaMaLib



#endif

















