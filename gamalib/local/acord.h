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
 *  $Id: acord.h,v 1.3 2003/03/13 20:22:42 cepek Exp $
 */

 
#ifndef GaMaLib_Acord___accord___header___h
#define GaMaLib_Acord___accord___header___h

#include <gamalib/local/gamadata.h>
#include <gamalib/local/acord/reduce_observations.h>
#include <fstream>
#include <algorithm>
#include <list>
#include <set>

namespace GaMaLib {


  class Acord 
    {
    public: 
      
      PointData&          PD;
      ObservationData&    OD;
      ReducedObservations RO;
	
      Acord(PointData& b, ObservationData& m);
      void execute();

      int  observations;
      int  given_xy, given_z, given_xyz;
      int  computed_xy, computed_z, computed_xyz;
      int  total_xy, total_z, total_xyz;
      bool missing_coordinates;

    private:

      class CountObs {
        mutable int& count;
      public:
        CountObs(int& c) : count(c) {}
        void operator()(const Observation*) const { count++; }
      };

      std::set<PointID> set_xyz, set_xy, set_z;

      class SlopeToHorizontal {  // s-dist & z-angle => horizontal distance
        mutable ObservationList& OL;
      public:
        SlopeToHorizontal(ObservationList& ol) : OL(ol) {}
        void operator()(const Observation*) const;
      }; 
      
    };

}   // namespace GaMaLib

#endif





