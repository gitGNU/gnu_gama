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
 *  $Id: acord.cpp,v 1.2 2002/04/04 13:19:37 cepek Exp $
 */

 
#include <gamalib/local/acord.h>
#include <gamalib/local/orientation.h>
#include <gamalib/local/median/g2d_coordinates.h>
#include <gamalib/local/acord/approx_heights.h>

#include <iomanip>
#include <cmath>

using namespace std;
using namespace GaMaLib;

Acord::Acord(PointData& b, ObservationData& m) 
  : PD(b), OD(m)
{
  missing_coordinates = false;
  observations = 0;
  given_xy = given_z = given_xyz = 0;
  computed_xy = computed_z = computed_xyz = 0;
  total_xy = total_z = total_xyz = 0;

  for (PointData::const_iterator i=PD.begin(); i!=PD.end(); ++i)
    {
      const PointID& c = (*i).first;
      const Point&   p = (*i).second;
      bool cp = p.test_xy();
      bool hp = p.test_z();

      if (cp && hp) given_xyz++, set_xyz.insert(c);
      else if (cp)  given_xy++,  set_xy .insert(c);
      else if (hp)  given_z++,   set_z  .insert(c);
    }

  OD.for_each(CountObs(observations));

  if (Consistent(PD, OD)) return;

  // switch (PD.local_coordinate_system)
  //   {
  //   case LocalCoordinateSystem::NE: 
  //     PD.local_coordinate_system = LocalCoordinateSystem::NW;
  //     break;
  //   case LocalCoordinateSystem::SW: 
  //     PD.local_coordinate_system = LocalCoordinateSystem::SE;
  //     break;
  //   case LocalCoordinateSystem::ES: 
  //     PD.local_coordinate_system = LocalCoordinateSystem::EN;
  //     break;
  //   case LocalCoordinateSystem::WN: 
  //     PD.local_coordinate_system = LocalCoordinateSystem::WS;
  //     break;
  //   case LocalCoordinateSystem::EN: 
  //     PD.local_coordinate_system = LocalCoordinateSystem::ES;
  //     break;
  //   case LocalCoordinateSystem::NW: 
  //     PD.local_coordinate_system = LocalCoordinateSystem::NE;
  //     break;
  //   case LocalCoordinateSystem::SE: 
  //     PD.local_coordinate_system = LocalCoordinateSystem::SW;
  //     break;
  //   case LocalCoordinateSystem::WS: 
  //     PD.local_coordinate_system = LocalCoordinateSystem::WN;
  //     break;
  //   default:
  //     return;
  //   }
  for (PointData::iterator ii=PD.begin(); ii!=PD.end(); ++ii)
    {
      Point& p = (*ii).second;

      if (p.test_xy()) p.set_xy(p.x(), -p.y());
    }
}

void Acord::execute()
{
  try
    {
      int all;

      do {
        all = total_z + total_xy + total_xyz;

        ApproximateHeights ah(PD, OD);
        ah.execute();
        
        {
          // all transformed slope distances go to a single standpoint
          StandPoint* standpoint = new StandPoint(&OD);
          OD.for_each(Acord::SlopeToHorizontal(standpoint->observation_list));
          // bind observations to the cluster
          standpoint->update();
          // insert standpoint into `observation data'
          OD.CL.push_back(standpoint);
          
          ApproximateCoordinates ps(PD, OD);
          ps.Calculation();
          
          OD.CL.pop_back();
          delete standpoint;
        }
        
        ObservationList local;
        OD.for_each(Observation::CopyTo(local));
        
        Orientation orp(PD, local);
        orp.add_all();
        
        for (PointData::const_iterator i=PD.begin(); i!=PD.end(); ++i)
          {
            const PointID& c = (*i).first;
            const Point&   p = (*i).second;
            bool cp = p.test_xy();
            bool hp = p.test_z();
            
            if (cp && hp) total_xyz++;
            else if (cp)  total_xy++;
            else if (hp)  total_z++;
            
            if (p.active_xy() && !cp) missing_coordinates = true;
            if (p.active_z()  && !hp) missing_coordinates = true;
            
            int t = 0;
            if      (set_xyz.find(c) != set_xyz.end()) t = 1;
            else if (set_xy .find(c) != set_xy .end()) t = 2;
            else if (set_z  .find(c) != set_z  .end()) t = 3;
            
            if (cp && hp)
              {
                switch (t) 
                  {
                  case 0 : computed_xyz++; break;
                  case 1 :                 break;
                  case 2 : computed_z++;   break;
                  case 3 : computed_xy++;  break;
                  default:                 break;
                  }
              }
            else if (cp)
              {
                if (t!=2)  computed_xy++; 
              }
            else if (hp)
              {
                if (t!=3)  computed_z++;
              }
          }

      } while (all != (total_z + total_xy + total_xyz));
    }
  catch(...)
    {
      // we must handle the case when there are no observations here
      throw;
    }
}

void Acord::SlopeToHorizontal::operator()(Observation* obs)
{
  S_Distance* s = dynamic_cast<S_Distance*>(obs);
  if (s == 0) return;

  StandPoint* c = dynamic_cast<StandPoint*>(s->ptr_cluster());
  if (c == 0) return;   // this should newer happen

  PointID from = s->from();
  PointID to   = s->to();

  // look for a zenith angle corresponding to the given slope distance
  ObservationList::iterator i   = c->observation_list.begin();
  ObservationList::iterator end = c->observation_list.end();
  for ( ;i!=end; ++i)
    if (Z_Angle* z = dynamic_cast<Z_Angle*>(*i))
      if (from == z->from() && to == z->to())
        {
          // ... and fake a horizontal distance
          OL.push_back(new Distance(from, to, s->value()*sin(z->value())));
          return;
        }

}







