/*  
    Geodesy and Mapping C++ Library (GNU GaMa / GaMaLib)
    Copyright (C) 1999  Jiri Vesely <vesely@gama.fsv.cvut.cz>
                  2001  Ales Cepek  <cepek@fsv.cvut.cz>

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
 *  $Id: g2d_coordinates.h,v 1.4 2004/09/01 21:59:29 cepek Exp $
 */

/*************************************************************
 * computation of approximate coordinates:                   *
 * - single point                                            *
 * - list of points                                          *
 * - all points with unknown coordinates                     *
 *************************************************************/
 
#ifndef GaMaLib_g2d_coordinates_h__GaMaLib_Median_Pribl_s_H
#define GaMaLib_g2d_coordinates_h__GaMaLib_Median_Pribl_s_H

#include <algorithm>
#include <gamalib/local/gamadata.h>
#include <gamalib/local/median/g2d_helper.h>


namespace GaMaLib {

  class ApproximateCoordinates
    {
    private:
      
      // point list used to return results
      PointData&       SB;
      ObservationData& OD;
      ObservationList  SM;

      // Solution_state_tag -> see g2d_helper.h
      Solution_state_tag state;

      // list of selected points (with unknown coordinates)
      PointIDList selected;
      
      // list of solved points
      PointData solved;
      int depth;

      // number of points with known coordinates
      int known_coordinates_;
      
      bool Absent(PointID cb)
        {
          return 
            (std::find(selected.begin(),selected.end(),cb) == selected.end());
        }

      bool Local_observations(ObservationList::iterator sm, PointIDList sb)
        {
          bool pom = false;
          pom = (std::find(sb.begin(), sb.end(), (*sm)->from()) != sb.end()) &&
            (std::find(sb.begin(), sb.end(), (*sm)->to()) != sb.end());
          if(Angle* u = dynamic_cast<Angle*>(*sm))
            {
              pom = pom && (std::find(sb.begin(), sb.end(), u->fs()) 
                            != sb.end());
            }
          return pom;
        }

      void Reset();

      bool Observation_hasID(ObservationList::iterator m, 
                             PointIDList::iterator cb)
        {
          bool pom = (((*m)->from() == (*cb)) || ((*m)->to() == (*cb)));

          if (Angle* u = dynamic_cast<Angle*>(*m))
            {
              pom = pom || (u->fs() == (*cb));
            };
          return pom;
        }

      // true - at least two points exist with known coordinates
      bool Solvable_data(PointData& b);

      // true - at least 2 observations with point ID exist
      bool Necessary_observations(PointID ID);

      // in lists SM and SB finds points without coordinates and
      // stores their IDs in "selected" list
      void Find_missing_coordinates();

      // move point with ID from list From to list To
      void Move_point(PointData& From, PointData& To, PointID& ID);

      // solution by simple intersection; true - at least one point solved
      bool Solve_intersection(PointData& body, PointIDList& co);

      // computation of points that cannot be solved by a simple
      // intersection, eg polygonal traverse inserted between two
      // known points; true - coordinates of at least one point solved
      bool Solve_insertion();

      // combines both previous methods
      void Computational_loop();

      void copy_horizontal(const ObservationData& from, ObservationList& to);

    public:
      
      ApproximateCoordinates(PointData& b, ObservationData& m) 
        : SB(b), OD(m),  depth(0)
        {
          copy_horizontal(OD, SM);
          Reset();
        }
      
      ApproximateCoordinates(PointData& b, ObservationData& m, int vn)   
        : SB(b), OD(m), depth(vn)
        {
          copy_horizontal(OD, SM);
          Reset();
        }
      
      void Reset(PointData& b, ObservationList& m)
        {
          SB = b;
          SM = m;
          depth = 0;
          Reset();
        }
      
      // one point (even if already solved); true - succeeded to get
      // coordinates
      bool Calculation(PointID cb);
      
      // poinst in PointIDList (even those already solved); true -
      // succeded to get all coordinates
      bool Calculation(PointIDList cb);
      
      // all points without coordinates; true - succeeded to get all
      // coordinates
      bool Calculation()
        {
          // looking for all points not placed in the point list (ie
          // found only in observation list) and points without
          // coordinates
          state = calculation_done;

          // gamalib-0.9.51 (AC) ... added test on empty lists
          // gamalib-1.1.51 (AC) ... the test moved where it belongs ;-)
          if (SB.empty() || SM.empty()) return true;

          Find_missing_coordinates();
          if(!Solvable_data(SB))
            return false;
          Computational_loop();
          return All_is_solved();
        }

      bool All_is_solved() const
        {
          if(state <= calculation_not_done)
            throw g2d_exc("ApproximateCoordinates::"
                          "All_is_solved - nothing to do");
          return selected.empty();
        }

      PointIDList Unsolved() const
        {
          if(state <= calculation_not_done)
            throw g2d_exc("ApproximateCoordinates::"
                          "Unsolved - calculation not done");
          return selected;
        }

      PointData Solved() const
        {
          if(state <= calculation_not_done)
            throw g2d_exc("ApproximateCoordinates::"
                          "Solved - calculation not done");
          return solved;
        }  

      int Total_points () const
        {
          return SB.size();
        }

      int Total_observations () const
        {
          return SM.size();
        }

      // number of points with known coordinates (see Reset)
      int Known_coordinates() const
        { 
          return known_coordinates_; 
        }

    };


} // namespace GaMaLib


#endif




