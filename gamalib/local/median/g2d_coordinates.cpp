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
    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/

/*
 *  $Id: g2d_coordinates.cpp,v 1.10 2005/05/07 18:06:20 cepek Exp $
 */

/*************************************************************
 * computation of approximate coordinates:                   *
 *************************************************************/
 
#include <gamalib/local/median/g2d_coordinates.h>
#include <gamalib/local/median/g2d_point.h>

using namespace std;
using namespace GaMaLib;

// private

void ApproximateCoordinates::
copy_horizontal(const ObservationData& from, ObservationList& to)
{
  for (ObservationData::
         const_iterator i=from.begin(), e=from.end(); i!=e; ++i)
    {
        Observation* obs = const_cast<Observation*>(*i);

        if      (dynamic_cast<Direction*>(obs))  to.push_back(obs);
        else if (dynamic_cast<Angle*    >(obs))  to.push_back(obs);
        else if (dynamic_cast<Distance *>(obs))  to.push_back(obs);
    }
}

void ApproximateCoordinates::Reset()
{

  state = calculation_not_done;
  selected.erase(selected.begin(), selected.end());
  solved.erase(solved.begin(), solved.end());
  
  known_coordinates_ = 0;
  for (PointData::const_iterator i=SB.begin(); i!=SB.end(); ++i)
    if ((*i).second.test_xy())
      ++known_coordinates_;

}


// to be solvable, at least 2 points with known coordinates must exist

bool ApproximateCoordinates::Solvable_data(PointData& b)
{

  bool first = false, second = false;
  bool tmp;
  PointData::iterator i = b.begin();
  do
    {
      tmp = (*i).second.test_xy() && Absent((*i).first);
      if(first)
        second = tmp;
      else
        first = tmp;
      i++;
    }
  while(!(second || i == SB.end()));
  return second;

}  // bool ApproximateCoordinates::Solvable_data(PointData& b)


bool ApproximateCoordinates::Necessary_observations(PointID id)
{

  // for coordinates to be solvable, at least 2 observations are needed
  // for details see ApproxPoint

  bool first = false, second = false;
  bool tmp;
  ObservationList::iterator i = SM.begin();
  do
    {
      tmp = (((*i)->from() == id) || ((*i)->to() == id));
      if(first)
        second = tmp;
      else
        first = tmp;
      // is second target available?
      Angle* u = dynamic_cast<Angle*>(*i);
      if( u && (u->fs() == id))
        if(first)
          second = true;
        else
          first = true;
      i++;
    }
  while(!(second || i == SM.end()));
  return second;

}	// bool ApproximateCoordinates::Necessary_observations()


void ApproximateCoordinates::Find_missing_coordinates()
{

  selected.erase(selected.begin(), selected.end());

  // from observation list points we fetch points that are not in SB
  for(ObservationList::iterator i = SM.begin(); i != SM.end(); i++)
    {
      if((SB.find((*i)->from()) == SB.end()) && Absent((*i)->from()))
        selected.push_back((*i)->from());
      if((SB.find((*i)->to()) == SB.end()) && Absent((*i)->to()))
        selected.push_back((*i)->to());
      // is second target available?
      Angle*u = dynamic_cast<Angle*>(*i);
      if(u && (SB.find(u->fs()) == SB.end()) && Absent(u->fs()))
        selected.push_back(u->fs());
    }

  // from point list we fetch the points with test_xy() == false
  {  // VC++ {}
    for(PointData::iterator i = SB.begin(); i != SB.end(); i++)
      if((!(*i).second.test_xy()) && Absent((*i).first))
        selected.push_back((*i).first);
  }

  // final sort and removal of duplicities in the list of seleted points
  selected.sort();
  selected.unique();

}	// void ApproximateCoordinates::Find_missing_coordinates()


void ApproximateCoordinates::Move_point(PointData& from, PointData& to, 
                                        PointID& what)
{

  PointData::iterator i = from.find(what);
  PointData::iterator j = to.find(what);
  if(i != from.end())
    if(j != to.end())
      (*j).second.set_xy((*i).second.x(), (*i).second.y());
    else
      to[what] = LocalPoint::XY((*i).second.x(), (*i).second.y());

}    // ApproximateCoordinates::Move_point(PointData&, PointData&, PointID&)


bool ApproximateCoordinates::Solve_intersection(PointData& points, 
                                                PointIDList& what)
{
  if(what.empty()) return false;

  ApproxPoint PB(&points,&SM);
  bool finished = false;
  bool success = false;
  LocalPoint bb;
  PointIDList::iterator i;
  PointData::iterator j;
  while (!finished)
    {
      finished = true;
      i = what.begin();
      while(i != what.end())
        {
          PB.Calculation(*i);
          if(PB.State() == unique_solution)
            {
              finished = false;
              success = true;
              bb = PB.Solution();
              j = points.find(*i);
              if(j != points.end())
                (*j).second.set_xy(bb.x(), bb.y());
              else
                points[*i] = LocalPoint::XY(bb.x(), bb.y());
              solved[*i] = bb;
              i = what.erase(i);
            }
          else
            i++;
        }
    }

  return success;

}  // ApproximateCoordinates::Solve_intersection(PointData&, PointIDList&)


// points that cannot be solved by intersections - eg inserted traverse

bool ApproximateCoordinates::Solve_insertion()
{
  
  const int max_depth = 100;
  if(selected.empty() || (depth >= max_depth)) return false;

  // building a point list in local coordinate system.  during loop
  // through all observations connecting selected points the involved
  // points are fetched

  PointIDList obs_points;
  bool prv_distance = true;
  bool prv_observation = true;
  ObservationList::iterator first_distance = SM.end();
  ObservationList::iterator first_observation = SM.end();
  {  // VC++ {}
    for(PointIDList::iterator i = selected.begin(); i != selected.end(); i++)
      for(ObservationList::iterator j = SM.begin(); j != SM.end(); j++)
        if(Observation_hasID(j,i))
          // all point IDs are stored and then duplicities are removed
          // - it's faster
          {
            if(Angle *u = dynamic_cast<Angle*>(*j))
              {
                obs_points.push_back((*j)->from());
                obs_points.push_back((*j)->to());
                obs_points.push_back(u->fs());
                if(prv_observation)
                  {
                    first_observation = j;
                    prv_observation = false;
                  }
              }
            else
              {
                if(dynamic_cast<Distance*>(*j) && prv_distance)
                  {
                    first_distance = j;
                    prv_distance = false;
                  }
                else
                  if(prv_observation)
                    {
                      first_observation = j;
                      prv_observation = false;
                    }
                obs_points.push_back((*j)->from());
                obs_points.push_back((*j)->to());
              }
          }
  }  // VC++ {}
  obs_points.sort();
  obs_points.unique();
  /*
   * // removing selected with less the 2 observations - is it correct ???
   * {  // VC++ {}
   * for(PointIDList::iterator i = selected.begin(); i != selected.end(); i++)
   * if(!Necessary_observations(*i))
   * {
   * PointIDList::iterator pozice = find(obs_points.begin(), 
   *                                     obs_points.end(), *i);
   * obs_points.erase(pozice);    
   * }; 
   * }  // VC++ {}
   */
  
  int number_of_given = 0;
  {  // VC++ {}
    for(PointIDList::iterator i = obs_points.begin(); i != obs_points.end(); i++)
      if(Absent(*i))
        number_of_given++;
  }  // VC++ {}

  // not enough given points needed for transformation from local
  // coordinate system (cs)
  if(number_of_given < 2)
    return false;

  PointData local_s;
  {  // VC++ {}
    for(PointIDList::iterator i = obs_points.begin(); 
        i != obs_points.end(); i++)
      local_s[*i] = LocalPoint();
  }  // VC++ {}

  ObservationData OD_local = OD;                // ... deep copy
  /* ObservationList& SM_local(OD_local.OL);   // gamalib-1.1.13 (AC)
   * {  // VC++ {}
   *   for(ObservationList::iterator i = SM.begin(); i != SM.end(); i++)
   *     if(Local_observation(i, obs_points))
   *       SM_local.push_back(*i);  
   * }  // VC++ {}
   */

  // local coordinate system (cs) must be defined now
  const Double pom_Y = 1000;
  const Double pom_X = 5000;
  const Double const_distance = 1000;
  PointID local_cs_1, local_cs_2;
  if(first_distance != SM.end())
    {
      local_cs_1 = (*first_distance)->from();
      local_cs_2 = (*first_distance)->to();
      local_s[local_cs_1].set_xy(pom_X, pom_Y);
      local_s[local_cs_2].set_xy(pom_X, pom_Y + (*first_distance)->value());
    }
  else
    {
      local_cs_1 = (*first_observation)->from();
      local_cs_2 = (*first_observation)->to();
      local_s[local_cs_1].set_xy(pom_X, pom_Y);
      local_s[local_cs_2].set_xy(pom_X, pom_Y + const_distance);
    }
  obs_points.remove(local_cs_1);
  obs_points.remove(local_cs_2);

  // Calculation of coordinates in local coordinate system
  ApproximateCoordinates local_solution(local_s, OD_local, depth + 1);
  local_solution.Calculation();
  if(local_solution.Solved().empty())
    return false;

  // computing in local coordinate system is done, now follows
  // transformation into the current coordinate system
  PointIDList vypoctene_urc;
  PointData::iterator i;
  for(PointIDList::iterator id = selected.begin(); id != selected.end(); id++)
    {
      i = local_s.find(*id);
      if((i != local_s.end()) && (*i).second.test_xy())
        vypoctene_urc.push_back(*id);
    }
  PointData local_s_yx;
  for(PointData::iterator sb_it = local_s.begin(); 
      sb_it != local_s.end(); sb_it++)
    if((*sb_it).second.test_xy())
      local_s_yx[(*sb_it).first] = (*sb_it).second;
  
  SimilarityTr2D transf(SB, local_s_yx, vypoctene_urc);
  transf.Calculation();
  if(transf.State() < unique_solution)
    return false;
  PointData vysledek = transf.Transf_points();
  if(vysledek.empty())
    return false;

  // moving solved points to the point list
  for(PointData::iterator bod = vysledek.begin(); bod != vysledek.end(); bod++)
    {
      if((i = SB.find((*bod).first)) == SB.end())
        SB[(*bod).first] = (*bod).second;
      else
        (*i).second.set_xy((*bod).second.x(), (*bod).second.y());
      solved[(*bod).first] = (*bod).second;
    }

  // removing solved points from "selected" list
  {  // VC++ {}
    for(PointIDList::iterator id = vypoctene_urc.begin(); 
        id != vypoctene_urc.end(); id++)
      selected.remove(*id);
  }  // VC++ {}

  return true;

}  // ApproximateCoordinates::Solve_insertion()


void ApproximateCoordinates::Computational_loop()
{

  PointIDList unsolvable;
  PointIDList::iterator i = selected.begin();
  if (i==selected.end()) return;
  do
    {
      if(!Necessary_observations(*i))
        unsolvable.push_back(*i);
      i++;
    }
  while(i != selected.end());
  for(i = unsolvable.begin(); i != unsolvable.end(); i++)
    {
      PointIDList::iterator pozice = std::find(selected.begin(), 
                                               selected.end(), *i);
      selected.erase(pozice);    
    } 
  
  bool finished = false;
  while(!finished)
    {
      finished = !Solve_intersection(SB, selected);
      // compute insertion
      if(!selected.empty())
        finished = finished && (!Solve_insertion());
    }

  // return unsolvable points back
  selected.insert(selected.end(), unsolvable.begin(), unsolvable.end());
  unsolvable.erase(unsolvable.begin(), unsolvable.end());

}	// ApproximateCoordinates::Computational_loop()


// public


// single point (even if it has coordinates already)

bool ApproximateCoordinates::Calculation(PointID id)
{

  state = calculation_done;
  selected.erase(selected.begin(), selected.end());
  selected.push_back(id);
  if(!Necessary_observations(id) || !Solvable_data(SB))
    return false;
  Solve_intersection(SB, selected);
  Solve_insertion();
  if(!selected.empty())
    {
      Find_missing_coordinates();
      // id already has coordinates
      if(std::find(selected.begin(), selected.end(), id) == selected.end())
        selected.push_back(id);
      PointData SB_puv = SB;
      Computational_loop();
      if(std::find(selected.begin(), selected.end(), id) == selected.end())
        {
          Move_point(SB, SB_puv, id);
          selected.erase(selected.begin(), selected.end());
        }
      else
        {
          selected.erase(selected.begin(), selected.end());
          selected.push_back(id);
        }
      SB = SB_puv;
    }

  return All_is_solved();

}	// bool ApproximateCoordinates::Calculation(PointID id)


// point list (even the points with coordinates)

bool ApproximateCoordinates::Calculation(PointIDList id)
{

  state = calculation_done;
  selected = id;
  // removing duplicities
  selected.sort();
  selected.unique();
  if(selected.empty())
    return true;
  if(!Solvable_data(SB))
    return false;
  PointIDList unsolvable;
  PointIDList::iterator i = selected.begin();
  do
    {
      if(!Necessary_observations(*i))
        {
          unsolvable.push_back(*i);
          i = selected.erase(i);
        }
      else
        i++;
    }
  while(i != selected.end());
  Computational_loop();
  if(!selected.empty())
    {
      PointIDList selected_puv = selected;
      Find_missing_coordinates();
      // id already has coordinates
      for(i = selected_puv.begin(); i != selected_puv.end(); i++)
        if(std::find(selected.begin(), selected.end(), (*i)) == selected.end())
          selected.push_back(*i);
      PointData SB_puv = SB;
      Computational_loop();
      i = selected_puv.begin();
      do
        if(std::find(selected.begin(), selected.end(), (*i)) == selected.end())
          {
            Move_point(SB, SB_puv, *i);
            i = selected_puv.erase(i);
          }
        else
          i++;
      while(i != selected_puv.end());
      selected = selected_puv;
      SB = SB_puv;
    }

  // return unsolvable points back
  selected.insert(selected.end(), unsolvable.begin(), unsolvable.end());

  return All_is_solved();

}	// bool ProblSour::Calculation(PointIDList id)






