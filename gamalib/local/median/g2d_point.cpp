/*  
    Geodesy and Mapping C++ Library (GNU GaMa / GaMaLib)
    Copyright (C) 1999  Jiri Vesely <vesely@gama.fsv.cvut.cz>
                  2001  Ales Cepek  <cepek@gama.fsv.cvut.cz>

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
 *  $Id: g2d_point.cpp,v 1.2 2002/05/24 19:30:51 cepek Exp $
 */

/*************************************************************
 * approximate coordinates of a single point                 *
 *************************************************************/


#include <gamalib/local/median/g2d_point.h>
#include <gamalib/local/median/g2d_cogo.h>
#include <gamalib/local/orientation.h>

namespace GaMaLib {

  // private

  void ApproxPoint::ClearLists()
  {
    for(ObservationList::iterator i = SM.begin(); i < SM.end(); i++)
      delete(*i);
    SM.erase(SM.begin(),SM.end());
    Solved_points.erase(Solved_points.begin(),Solved_points.end());
  }

  
  void ApproxPoint::ArangeObservations(ObservationList& sm_pom)
  {
    ObservationList::iterator i, j;
    // distances
    Double med;
    std::vector<Double> pom_sez; 
    i = sm_pom.begin();
    while(i < sm_pom.end())
      {
        if(Distance* d1 = dynamic_cast<Distance*>(*i))
          {
            pom_sez.push_back(d1->value());
            j = i+1;
            while(j < sm_pom.end())
              {
                Distance* d2 = dynamic_cast<Distance*>(*j);
                if(d2 && ((d1->from() == d2->from() && 
                           d1->to() == d2->to()) ||
                          (d1->to() == d2->from() && 
                           d1->from() == d2->to())))
                  {
                    pom_sez.push_back(d2->value());
                    sm_pom.erase(j);
                  }
                else
                  j++;
              };
            std::sort(pom_sez.begin(),pom_sez.end());
            std::vector<Double>::size_type size = pom_sez.size();
            med = (g2d_even(size) ? (pom_sez[size/2-1] + pom_sez[size/2])/2 : 
                   pom_sez[(size+1)/2-1]);
            Distance* DD 
              = new Distance(CB, (d1->from()==CB?d1->to():d1->from()), med);
            SM.push_back(DD);
            sm_pom.erase(i);
            pom_sez.erase(pom_sez.begin(), pom_sez.end());
          }
        else
          i++;
      };
    
    // outer bearings
    i = SM_S.begin();
    while(i < SM_S.end())
      {
        Direction* s1 = static_cast<Direction*>(*i);
        pom_sez.push_back(s1->value());
        j = i+1;
        while(j < SM_S.end())
          {
            Direction* s2 = static_cast<Direction*>(*j);
            if(s1->from() == s2->from())
              {
                pom_sez.push_back(s2->value());
                SM_S.erase(j);
                delete(s2);
              }
            else
              j++;
          };
        std::sort(pom_sez.begin(),pom_sez.end());
        std::vector<Double>::size_type size = pom_sez.size();
        med = (g2d_even(size) ? (pom_sez[size/2-1] + pom_sez[size/2])/2 : 
               pom_sez[(size+1)/2-1]);
        Direction* SS = new Direction(s1->from(),s1->to(),med);
        SM.push_back(SS);
        SM_S.erase(i);
        delete(s1);
        pom_sez.erase(pom_sez.begin(), pom_sez.end());
      };

    // inner angles
    // in sm_pom are only inner angels now; dists. have been already removed 
    // in SM_U are only angles
    i = sm_pom.begin();
    Double u_mer;
    while(i < sm_pom.end())
      {
        Angle* u1 = static_cast<Angle*>(*i);
        j = SM_U.begin();
        pom_sez.push_back(u1->value());
        while(j < SM_U.end())
          {
            Angle* u2 = static_cast<Angle*>(*j);
            if(((u1->to()==u2->to())&&(u1->fs()==u2->fs()))||
               ((u1->to()==u2->fs())&&(u1->fs()==u2->to())))
              {
                u_mer = (u1->to() == u2->to() ? u2->value() : 2*M_PI-u2->value());
                pom_sez.push_back(u_mer);
                SM_U.erase(j);
                delete(u2);
              }
            else
              j++;
          };
        std::sort(pom_sez.begin(),pom_sez.end());
        std::vector<Double>::size_type size = pom_sez.size();
        med = (g2d_even(size) ? (pom_sez[size/2-1] + pom_sez[size/2])/2 : 
               pom_sez[(size+1)/2-1]);
        Angle* UU;
        if(med >= M_PI)
          UU = new Angle(u1->from(),u1->fs(),u1->to(),
                         M_PI*2-med);
        else
          UU = new Angle(u1->from(),u1->to(),u1->fs(),med);
        SM.push_back(UU);
        sm_pom.erase(i);
        pom_sez.erase(pom_sez.begin(), pom_sez.end());
      };
    i = SM_U.begin();
    // finishing remaining angles
    while(i < SM_U.end())
      {
        Angle* u1 = static_cast<Angle*>(*i);
        j = i+1;
        pom_sez.push_back(u1->value());
        while(j < SM_U.end())
          {
            Angle* u2 = static_cast<Angle*>(*j);
            if(((u1->to()==u2->to())&&(u1->fs()==u2->fs()))||
               ((u1->to()==u2->fs())&&(u1->fs()==u2->to())))
              {
                u_mer = (u1->to() == u2->to() ? u2->value() : 2*M_PI-u2->value());
                pom_sez.push_back(u_mer);
                SM_U.erase(j);
                delete(u2);
              }
            else
              j++;
          };
        std::sort(pom_sez.begin(),pom_sez.end());
        std::vector<Double>::size_type size = pom_sez.size();
        med = (g2d_even(size) ? (pom_sez[size/2-1] + pom_sez[size/2])/2 : 
               pom_sez[(size+1)/2-1]);
        Angle* UU;
        if(med >= M_PI)
          UU = new Angle(u1->from(),u1->fs(),u1->to(),
                         M_PI*2-med);
        else
          UU = new Angle(u1->from(),u1->to(),u1->fs(),med);
        SM.push_back(UU);
        SM_U.erase(i);
        delete(u1);
        pom_sez.erase(pom_sez.begin(), pom_sez.end());
      };
  };  // void ApproxPoint::ArangeObservations()


  void ApproxPoint::Reset(PointData* sb, ObservationList* sm, 
                          const PointID& cb)
  {
    state = calculation_not_done;
    ClearLists();  // clearin lists (if a distance is used more then once)
    SB = *sb;
    CB = cb;
    ObservationList sm_s;
    ObservationList sm_pom;

    // computing orientation shift again - solved points are considered as well
    Orientation ors(SB,*sm);
    ors.add_all();

    // selecting observations related to the computed point
    for(ObservationList::const_iterator i = sm->begin(); i < sm->end(); i++)
      {
        if(((*i)->from() == CB) && KnownTarget(i))
          if(Direction* s = dynamic_cast<Direction*>(*i))
            sm_s.push_back(s);
          else
            sm_pom.push_back(*i);
        else
          if(KnownStandpoint(i))
            if(Angle* u = dynamic_cast<Angle*>(*i))
              {
                if((u->to() == CB && KnownTarget2(u)) ||
                   (u->fs() == CB && KnownTarget1(i)))
                  SM_S.push_back(MakeBearing(u,CB));
              }
            else
              if((*i)->to() == CB)
                if(Direction* s = dynamic_cast<Direction*>(*i))
                  SM_S.push_back(MakeBearing(s,CB));
                else
                  sm_pom.push_back(*i);
      };

    // transforming directions on the computed standpoint to inner angels
    { // VC++ {} ...... here and elsewhere curly braces are added to
      // enable processing of standard C++ code with MS compiler
      for(ObservationList::iterator i = sm_s.begin(); i < sm_s.end(); i++)
        for(ObservationList::iterator j = i+1; j < sm_s.end(); j++)
          if ((*i)->ptr_cluster() == (*j)->ptr_cluster())
            SM_U.push_back(MakeAngle(i,j));
      sm_s.erase(sm_s.begin(), sm_s.end());
    }  // VC++ {}

    // now putting selected observations in good order - both distances
    // from and to the stanpoint, identical angles, ... etc.
    ArangeObservations(sm_pom);

#ifdef PB_Debug
    std::cout << "\n***** reset done *****\nfor poind ID " << CB << '\n';
    std::cout << "***** directions on the standpoint *****\n" << sm_s << '\n';
    std::cout << "***** inner angles from directions  *****\n" << SM_U << '\n';
    std::cout << "***** outer bearing *****\n"  << SM_S << '\n';
    std::cout << "***** other observations *****\n" << sm_pom << '\n';
    std::cout << "+++++ observations with a point +++++\n" << SM << '\n';
#endif
  };  /* ApproxPoint::reset(PointData& sb, const ObservationList& sm, 
         const PointID& cb) */


  // public

  void ApproxPoint::Calculation()
  {
    CoordinateGeometry2D* GU;
    Select_solution_g2d* VR = new Select_solution_g2d(&SB,&SM);
    bool two_solutions = false;
    Point prv, dru;
    for(ObservationList::const_iterator i = SM.begin(); i != SM.end(); i++)
      for(ObservationList::const_iterator j = i+1; j != SM.end(); j++)
        {
          switch (ObservationType(*i) + ObservationType(*j))
            {
            case 2*is_Distance :
              {
                Distance_distance* V = new Distance_distance(*i, *j, &SB, CB);
                V->Calculation();
                GU = V;
              }
              break;
            case is_Distance + is_Direction :
              {
                Direction_distance* V = new Direction_distance(*i, *j, &SB);
                V->Calculation();
                GU = V;
              }
              break;
            case is_Distance + is_Angle :
              {
                Distance_angle* V = new Distance_angle(*i, *j, &SB);
                V->Calculation();
                GU = V;
              }
              break;
            case 2*is_Direction :
              {
                Direction_direction* V = new Direction_direction(*i, *j, &SB);
                V->Calculation();
                GU = V;
              }
              break;
            case is_Direction + is_Angle :
              {
                Direction_angle* V = new Direction_angle(*i, *j, &SB);
                V->Calculation();
                GU = V;
              }
              break;
            case 2*is_Angle :
              {
                Angle_angle* V = new Angle_angle(*i, *j, &SB);
                V->Calculation();
                GU = V;
              }
              break;
            };
          switch (GU->Number_of_solutions())
            {
            case 2 :
              {
                VR->Calculation(GU->Solution_1(),GU->Solution_2());
                if(VR->State() == 1)         // unique solution selection
                  {
                    Solved_points.push_back(VR->Solution());
#ifdef PB_Debug
                    std::cout.precision(10);
                    std::cout << "-> " << VR->Solution().y() << ' ' 
                              << VR->Solution().x() << '\n';
#endif
                  };
                if(Solved_points.empty() && 
                   (VR -> State() == 0)  &&  (!two_solutions))
                  {
                    prv = GU->Solution_1();
                    dru = GU->Solution_2();
                    two_solutions = true;
                  }
              }
              break;
            case 1 :
              {
                Solved_points.push_back(GU->Solution_1());
#ifdef PB_Debug
                std::cout << "-> " << GU->Solution_1().y() << ' ' 
                          << GU->Solution_1().x() << '\n';
#endif
              }
              break;
            default : // void solution remains - nothing else could be done
              break;
            };
          delete GU;
        };
    // making final coordinates
    if(Solved_points.size() > 0)
      {
        state = unique_solution;
        Statistics_g2d* ST = new Statistics_g2d(&Solved_points);
        ST->Calculation();
        v_point = ST->Median();
#ifdef PB_Debug
        std::cout.precision(10);
        std::cout << '*' << ST->Median().y() << ' ' 
                  << ST->Median().x() << '\n';
#endif
      }
    else
      state = no_solution;
    if(two_solutions && Solved_points.empty())
      {
        v_point = prv;
        v_point2 = dru;
        state = ambiguous_solution;   // two solutions
      }

#ifdef PB_Debug
    for(std::vector<Point>::const_iterator i = Solved_points.begin(); 
        i < Solved_points.end(); i++)
      std::cout << "-> " << i->y() << ' ' << i->x() << '\n';
#endif
  };      // void ApproxPoint::Calculation()

}       // namespace GaMaLib



