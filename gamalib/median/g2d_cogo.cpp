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
 *  $Id: g2d_cogo.cpp,v 1.1 2001/12/07 12:46:44 cepek Exp $
 */

/**************************************************************
 * 2d coordinate geometry                                     *
 **************************************************************/
 
/*
 * In Median in all methods "Calculation" was added try / catch
 * block. Exceptions g2d_exc are sent by command throw for next
 * processing, for all other exceptions is returned state =
 * no_solution or number_of_solutions = -1. This change was motivated
 * by a bug, when in Median computaion of bearing for two identical
 * points failed.
 *
 * AC 2000.04.26 change median-0.7.5 / gamalib-0.9.57 
 * AC 2001.04.20 gamalib-1.1.61 */
 
#include <gamalib/local/median/g2d_cogo.h>
#include <gamalib/local/median/g2d_exception.h>
#include <gamalib/local/median/g2d_helper.h>
#include <gamalib/local/pobs/bearing.h>


namespace GaMaLib {
  

  // ************** Distance_distance *********************

  void Distance_distance::Observation_check(Observation* m1, Observation* m2)
  {
    h1 = dynamic_cast<Distance*>(m1);
    h2 = dynamic_cast<Distance*>(m2);
    if(!(h1 && h2))
      throw g2d_exc("Distance_distance: wrong observation type");
  };  

  // ** computation: CHARAMZA, pp.127-130 **

  void Distance_distance::Calculation()
  {
    try {
#ifdef PB_Debug
      std::cout << "computing distance_distance....\n ";
#endif

      number_of_solutions = 0;             // -1 when computation not done
      if(r1 == -1)
        {
          PointID CB1 = h1->to();
          PointID CB2 = h2->to();
          B1 = (*(SB->find(CB1))).second;
          B2 = (*(SB->find(CB2))).second;
          r1 = h1->value();
          r2 = h2->value();
        };
      Double dy = B2.y() - B1.y();
      Double dx = B2.x() - B1.x();
      Double s12 = sqrt(g2d_sqr(dy)+g2d_sqr(dx));
      if(s12 == 0)                  // given identical points; no solution
        return;
      Double s1 = r1 / s12;
      Double s2 = r2 / s12;
      Double f = ((s1+s2)*(s1-s2)+1)/2;
      Double g = (s1+f)*(s1-f);
      if(g < 0)                     // intersection doesn't exist
        return;
      if(sqrt(g) < (0.15*s1*s2))    // intersection angle < 10 gon
        return;
      point1->set_xy(B1.x()+dx*f-dy*sqrt(g), B1.y()+dy*f+dx*sqrt(g));
      number_of_solutions = 1;
      if(g > 0)
        {
          point2->set_xy(B1.x()+dx*f+dy*sqrt(g), B1.y()+dy*f-dx*sqrt(g));
          number_of_solutions = 2;
        };
      return;

    }
    catch (g2d_exc& exc) 
      {
        throw exc;
      }
    catch (...) 
      {
        number_of_solutions = -1;
        return;
      }
  
  };  // void Distance_distance::Calculation


  //----------------------------------------------------------------
  // ************** Direction_direction *********************

  void Direction_direction::Observation_check(Observation* m1, Observation* m2)
  {
    h1 = dynamic_cast<Direction*>(m1);
    h2 = dynamic_cast<Direction*>(m2);
    if(!(h1 && h2))
      throw g2d_exc("Direction_direction: wrong observation type");
  };


  // ** computation: CHARAMZA, pp.131-133 **

  void Direction_direction::Calculation()
  {
    try {
#ifdef PB_Debug
      std::cout << "computing direction_direction...." << h1->from() 
                << " and " << h2->from() << '\n';
#endif

      number_of_solutions = 0;           // -1 when computation not done
      if(fabs(sin(h2->value())) < fabs(sin(h1->value())))
        {
          Direction* h = h1;
          h1 = h2;
          h2 = h;
        };
      const Point B1 = (*(SB->find(h1->from()))).second;
      const Point B2 = (*(SB->find(h2->from()))).second;
      Double jmen = cos(h1->value())*sin(h2->value())
        -sin(h1->value())*cos(h2->value());
      if(fabs(jmen) < 0.15)       // unreliable intersection
        return;
      Double dy = (sin(h1->value())*sin(h2->value())*(B2.x()-B1.x())-
                   cos(h1->value())*sin(h2->value())*(B2.y()-B1.y()))/jmen;
      /*
       * if((signum(dy) != signum(sin(h2->value()))) ||
       *   (signum(B2.x()+(dy*cos(h2->value()))/sin(h2->value())-B1.x()) != 
       *    signum(cos(h1->value()))))
       *  return; //intersection doesn't exist; intersection in the given point
       *
       * AC gamalib-1.1.61 : problems with MSC
       */
      const int s_1=signum(dy);
      const int s_2=signum(sin(h2->value()));
      const int s_3=signum(B2.x()+(dy*cos(h2->value()))
                           /sin(h2->value())-B1.x());
      const int s_4=signum(cos(h1->value()));
      if ((s_1 != s_2) || (s_3 != s_4)) return;

      point1->set_xy(B2.x()+(dy*cos(h2->value()))/sin(h2->value()), B2.y()+dy);
      number_of_solutions = 1;
      return;

    }
    catch (g2d_exc& exc) 
      {
        throw exc;
      }
    catch (...) 
      {
        number_of_solutions = -1;
        return;
      }
  
  };  // void Direction_direction::Calculation


  // -----------------------------------------------------------------
  // ************** Direction_distance *********************

  void Direction_distance::Observation_check(Observation* m1, Observation* m2)
  {
    if((h1 = dynamic_cast<Direction*>(m1)) != 0)
      h2 = dynamic_cast<Distance*>(m2);
    else
      {
        h1 = dynamic_cast<Direction *>(m2);
        h2 = dynamic_cast<Distance*>(m1);
      };
    if(!(h1 && h2))
      throw g2d_exc("Direction_distance: wrong observation type");
  }; 


  // ** computation: CHARAMZA, pp.115-118 **
  void Direction_distance::Calculation()
  {
    try {
#ifdef PB_Debug
      std::cout << "computing direction_distance.... \n";
#endif

      number_of_solutions = 0;    // -1 when computation not done
      const Point B1 = (*(SB->find(h1->from()))).second;
      Point B2;
      if(r == -1)                 // input Direction*, Distance*
        {
          if(SB->find(h1->to()) == SB->find(h2->to()))
            B2 = (*(SB->find(h2->from()))).second;
          else
            B2 = (*(SB->find(h2->to()))).second;
          r = h2->value();
        }
      else                        // input Direction*, Double, Point
        B2 = B;
      if(r <= 0)                  // radius <= 0
        return;
      Double yp = (B1.y()-B2.y())*cos(h1->value())
        -(B1.x()-B2.x())*sin(h1->value());
      if(fabs(yp) > r) // semi-line outside circle; intersection doesn't exist
        return;
      Double xp = (B1.x()-B2.x())*cos(h1->value())
        +(B1.y()-B2.y())*sin(h1->value());
      Double x1 = sqrt(g2d_sqr(r)-g2d_sqr(yp));
      if(x1 <= xp)     // semi-line outside circle; intersection doesn't exist
        return;
      if(x1 < (0.15*r))      // intersection angle < 10 gon
        return;
      point1->set_xy(B2.x()+x1*cos(h1->value())-yp*sin(h1->value()),
                     B2.y()+yp*cos(h1->value())+x1*sin(h1->value()));
      number_of_solutions = 1;
      if((-x1) <= (xp+1e-6)) // only one solution; 1e-6 ==> roundoff
        return;
      point2->set_xy(B2.x()-x1*cos(h1->value())-yp*sin(h1->value()),
                     B2.y()+yp*cos(h1->value())-x1*sin(h1->value()));
      number_of_solutions = 2;
      return;

    }
    catch (g2d_exc& exc) 
      {
        throw exc;
      }
    catch (...) 
      {
        number_of_solutions = -1;
        return;
      }
  
  };  // void Direction_distance::Calculation()


  //----------------------------------------------------------------
  // ************** Direction_angle *********************

  void Direction_angle::Observation_check(Observation* m1, Observation* m2)
  {
    if((h1 = dynamic_cast<Direction*>(m1)) != 0)
      h2 = dynamic_cast<Angle*>(m2);
    else
      {
        h1 = dynamic_cast<Direction*>(m2);
        h2 = dynamic_cast<Angle*>(m1);
      };
    if(!(h1 && h2))
      throw g2d_exc("Direction_angle: wrong observation type");
  };


  // ** computation: CHARAMZA, p.142 **

  void Direction_angle::Calculation()
  {
    try {
#ifdef PB_Debug
      std::cout << "computing direction_angle....\n";
#endif

      number_of_solutions = 0;           // -1 when computation not done
      Circle K(h2,SB);
      K.Calculation();
      if(K.Number_of_solutions() < 1)    // circle parameters not solved
        return;
      Direction_distance SD(h1,K.radius(),K.Solution_1(),SB);
      SD.Calculation();
      if(SD.Number_of_solutions() < 1)   // no intersection exist
        return;
      Point B1 = (*(SB->find(h2->to()))).second;
      Point B2 = (*(SB->find(h2->rs()))).second;
      Double uu = bearing(SD.Solution_1(),B2) - bearing(SD.Solution_1(),B1);
      uu += (uu < 0 ? 2*M_PI : 0);
      // uu should be equal to h2->value(), but ...  uu is either
      // value() or value()+-PI
      if((uu < (h2->value()+0.1)) && (uu > (h2->value()-0.1)))
        {
          point1->set_xy(SD.Solution_1().x(), SD.Solution_1().y());
          number_of_solutions = 1;
        };
      if(SD.Number_of_solutions() > 1)
        {
          uu = bearing(SD.Solution_2(),B2) - bearing(SD.Solution_2(),B1);
          uu += (uu < 0 ? 2*M_PI : 0);
          if((uu < (h2->value()+0.1)) && (uu > (h2->value()-0.1)))
            if(number_of_solutions == 1)
              {
                point2->set_xy(SD.Solution_2().x(), SD.Solution_2().y());
                number_of_solutions = 2;
              }
            else
              {
                point1->set_xy(SD.Solution_2().x(), SD.Solution_2().y());
                number_of_solutions = 1;
              };
        };
      return;

    }
    catch (g2d_exc& exc) 
      {
        throw exc;
      }
    catch (...) 
      {
        number_of_solutions = -1;
        return;
      }
  
  };  // void Direction_angle::Calculation()


  //----------------------------------------------------------------
  // ************** Distance_angle *********************

  void Distance_angle::Observation_check(Observation* m1, Observation* m2)
  {
    if((h1 = dynamic_cast<Distance*>(m1)) != 0)
      h2 = dynamic_cast<Angle*>(m2);
    else
      {
        h1 = dynamic_cast<Distance*>(m2);
        h2 = dynamic_cast<Angle *>(m1);
      };
    if(!(h1 && h2))
      throw g2d_exc("Distance_angle: wrong observation type");
  };


  // ** computiong: CHARAMZA, p.143 **

  void Distance_angle::Calculation()
  {
    try {
#ifdef PB_Debug
      std::cout << "computing distance_angle....\n";
#endif

      number_of_solutions = 0;           // -1 when computation not done
      Circle K(h2,SB);
      K.Calculation();
      if(K.Number_of_solutions() < 1)    // computation of circle falied
        return;
      PointID CBB = (h1->from() == h2->from() ? h1->to() : h1->from());
      Point  BB = (*(SB->find(CBB))).second;
      Double dd1 = h1->value();
      Double dd2 = K.radius();
      Distance_distance DD(dd1,dd2,BB,K.Solution_1(),SB);
      DD.Calculation();
      if(DD.Number_of_solutions() < 1)   // intersection doesn't exist
        return;
      Point B1 = (*(SB->find(h2->to()))).second;
      Point B2 = (*(SB->find(h2->rs()))).second;
      Double uu = bearing(DD.Solution_1(),B2) - bearing(DD.Solution_1(),B1);
      uu += (uu < 0 ? 2*M_PI : 0);
      // uu should be equalto h2->value(), but ...  uu is either
      // value() or value()+-PI
      if((uu < (h2->value()+0.1)) && (uu > (h2->value()-0.1)))
        {
          point1->set_xy(DD.Solution_1().x(), DD.Solution_1().y());
          number_of_solutions = 1;
        };
      if(DD.Number_of_solutions() > 1)
        {
          uu = bearing(DD.Solution_2(),B2) - bearing(DD.Solution_2(),B1);
          uu += (uu < 0 ? 2*M_PI : 0);
          if((uu < (h2->value()+0.1)) && (uu > (h2->value()-0.1)))
            if(number_of_solutions == 1)
              {
                point2->set_xy(DD.Solution_2().x(), DD.Solution_2().y());
                number_of_solutions = 2;
              }
            else
              {
                point1->set_xy(DD.Solution_2().x(), DD.Solution_2().y());
                number_of_solutions = 1;
              };
        };
      return;  

    }
    catch (g2d_exc& exc) 
      {
        throw exc;
      }
    catch (...) 
      {
        number_of_solutions = -1;
        return;
      }
  
  };  // void Distance_angle::Calculation()


  //-----------------------------------------------------------
  // ************** Angle_angle *********************

  void Angle_angle::Observation_check(Observation* m1, Observation* m2)
  {
    h1 = dynamic_cast<Angle*>(m1);
    h2 = dynamic_cast<Angle*>(m2);
    if(!(h1 && h2))
      throw g2d_exc("Angle_angle: wrong observation type");
  };


  // ** computation: CHARAMZA, p.143 **

  void Angle_angle::Calculation()
  {
    try {
#ifdef PB_Debug
      std::cout << "computing angle_angle....\n";
#endif

      number_of_solutions = 0;           // -1 when computation not done
      Circle K1(h1,SB);
      K1.Calculation();
      Circle K2(h2,SB);
      K2.Calculation();
      // circle parameters were not solved ?
      if((K1.Number_of_solutions() < 1) || (K2.Number_of_solutions() < 1))    
        return;
      Double rr1 = K1.radius();
      Double rr2 = K2.radius();
      Distance_distance DD(rr1,rr2,K1.Solution_1(),K2.Solution_1(),SB);
      DD.Calculation();
      if(DD.Number_of_solutions() < 1)   // intersection doesn't exist
        return;
      Point B1 = (*(SB->find(h1->to()))).second;
      Point B2 = (*(SB->find(h1->rs()))).second;
      Point B3 = (*(SB->find(h2->to()))).second;
      Point B4 = (*(SB->find(h2->rs()))).second;
      bool Vyhovuje1, Vyhovuje2;
      Double uu1, uu2;
      // in the case of common point at both angles is one of
      // intersections this point
      if(!(((B1.x()==DD.Solution_1().x()) && (B1.y()==DD.Solution_1().y())) ||
           ((B2.x()==DD.Solution_1().x()) && (B2.y()==DD.Solution_1().y()))))
        {
          uu1 = bearing(DD.Solution_1(),B2) - bearing(DD.Solution_1(),B1);
          uu2 = bearing(DD.Solution_1(),B4) - bearing(DD.Solution_1(),B3);
          uu1 += (uu1 < 0 ? 2*M_PI : 0);
          uu2 += (uu2 < 0 ? 2*M_PI : 0);
          // uu should be equal to h2->value(), but ...  uu is either
          // value() or value()+-PI
          Vyhovuje1 = (uu1 < (h1->value()+0.1)) && (uu1 > (h1->value()-0.1));
          Vyhovuje2 = (uu2 < (h2->value()+0.1)) && (uu2 > (h2->value()-0.1));
          if(Vyhovuje1 && Vyhovuje2)
            {
              point1->set_xy(DD.Solution_1().x(), DD.Solution_1().y());
              number_of_solutions = 1;
            };
        };
      if(DD.Number_of_solutions() > 1)
        if(!(((B1.x()==DD.Solution_2().x()) && 
              (B1.y()==DD.Solution_2().y())) ||
             ((B2.x()==DD.Solution_2().x()) && 
              (B2.y()==DD.Solution_2().y()))))
          {
            uu1 = bearing(DD.Solution_2(),B2) - bearing(DD.Solution_2(),B1);
            uu2 = bearing(DD.Solution_2(),B4) - bearing(DD.Solution_2(),B3);
            uu1 += (uu1 < 0 ? 2*M_PI : 0);
            uu2 += (uu2 < 0 ? 2*M_PI : 0);
            Vyhovuje1 = (uu1 < (h1->value()+0.1)) && (uu1 > (h1->value()-0.1));
            Vyhovuje2 = (uu2 < (h2->value()+0.1)) && (uu2 > (h2->value()-0.1));
            if(Vyhovuje1 && Vyhovuje2)
              if(number_of_solutions == 1)
                {
                  point2->set_xy(DD.Solution_2().x(), DD.Solution_2().y());
                  number_of_solutions = 2;
                }
              else
                {
                  point1->set_xy(DD.Solution_2().x(), DD.Solution_2().y());
                  number_of_solutions = 1;
                };
          };
      return;

    }
    catch (g2d_exc& exc) 
      {
        throw exc;
      }
    catch (...) 
      {
        number_of_solutions = -1;
        return;
      }
  
  };  // void Angle_angle::Calculation()


  //------------------------------------------------------------
  // ************** Circle *********************

  // ** computation: CHARAMZA, pp.107-109 **

  void Circle::Calculation()
  {
    try {
#ifdef PB_Debug
      std::cout << "computing parameters of circle....\n";
#endif

      number_of_solutions = 0;       // -1 when computation not done
      Double u = h1->value();	     // just to spare typing h1->value()
      if(fabs(sin(u)) < 0.15)        // small angle
        return;
      B1 = (*(SB->find(h1->to()))).second;
      B2 = (*(SB->find(h1->rs()))).second;
      Double sm, d;
      bearing_distance(B1,B2,sm,d);
      if(d == 0)                     // identical points
        return;
      Double rr = d/sin(u)/2;
      R = fabs(rr);
      point1->set_xy(B1.x()-rr*sin(sm - u), B1.y()+rr*cos(sm - u));
      number_of_solutions = 1;
      return;

    }
    catch (g2d_exc& exc) 
      {
        throw exc;
      }
    catch (...) 
      {
        number_of_solutions = -1;
        return;
      }
  
  };  // void Circle::Calculation()

} // namespace GaMaLib



