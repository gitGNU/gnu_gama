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
 *  $Id: g2d_cogo.h,v 1.1 2001/12/07 12:46:44 cepek Exp $
 */

/*************************************************************
 * 2d coordinate geometry                                    *
 * --------------------------------------------------------- *
 * - class CoordinateGeometry2D                              *
 * - class Distance_distance : public CoordinateGeometry2D   *
 * - class Direction_direction : public CoordinateGeometry2D *
 * - class Direction_distance : public CoordinateGeometry2D  *
 * - class Direction_angle : public CoordinateGeometry2D     *
 * - class Distance_angle : public CoordinateGeometry2D      *
 * - class Angle_angle : public CoordinateGeometry2D         *
 * - class Circle : public CoordinateGeometry2D              *
 *************************************************************/
 
#ifndef GaMaLib_g2d_cogo_h__GaMaLib_Median_G_ulohy_H
#define GaMaLib_g2d_cogo_h__GaMaLib_Median_G_ulohy_H

#include <gamalib/local/gamadata.h>
#include <gamalib/local/median/g2d_exception.h>

namespace GaMaLib {
  
  class CoordinateGeometry2D
    {
      
    protected:
      int number_of_solutions;
      Point*  point1;
      Point*  point2;
      PointData* SB;
      virtual void Observation_check(Observation*, Observation*) = 0;
      
    public:
      CoordinateGeometry2D(PointData* sb) : number_of_solutions(-1), SB(sb)
        { 
          point1 = new Point; 
          point2 = new Point; 
        }
      virtual ~CoordinateGeometry2D() 
        { 
          delete point1; 
          delete point2; 
        }
      virtual void Calculation() = 0;
      int Number_of_solutions() const 
        { 
          return number_of_solutions; 
        }
      Point Solution_1() const
        {
          if(number_of_solutions == -1)
            throw g2d_exc("CoordinateGeometry2D: calculation not done");
          if(number_of_solutions < 1)
            throw g2d_exc("CoordinateGeometry2D: no solution");
          return *point1;
        }
      Point Solution_2() const
        {
          if(number_of_solutions == -1)
            throw g2d_exc("CoordinateGeometry2D: calculation not done");
          if(number_of_solutions < 2)
            throw g2d_exc("CoordinateGeometry2D: two solutions");
          return *point2;
        }
    };
  
  
  //---------------------------------------------------------------
  
  class Distance_distance : public CoordinateGeometry2D
    {
      
    private:
      Distance* h1;
      Distance* h2;
      PointID   CB;     // for easier search in computed distances
      Double    r1, r2;
      Point     B1;
      Point     B2;
      void Observation_check(Observation*, Observation*);
      
    public:
      Distance_distance() : CoordinateGeometry2D(0), h1(0), h2(0), r1(-1) 
        {
        }
      Distance_distance(Observation* m1, Observation* m2, 
                        PointData* sb, PointID cb)
        : CoordinateGeometry2D(sb), CB(cb), r1(-1)
        {
          Observation_check(m1, m2);
        }
      Distance_distance(Double& m1, Double& m2, 
                        Point b1, Point b2, PointData* sb)
        : CoordinateGeometry2D(sb), r1(m1), r2(m2), B1(b1), B2(b2) 
        {
        }
      ~Distance_distance() 
        {
        }
      void Calculation();
      void New_calculation(Observation* m1, Observation* m2, 
                           PointData* sb, PointID cb)
        {
          point1 = point2 = 0;
          SB = sb;
          CB = cb;
          r1 = -1;
          Observation_check(m1, m2);
          Calculation();
        }
    };
  
  
  //---------------------------------------------------------------

  class Direction_direction : public CoordinateGeometry2D
    {
      
    private:
      Direction* h1;
      Direction* h2;
      void Observation_check(Observation*, Observation*);
      
    public:
      Direction_direction() : CoordinateGeometry2D(0), h1(0), h2(0) 
        {
        }
      Direction_direction(Observation* m1, Observation* m2, PointData* sb)
        : CoordinateGeometry2D(sb)
        {
          Observation_check(m1, m2);
        }
      ~Direction_direction() 
        {
        }
      void Calculation();
      void New_calculation(Observation* m1, Observation* m2, PointData* sb)
        {
          point1 = point2 = 0;
          SB = sb;
          Observation_check(m1, m2);
          Calculation();
        }
    };
  
  
  //---------------------------------------------------------------
  
  class Direction_distance : public CoordinateGeometry2D
    {
      
    private:
      
      Direction*  h1;
      Distance* h2;
      Point  B;
      Double r;
      void Observation_check(Observation*, Observation*);
      
    public:
      
      Direction_distance() : CoordinateGeometry2D(0), h1(0), h2(0), r(-1) 
        {
        }
      Direction_distance(Observation* m1, Observation* m2, PointData* sb) 
        : CoordinateGeometry2D(sb), r(-1)
        {
          Observation_check(m1, m2);
        }
      Direction_distance(Direction* m1, Double m2, Point b, PointData* sb)
        : CoordinateGeometry2D(sb), h1(m1), B(b), r(m2) 
        {
        }
      ~Direction_distance() 
        {
        }
      void Calculation();
      void New_calculation(Observation* m1, Observation* m2, PointData* sb)
        {
          point1 = point2 = 0;
          SB = sb;
          r = -1;
          Observation_check(m1, m2);
          Calculation();
        }
    };


  //---------------------------------------------------------------
  
  class Direction_angle : public CoordinateGeometry2D
    {

    private:
      Direction* h1;
      Angle* h2;
      void Observation_check(Observation*, Observation*);
      
    public:
      Direction_angle() : CoordinateGeometry2D(0), h1(0), h2(0) 
        {
        }
      Direction_angle(Observation* m1, Observation* m2, PointData* sb)
        : CoordinateGeometry2D(sb)
        {
          Observation_check(m1, m2);
        }
      ~Direction_angle() 
        {
        }
      void Calculation();
      void New_calculation(Observation* m1, Observation* m2, PointData* sb)
        {
          point1 = point2 = 0;
          SB = sb;
          Observation_check(m1, m2);
          Calculation();
        }
    };
  
  
  //---------------------------------------------------------------
  
  class Distance_angle : public CoordinateGeometry2D
    {

    private:
      Distance* h1;
      Angle* h2;
      void Observation_check(Observation*, Observation*);

    public:
      Distance_angle() : CoordinateGeometry2D(0), h1(0), h2(0) 
        {
        }
      Distance_angle(Observation* m1, Observation* m2, PointData* sb)
        : CoordinateGeometry2D(sb)
        {
          Observation_check(m1, m2);
        }
      ~Distance_angle() 
        {
        }
      void Calculation();
      void New_calculation(Observation* m1, Observation* m2, PointData* sb)
        {
          point1 = point2 = 0;
          SB = sb;
          Observation_check(m1, m2);
          Calculation();
        }
    };
  
  
  //---------------------------------------------------------------
  
  class Angle_angle : public CoordinateGeometry2D
    {

    private:
      Angle* h1;
      Angle* h2;
      void Observation_check(Observation*, Observation*);

    public:
      Angle_angle() : CoordinateGeometry2D(0), h1(0), h2(0) 
        {
        }
      Angle_angle(Observation* m1, Observation* m2, PointData* sb)
        : CoordinateGeometry2D(sb)
        {
          Observation_check(m1, m2);
        }
      ~Angle_angle() 
        {
        }
      void Calculation();
      void New_calculation(Observation* m1, Observation* m2, PointData* sb)
        {
          point1 = point2 = 0;
          SB = sb;
          Observation_check(m1, m2);
          Calculation();
        }
    };
  
  
  //---------------------------------------------------------------
  // ** circle parameters from two points and perimeter angle **
  
  class Circle : public CoordinateGeometry2D
    {

    private:
      Angle*  h1;
      Point   B1, B2;
      Double  R;
      void Observation_check(Observation*, Observation*) {}

    public:
      Circle() : CoordinateGeometry2D(0), h1(0) 
        {
        }
      Circle(Angle* u,PointData* sb) : CoordinateGeometry2D(sb), h1(u) 
        {
        }
      ~Circle() 
        {
        }
      void Calculation();
      void New_calculation(Angle* u,PointData* sb)
        {
          point1 = point2 = 0;
          SB = sb;
          h1 = u;
          Calculation();
        }
      Double radius() const
        {
          if(number_of_solutions == -1)
            throw g2d_exc("Circle: computation not done");
          if(number_of_solutions < 1)
            throw g2d_exc("Circle: two solutions");
          return R;
        }
    };
  
} // namespace GaMaLib

#endif

