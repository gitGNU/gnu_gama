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
 *  $Id: g2d_helper.h,v 1.1 2001/12/07 12:46:44 cepek Exp $
 */

/****************************************************************
 * helper functions and classes:                                *
 * ------------------------------------------------------------ *
 * - typedef std::vector<Point> Helper_list;                    *
 * - inline int signum(const Double& d1)                        *
 * - inline Double g2d_sqr(const Double& d)                     *
 * - enum Solution_state_tag                                    *
 * - enum Observation_types                                     *
 * - inline Observation_types ObservationType(Observation* m)   *
 * - inline Double g2d_distance(const Point& b1, const Bod& b2) *
 * - inline bool g2d_even(std::vector<Double>::size_type& x)    *
 * - class Select_solution_g2d                                  *
 * - class Statistics_g2d                                       *
 * - class Transformation_g2d                                   *
 * - class SimilarityTr2D : public Transformation_g2d           *
 ****************************************************************/
 
#ifndef GaMaLib_g2d_helper_h__GaMaLib_Median_G_fce_H
#define GaMaLib_g2d_helper_h__GaMaLib_Median_G_fce_H

#include <algorithm>
#include <gamalib/local/gamadata.h>
#include <gamalib/local/median/g2d_exception.h>

namespace GaMaLib {

  // --------------------------------------------------------------

  typedef std::vector<Point> Helper_list;
  
  // --------------------------------------------------------------
  // inline int signum(const Double& d1) wass added into
  // gamalib-1.1.61 due to problems with MSC (AC)

  inline int signum(Double d1)
    {
      if(d1 < 0) return -1;
      if(d1 > 0) return  1;
      return 0;
    };
  
  // --------------------------------------------------------------
  
  inline Double g2d_sqr(const Double& d)
    {
      return (d*d);
    };

  // --------------------------------------------------------------

  enum Solution_state_tag
  {
    missing_init = -2,
    calculation_not_done = -1,
    no_solution = 0,
    unique_solution = 1,
    ambiguous_solution = 2,
    calculation_done = 3
  };
  
  // --------------------------------------------------------------

  enum Observation_types
  {
    is_Distance = 1,
    is_Direction = 10,
    is_Angle = 100
  };
  
  // --------------------------------------------------------------

  inline Observation_types ObservationType(Observation* m)
    {
      Distance  *d = dynamic_cast<Distance* >(m);
      Direction *s = dynamic_cast<Direction*>(m);

      if(d) return is_Distance;
      if(s) return is_Direction;
      return is_Angle;
    };

  // --------------------------------------------------------------
  inline Double g2d_distance(const Point& b1, const Point& b2)
    {
      return sqrt(g2d_sqr(b1.x()-b2.x())+g2d_sqr(b1.y()-b2.y()));
    };

  // --------------------------------------------------------------
  inline bool g2d_even(std::vector<Double>::size_type& x)
    {
      return (fmod(x,2) == 0);
    };

  
  // -------------------------------------------------------------- 
  // in the case of ambiguous (equivocal) solution, the one is chosen
  // to be in accordance with the others

  class Select_solution_g2d
    {
    private:
      
      Solution_state_tag state;
      Point   B1, B2;
      PointData* SB;
      ObservationList* SM;
      
    public:    
      Select_solution_g2d(PointData* sb, ObservationList* sm) 
        : state(calculation_not_done), SB(sb), SM(sm) 
        {
        }
      Select_solution_g2d(Point& b1, Point& b2, PointData* sb, 
                          ObservationList* sm) :
        state(calculation_not_done), B1(b1), B2(b2), SB(sb), SM(sm) 
        {
        }
      void Calculation();
      void Calculation(Point b1, Point b2)
        {
          state = calculation_not_done;
          B1 = b1; B2 = b2;
          Calculation();
        }
      Point Solution()
        {
          if(state == calculation_not_done)
            throw g2d_exc("Select_solution_g2d: calculation not done");
          if(state == no_solution)
            throw g2d_exc("Select_solution_g2d: ambiguous solution");
          return B1;
        }
      int State() const 
        { 
          return (state > no_solution ? unique_solution : state); 
        }
      
    };  // class Select_solution_g2d
  
  
  // --------------------------------------------------------------
  // Statistics_g2d - calculation of medina from coordinate list

  class Statistics_g2d
    {
    private:

      Helper_list* PS;
      Point median;
      Solution_state_tag state;
      
    public:
      Statistics_g2d() : state(missing_init) 
        {
        }
      Statistics_g2d(Helper_list* ps) 
        : PS(ps), state(calculation_not_done) 
        {
        }
      void Calculation();
      void Calculation(Helper_list* ps)
        {
          state = calculation_not_done;
          PS = ps;
          Calculation();
        }
      Solution_state_tag State() const 
        { 
          return state; 
        }
      Point Median()       // resulting coordinate
        {
          if(state < no_solution)
            throw g2d_exc("Statistics_g2d: calculation not done");
          
          return median;
        }

    };	// class Statistics_g2d


  // --------------------------------------------------------------
  // base transfromation class

  class Transformation_g2d
    {
    protected:
      
      PointData& SB;              // point list in target syste
      PointData& local;           // point list in local system
      PointIDList& computed;      // list of computed points
      PointData transf_points;    // points transformed into target system
      virtual void Reset() = 0;
      
      Solution_state_tag state;
      
    public:
      Transformation_g2d(PointData& sb, PointData& locl, 
                         PointIDList& comp)
        : SB(sb), local(locl), computed(comp), state(calculation_not_done) 
        {
        }
      virtual ~Transformation_g2d()
        {
          transf_points.erase(transf_points.begin(), transf_points.end());
        }
      void Reset(PointData& sb, PointData& locl, PointIDList& comp)
        {
          SB = sb;
          local = locl;
          computed = comp;
          state = calculation_not_done;
          transf_points.erase(transf_points.begin(), transf_points.end());
          Reset();
        }
      virtual void Calculation() = 0;
      Solution_state_tag State() const 
        { 
          return state; 
        }
      PointData Transf_points() const
        {
          if(state == calculation_not_done)
            throw g2d_exc("Transformation_g2d: calculation not done");
          if(state == no_solution)
            throw 
              g2d_exc("Transformation_g2d: not enough of identical points");
          
          return transf_points;
        }

    };  // class Transformation_g2d

  // --------------------------------------------------------------
  // Transformation in plane
  // Transformation from minimal number of point - selects the best tuple
  // - similarity transformation
  /*
   * transformation key for similarity transformation:
   *
   *  | x |   |tx | |a1  -a2| |x'|
   *  |   | = |   |+|       |*|  |
   *  | y |   |ty | |a2   a1| |y'|
   *
   *  transf_key[0] = a2
   *  transf_key[1] = a1
   *  transf_key[2] = ty
   *  transf_key[3] = tx
   *
   */

  class SimilarityTr2D : public Transformation_g2d
    {
    private:

      std::vector<Double> transf_key;
      void Reset();
      bool Given_point(const PointID& cb)
        {
          return (std::find(computed.begin(), computed.end(), cb) == 
                  computed.end());
        }
      void Identical_points(PointData::iterator& b1, 
                            PointData::iterator& b2);
      void Transformation_key(PointData::iterator& b1, 
                              PointData::iterator& b2);

    public:
      SimilarityTr2D(PointData& sb, PointData& locl, PointIDList& comp)
        : Transformation_g2d(sb, locl, comp) 
        { 
          Reset(); 
        }
      void Calculation();
      std::vector<Double> Transf_key() const
        {
          if(state == calculation_not_done)
            throw g2d_exc("SimilarityTr2D: calculation not done");
          if(state == no_solution)
            throw g2d_exc("SimilarityTr2D: not enough of identical poinst");
          return transf_key;
        }

    };  // class SimilarityTr2D

} // namespace GaMaLib

#endif 








