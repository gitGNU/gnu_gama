/*  
    GNU Gama -- adjustment of geodetic networks
    Copyright (C) 2003  Ales Cepek <cepek@fsv.cvut.cz>

    This file is part of the GNU Gama C++ library.
    
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
 *  $Id: g3_point.h,v 1.17 2003/05/29 16:04:14 cepek Exp $
 */

#include <gamalib/pointid.h>
#include <gnu_gama/g3/g3_parameter.h>

#ifndef GNU_gama__g3_point_h_gnugamag3pointh___gnu_gama_g3point
#define GNU_gama__g3_point_h_gnugamag3pointh___gnu_gama_g3point


namespace GNU_gama {  namespace g3 {

  class Model;
  class Point;
  class Distance;
  class HeightDiff;

  class Parameter_LocalPosition : public Parameter
    {
      Point* pt;

    public:

      Parameter_LocalPosition(Point* p) : pt(p) {}

      Point* point() const { return pt; }
      void   set_point(Point* p) { pt = p; }
    };

  class Parameter_N : 
    public Parameter_LocalPosition,
    public Derivative<Distance>
    {
    public:
      Parameter_N(Point* p) : Parameter_LocalPosition(p) {}

      double derivative(Distance*);
    };
  
  class Parameter_E : 
    public Parameter_LocalPosition,
    public Derivative<Distance>
    {
    public:
      Parameter_E(Point* p) : Parameter_LocalPosition(p) {}

      double derivative(Distance*);
    };
  
  class Parameter_U : 
    public Parameter_LocalPosition,
    public Derivative<Distance>,
    public Derivative<HeightDiff>
    {
    public:
      Parameter_U(Point* p) : Parameter_LocalPosition(p) {}

      double derivative(Distance*);
      double derivative(HeightDiff*);
    };



  class Point {
  public:
  
    typedef GaMaLib::PointID Name;
    typedef Model            Common;

    Name    name;
    Common* common; 

    Point();
    Point(const Point&);
    Point& operator=(const Point&);
    ~Point();

    // -----------------------------

    Parameter  B, L, H, X, Y, Z, height;

    Parameter_N*  N;
    Parameter_E*  E;
    Parameter_U*  U; 

    ParameterList        parlist;


    void set_unused();
    void set_fixed_horizontal_position();
    void set_fixed_height();
    void set_fixed_position();
    void set_free_horizontal_position();
    void set_free_height();
    void set_free_position();
    void set_constr_horizontal_position();
    void set_constr_height();
    void set_constr_position();

    bool unused() const;
    bool fixed_horizontal_position() const;
    bool fixed_height() const;
    bool fixed_position() const;
    bool free_horizontal_position() const;
    bool free_height() const;
    bool free_position() const;
    bool constr_horizontal_position() const;
    bool constr_height() const;
    bool constr_position() const;

    void set_blh(double, double, double);
    void set_xyz(double, double, double);
    void set_height(double);

    bool has_xyz()    const { return has_xyz_;    } 
    bool has_blh()    const { return has_blh_;    } 
    bool has_height() const { return has_height_; }


    double x_transform(double n, double e, double u);
    double y_transform(double n, double e, double u);
    double z_transform(double n, double e, double u);

    void   set_diff_XYZ(double dx, double dy, double dz);
    double diff_N() const;
    double diff_E() const;
    double diff_U() const;

  private:

    enum {
      unused_          = 0,
      fixed_h_pos_     = 1,
      fixed_height_    = 2,
      fixed_position_  = fixed_h_pos_ + fixed_height_,
      free_h_pos_      = 4,
      free_height_     = 8,
      free_position_   = free_h_pos_  + free_height_,
      constr_h_pos_    = 16 + free_position_,
      constr_height_   = 32 + free_height_,
      constr_position_ = constr_h_pos_ + constr_height_,
      h_pos_           = fixed_h_pos_  + free_h_pos_, 
      height_          = fixed_height_ + free_height_,
      position_        = h_pos_ + height_  
    };


    void     transformation_matrix(double b, double l);

    // rotation matrix of transformation from local to global
    // Cartesian coordinates (NEU --> XYZ)

    double   r11, r12, r13,   r21, r22, r23,   r31, r32, r33;

  private:
    double   dX, dY, dZ ; 
    bool     has_xyz_, has_blh_, has_height_;

  };


}}

#endif
