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
 *  $Id: g3_observation.h,v 1.12 2005/03/27 17:43:26 cepek Exp $
 */


#include <gnu_gama/model.h>
#include <gnu_gama/g3/g3_point.h>
#include <gnu_gama/matvec.h>


#ifndef GNU_gama__g3_observation_h_gnugamag3obs_baseh___gnu_gama_g3obs
#define GNU_gama__g3_observation_h_gnugamag3obs_baseh___gnu_gama_g3obs


namespace GNU_gama {  namespace g3 {


  class Observation :
    public GNU_gama::Observation<Cluster<Observation>, GNU_gama::Cov<> >
  {
  };


  class FromTo {
  public:

    FromTo() : from_dh(0), to_dh(0) {}

    Point::Name  from;
    Point::Name  to;

    double from_dh, to_dh;
  };


  class Value {
    double value;
  public:

    Value() {}
    Value(double d) : value(d) {}

    double obs() const   { return value; }
    void   set(double d) { value = d;    }
  };


  class Distance : public Observation, public FromTo, public Value {
  public:  

    Distance() {}
    Distance(double d) : Value(d) {}

    int dimension() const { return 1; }

    bool revision_accept(ObservationVisitor* visitor)
    {
      if (Revision<Distance>* rv = dynamic_cast<Revision<Distance>*>(visitor))
        {          
          return rv->revision_visit(this);
        }
      else
        return  set_active(false);
    }
  
    void linearization_accept(ObservationVisitor* visitor)
    {
      if (Linearization<Distance>* 
          lv = dynamic_cast<Linearization<Distance>*>(visitor))
        {
          lv->linearization_visit(this);
        }
    }
  };


  class ZenithAngle : public Observation, public FromTo, public Value {
  public:  

    ZenithAngle() {}
    ZenithAngle(double d) : Value(d) {}

    int dimension() const { return 1; }

    bool revision_accept(ObservationVisitor* visitor)
    {
      if (Revision<ZenithAngle>* 
          rv = dynamic_cast<Revision<ZenithAngle>*>(visitor))
        {          
          return rv->revision_visit(this);
        }
      else
        return  set_active(false);
    }
  
    void linearization_accept(ObservationVisitor* visitor)
    {
      if (Linearization<ZenithAngle>* 
          lv = dynamic_cast<Linearization<ZenithAngle>*>(visitor))
        {
          lv->linearization_visit(this);
        }
    }
  };


  class Azimuth : public Observation, public FromTo, public Value {
  public:  

    Azimuth() {}
    Azimuth(double d) : Value(d) {}

    int dimension() const { return 1; }

    bool revision_accept(ObservationVisitor* visitor)
    {
      if (Revision<Azimuth>* 
          rv = dynamic_cast<Revision<Azimuth>*>(visitor))
        {          
          return rv->revision_visit(this);
        }
      else
        return  set_active(false);
    }
  
    void linearization_accept(ObservationVisitor* visitor)
    {
      if (Linearization<Azimuth>* 
          lv = dynamic_cast<Linearization<Azimuth>*>(visitor))
        {
          lv->linearization_visit(this);
        }
    }
  };


  class Vector : public Observation, public FromTo {
  public:  
        
    Vector()
    {
    }
    Vector(double x, double y, double z)
    {
      dx_ = x; dy_ = y; dz_ = z;
    }
      
    int  dimension() const { return 3; }
    void set_dxyz(double x, double y, double z)
    {
      dx_ = x; dy_ = y; dz_ = z;
    }
    
    double dx() const { return dx_; }
    double dy() const { return dy_; }
    double dz() const { return dz_; }
    
    bool revision_accept(ObservationVisitor* visitor)
    {
      if (Revision<Vector>* rv = dynamic_cast<Revision<Vector>*>(visitor))
        {
          return rv->revision_visit(this);
        }
      else
        return  set_active(false);
    }
  
    void linearization_accept(ObservationVisitor* visitor)
    {
      if (Linearization<Vector>* 
          lv = dynamic_cast<Linearization<Vector>*>(visitor))
        {
          lv->linearization_visit(this);
        }
    }
    
  private:    
    double dx_, dy_, dz_;
  };


  class XYZ : public Observation {
  public:  

    Point::Name id;
        
    XYZ()
    {
    }
    XYZ(double xp, double yp, double zp)
    {
      x_ = xp; y_ = yp; z_ = zp;
    }
      
    int  dimension() const { return 3; }
    void set_xyz(double xp, double yp, double zp)
    {
      x_ = xp; y_ = yp; z_ = zp;
    }
    
    double x() const { return x_; }
    double y() const { return y_; }
    double z() const { return z_; }
    
    bool revision_accept(ObservationVisitor* visitor)
    {
      if (Revision<XYZ>* rv = dynamic_cast<Revision<XYZ>*>(visitor))
        {
          return rv->revision_visit(this);
        }
      else
        return set_active(false);
    }
  
    void linearization_accept(ObservationVisitor* visitor)
    {
      if (Linearization<XYZ>* 
          lv = dynamic_cast<Linearization<XYZ>*>(visitor))
        {
          lv->linearization_visit(this);
        }
    }
    
  private:    
    double x_, y_, z_;
  };

}}

#endif



