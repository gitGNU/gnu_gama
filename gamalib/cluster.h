/*  
    Geodesy and Mapping C++ Library (GNU GaMa / GaMaLib)
    Copyright (C) 2000  Ales Cepek <cepek@fsv.cvut.cz>

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
 *  $Id: cluster.h,v 1.2 2002/09/06 16:08:49 cepek Exp $
 */

#ifndef GaMaLib_Cluster_of_observations__h
#define GaMaLib_Cluster_of_observations__h

#ifdef _MSC_VER
      #define C_CLONE(ptr) Cluster*
#else
      #define C_CLONE(ptr) ptr
#endif

#include <gamalib/observation.h>
#include <gamalib/matvec.h>
#include <vector>
#include <cmath>

namespace GaMaLib {

  class ObservationData;

  class Cluster {
  public:

    const ObservationData*  observation_data;
    ObservationList         observation_list;
    Cov                     covariance_matrix;  

  
    Cluster(const ObservationData* od) : observation_data(od), act_count(0) {}
    virtual ~Cluster();

    virtual Cluster* clone(const ObservationData*) const = 0;

    Double stdDev(int i) const 
      { 
        i++; using namespace std; return sqrt(covariance_matrix(i,i)); 
      }
    int size() const 
      { 
        return observation_list.size(); 
      }

    void update();

    int  activeCount() const { return act_count; }
    Cov  activeCov() const; 

  private:    // no copy ctor and no assignment

    Cluster(const Cluster&); 
    void operator=(const Cluster&);

    friend class Observation;
    int act_count;
  };



  // simple observation types: horizontal directions, distances and
  // angles

  class StandPoint : public Cluster {
  public:

    PointID  station;
    
    StandPoint(const ObservationData* od) 
      : Cluster(od), test_or(false), indx_or(0) 
      {
      }

    C_CLONE(StandPoint*) clone(const ObservationData*p) const 
      {
        return new StandPoint(p); 
      } 

    double orientation() const       
      {
        if (!test_or) throw GaMaLib::Exception(T_POBS_bad_data);
        return attr_or;                  
      }
    void   set_orientation(Double p) { attr_or = p; test_or = true; }
    bool   test_orientation() const  { return test_or;              }
    void   delete_orientation()      { test_or = false;             }
    void   index_orientation(int n)  { indx_or = n;                 }
    int    index_orientation() const { return indx_or;              }

  private:

    Double   attr_or;
    bool     test_or;
    int      indx_or;

  };


  // coordinate observations (observed coordinates) x, y, z

  class Coordinates : public Cluster {
  public:

    Coordinates(const ObservationData* od) : Cluster(od) {}
    C_CLONE(Coordinates*) clone(const ObservationData*p) const 
      { 
        return new Coordinates(p); 
      } 
  };


  // height differences (levelling)

  class HeightDifferences : public Cluster {
  public:

    HeightDifferences(const ObservationData* od) : Cluster(od) {}
    C_CLONE(HeightDifferences*) clone(const ObservationData*p) const 
      { 
        return new HeightDifferences(p); 
      } 
  };


  // vectors (coordinate differences  diff_x, diff_y, diff_z)

  class Vectors : public Cluster {
  public:

    Vectors(const ObservationData* od) : Cluster(od) {}
    C_CLONE(Vectors*) clone(const ObservationData*p) const 
      { 
        return new Vectors(p); 
      } 
  };


}   // namespace GaMaLib

#undef C_CLONE

#endif









