/*  
    GNU Gama -- adjustment of geodetic networks
    Copyright (C) 2003  Ales Cepek <cepek@fsv.cvut.cz>

    This file is part of the GNU Gama library.
    
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
 *  $Id: pointbase.h,v 1.1 2003/03/04 21:46:26 cepek Exp $
 */

#include <map>


#ifndef GNU_gama__point_data_h_gnugamapointbase___gnu_gama_pointbase
#define GNU_gama__point_data_h_gnugamapointbase___gnu_gama_pointbase


namespace GNU_gama {



  // template <class Point> 
  //   class PointList : public std::vector<Point*>
  //   {
  //   };



  template <class Point>
    class PointBase
    {
    private:

      typedef std::map<typename Point::Name, Point*>  Points;
      Points  points;

    public:    
      
      PointBase() {}
      PointBase(const PointBase& cod);
      ~PointBase();
      
      PointBase& operator=(const PointBase& cod);

      void put(const Point&);

      
      class const_iterator {
      public:

        const_iterator(typename Points::const_iterator p) : pit(p) 
          {
          }
        bool operator==(const const_iterator& x) const 
          { 
            return pit==x.pit; 
          }
        bool operator!=(const const_iterator& x) const 
          { 
            return pit!=x.pit; 
          }
        const_iterator& operator++()
          {
            pit++;
            return *this;
          }
        const_iterator operator++(int)
          {
            const_iterator tmp(pit);
            return ++pit;
          }
        const Point* operator*()
          {
            return (*pit).second;
          }

      private:
        typename Points::const_iterator pit;

      };

      const_iterator  begin() const { return points.begin(); }
      const_iterator  end  () const { return points.end  (); }

      class iterator {
      public:

        iterator(typename Points::iterator p) : pit(p) 
          {
          }
        operator const_iterator() const
          {
            return const_iterator(pit);
          }
        bool operator==(const iterator& x) const 
          { 
            return pit==x.pit; 
          }
        bool operator!=(const iterator& x) const 
          { 
            return pit!=x.pit; 
          }
        iterator& operator++()
          {
            pit++;
            return *this;
          }
        iterator operator++(int)
          {
            iterator tmp(pit);
            return ++pit;
          }
        Point* operator*()
          {
            return (*pit).second;
          }

      private:
        typename Points::iterator pit;

      };

      iterator  begin() { return points.begin(); }
      iterator  end  () { return points.end  (); }

      
    };
  

  template <class Point>
    PointBase<Point>::~PointBase()
    {

    }

  
  template <class Point>
    PointBase<Point>::PointBase(const PointBase& cpd)
    {
      if (this != &cpd)
        {

        }

      return *this;
    }


  template <class Point>
    PointBase<Point>& PointBase<Point>::operator=(const PointBase& cpd)
    {
      if (this != &cpd)
        {

        }

      return *this;
    }


  template <class Point>
    void PointBase<Point>::put(const Point& point)
    {
      points[point.name] = new Point(point);
    }


}

#endif
