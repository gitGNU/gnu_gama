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
 *  $Id: pointbase.h,v 1.3 2003/03/08 20:39:31 cepek Exp $
 */

#include <map>


#ifndef GNU_gama__point_data_h_gnugamapointbase___gnu_gama_pointbase
#define GNU_gama__point_data_h_gnugamapointbase___gnu_gama_pointbase


namespace GNU_gama {


  template <class Point>
    class PointBase
    {
    private:

      typedef std::map<typename Point::Name, Point*>  Points;
      Points  points;

    public:    

      typename Point::Common* common;
      
      PointBase() : common(0) {}
      PointBase(const PointBase& cod);
      ~PointBase();
      
      PointBase& operator=(const PointBase& cod);

      void put(const Point&);
      void put(Point*);

      Point*       find(const typename Point::Name&);
      const Point* find(const typename Point::Name&) const;
      
      void erase(const typename Point::Name&);
      void erase();
      
      class const_iterator {
      public:

        const_iterator(const typename Points::const_iterator& p) : pit(p) 
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
            ++pit;
            return *this;
          }
        const_iterator operator++(int)
          {
            const_iterator tmp(pit);
            ++pit;
            return tmp;
          }
        const Point* operator*() const
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

        iterator(const typename Points::iterator& p) : pit(p) 
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
            ++pit;
            return *this;
          }
        iterator operator++(int)
          {
            iterator tmp(pit);
            ++pit;
            return tmp;
          }
        Point* operator*() const
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
      erase();
    }

  
  template <class Point>
    PointBase<Point>::PointBase(const PointBase& cpd)
    {
      common = cpd.common;
      for (const_iterator p=cpd.begin(), e=cpd.end(); p!=e; ++p)
        {
          put( **p );
        }
    }


  template <class Point>
    PointBase<Point>& PointBase<Point>::operator=(const PointBase& cpd)
    {
      if (this != &cpd)
        {
          erase();
          common = cpd.common;
          for (const_iterator p=cpd.begin(), e=cpd.end(); p!=e; ++p)
            {
              put( **p );
            }
        }

      return *this;
    }


  template <class Point>
    void PointBase<Point>::put(const Point& point)
    {
      typename Points::iterator t = points.find(point.name);
      if (t != points.end())
        {
          *((*t).second) = point;
        }
      else
        {
          Point* ptr  = new Point(point);
          ptr->common = common;
          points[point.name] = ptr;
        }
    }


  template <class Point>
    void PointBase<Point>::put(Point* point_ptr)
    {
      typename Points::iterator t = points.find(point_ptr->name);
      if (t != points.end())
        {
          Point* p = (*t).second;
          if (p != point_ptr)
            {
              delete p;
            }
          else
            {
              return;
            }
        }

      point_ptr->common = common;
      points[point_ptr->name] = point_ptr;
    }


  template <class Point>
    Point* PointBase<Point>::find(const typename Point::Name& name)
    {
      typename Points::iterator t = points.find(name);
      if (t != points.end())
        {
          return (*t).second;
        }

      return 0;
    }


  template <class Point>
    const Point* PointBase<Point>::find(const typename Point::Name& name) const
    {
      typename Points::const_iterator t = points.find(name);
      if (t != points.end())
        {
          return (*t).second;
        }

      return 0;
    }
  

  template <class Point>
    void PointBase<Point>::erase(const typename Point::Name& name)
    {
      typename Points::iterator t = points.find(name);
      if (t != points.end())
        {
          delete (*t).second;
          points.erase(t);
        }
    }
  
  
  template <class Point>
    void PointBase<Point>::erase()
    {
      typename Points::iterator t = points.begin();
      typename Points::iterator e = points.end();
      while (t != e)
        {
          delete (*t).second;
          ++t;
        }

      points.erase(points.begin(), points.end());
    }  

}

#endif
