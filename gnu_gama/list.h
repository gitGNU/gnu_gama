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
 *  $Id: list.h,v 1.7 2004/01/25 11:07:13 cepek Exp $
 */

#include <vector>
#include <iterator>


#ifndef GNU_gama__pointer_list_h_gnugamapointerlist___gnu_gama_pointerlist
#define GNU_gama__pointer_list_h_gnugamapointerlist___gnu_gama_pointerlist


namespace GNU_gama {


  template <class T> class List;

  template <class T> class List<T*>
    {
      typedef std::vector<void*> Vector;
      Vector  vec;

      List(const List& cod);
      // List& operator=(const List& cod);
      
    public:    
      
      class iterator;
      class const_iterator;
      
      List() {}
      
      std::size_t size()  const { return vec.size();  }
      bool        empty() const { return vec.empty(); }
      
      void push_back(T* t)                        { vec.push_back(t);     }
      void pop_back()                             { vec.pop_back();       }
      void clear()                                { vec.clear();          }
      void erase (typename List<T*>::iterator i)  { vec.erase(i.vit);     }
      void insert(typename List<T*>::iterator i,  T* value)  
        { 
          vec.insert(i.vit, value); 
        }


      T* operator[](std::size_t n)
        { 
          return static_cast<T*>(vec[n]); 
        }
      const T* operator[](std::size_t n) const 
        { 
          return static_cast<T*>(vec[n]); 
        }

    
      class const_iterator
        // : public std::iterator <std::forward_iterator_tag, T> 
        {
        public:
          
          const_iterator()
            {
            }
          const_iterator(const typename Vector::const_iterator& p) : vit(p) 
            {
            }
          bool operator==(const const_iterator& x) const 
            { 
              return vit==x.vit; 
            }
          bool operator!=(const const_iterator& x) const 
            { 
              return vit!=x.vit; 
            }
          const_iterator& operator++()
            {
              ++vit;
              return *this;
            }
          const_iterator operator++(int)
            {
              const_iterator tmp(vit);
              ++vit;
              return tmp;
            }
          const T* operator*() const
            {
              return static_cast<T*>(*vit);
            }
          
        private:
          typename Vector::const_iterator vit;
          
        };
      
      const_iterator  begin() const { return vec.begin(); }
      const_iterator  end  () const { return vec.end  (); }


      class iterator 
        // : public std::iterator <std::forward_iterator_tag, T> 
        {
        public:

          iterator()
            {
            }
          iterator(const typename Vector::iterator& p) : vit(p)
            {
            }
          operator const_iterator() const
            {
              return const_iterator(vit);
            }
          bool operator==(const iterator& x) const 
            { 
              return vit==x.vit; 
            }
          bool operator!=(const iterator& x) const 
            { 
              return vit!=x.vit; 
            }
          iterator& operator++()
            {
              ++vit;
              return *this;
            }
          iterator operator++(int)
            {
              iterator tmp(vit);
              ++vit;
              return tmp;
            }
          T* operator*() const
            {
              return static_cast<T*>(*vit);
            }

        private:
          typename Vector::iterator vit;
          #ifdef __BORLANDC__
          friend class List<T*>;
          #else
          friend void List<T*>::erase (typename List<T*>::iterator i);
          friend void List<T*>::insert(typename List<T*>::iterator i, T* val);
          #endif

        };

      iterator  begin() { return vec.begin(); }
      iterator  end  () { return vec.end  (); }
    };

}

#endif
