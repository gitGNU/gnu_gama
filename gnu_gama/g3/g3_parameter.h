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

/* $Id: g3_parameter.h,v 1.6 2003/03/26 17:33:47 cepek Exp $  */

#include <cstddef>
#include <gnu_gama/list.h>


#ifndef GNU_gama_____g3______parameter______h__________GNUgamag3parameterh
#define GNU_gama_____g3______parameter______h__________GNUgamag3parameterh


namespace GNU_gama { namespace g3 {

  using std::size_t;
  class Parameter;


  class ParameterList {
  public:

    ParameterList() : begin_(0), end_(0)   {}
    ParameterList(int n);
    ~ParameterList();

    Parameter** begin() const { return begin_; }
    Parameter** end  () const { return end_;   }


  private:

    Parameter** begin_;
    Parameter** end_;
    
    ParameterList(const ParameterList& pl);
    ParameterList& operator=(const ParameterList&);

  };



  class Parameter {
  public:
    
    Parameter() : cor(0) {}
    Parameter(const Parameter&);
    virtual ~Parameter() {}
    
    virtual Parameter* clone() = 0;

    double value     () const { return val + cor; }
    double init_value() const { return val; }
    double correction() const { return cor; }
    size_t index     () const { return ind; }

    void set_init_value(double p) { val = p; cor = 0; }
    void set_correction(double p) { cor = p; }
    void set_index     (size_t t) { ind = t; }


    ParameterList  parlist;

    
    void set_unused() { state_ = unused_; }
    void set_fixed () { state_ = fixed_;  }
    void set_free  () { state_ = free_;   }
    void set_constr() { state_ = constr_; }

    bool active() const { return state_ != unused_; }
    bool unused() const { return state_ == unused_; }
    bool fixed () const { return state_ == fixed_;  }
    bool free  () const { return state_ &  free_;   }
    bool constr() const { return state_ == constr_; }

  private:
    
    Parameter& operator=(const Parameter&);

    double val;
    double cor;
    size_t ind;

    enum 
      {
        unused_ = 0,
        fixed_  = 1,
        free_   = 2,
        constr_ = 4 + free_
      } state_;

  };


  inline ParameterList::ParameterList(int n) 
    : begin_(new Parameter*[n]), end_(begin_+n) 
    {
    }

  inline Parameter::Parameter(const Parameter& par)
    {
      val = par.val;
      cor = par.cor;
      ind = par.ind;
      state_ = par.state_;
    }

  inline ParameterList::~ParameterList() 
    { 
      delete[] begin_; 
    }



  class ParameterTree {
  private:

    typedef GNU_gama::List<Parameter*> Tree;
    Tree tree;

    void add_parlist(const ParameterList& parlist);

  public:

    ParameterTree(const ParameterList& parlist) { add_parlist(parlist); }

    class const_iterator
      // : public std::iterator <std::forward_iterator_tag, Parameter*> 
      {
      public:
        
        const_iterator()
          {
          }
        const_iterator(const Tree::const_iterator& p) : tree(p) 
          {
          }
        bool operator==(const const_iterator& x) const 
          { 
            return tree==x.tree; 
          }
        bool operator!=(const const_iterator& x) const 
          { 
            return tree!=x.tree; 
          }
        const_iterator& operator++()
          {
            ++tree;
            return *this;
          }
        const_iterator operator++(int)
          {
            const_iterator tmp(tree);
            ++tree;
            return tmp;
          }
        const Parameter* operator*() const
          {
            return *tree;
          }
        
      private:
        Tree::const_iterator tree;
        
      };
    
    const_iterator  begin() const { return tree.begin(); }
    const_iterator  end  () const { return tree.end  (); }
    
    
    class iterator 
      // : public std::iterator <std::forward_iterator_tag, Parameter*> 
      {
      public:
        
        iterator()
          {
          }
        iterator(const Tree::iterator& p) : tree(p)
          {
          }
        operator const_iterator() const
          {
            return const_iterator(tree);
          }
        bool operator==(const iterator& x) const 
          { 
            return tree==x.tree; 
          }
        bool operator!=(const iterator& x) const 
          { 
            return tree!=x.tree; 
          }
        iterator& operator++()
          {
            ++tree;
            return *this;
          }
        iterator operator++(int)
          {
            iterator tmp(tree);
            ++tree;
            return tmp;
          }
        Parameter* operator*() const
          {
            return *tree;
          }
        
      private:
        Tree::iterator tree;
        
      };
    
    iterator  begin() { return tree.begin(); }
    iterator  end  () { return tree.end  (); }


  };


}}
  
#endif
