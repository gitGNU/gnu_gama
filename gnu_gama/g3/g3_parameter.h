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

/* $Id: g3_parameter.h,v 1.4 2003/03/23 18:39:53 cepek Exp $  */

#include <cstddef>

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

    enum {
      unused      = 0,
      fixed       = 1,
      free        = 2,
      constrained = 4+free
    };

    int  state(int s) const { return s & state_; }
    void set_state(int s)   { state_ = s;        }

  private:
    
    Parameter& operator=(const Parameter&);

    double val;
    double cor;
    size_t ind;
    int    state_;

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


}}
  
#endif
