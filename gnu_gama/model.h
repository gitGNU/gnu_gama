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
 *  $Id: model.h,v 1.6 2003/12/27 21:00:58 uid66336 Exp $
 */


#include <gnu_gama/obsdata.h>

#ifndef GNU_gama__mathematical_model_h_gnugamamodel___gnu_gama_gmodel___h
#define GNU_gama__mathematical_model_h_gnugamamodel___gnu_gama_gmodel___h


namespace GNU_gama {

  // Three basic components of mathematical model (of geodetic
  // adjustment) are functional relations (class Model), unknown
  // parameters and constants (class Parameter) and observables (class
  // Observation). 

  // Model, Parameter and Observation classes are logically
  // related. To brake the source code dependency we use the 'acyclic
  // visitor' pattern, where Model objects are visiting Observation
  // objects.

  // ObservationVisitor is a completely degenerated class having only
  // the virtual destructor.

  class ObservationVisitor 
  {
  public: 
    virtual ~ObservationVisitor() {}
  };
  
  class Observation
  {
  public:

    Observation() : active_(true) {}

    virtual ~Observation() {}
    virtual int  dimension() const { return 1; }
    virtual bool revision_accept(ObservationVisitor* visitor) = 0;
    virtual void linearization_accept(ObservationVisitor* visitor) = 0;

    bool active() const     { return  active_;      }
    bool set_active(bool b) { return (active_ = b); }


  private:

    bool active_;
  };


  // .....................................................................
    
  template <class Observation> class Revision 
  {
  public:
    virtual bool revision_visit(Observation* observation) = 0;
  };
  
  template <class Observation> class Linearization 
  {
  public:
    virtual void linearization_visit(Observation* observation) = 0;
  };
    
  template <class Observation> 
  class Model : public ObservationVisitor
  {
  public:
    typedef Observation                           ObservationType;
    typedef ObservationData<Observation>          ObsData;

    ObsData  obsdata; 
  };
  
}

#endif
