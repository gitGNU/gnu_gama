/*  
    Geodesy and Mapping C++ Library (GNU GaMa / GaMaLib)
    Copyright (C) 2002,2003  Jan Pytel  <pytel@gama.fsv.cvut.cz>

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
 *  $Id: reduce_observations.h,v 1.3 2003/01/20 17:57:17 cepek Exp $
 */

 
#ifndef GaMaLib_acord_reduce_observations_h
#define GaMaLib_acord_reduce_observations_h

#include <gamalib/local/gamadata.h>
#include <fstream>
#include <list>
#include <cstddef>

namespace GaMaLib {

  struct TypeOfZAngles;
    
  class ReducedObservations
    {
    public:

	enum TypeOfReduction {
	    none_     = 1,
	    approx_   = 2,
	    precise_  = 4,
	    nonexist_ = 8
	};
		
	class ReducedObs
	{
	public:
	    
	    ReducedObs(Observation* _obs, TypeOfReduction _type = none_)
		:ptr_obs(_obs),type_of_reduction(_type)
	    { 
		orig_value_ = ptr_obs->value();
	    }
	    
	    Observation*	ptr_obs;
	    TypeOfReduction	type_of_reduction;
	    
	    Double  orig_value() const
	    {
		return orig_value_;
	    }
	    
	    bool reduced() const
		{
		    return !( type_of_reduction & (none_ | nonexist_) );
		}
	    
	private:
	    friend class ReducedObservations;
	    Double orig_value_;
	};

	typedef std::list<ReducedObs> ListReducedObs;
	typedef std::list<ReducedObs>::iterator ListReducedObs_iter;
	typedef std::list<ReducedObs>::const_iterator ListReducedObs_c_iter;

	ListReducedObs_iter begin()
	{
	    return list_reduced_obs.begin();
	}

	ListReducedObs_iter end()
	{
	    return list_reduced_obs.end();
	}

	ListReducedObs_c_iter begin() const
	{
	    return list_reduced_obs.begin();
	}

	ListReducedObs_c_iter end()const
	{
	    return list_reduced_obs.end();
	}
	
    private:
	
	ReducedObs* giveReducedObs(const Observation* _obs)
	{
	    for (ListReducedObs_iter i  = list_reduced_obs.begin();
		                     i != list_reduced_obs.end(); ++i)
		if (i->ptr_obs == _obs)
 		    return &(*i);
	    return 0;
	} 
	
	
	struct CopyReducedObservation 
	{
	    ListReducedObs&  LRO;
	    ObservationList& OL;
	    
	    CopyReducedObservation(ListReducedObs& lro,ObservationList& ol) :
		LRO(lro),OL(ol) {}
	    
	    void operator()(Observation* obs) const
	    {
		
		if ( !obs->active() )
		    return;
		
		OL.push_back(obs);
		
		if ( (obs->from_dh() == 0) && (obs->to_dh() == 0 ) )
		    return;
		
		if (	dynamic_cast<S_Distance*>(obs) || 
			dynamic_cast<Z_Angle*   >(obs) ||
			dynamic_cast<Ydiff*     >(obs)  )  
		    LRO.push_back(obs);
	    }
	};
	
	struct RemoveNonActiveObs
	{
	    bool operator()(const ReducedObs& red_obs)
	    {
		return !red_obs.ptr_obs->active();
	    }
	};
	
	
	PointData&          PD;
	ObservationData&    OD;
	
	ListReducedObs  list_reduced_obs;
	ObservationList list_obs;
	
    protected:
	
	void reduce(ReducedObs&);
	
	void reduce_sdistance( ReducedObs* );
	void reduce_zangle   ( ReducedObs* ); 
	void reduce_ydiff    ( ReducedObs* );

	size_t number_of_reduced_observations_with_attribute(size_t attrib) const
	{
	    if ( ! list_reduced_obs.size() )
		return 0;
	    
	    size_t number = 0;
	    
	    for (ListReducedObs_c_iter ci  = list_reduced_obs.begin();
		 ci != list_reduced_obs.end(); ++ci)
		if ( ci->type_of_reduction & attrib )
		    number++;
	    
	    return number;
	}
	    
	size_t number_of_not_reduced_observations() const
	{	 
	    return number_of_reduced_observations_with_attribute(none_ | approx_);
	}
	
    public:
	
	ReducedObservations(PointData& b, ObservationData& m);

	size_t size() const
	{
	    return list_reduced_obs.size();
	}

	size_t size_nonexist() const
	{
	    return number_of_reduced_observations_with_attribute( nonexist_ );
	}
	
	void execute();
	void print(std::ostream&);
  };
    
}   // namespace GaMaLib

#endif


