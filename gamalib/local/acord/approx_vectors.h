/*  
    Geodesy and Mapping C++ Library (GNU GaMa / GaMaLib)
    Copyright (C) 2003  Jan Pytel  <pytel@gama.fsv.cvut.cz>

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
 *  $Id: approx_vectors.h,v 1.1 2003/08/16 16:30:35 cepek Exp $
 */

 
#ifndef GaMaLib_acord___ApproximateVectors__header___h
#define GaMaLib_acord___ApproximateVectors__header___h

#include <gamalib/local/gamadata.h>
#include <algorithm>
#include <vector>

namespace GaMaLib {


    class ApproximateVectors
	{
	private:
	    
	    struct ObservedVData {
		std::vector<const Xdiff*>  XD;
		std::vector<const Ydiff*>  YD;        
		std::vector<const Zdiff*>  ZD;
	    };

	    typedef std::vector<const Xdiff*>::const_iterator X_const_iterator;
	    typedef std::vector<const Ydiff*>::const_iterator Y_const_iterator;
	    typedef std::vector<const Zdiff*>::const_iterator Z_const_iterator;
	    
	    struct CopyVectors 
	    {
		mutable ObservedVData& OVD;
		PointData&       PD;
	  
		CopyVectors(ObservedVData& ovd, PointData& pd): OVD(ovd),PD(pd) {}
		
		void operator()(const Observation* obs) const
		{
		    LocalPoint& from = PD[obs->from()];
		    LocalPoint& to   = PD[obs->to  ()];
		    
		    if ( ApproximateVectors::unknown_xy(from) ||
			 ApproximateVectors::unknown_xy(to) )
		    {
			if (const Xdiff* v = dynamic_cast<const Xdiff*>(obs))
			    OVD.XD.push_back(v);
			else
			    if (const Ydiff* v = dynamic_cast<const Ydiff*>(obs))
				OVD.YD.push_back(v);
		    }
		    if ( ApproximateVectors::unknown_z(from) ||
			 ApproximateVectors::unknown_z(to) )
			if (const Zdiff* v = dynamic_cast<const Zdiff*>(obs))
			    OVD.ZD.push_back(v);
		}
	    };

	    bool obs_list_empty() const
		{
		    return OVD.XD.empty() && OVD.YD.empty() && OVD.ZD.empty();
		}

	    Y_const_iterator find_ydiff(const PointID& from_id, const PointID& to_id )
		{
		    for (Y_const_iterator ci = OVD.YD.begin(); ci != OVD.YD.end(); ++ci)
			if ( ( (*ci)->from() == from_id ) && ( (*ci)->to() == to_id) )
			    return ci;
		    return OVD.YD.end();
		}

	    static bool known_xy(const LocalPoint& lp)
		{
		    return lp.active_xy() && lp.test_xy();
		}
	    
	    static bool unknown_xy(const LocalPoint& lp)
		{
		    return lp.active_xy() && !lp.test_xy();
		}
	    
	    static bool known_z(const LocalPoint& lp)
		{
		    return lp.active_z() && lp.test_z();
		}
	    
	    static bool unknown_z(const LocalPoint& lp)
		{
		    return lp.active_z() && !lp.test_z();
		}

	    bool                missing_coords;
	    PointData&          PD;
	    ObservationData&    OD;

	    ObservedVData      OVD;

	public:
	    
	    ApproximateVectors(PointData& b, ObservationData& m);
	    void execute();
	    void print(std::ostream&);

	};

}   // namespace GaMaLib

#endif

