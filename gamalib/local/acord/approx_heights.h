/*  
    Geodesy and Mapping C++ Library (GNU GaMa / GaMaLib)
    Copyright (C) 2001  Ales Cepek <cepek@fsv.cvut.cz>,
                        Jan Pytel  <pytel@gama.fsv.cvut.cz>

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
 *  $Id: approx_heights.h,v 1.1 2001/12/07 12:40:58 cepek Exp $
 */

 
#ifndef GaMaLib_acord___ApproximateHeights_ApproxHeights___header___h
#define GaMaLib_acord___ApproximateHeights_ApproxHeights___header___h

#include <gamalib/local/gamadata.h>
#include <fstream>
#include <algorithm>
#include <list>

namespace GaMaLib {


  class ApproximateHeights 
    {
    private:

      struct ObservedHData {
        std::list<H_Diff*>     HD;        
        std::list<H_Diff>      tmpHD;        
        std::list<Distance*>   DI;
        std::list<S_Distance*> SD;
        std::list<Z_Angle*>    ZA;
      };

      struct CopyHeights 
      {
        ObservedHData& OHD;
        CopyHeights(ObservedHData& ohd) : OHD(ohd) {}
        void operator()(Observation* obs) const
        {
          if (H_Diff* h = dynamic_cast<H_Diff*>(obs))  
            OHD.HD.push_back(h);
          if (Distance* d = dynamic_cast<Distance*>(obs))  
            OHD.DI.push_back(d);
          if (S_Distance* s = dynamic_cast<S_Distance*>(obs))  
            OHD.SD.push_back(s);
          if (Z_Angle* z = dynamic_cast<Z_Angle*>(obs))  
            OHD.ZA.push_back(z);
        }
      };

      void make_heights();

      int                 missing_heights;
      PointData&          PD;
      ObservationData&    OD;

      ObservedHData      OHD;

    public: 
      
      ApproximateHeights(PointData& b, ObservationData& m);
      void execute();
      void print(std::ostream&);
    };

}   // namespace GaMaLib

#endif





