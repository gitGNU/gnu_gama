/*  
    GNU Gama --- Geodesy and Mapping C++ library 
    Copyright (C) 2004  Ales Cepek <cepek@fsv.cvut.cz>

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
 *  $Id: radian.h,v 1.1 2004/04/21 16:57:44 cepek Exp $
 */


#ifndef GNU_gama___gnu_gama_radian_h________________gnugamaradian_h
#define GNU_gama___gnu_gama_radian_h________________gnugamaradian_h


namespace GNU_gama {

#ifndef M_PI
#define M_PI          3.14159265358979323846264338328
#endif

#ifndef RAD_TO_GON
#define RAD_TO_GON    200.0/M_PI
#endif

#ifndef GON_TO_RAD
#define GON_TO_RAD    M_PI/200.0
#endif

#ifndef RAD_TO_CC
#define RAD_TO_CC     200.0e4/M_PI
#endif

#ifndef CC_TO_RAD
#define CC_TO_RAD     M_PI/200.0e4
#endif

}


#endif
