/*  
    Geodesy and Mapping C++ Library (GNU GaMa / GaMaLib)
    Copyright (C) 2001  Ales Cepek <cepek@fsv.cvut.cz>

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
 *  $Id: angle.cpp,v 1.2 2002/05/24 19:30:51 cepek Exp $
 */

#include <iostream>
#include <iomanip>
#include <gamalib/observation.h>
#include <gamalib/local/pobs/bearing.h>
#include <gamalib/local/pobs/format.h>

using namespace GaMaLib;
using namespace std;


void Angle::write(std::ostream& out, bool print_at) const
{
  using namespace std;
  out << "<angle";
  if (print_at)
    out << " from=\"" << from() << '"';
  out << " bs=\"" << bs()  << '"'
      << " fs=\"" << fs() << '"'
      << " val=\""   << setprecision(Format::gon_p()  ) << value()*R2G << '"'
      << " stdev=\"" << setprecision(Format::stdev_p()) << stdDev()    << '"'
      << " />";
}




