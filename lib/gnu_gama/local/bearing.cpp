/* GNU Gama -- adjustment of geodetic networks
   Copyright (C) 1999, 2015  Ales Cepek <cepek@fsv.cvut.cz>

   This file is part of the GNU Gama C++ library.

   This library is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License as
   published by the Free Software Foundation; either version 3 of the
   License, or (at your option) any later version.

   This library is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this library; if not, write to the Free Software
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
   02110-1301 USA */

#include <gnu_gama/local/bearing.h>
#include <gnu_gama/local/language.h>

namespace GNU_gama { namespace local {

    double
    bearing(double dx, double dy, bool consistent)
    {
       if (dy == 0 && dx == 0)
          throw Exception(T_POBS_computation_of_bearing_for_identical_points);

       double b = std::atan2(consistent ? dy : -dy, dx);
       return b >= 0 ? b : b + 2 * M_PI;
    }

    double
    bearing(const LocalPoint& from, const LocalPoint& to, bool consistent)
    {
       return bearing(to.x()- from.x(), to.y()- from.y(), consistent);
    }

    void
    bearing_distance(const LocalPoint& from, const LocalPoint& to,
                     bool consistent, double& br, double& d)
    {
        double dx = to.x() - from.x();
        double dy = to.y() - from.y();
        br = bearing(dx, dy, consistent);
        d  = std::sqrt(dx*dx + dy*dy);
    }

}}  // GNU_gama::local