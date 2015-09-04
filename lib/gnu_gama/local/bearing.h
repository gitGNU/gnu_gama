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

#ifndef GAMA_BEARING_DISTANCE_H
#define GAMA_BEARING_DISTANCE_H

#include <tuple>
#include <gnu_gama/local/gamadata.h>

namespace GNU_gama { namespace local {

double bearing(double dx, double dy);
double bearing(const LocalPoint &from, const LocalPoint &to);

void bearing_distance(const LocalPoint &from, const LocalPoint &to,
                      double &br, double &d);

}}

#endif //GAMA_BEARING_DISTANCE_H
