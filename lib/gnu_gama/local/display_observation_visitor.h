/*
    GNU Gama -- adjustment of geodetic networks
    Copyright (C) 2000  Ales Cepek <cepek@fsv.cvut.cz>
                  2011  Vaclav Petras <wenzeslaus@gmail.com>
                  2013, 2014, 2015  Ales Cepek <cepek@gnu.org>

    This file is part of the GNU Gama C++ library.

    This library is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this library; if not, write to the Free Software Foundation,
    Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/

/** \file display_observation_visitor.h
 * \brief Display Observation Visitor
 *
 * \author Ales Cepek
 * \author Vaclav Petras (acyclic visitor pattern)
 */

#ifndef DISPLAY_OBSERVATION_VISITOR_H
#define	DISPLAY_OBSERVATION_VISITOR_H

#include <gnu_gama/local/network.h>

namespace GNU_gama { namespace local {

/** \brief Helper class for printing observational data.
   */

class LocalNetwork;

class DisplayObservationVisitor final : public AllObservationsVisitor
  {
  public:

    DisplayObservationVisitor(LocalNetwork* ln);

    std::string xml_name;
    std::string str_val;
    std::string str_stdev;
    std::string str_from;
    std::string str_to;
    std::string str_bs;
    std::string str_fs;

    void visit(Distance* obs);
    void visit(Direction* obs);
    void visit(Angle* obs);
    void visit(H_Diff* obs);
    void visit(S_Distance* obs);
    void visit(Z_Angle* obs);
    void visit(X* obs);
    void visit(Y* obs);
    void visit(Z* obs);
    void visit(Xdiff* obs);
    void visit(Ydiff* obs);
    void visit(Zdiff* obs);
    void visit(Azimuth* obs);

  private:

    LocalNetwork* lnet;
    const double  scale;
    const bool    consistent;
    void clear();
  };

}}

#endif	/* DISPLAY_OBSERVATION_VISITOR_H */

