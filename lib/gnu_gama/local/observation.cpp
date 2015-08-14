/*
    GNU Gama --- adjustment of geoC++ library
    Copyright (C) 2000, 2010, 2014  Ales Cepek <cepek@gnu.org>

    This file is part of the GNU Gama C++ library

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

#include <gnu_gama/local/observation.h>
#include <gnu_gama/local/cluster.h>
#include <gnu_gama/local/network.h>
#include <gnu_gama/gon2deg.h>
#include "observation.h"
#include "network.h"

namespace GNU_gama { namespace local
{
  // the static member Observation::gons was introduced in 1.7.09 to
  // enable selection between grades and degrees in virtual functions
  // Observation::write() with minimal changes in existing old code of
  // GNU_gama::local

  bool Observation::gons = true;
}}

using namespace GNU_gama::local;
using namespace std;

bool  Observation::check_std_dev() const
{
  return cluster && cluster->covariance_matrix.dim();
}


Double Observation::stdDev() const
{
   return cluster->stdDev(cluster_index);
}

Double Direction::orientation() const
{
  StandPoint* sp = static_cast<StandPoint*>(cluster);
  return sp->orientation();
}

void Direction::set_orientation(Double p)
{
  StandPoint* sp = static_cast<StandPoint*>(cluster);
  sp->set_orientation(p);
}

bool Direction::test_orientation() const
{
  StandPoint* sp = static_cast<StandPoint*>(cluster);
  return sp->test_orientation();
}

void Direction::delete_orientation()
{
  StandPoint* sp = static_cast<StandPoint*>(cluster);
  sp->delete_orientation();
}

void Direction::index_orientation(int n)
{
  StandPoint* sp = static_cast<StandPoint*>(cluster);
  sp->index_orientation(n);
}

int Direction::index_orientation() const
{
  StandPoint* sp = static_cast<StandPoint*>(cluster);
  return sp->index_orientation();
}

// -------------------------------------------------------------------

DisplayObservationVisitor::DisplayObservationVisitor(LocalNetwork* ln)
  : lnet(ln), scale(ln->gons() ? 1.0 : 0.324),
    consistent(ln->PD.consistent())
{
}

void DisplayObservationVisitor::visit(Distance* obs)
{
  clear();
  xml_name  = "distance";
  str_from  = obs->from().str();
  str_to    = obs->to().str();
  str_val   = std::to_string(obs->value());
  str_stdev = std::to_string(obs->stdDev());
}

void DisplayObservationVisitor::visit(Direction* obs)
{
  clear();
  xml_name = "direction";
  str_from = obs->from().str();
  str_to   = obs->to().str();
  double obsval = consistent ? obs->value() : -obs->value();
  double m = R2G*obsval;
  while (m >= 400) m -= 400;
  while (m <   0 ) m += 400;
  if (lnet->gons())
    str_val =  std::to_string(m);
  else
    str_val =  GNU_gama::gon2deg(m, 0, 2);
  str_stdev = std::to_string(obs->stdDev()*scale);
}

void DisplayObservationVisitor::visit(Angle* obs)
{
  clear();
  xml_name = "angle";
  str_from = obs->from().str();
  str_bs   = obs->bs().str();
  str_fs   = obs->fs().str();
  double m = R2G*(obs->value());
  if (lnet->gons())
    str_val =  std::to_string(m);
  else
    str_val =  GNU_gama::gon2deg(m, 0, 2);
  str_stdev = std::to_string(obs->stdDev()*scale);
}

void DisplayObservationVisitor::visit(H_Diff* obs)
{
  clear();
  xml_name = "dh";
  str_from = obs->from().str();
  str_to   = obs->to().str();
  str_val   = std::to_string(obs->value());
  str_stdev = std::to_string(obs->stdDev());
}

void DisplayObservationVisitor::visit(S_Distance* obs)
{
  clear();
  xml_name = "s-distance";
  str_from = obs->from().str();
  str_to   = obs->to().str();
  str_val   = std::to_string(obs->value());
  str_stdev = std::to_string(obs->stdDev());
}

void DisplayObservationVisitor::visit(Z_Angle* obs)
{
  clear();
  xml_name = "z-angle";
  str_from = obs->from().str();
  str_to   = obs->to().str();
  double m = R2G*(obs->value());
  if (lnet->gons())
    str_val =  std::to_string(m);
  else
    str_val =  GNU_gama::gon2deg(m, 0, 2);
  str_stdev = std::to_string(obs->stdDev()*scale);
}

void DisplayObservationVisitor::visit(X* obs)
{
  clear();
  xml_name = "x";
  str_from = obs->from().str();
  str_val   = std::to_string(obs->value());
  str_stdev = std::to_string(obs->stdDev());
}

void DisplayObservationVisitor::visit(Y* obs)
{
  clear();
  xml_name = "y";
  str_from = obs->from().str();
  str_val   = std::to_string(obs->value());
  str_stdev = std::to_string(obs->stdDev());
}

void DisplayObservationVisitor::visit(Z* obs)
{
  clear();
  xml_name = "z";
  str_from = obs->from().str();
  str_val   = std::to_string(obs->value());
  str_stdev = std::to_string(obs->stdDev());
}

void DisplayObservationVisitor::visit(Xdiff* obs)
{
  clear();
  xml_name = "dx";
  str_from = obs->from().str();
  str_to   = obs->to().str();
  str_val   = std::to_string(obs->value());
  str_stdev = std::to_string(obs->stdDev());
}

void DisplayObservationVisitor::visit(Ydiff* obs)
{
  clear();
  xml_name = "dy";
  str_from = obs->from().str();
  str_to   = obs->to().str();
  str_val   = std::to_string(obs->value());
  str_stdev = std::to_string(obs->stdDev());
}

void DisplayObservationVisitor::visit(Zdiff* obs)
{
  clear();
  xml_name = "dz";
  str_from = obs->from().str();
  str_to   = obs->to().str();
  str_val   = std::to_string(obs->value());
  str_stdev = std::to_string(obs->stdDev());
}

void DisplayObservationVisitor::visit(Azimuth* obs)
{
  clear();
  xml_name = "azimuth";
  str_from = obs->from().str();
  str_to   = obs->to().str();
  double m = R2G*(obs->value());
  if (lnet->gons())
    str_val =  std::to_string(m);
  else
    str_val =  GNU_gama::gon2deg(m, 0, 2);
  str_stdev = std::to_string(obs->stdDev()*scale);
}

void DisplayObservationVisitor::clear()
{
  str_to  .clear();
  str_bs  .clear();
  str_fs  .clear();
}
