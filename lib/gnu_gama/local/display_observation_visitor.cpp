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

#include <gnu_gama/local/display_observation_visitor.h>
#include <gnu_gama/local/network.h>
#include <gnu_gama/gon2deg.h>

using namespace GNU_gama::local;

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
