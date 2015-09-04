/*
    GNU Gama -- adjustment of geodetic networks
    Copyright (C) 2001, 2011  Ales Cepek <cepek@fsv.cvut.cz>
                  2011  Vaclav Petras <wenzeslaus@gmail.com>
                  2013, 2015  Ales Cepek <cepek@gnu.org>

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
    along with this library; if not, write to the Free Software
    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/

#include <gnu_gama/local/local_linearization_visitor.h>
#include <gnu_gama/local/bearing.h>

using namespace GNU_gama::local;
using namespace std;

void LocalLinearizationVisitor::direction(const Direction* obs) const
{
   LocalPoint& sbod = PD[obs->from()];
   LocalPoint& cbod = PD[obs->to()];
   Double s, d;
   bearing_distance(sbod, cbod, s, d);
   // const Double p = m0 / obs->stdDev();
   const Double K = 10*R2G/d;
   const Double ps = K*sin(s);
   const Double pc = K*cos(s);

   const StandPoint*  csp = static_cast<const StandPoint*>(obs->ptr_cluster());
   StandPoint* sp = const_cast<StandPoint*>(csp);

   // Double w = p*p;                                          // weight
   double obsval = consistent ? obs->value() : 2*M_PI-obs->value();
   Double a = (obsval + sp->orientation() - s)*R2CC;
   while (a >  200e4) a -= 400e4;
   while (a < -200e4) a += 400e4;
   rhs = a;                                                    // rhs in cc

   size = 0;
   if (!sp->index_orientation()) sp->index_orientation(++maxn);
   index[ size ] = sp->index_orientation();
   coeff[ size ] = -1;
   size++;
   if (sbod.free_xy())
   {
      if (!sbod.index_x()) sbod.index_x() = ++maxn;
      if (!sbod.index_y()) sbod.index_y() = ++maxn;
      index[ size ] = sbod.index_y();
      coeff[ size ] = -pc;
      size++;
      index[ size ] = sbod.index_x();
      coeff[ size ] = ps;
      size++;
   }
   if (cbod.free_xy())
   {
      if (!cbod.index_x()) cbod.index_x() = ++maxn;
      if (!cbod.index_y()) cbod.index_y() = ++maxn;
      index[ size ] = cbod.index_y();
      coeff[ size ] = pc;
      size++;
      index[ size ] = cbod.index_x();
      coeff[ size ] = -ps;
      size++;
   }
}


void LocalLinearizationVisitor::distance(const Distance* obs) const
{
   LocalPoint& sbod = PD[obs->from()];
   LocalPoint& cbod = PD[obs->to()];
   Double s, d;
   bearing_distance(PD[obs->from()], PD[obs->to()], s, d);
   // Double p = M_0 / stdDev();
   Double ps = sin(s);
   Double pc = cos(s);

   // Double w = p*p;                  // weight
   rhs = (obs->value() - d)*1e3;       // abs. term in millimetres

   size = 0;
   if (sbod.free_xy())
   {
      if (!sbod.index_x()) sbod.index_x() = ++maxn;
      if (!sbod.index_y()) sbod.index_y() = ++maxn;
      index[ size ] = sbod.index_y();
      coeff[ size ] = -ps;
      size++;
      index[ size ] = sbod.index_x();
      coeff[ size ] = -pc;
      size++;
   }
   if (cbod.free_xy())
   {
      if (!cbod.index_x()) cbod.index_x() = ++maxn;
      if (!cbod.index_y()) cbod.index_y() = ++maxn;
      index[ size ] = cbod.index_y();
      coeff[ size ] = ps;
      size++;
      index[ size ] = cbod.index_x();
      coeff[ size ] = pc;
      size++;
   }
}


void LocalLinearizationVisitor::h_diff(const H_Diff* obs) const
{
   LocalPoint& sbod = PD[obs->from()];
   LocalPoint& cbod = PD[obs->to()];
   Double h = cbod.z() - sbod.z();
   // Double p = M_0 / stdDev();

   // Double w = p*p;                  // weight
   rhs = (obs->value() - h)*1e3;       // abs. term in millimetres

   size = 0;
   if (sbod.free_z())
   {
      if (!sbod.index_z())  sbod.index_z() = ++maxn;
      index[ size ] = sbod.index_z();
      coeff[ size ] = -1;
      size++;
   }
   if (cbod.free_z())
   {
      if (!cbod.index_z())  cbod.index_z() = ++maxn;
      index[ size ] = cbod.index_z();
      coeff[ size ] = 1;
      size++;
   }
}


void LocalLinearizationVisitor::s_distance(const S_Distance* obs) const
{
   LocalPoint& sbod = PD[obs->from()];
   LocalPoint& cbod = PD[obs->to()];
   // Double s, sd;
   // bearing_sdistance(PD[obs->from()], PD[obs->to()], s, sd);
   // Double p = M_0 / stdDev();
   Double dx = cbod.x() - sbod.x();
   Double dy = cbod.y() - sbod.y();
   Double dz = cbod.z() - sbod.z();
   Double sd = sqrt(dx*dx + dy*dy + dz*dz);
   if (sd == 0)
     throw GNU_gama::local::Exception(T_POBS_zero_or_negative_slope_distance);

   Double px = dx / sd;
   Double py = dy / sd;
   Double pz = dz / sd;

   // Double w = p*p;                // weight
   rhs = (obs->value() - sd)*1e3;    // abs. term in millimetres

   size = 0;
   if (sbod.free_xy())
   {
      if (!sbod.index_x()) sbod.index_x() = ++maxn;
      if (!sbod.index_y()) sbod.index_y() = ++maxn;
      index[ size ] = sbod.index_y();
      coeff[ size ] = -py;
      size++;
      index[ size ] = sbod.index_x();
      coeff[ size ] = -px;
      size++;
   }
   if (sbod.free_z())
   {
      if (!sbod.index_z()) sbod.index_z() = ++maxn;
      index[ size ] = sbod.index_z();
      coeff[ size ] = -pz;
      size++;
   }
   if (cbod.free_xy())
   {
      if (!cbod.index_x()) cbod.index_x() = ++maxn;
      if (!cbod.index_y()) cbod.index_y() = ++maxn;
      index[ size ] = cbod.index_y();
      coeff[ size ] = py;
      size++;
      index[ size ] = cbod.index_x();
      coeff[ size ] = px;
      size++;
   }
   if (cbod.free_z())
   {
      if (!cbod.index_z()) cbod.index_z() = ++maxn;
      index[ size ] = cbod.index_z();
      coeff[ size ] = pz;
      size++;
   }
}


void LocalLinearizationVisitor::x(const X* obs) const
{
   LocalPoint& point = PD[obs->from()];
   // Double p = M_0 / stdDev();

   // Double w = p*p;                          // weight
   rhs = (obs->value() - point.x())*1e3;       // abs. term in millimetres

   size = 0;
   if (point.free_xy())
   {
      if (!point.index_x()) point.index_x() = ++maxn;
      index[ size ] = point.index_x();
      coeff[ size ] = 1;
      size++;
   }
}


void LocalLinearizationVisitor::y(const Y* obs) const
{
   LocalPoint& point = PD[obs->from()];
   // Double p = M_0 / stdDev();

   // Double w = p*p;                          // weight
   rhs = (obs->value() - point.y())*1e3;       // abs. term in millimetres

   size = 0;
   if (point.free_xy())
   {
      if (!point.index_y()) point.index_y() = ++maxn;
      index[ size ] = point.index_y();
      coeff[ size ] = 1;
      size++;
   }
}


void LocalLinearizationVisitor::z(const Z* obs) const
{
   LocalPoint& point = PD[obs->from()];
   // Double p = M_0 / stdDev();

   // Double w = p*p;                          // weight
   rhs = (obs->value() - point.z())*1e3;       // abs. term in millimetres

   size = 0;
   if (point.free_z())
   {
      if (!point.index_z()) point.index_z() = ++maxn;
      index[ size ] = point.index_z();
      coeff[ size ] = 1;
      size++;
   }
}


void LocalLinearizationVisitor::xdiff(const Xdiff* obs) const
{
  LocalPoint& spoint = PD[obs->from()];               // stand point
  LocalPoint& tpoint = PD[obs-> to() ];               // target
  Double df = tpoint.x() - spoint.x();
  // Double p  = M_0 / stdDev();

  // Double w = p*p;                             // weight
  rhs = (obs->value() - df)*1e3;                 // abs. term in millimetres

  size = 0;
  if (spoint.free_xy())
    {
      if (!spoint.index_x()) spoint.index_x() = ++maxn;
      index[ size ] =  spoint.index_x();
      coeff[ size ] = -1;
      size++;
    }
  if (tpoint.free_xy())
    {
      if (!tpoint.index_x()) tpoint.index_x() = ++maxn;
      index[ size ] =  tpoint.index_x();
      coeff[ size ] = +1;
      size++;
    }
}


void LocalLinearizationVisitor::ydiff(const Ydiff* obs) const
{
  LocalPoint& spoint = PD[obs->from()];
  LocalPoint& tpoint = PD[obs-> to() ];
  Double df = tpoint.y() - spoint.y();
  // Double p = M_0 / stdDev();

  // Double w = p*p;
  rhs = (obs->value() - df)*1e3;

  size = 0;
  if (spoint.free_xy())
    {
      if (!spoint.index_y()) spoint.index_y() = ++maxn;
      index[ size ] =  spoint.index_y();
      coeff[ size ] = -1;
      size++;
    }
  if (tpoint.free_xy())
    {
      if (!tpoint.index_y()) tpoint.index_y() = ++maxn;
      index[ size ] =  tpoint.index_y();
      coeff[ size ] = +1;
      size++;
    }
}


void LocalLinearizationVisitor::zdiff(const Zdiff* obs) const
{
  LocalPoint& spoint = PD[obs->from()];
  LocalPoint& tpoint = PD[obs-> to() ];
  Double df = tpoint.z() - spoint.z();
  // Double p = M_0 / stdDev();

  // Double w = p*p;
  rhs = (obs->value() - df)*1e3;

  size = 0;
  if (spoint.free_z())
    {
      if (!spoint.index_z()) spoint.index_z() = ++maxn;
      index[ size ] =  spoint.index_z();
      coeff[ size ] = -1;
      size++;
    }
  if (tpoint.free_z())
    {
      if (!tpoint.index_z()) tpoint.index_z() = ++maxn;
      index[ size ] =  tpoint.index_z();
      coeff[ size ] = +1;
      size++;
    }
}


void LocalLinearizationVisitor::z_angle(const Z_Angle* obs) const
{
   LocalPoint& sbod = PD[obs->from()];
   LocalPoint& cbod = PD[obs->to()];
   // Double s, d, sd;
   // bearing_distance(PD[obs->from()], PD[obs->to()], s, d);
   // bearing_sdistance(PD[obs->from()], PD[obs->to()], s, sd);

   Double dx = cbod.x() - sbod.x();
   Double dy = cbod.y() - sbod.y();
   Double dz = cbod.z() - sbod.z();

   Double d2 = dx*dx + dy*dy;
   Double d  = sqrt(d2);
   Double sd = sqrt(d2 + dz*dz);
   if (d == 0 || sd == 0)
     throw GNU_gama::local::Exception(T_POBS_zero_or_negative_zenith_angle);

   Double k  = 10*R2G/(d*sd*sd);

   // Double p = M_0 / stdDev();
   Double px =  k * dz * dx;
   Double py =  k * dz * dy;
   Double pz = -k *  d *  d;
   // Double w = p*p;                  // weight

   Double za = acos(dz/sd);

   if (obs->value() > M_PI) za = 2*M_PI - za;
   Double a  = (obs->value() - za);

   rhs = a * R2CC; // abs. term in cc
   size = 0;
   if (sbod.free_xy())
   {
      if (!sbod.index_x()) sbod.index_x() = ++maxn;
      if (!sbod.index_y()) sbod.index_y() = ++maxn;
      index[ size ] = sbod.index_y();
      coeff[ size ] = -py;
      size++;
      index[ size ] = sbod.index_x();
      coeff[ size ] = -px;
      size++;
   }
   if (sbod.free_z())
   {
      if (!sbod.index_z()) sbod.index_z() = ++maxn;
      index[ size ] = sbod.index_z();
      coeff[ size ] = -pz;
      size++;
   }
   if (cbod.free_xy())
   {
      if (!cbod.index_x()) cbod.index_x() = ++maxn;
      if (!cbod.index_y()) cbod.index_y() = ++maxn;
      index[ size ] = cbod.index_y();
      coeff[ size ] = py;
      size++;
      index[ size ] = cbod.index_x();
      coeff[ size ] = px;
      size++;
   }
   if (cbod.free_z())
   {
      if (!cbod.index_z()) cbod.index_z() = ++maxn;
      index[ size ] = cbod.index_z();
      coeff[ size ] = pz;
      size++;
   }
}


void LocalLinearizationVisitor::angle(const Angle* obs) const
{
   LocalPoint& sbod  = PD[obs->from()];
   LocalPoint& cbod1 = PD[obs->bs()];
   LocalPoint& cbod2 = PD[obs->fs()];
   Double s1, d1, s2, d2;
   bearing_distance(PD[obs->from()], PD[obs->bs()], s1, d1);
   bearing_distance(PD[obs->from()], PD[obs->fs()], s2, d2);
   // Double p = m0 / obs->stdDev();
   const Double K1 = 10*R2G/d1;
   const Double K2 = 10*R2G/d2;
   const Double ps1 = K1*sin(s1);
   const Double pc1 = K1*cos(s1);
   const Double ps2 = K2*sin(s2);
   const Double pc2 = K2*cos(s2);

   // Double w = p*p;                          // weight
   Double ds = s2 - s1;
   if (ds < 0) ds += 2*M_PI;
   Double a = (obs->value() - ds)*R2CC;        // rhs
   // "big" positive/negative angle transformed to "lesser" positive/negative
   while (a > 200e4)
      a -= 400e4;
   while (a < -200e4)
      a += 400e4;
   rhs = a;          // rhs in cc

   size = 0;
   if (sbod.free_xy())
   {
      if (!sbod.index_x()) sbod.index_x() = ++maxn;
      if (!sbod.index_y()) sbod.index_y() = ++maxn;
      index[ size ] =  sbod.index_y();
      coeff[ size ] = -pc2 + pc1;
      size++;
      index[ size ] =  sbod.index_x();
      coeff[ size ] =  ps2 - ps1;
      size++;
   }
   if (cbod1.free_xy())
   {
      if (!cbod1.index_x()) cbod1.index_x() = ++maxn;
      if (!cbod1.index_y()) cbod1.index_y() = ++maxn;
      index[ size ] =  cbod1.index_y();
      coeff[ size ] = -pc1;
      size++;
      index[ size ] =  cbod1.index_x();
      coeff[ size ] =  ps1;
      size++;
   }
   if (cbod2.free_xy())
   {
      if (!cbod2.index_x()) cbod2.index_x() = ++maxn;
      if (!cbod2.index_y()) cbod2.index_y() = ++maxn;
      index[ size ] =  cbod2.index_y();
      coeff[ size ] =  pc2;
      size++;
      index[ size ] =  cbod2.index_x();
      coeff[ size ] = -ps2;
      size++;
   }
}


void LocalLinearizationVisitor::azimuth(const Azimuth* obs) const
{
   LocalPoint& sbod = PD[obs->from()];
   LocalPoint& cbod = PD[obs->to()];
   Double s, d;
   bearing_distance(sbod, cbod, s, d);
   const Double K = 10*R2G/d;
   const Double ps = K*sin(s);
   const Double pc = K*cos(s);

   // const StandPoint*  csp = static_cast<const StandPoint*>(obs->ptr_cluster()); ... unused
   // StandPoint* sp = const_cast<StandPoint*>(csp); ... unused

   Double a = (obs->value() + PD.xNorthAngle() - s)*R2CC;  // rhs

   while (a > 200e4)
      a -= 400e4;
   while (a < -200e4)
      a += 400e4;
   rhs = a;          // absolute termm is in cc

   size = 0;
   if (sbod.free_xy())
   {
      if (!sbod.index_x()) sbod.index_x() = ++maxn;
      if (!sbod.index_y()) sbod.index_y() = ++maxn;
      index[ size ] = sbod.index_y();
      coeff[ size ] = -pc;
      size++;
      index[ size ] = sbod.index_x();
      coeff[ size ] = ps;
      size++;
   }
   if (cbod.free_xy())
   {
      if (!cbod.index_x()) cbod.index_x() = ++maxn;
      if (!cbod.index_y()) cbod.index_y() = ++maxn;
      index[ size ] = cbod.index_y();
      coeff[ size ] = pc;
      size++;
      index[ size ] = cbod.index_x();
      coeff[ size ] = -ps;
      size++;
   }
}

