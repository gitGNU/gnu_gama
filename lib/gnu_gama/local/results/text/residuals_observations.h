/*
    GNU Gama C++ library
    Copyright (C) 1999, 2010  Ales Cepek <cepek@fsv.cvut.cz>
                  2011  Vaclav Petras <wenzeslaus@gmail.com>
                  2013  Ales Cepek <cepek@gnu.org>

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
    along with this library; if not, write to the Free Software
    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/

/** \file residuals_observations.h
 * \brief Function for writing residuals observations
 *
 * \author Ales Cepek
 * \author Vaclav Petras (acyclic visitor pattern)
 */

#ifndef GaMa_GaMaProg_Opravy_Pozorovani__h_
#define GaMa_GaMaProg_Opravy_Pozorovani__h_

#include <gnu_gama/local/network.h>
#include <gnu_gama/local/observation.h>
#include <gnu_gama/statan.h>
#include <gnu_gama/utf8.h>
#include <algorithm>

#if 1
static float ResidualsObservations_N01(float x)   // local helper function
{
   double D, f;
   GNU_gama::NormalDistribution(double(x), D, f);
   return D;
}
#endif


/* *******************************************************************
 * local helper class StOpSort for sorting outlying observations (sort
 * by "studentized" residuals)
 * ******************************************************************* */

class StOpSort {

  GNU_gama::local::LocalNetwork* IS;

public:

  StOpSort(GNU_gama::local::LocalNetwork* is) : IS(is) {}
  bool operator()(int a, int b)
    {
      using namespace std;
      GNU_gama::local::Double sa = fabs(IS->studentized_residual(a));
      GNU_gama::local::Double sb = fabs(IS->studentized_residual(b));
      return sa > sb;
    }
};

namespace GNU_gama { namespace local {

/** \brief Writes observation short name to output stream. */
template <typename OutStream>
class WriteShortObservationName : public GNU_gama::local::AllObservationsVisitor
{
private:
    OutStream& out;
public:
    WriteShortObservationName(OutStream& outStream) : out(outStream)
    {}

    void visit(Distance* obs) { out << T_GaMa_distance; }
    void visit(Direction* obs) { out << T_GaMa_direction; }
    void visit(Angle* obs) { out << T_GaMa_angle; }
    void visit(H_Diff* obs) { out << T_GaMa_levell; }
    void visit(S_Distance* obs) { out << T_GaMa_s_distance; }
    void visit(Z_Angle* obs) { out << T_GaMa_z_angle; }
    void visit(X* obs) { out << T_GaMa_x; }
    void visit(Y* obs) { out << T_GaMa_y; }
    void visit(Z* obs) { out << T_GaMa_z; }
    void visit(Xdiff* obs) { out << T_GaMa_xdiff; }
    void visit(Ydiff* obs) { out << T_GaMa_ydiff; }
    void visit(Zdiff* obs) { out << T_GaMa_zdiff; }
    void visit(Azimuth* obs) { out << T_GaMa_azimuth; }

};

template <typename OutStream>
void ResidualsObservations(GNU_gama::local::LocalNetwork* IS, OutStream& out)
{
  if(IS->degrees_of_freedom() <= 1) return;

  using namespace std;
  using namespace GNU_gama::local;
  using GNU_gama::local::Double;

  const Vec&    v      = IS->residuals();
  const int     pocmer = IS->sum_observations();
  const double  scale  = IS->gons() ? 1.0 : 0.324;

  vector<int> odlehla;

  Double kki = IS->conf_int_coef();
  int imax = 1;         // index of maximal studentized residual
  {
    using namespace std;
    Double maxno = 0;
    for (int i=1; i<=pocmer; i++)
      {
        if (IS->obs_control(i) < 0.1) continue;

        Double no = fabs(IS->studentized_residual(i));
        if (no > maxno) {
          maxno = no;
          imax = i;
        }
        if (no > kki) odlehla.push_back(i);
      }
    if (odlehla.size() > 0)
      sort(odlehla.begin(), odlehla.end(), StOpSort(IS));
  }

  /* *****************************************************************
   * Review of residuals is printed twice. Firstly all observations
   * and then only outlying (if any are apresent)
   * ***************************************************************** */
  int max_pruchod = odlehla.size()==0 ? 1 : 2;
  for (int pruchod=1; pruchod<=max_pruchod; pruchod++)
    {

      if (pruchod == 1)
        out << T_GaMa_resobs_Review_of_residuals_analysis_obs << "\n"
            << underline(T_GaMa_resobs_Review_of_residuals_analysis_obs, '*')
            << "\n\n";
      else
        out << "\n\n"
            << T_GaMa_resobs_Outlying_observations << "\n"
            << underline(T_GaMa_resobs_Outlying_observations, '*') << "\n\n";

      out.width(IS->maxw_obs());
      out << "i" << " ";
      out.width(IS->maxw_id());
      out << T_GaMa_standpoint << " ";
      out.width(IS->maxw_id());
      out << T_GaMa_target
          << T_GaMa_resobs_header1;
      {   // for ...
        for (int i=0; i < (IS->maxw_obs() + 2*(IS->maxw_id()) + 10); i++)
          out << "=";
      }   // for ...
      // removed in 1.7.09 : out << T_GaMa_resobs_header2;
      if (IS->gons())
        out << "======== [mm|cc] =========== [mm|cc] ===\n\n";
      else
        out << "======== [mm|ss] =========== [mm|ss] ===\n\n";
      out.flush();

      WriteShortObservationName<OutStream> nameVisitor(out);

      PointID predcs = "";   // previous standpoint ID
      int max_ii = pruchod==1 ? pocmer : odlehla.size();
      for (int ii=1; ii<=max_ii; ii++)
        {
          int i = pruchod==1 ? ii : odlehla[ii-1];
          Observation* pm = IS->ptr_obs(i);
          out.width(IS->maxw_obs());
          out << i << " ";
          PointID cs = pm->from();
          out.width(IS->maxw_id());
          if (cs != predcs)
            out << Utf8::leftPad(cs.str(), IS->maxw_id());
          else
            out << " ";
          out << " ";
          PointID cc = pm->to();
          out << Utf8::leftPad(cc.str(), IS->maxw_id());
          out.setf(ios_base::fixed, ios_base::floatfield);

          // special case for angle (before calling name visitor)
          if (Angle* u = dynamic_cast<Angle*>(pm))
            {
              out << '\n';
              const int w = IS->maxw_obs() + 2 + 2*(IS->maxw_id());
              out << Utf8::leftPad((u->fs()).str(), w);
            }

          pm->accept(&nameVisitor);

          Double f  = IS->obs_control(i);
          out.precision(1);
          out.width(5);
          out << f;
          if (f < 0.1)    out << T_GaMa_resobs_no_control;   // uncontrolled
          else if (f < 5) out << T_GaMa_resobs_weak_control; // weak control
          else            out << "  ";
          out << ' ';

          double sc = 1.0;
          if (dynamic_cast<Direction*>(pm))
            sc = IS->PD.consistent() ? scale : -scale;
          else if (dynamic_cast<Angle*>(pm))
            sc = scale;
          else if (dynamic_cast<Z_Angle*>(pm))
            sc = scale;

          out.precision(3);
          out.width(9);
          out << v(i)*sc << ' ';
          out.precision(1);
          out.width(4);

          if (f >= 0.1)
            {
              using namespace std;
              Double no = fabs(IS->studentized_residual(i));
              out << no;

              if (i == imax)
                {
                  if (no > kki)  out << T_GaMa_resobs_mc_max_critical;
                  else           out << T_GaMa_resobs_mc_max;
                }
              else if (no > kki) out << T_GaMa_resobs_mc_critical;
              else               out << "   ";


              if ( (pm->ptr_cluster())->covariance_matrix.bandWidth() == 0 &&
                  (f >=5 || (f >= 0.1 && no > kki)))
                {
                  Double em = v(i) / (IS->wcoef_res(i)*IS->weight_obs(i));
                  out.width(7);
                  out << em*sc;

                  Double ev = em - v(i);
                  out.width(7);
                  out << ev*sc;
                }
            }

          out << '\n';
          out.flush();

          predcs = cs;  // previous standpoint ID
        }
    }

#if 1
  if (pocmer >= 30)
    {
      using namespace GNU_gama::local;

      out << "\n\n"
          << T_GaMa_resobs_normality_test << "\n"
          << underline(T_GaMa_resobs_normality_test, '=') << "\n\n";

      { // ****** Kolmogorov-Smirnov

        Vec   pv(pocmer);
        Float pvvar = 0, pvstr = 0, p;
        {
          for (int i=1; i<=pocmer; i++)
            {
              p = sqrt(IS->weight_obs(i))*v(i);

              Observation* obs = IS->ptr_obs(i);
              if (!IS->PD.consistent() && (dynamic_cast<Direction*>(obs) != nullptr)) p = -p;

              pv(i)  = p;
              pvvar += p*p;
              pvstr += p;
            }
        }
        pvstr /= pocmer;
        pvvar  = pvvar/pocmer - pvstr*pvstr;
        if (pvvar > 0)
          pvvar = sqrt(pvvar);
        else
          pvvar = 0;   // random noise

        if (pvvar)
          for (int i=1; i<=pocmer; i++) pv(i) = (pv(i) - pvstr) / pvvar;

        float  ks, prob;
        GNU_gama::KStest(pv.begin(),
                         pocmer, ResidualsObservations_N01, ks, prob);


        out.setf(ios_base::fixed, ios_base::floatfield);
        out.precision(1);
        out.width(5);
        out << "Test Kolmogorov-Smirnov : " << 100*prob << " %\n";

      }
    }

  if (Double cond = IS->cond())
    {
      out.setf(ios_base::scientific, ios_base::floatfield);
      out.precision(1);
      out << "\n"
          << T_GaMa_resobs_condition_number << cond << "\n";
    }
#endif

  out << "\n\n";
  out.flush();
}

}}

#endif















