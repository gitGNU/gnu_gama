/*  
    Geodesy and Mapping C++ Library (GNU GaMa / GaMaLib)
    Copyright (C) 2000  Ales Cepek <cepek@fsv.cvut.cz>

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
 *  $Id: gkfparser.h,v 1.1 2001/12/07 13:02:30 cepek Exp $
 */

#ifndef GaMaLib_GKF__XML__parser__h_
#define GaMaLib_GKF__XML__parser__h_


// GKFparser is just a simple C++ wrapper for XML parser expat

#include <expat/xmlparse/xmlparse.h>
#include <gamalib/local/gamadata.h>
#include <string>

namespace GaMaLib {

  class GKFexception : public GaMaLib::Exception 
    {
    public:
      int line;
      GKFexception(std::string s, int r) : GaMaLib::Exception(s), line(r) {}
    };
  
  class GKFparser 
    {
    public:
   
      // public data members
      
      std::string description;          // network description
      std::string gama_xml_version;
      
      // adjustment parameters
      
      std::string   TXT, STX, OPR;      // file names from gkf specification
      double m0_apr, konf_pr, tol_abs;  // implicitly 10, 0.95, 1000
      bool   typ_m0_apriorni;           // implicitly false
      
      std::string errString;
      int         errLineNumber;  
      int         errCode;              // -1 bad data in gkf; 0 OK; >0 expat
      
      // constructor and destructor
      
      GKFparser(GaMaLib::PointData& sb, GaMaLib::ObservationData& od);
      ~GKFparser();
      
      // expat parser interface
      
      void xml_parse(const char *s, int len, int isFinal) 
        { 
          int err = XML_Parse(parser, s, len, isFinal);
          if (err == 0)
            {
              // fatal error
              
              errString=std::string(XML_ErrorString(XML_GetErrorCode(parser)));
              errCode  =XML_GetErrorCode(parser);
              errLineNumber = XML_GetCurrentLineNumber(parser);
              
              throw GKFexception(errString, errLineNumber);
            }
          
          if (state == state_error)
            {
              errCode = -1;   
              // errLineNumber is set by function  error("...");    
              throw GKFexception(errString, errLineNumber);
            }
        }
      
      int gkf_characterDataHandler(const char* s, int len);
      int gkf_startElement(const char *cname, const char **atts);
      int gkf_endElement(const char * name);
      
      double implicit_stdev_direction() const { return smer_str; }
      double implicit_stdev_angle() const { return uhel_str; }
      double implicit_stdev_distance(double d) const 
        { 
          return delka_str + delka_str_km * pow(d/1000, delka_str_exp);
        }
      double implicit_stdev_distance_a() const { return delka_str;     }
      double implicit_stdev_distance_b() const { return delka_str_km;  }
      double implicit_stdev_distance_c() const { return delka_str_exp; }
      
    private: 
      
      GaMaLib::PointData&       SB;        // point list
      GaMaLib::ObservationData& OD;        // observation list
      
      enum gkf_tag {
        tag_unknown,
        tag_gama_xml,
        tag_network,
        tag_description,
        tag_parameters,
        tag_points_observations,
        tag_point,
        tag_obs,
        tag_cov_mat,
        tag_direction,
        tag_distance,
        tag_angle,
        tag_s_distance,
        tag_z_angle,
        tag_height_differences,
        tag_dh,
        tag_coordinates,
        tag_vectors,
        tag_vec
      };
      
      gkf_tag tag(const char* cname);

      enum gkf_state {
        state_error,
        state_start,
        state_gama_xml,
        state_network,
        state_description,
        state_parameters,
        state_point_obs,
        state_point,
        state_obs,
        state_obs_direction,
        state_obs_distance,
        state_obs_angle,
        state_obs_sdistance,
        state_obs_zangle,
        state_obs_dh,
        state_obs_cov,
        state_obs_after_cov,
        state_coords,
        state_coords_point,
        state_coords_cov,
        state_coords_after_cov,
        state_hdiffs,
        state_hdiffs_dh,
        state_hdiffs_cov,
        state_hdiffs_after_cov,
        state_vectors,
        state_vectors_vec,
        state_vectors_cov,
        state_vectors_after_cov,
        state_stop
      } state;
      
      std::vector<Double> sigma;
      Index        idim, iband;            // covariance matrix dim. / band
      bool         pp_xydef, pp_zdef;      // process_point();
      Double       pp_x, pp_y, pp_z;
      PointID      pp_id;      
      std::string  cov_mat_data;
      XML_Parser   parser;

      int error(std::string s) { return error(s.c_str()); }
      int error(const char* text)
        {
          // store only the first detected error
          if(errCode) return 1;

          errString = std::string(text);
          errCode   = -1;
          errLineNumber = XML_GetCurrentLineNumber(parser);
          state = state_error;
          return 1;
        }
      
      bool toDouble(const std::string&, double&) const;
      bool toIndex (const std::string&, Index& ) const;
      
      // Implicit value of stanpoint ID is set for sets of
      // directions/distances and/or angles.
      
      std::string         standpoint_id;
      
      StandPoint        * standpoint;
      Coordinates       * coordinates;
      HeightDifferences * heightdifferences;
      Vectors           * vectors; 
      
      int process_gama_xml   (const char** atts);
      int process_network    (const char** atts);
      int process_parameters (const char** atts);
      int process_point_obs  (const char** atts);
      int process_point      (const char** atts);
      int process_distance   (const char** atts);
      int process_angle      (const char** atts);
      int process_direction  (const char** atts);
      int process_sdistance  (const char** atts);
      int process_zangle     (const char** atts);
      int process_obs_dh     (const char** atts);
      
      int process_obs(const char** atts);    
      int finish_obs();
      
      int process_coords(const char** atts);
      int finish_coords();
      int process_coords_point(const char** atts);
      
      int process_hdiffs(const char** atts);
      int finish_hdiffs();
      int process_dh(const char** atts);
      
      int process_cov(const char** atts);
      int finish_cov(Cov&);
      int process_obs_cov(const char** atts)
        {
          state = state_obs_cov;
          return process_cov(atts);
        }
      int process_coords_cov(const char** atts)
        {
          state = state_coords_cov;
          return process_cov(atts);
        }
      int process_hdiffs_cov(const char** atts)
        {
          state = state_hdiffs_cov;
          return process_cov(atts);
        }
      int process_vectors_cov(const char** atts)
        {
          state = state_vectors_cov;
          return process_cov(atts);
        }
      
      int process_vectors(const char** atts);
      int finish_vectors();
      int process_vec(const char** atts);
      
      // implicit values of standard deviations
      
      double delka_str;
      double delka_str_km;
      double delka_str_exp;
      double smer_str;
      double uhel_str;
      
    };  // class GKFparser
}       // namespace GaMaLib


#endif

















