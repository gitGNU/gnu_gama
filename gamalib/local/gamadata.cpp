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
 *  $Id: gamadata.cpp,v 1.1 2001/12/07 12:38:37 cepek Exp $
 */

#include <gamalib/local/gamadata.h>

using namespace GaMaLib;


ObservationData::~ObservationData()
{
  for (ClusterList::iterator c=CL.begin(); c!=CL.end(); ++c) delete *c;
}

ObservationData& ObservationData::operator=(const ObservationData& cod)
{
  if (this != &cod)
    {
      for (ClusterList::iterator     c=CL.begin(); c!=CL.end(); ++c) delete *c;
      deepCopy(cod);
    }
  return *this;
}


void ObservationData::deepCopy(const ObservationData& cod)
{
  for (ClusterList::const_iterator ci=cod.CL.begin(); ci!=cod.CL.end(); ++ci)
    {
      Cluster* current = (*ci)->clone(this);
      
      ObservationList::const_iterator begin = (*ci)->observation_list.begin();
      ObservationList::const_iterator end   = (*ci)->observation_list.end();
      for (ObservationList::const_iterator m=begin; m!=end; ++m)
        {
          current->observation_list.push_back( (*m)->clone() );
        }
      
      current->update();
      CL.push_back( current);
    }
  
}


