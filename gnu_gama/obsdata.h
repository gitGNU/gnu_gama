/*  
    GNU Gama -- adjustment of geodetic networks
    Copyright (C) 2003  Ales Cepek <cepek@fsv.cvut.cz>

    This file is part of the GNU Gama C++ library.
    
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
 *  $Id: obsdata.h,v 1.6 2003/05/10 19:35:17 cepek Exp $
 */


#include <gnu_gama/list.h>
#include <vector>

#ifndef GNU_gama__obsdata_h_gnugamaobsdata_observation_data__gnu_gama_obsdata
#define GNU_gama__obsdata_h_gnugamaobsdata_observation_data__gnu_gama_obsdata


namespace GNU_gama {


  template <class Observation> 
    class ObservationData;


  template <class Observation>  
    class Cluster 
    {
    public:
      
      const ObservationData<Observation>*  observation_data;
      List<Observation*>                   observation_list;
      typename Observation::Cov            covariance_matrix;  
      
      
      Cluster(const ObservationData<Observation>* od) 
        : observation_data(od), act_count(0) 
        {
        }
      virtual ~Cluster();
      
      virtual Cluster* clone(const ObservationData<Observation>*) const = 0;
      
      double stdDev(int i) const 
        { 
          i++; using namespace std; return sqrt(covariance_matrix(i,i)); 
        }
      int size() const 
        { 
          return observation_list.size(); 
        }
      
      void update();
      
      int  activeCount() const { return act_count; }
      typename Observation::Cov  activeCov() const; 
      
    private:    // no copy ctor and no assignment
      
      Cluster(const Cluster&); 
      void operator=(const Cluster&);
      
      int act_count;
    };



  template <class Observation>
    class ObservationData
    {
    public:    

      List<Cluster<Observation>*>  CL;
      
      ObservationData() {}
      ObservationData(const ObservationData& cod) { deepCopy(cod); }
      ~ObservationData();
      
      ObservationData& operator=(const ObservationData& cod);
      
      template <class P> void for_each(const P& p) const
        {
          for (typename List<Cluster<Observation>*>::const_iterator 
                 c=CL.begin(); c!=CL.end(); ++c)
            {
              const Cluster<Observation>* cluster = (*c);
              std::for_each(cluster->observation_list.begin(),
                            cluster->observation_list.end(),  p);
            }
        }

      template <class P> void transform(P& p)
        {
          for (typename List<Cluster<Observation>*>::iterator 
                 c=CL.begin(); c!=CL.end(); ++c)
            {
              Cluster<Observation>* cluster = (*c);
              std::transform(cluster->observation_list.begin(),
                             cluster->observation_list.end(), 
                             cluster->observation_list.begin(),  p);
            }
        }
      
    private:
      void deepCopy(const ObservationData& at);
      
    };
  


  template <class Observation>
    Cluster<Observation>::~Cluster()
    {
      for (typename List<Observation*>::iterator 
             i=observation_list.begin(); i!=observation_list.end() ;  ++i)
        {
          delete *i;
        }
    }


 
  template <class Observation>
    void Cluster<Observation>::update()
    {
      act_count = 0;
      int index = 0;
      Observation* p;
      for (typename List<Observation*>::iterator 
             i=observation_list.begin(); i!=observation_list.end(); ++i)
        {
          p = (*i);
          p->cluster = this;
          p->cluster_index = index++;
          if (p->active()) act_count++;
        }
    }



  template <class Observation>
    typename Observation::Cov Cluster<Observation>::activeCov() const
    {
      typedef std::size_t Index;
      const Index M = covariance_matrix.rows();
      const Index B = covariance_matrix.bandWidth();
      const Index N = activeCount();
      Index temp = B;
      if (N-1 < B) temp = N-1;
      typename Observation::Cov C(N, temp);
      
      Index row = 1;
      Index col = row;
      for (Index i=0; i<M; ++i)
        if (observation_list[i]->active())
          {
            for (Index j=0; j<=B && i+j<M; ++j)
              {
                if (observation_list[i+j]->active())
                  {
                    C(row, col++) = covariance_matrix(i+1, i+j+1);
                  }
              }
            col = ++row;
          }
      
      return C;
    }



  template <class Observation>
    ObservationData<Observation>::~ObservationData()
    {
      for (typename List<Cluster<Observation>*>::iterator 
             c=CL.begin(); c!=CL.end(); ++c) 
        {
          delete *c;
        }
    }


  
  template <class Observation>
    ObservationData<Observation>& 
    ObservationData<Observation>::operator=(const ObservationData& cod)
    {
      if (this != &cod)
        {
          for (typename List<Cluster<Observation>*>::iterator 
                 c=CL.begin(); c!=CL.end(); ++c) delete *c;
          {
            deepCopy(cod);
          }
        }

      return *this;
    }
  

  
  template <class Observation>
    void ObservationData<Observation>::deepCopy(const ObservationData& cod)
    {
      for (typename List<Cluster<Observation>*>::const_iterator 
             ci=cod.CL.begin(); ci!=cod.CL.end(); ++ci)
        {
          Cluster<Observation>* current = (*ci)->clone(this);
          
          typename List<Observation*>::const_iterator 
            begin = (*ci)->observation_list.begin();
          typename List<Observation*>::const_iterator 
            end   = (*ci)->observation_list.end();
          for (typename List<Observation*>::const_iterator 
                 m=begin; m!=end; ++m)
            {
              current->observation_list.push_back( (*m)->clone() );
            }
          
          current->covariance_matrix = (*ci)->covariance_matrix;
          current->update();
          CL.push_back( current );
        }      
    }
  

}

#endif
