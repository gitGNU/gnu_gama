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
 *  $Id: obsdata.h,v 1.12 2004/01/17 21:36:14 cepek Exp $
 */


#include <gnu_gama/list.h>

#ifndef GNU_gama__obsdata_h_gnugamaobsdata_observation_data__gnu_gama_obsdata
#define GNU_gama__obsdata_h_gnugamaobsdata_observation_data__gnu_gama_obsdata


namespace GNU_gama {


  template <class Observation> 
    class ObservationData;


  template <class Observation>  
    class Cluster 
    {
    public:
      
      typedef Observation                     ObservationType;  
      const ObservationData<Observation>*     observation_data;
      List<Observation*>                      observation_list;
      typename Observation::CovarianceMatrix  covariance_matrix;  
      
      
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
      typename Observation::CovarianceMatrix  
           activeCov() const; 
      
    private:    // no copy ctor and no assignment
      
      Cluster(const Cluster&); 
      void operator=(const Cluster&);
      
      int act_count;
    };



  template <class Observation>
    class ObservationData
    {
    public:    

      typedef Observation           ObservationType;  
      typedef Cluster<Observation>  ClusterType;
      typedef List<ClusterType*>    ClusterList;
      ClusterList                   CL;
      
      ObservationData() {}
      ObservationData(const ObservationData& cod) { deepCopy(cod); }
      ~ObservationData();
      
      ObservationData& operator=(const ObservationData& cod);
      

      // ------------------------------------------------------------------

      class iterator;

      class const_iterator
        // : public std::iterator <std::forward_iterator_tag, Observation*>
        {
        public:
   
          bool operator==(const const_iterator& x) const 
            {
              const bool t1 =  cluster   ==   OD->CL.end();
              const bool t2 =  x.cluster == x.OD->CL.end();
              if (t1 || t2) return t1 == t2;

              return cluster == x.cluster && obs == x.obs; 
            }
          bool operator!=(const const_iterator& x) const 
            { 
              return !this->operator==(x); 
            }
          const_iterator& operator++()
            {
                goto next_cycle;
              
                while (cluster != OD->CL.end())
                  {
                    obs = (*cluster)->observation_list.begin();
                    while (obs != (*cluster)->observation_list.end())
                      {
                        return *this;
                        
                      next_cycle: 
                        ++obs;
                      }                    
                    ++cluster;
                  }
              
                return *this;
            }
          const_iterator operator++(int)
            {
              const_iterator tmp = *this;
              operator++();
              return tmp;
            }
          const Observation* operator*() const
            {
              return *obs;
            }

        private:
          friend class ObservationData;
          friend class iterator;

          const ObservationData*  OD;
          typename ObservationData::ClusterList::const_iterator  cluster;
          typename List<Observation*>::const_iterator            obs;
        };
 
      const_iterator  begin() const 
        {     
          const_iterator  iter;
          iter.OD       = this;
          iter.cluster  = iter.OD->CL.begin();

          while (iter.cluster != iter.OD->CL.end())
          {
            iter.obs     = (*iter.cluster)->observation_list.begin();

            if (iter.obs != (*iter.cluster)->observation_list.end()) 
              return iter;

            ++iter.cluster;
          }

          return iter;
        }

      const_iterator  end() const 
        {
          const_iterator  iter;
          iter.OD       = this;
          iter.cluster  = iter.OD->CL.end();

          return iter;
        }
      

      // ------------------------------------------------------------------
      
      class iterator
        // : public std::iterator <std::forward_iterator_tag, Observation*>
        {
        public:
   
          operator const_iterator() const
            {
              const_iterator ci;
              ci.OD      = OD;
              ci.cluster = cluster;
              ci.obs     = obs;
              return ci;
            }
          bool operator==(const iterator& x) const 
            {
              const bool t1 =  cluster   ==   OD->CL.end();
              const bool t2 =  x.cluster == x.OD->CL.end();
              if (t1 || t2) return t1 == t2;

              return cluster == x.cluster && obs == x.obs; 
            }
          bool operator!=(const iterator& x) const 
            { 
              return !this->operator==(x); 
            }
          iterator& operator++()
            {
                goto next_cycle;
              
                while (cluster != OD->CL.end())
                  {
                    obs = (*cluster)->observation_list.begin();
                    while (obs != (*cluster)->observation_list.end())
                      {
                        return *this;
                        
                      next_cycle: 
                        ++obs;
                      }                    
                    ++cluster;
                  }
              
                return *this;
            }
          iterator operator++(int)
            {
              iterator tmp = *this;
              operator++();
              return tmp;
            }
          Observation* operator*() const
            {
              return *obs;
            }

        private:
          friend class ObservationData;

          ObservationData*  OD;
          typename ObservationData::ClusterList::iterator  cluster;
          typename List<Observation*>::iterator            obs;
        };
 
      iterator  begin() 
        {     
          iterator  iter;
          iter.OD       = this;
          iter.cluster  = iter.OD->CL.begin();

          while (iter.cluster != iter.OD->CL.end())
          {
            iter.obs     = (*iter.cluster)->observation_list.begin();

            if (iter.obs != (*iter.cluster)->observation_list.end()) 
              return iter;

            ++iter.cluster;
          }

          return iter;
        }

      iterator  end() 
        {
          iterator  iter;
          iter.OD       = this;
          iter.cluster  = iter.OD->CL.end();

          return iter;
        }
      

    private:
      void deepCopy(const ObservationData& at);
      
    };
  

  // =====================================================================


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
    typename Observation::CovarianceMatrix Cluster<Observation>::activeCov() const
    {
      typedef std::size_t Index;
      const Index M = covariance_matrix.rows();
      const Index B = covariance_matrix.bandWidth();
      const Index N = activeCount();
      Index temp = B;
      if (N-1 < B) temp = N-1;
      typename Observation::CovarianceMatrix C(N, temp);
      
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
