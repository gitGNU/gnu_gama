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
 *  $Id: obsdata.h,v 1.17 2004/03/18 17:07:01 cepek Exp $
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
        : observation_data(od), act_obs(0), act_dim(0), act_nonz(0) 
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
      
      int  activeObs()  const { return act_obs;  }
      int  activeDim()  const { return act_dim;  }
      int  activeNonz() const { return act_nonz; }
      typename Observation::CovarianceMatrix  
           activeCov() const; 
      void scaleCov(int i, double sc);
      
    private:    // no copy ctor and no assignment
      
      Cluster(const Cluster&); 
      void operator=(const Cluster&);
      
      int act_obs, act_dim, act_nonz;
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
      act_obs   = 0;
      act_dim   = 0;
      act_nonz  = 0;
      int index = 0;
      Observation* p;
      for (typename List<Observation*>::iterator 
             i=observation_list.begin(); i!=observation_list.end(); ++i)
        {
          p = (*i);
          p->cluster = this;
          p->cluster_index = index++;
          if (p->active())
            {
              act_obs++;
              act_dim += p->dimension();
            }
        }

      if (act_dim)
        {
          int b = covariance_matrix.bandWidth();
          if (act_dim - 1 < b) b = act_dim - 1;
          act_nonz = act_dim*(b+1) - b*(b+1)/2;
        }
    }



  template <class Observation>
    typename Observation::CovarianceMatrix 
       Cluster<Observation>::activeCov() const
    {
      typedef std::size_t Index;
      const Index M      = covariance_matrix.rows();
      const Index N      = activeDim();
      const Index i_size = observation_list.size();
      Index active_band  = covariance_matrix.bandWidth();

      if (N)
        {
          if (N-1 < active_band) active_band = N-1;
        }
      else
        active_band = 0;

      typename Observation::CovarianceMatrix C(N, active_band);

      Index* ind = new Index[act_dim + 1];
      for (Index k=1, n=1, i=0; i<i_size; i++)
        {
          const Observation* obs = observation_list[i];

          if (obs->active())
            for (Index d=0; d<obs->dimension(); d++)
              {
                ind[k++] = n + d;
              }
          
          n += obs->dimension();
        }

      for (Index i=1; i<=N; i++)
        for (Index j=0; j<=active_band && i+j<=N; j++)
          C(i, i+j) = covariance_matrix(ind[i], ind[i+j]);

      delete[] ind;
      return C;
    }



  template <class Observation>
    void Cluster<Observation>::scaleCov(int p, double sc)
    {
      const int N = covariance_matrix.dim();
      int k = p + covariance_matrix.bandWidth();
      if (k > N) k = N;
      for (int i=p; i<=k; i++)
        {
          covariance_matrix(p, i) *= sc;
          covariance_matrix(i, p) *= sc;
        }
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
