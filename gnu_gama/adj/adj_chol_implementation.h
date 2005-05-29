/*
  GNU Gama is a package for adjustment and analysis of geodetic observations
  Copyright (C) 2005  Ales Cepek

  This file is part of the GNU Gama C++ library.
  
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/

/*
 * $Id: adj_chol_implementation.h,v 1.12 2005/05/29 11:12:32 cepek Exp $
 */

#ifndef GNU_gama_adjustment_cholesky_decomposition_implementation__h
#define GNU_gama_adjustment_cholesky_decomposition_implementation__h

#include <algorithm>
#include <limits>

namespace GNU_gama {

  template <typename Float, typename Exc> 
  Index 
  AdjCholDec<Float, Exc>::defect()
  {
    return nullity;
  }



  template <typename Float, typename Exc> 
  Float AdjCholDec<Float, Exc>::q_xx(Index, Index)
  {
    throw Exception::adjustment("AdjCholDec::q_xx() NOT implemented"); 
    return 0;
  }



  template <typename Float, typename Exc> 
  Float 
  AdjCholDec<Float, Exc>::q_bb(Index i, Index j)
  {
    if (!this->is_solved) solve_me();

    const Mat<Float, Exc>& A = *pA;
    Vec<Float, Exc> aq(N);
    Float s;


    // aq = A_row(i) * Q0,  linearly dependent columns are ignored    
    
    for (Index k=1; k<=N; k++)
      {
        s = Float();        
        for (Index l=1; l<=N; l++)
          {
            s += A(i,l)*Q0(l,k);
          }
        aq(k) = s;
      }
    
    // s = aq * trans(A)_column(i) 

    s = Float();
    for (Index c=1; c<=N; c++)
      {
        s += aq(c)*A(j,c);
      }
    
    return s;
  }



  template <typename Float, typename Exc> 
  Float 
  AdjCholDec<Float, Exc>::q_bx(Index, Index)
  {
    throw Exception::adjustment("AdjCholDec::q_bx() NOT implemented");
    return 0;
  }



  template <typename Float, typename Exc> 
  bool 
  AdjCholDec<Float, Exc>::lindep(Index n)
  {
    solve_me();
    return nullity && perm(n) > N0;
  }



  template <typename Float, typename Exc> 
  void 
  AdjCholDec<Float, Exc>::min_x()
  {
    delete[] minx_i;
    minx_t = ALL;
    minx_n = 0;
  }



  template <typename Float, typename Exc> 
  void 
  AdjCholDec<Float, Exc>::min_x(Index n, Index minx[])
  {
    delete[] minx_i;
    minx_t = SUBSET;
    minx_n = n;
    minx_i = new Index[n];
    for (Index i=0; i<n; i++) minx_i[i] = minx[i];
  }



  template <typename Float, typename Exc> 
  Float 
  AdjCholDec<Float, Exc>::dot(const Mat<Float,Exc>& M, Index i, Index j)
    {
      Float s = Float();

      for (Index r, k=0; k<minx_n; k++) 
        {
          r = minx_i[k];
          s += M(r,i)*M(r,j);
        }

      return s;
    }



  template <typename Float, typename Exc> 
  void 
  AdjCholDec<Float, Exc>::solve_me()
  {
    if (this->is_solved) return;

    if (this->pw)
      throw 
        Exception::adjustment("AdjCholDec::reset(Mat<>, Vec<>, Vec<>) "
                              " NOT implemented");


    // project equations Ax = b 

    const Mat<Float, Exc>& A = *pA;
    const Vec<Float, Exc>& b = *pb;
  
    M = A.rows();
    N = A.cols();    

    // permutation vector (used in pivoting during cholesky decomposition)
    //
    // perm(i) = k       means that the original node 'k' is the i-th
    //                   node in the new ordering
    //
    // inv(perm(i)) = i  inverse permutation; invp(k) gives the position
    //                   in perm where the originnally numberd 'k' resides
    //
    // see George & Liu: Computer Solution of Large Sparse Positive
    // Definite Systems, Prentice-Hall, Inc., Englewood Cliffs, 1981

    perm.reset(N);
    for (Index i=1; i<=N; i++)
      {
        perm(i) = i;
      }


    // normal equations mat*x = rhs

    mat.reset(N);
    rhs.reset(N);
    for (Index i=1; i<=N; i++)
      {
        Float s = 0;
        for (Index k=1; k<=M; k++)
          {
            s += A(k,i)*b(k);
          }
        rhs(i) = s;
    
        for (Index j=i; j<=N; j++)
          {
            Float s = 0;
            for (Index k=1; k<=M; k++)
              {
                s += A(k,i)*A(k,j);
              }
            mat(i,j) = s;
            
          }
      }


    if (s_tol <= Float()) 
      {
        s_tol = 1000*std::numeric_limits<Float>::epsilon();
      }
    nullity = Index();


    for (Index column=1; column<=N; column++)
      {
        // pivoting

        Float pivot = mat(perm(column), perm(column));
        Index ipvt  = 0;
        for (Index i=column+1; i<=N; i++)
          {
            const Float t = mat(perm(i),perm(i));

            if (t > pivot)
              {
                pivot = t;
                ipvt  = i;
              }
          }
        if (ipvt) std::swap(perm(column),perm(ipvt));


        // with first linearly dependent column we are done ...

        if (pivot <= s_tol)
          {
            for (Index i=column; i<=N; i++)   // remove junk
              for (Index j=i; j<=N; j++)  
                mat(perm(j),perm(i)) = 0;


            nullity = N - column + 1;
            break;
          }


        // update pivot's submatrix [pivot, v'; v S]; S -= v*v'

        for (Index j=column+1; j<=N; j++)
          {
            const Float t = mat(perm(j),perm(column))/pivot; 

            for (Index i=j; i<=N; i++)
              {
                mat(perm(i), perm(j)) -= t*mat(perm(i),perm(column));
              }
          }


        // update column elements under the pivot

        for (Index pivot_row=column+1; pivot_row<=N; pivot_row++)
          {
            mat(perm(pivot_row), perm(column)) /= pivot;
          }
      }

    
    // inverse permutation

    invp.reset(N);
    for (Index i=1; i<=N; i++) invp(perm(i)) = i;



    // the particular solution 'x0' with all parameters corresponding
    // to linearly dependent colunms set to zero
    // **************************************************************

    N0 = N-nullity;                     // last independent column
    x0 = rhs;
    for (Index i=N0+1; i<=N; i++)
      {
        x0(perm(i)) = Float();
      }


    // forward substitution

    for (Index ii=2; ii<=N0; ii++)      // mat(1,1) == 1
      {
        const Index i = perm(ii);
        for (Index jj=1; jj<ii; jj++)
          {
            const Index j=perm(jj);
            x0(i) -= mat(i,j)*x0(j);
          }
      }


    for (Index ii=1; ii<=N0; ii++)
      {
        const Index i = perm(ii);
        x0(i) /= mat(i,i);
      }



    // backward substitution
    
    for (Index ii=N0-1; ii>=1; ii--)
      {
        const Index i=perm(ii);
        for (Index jj=ii+1; jj<=N0; jj++)
          {
            const Index j=perm(jj);
            x0(i) -= mat(i,j)*x0(j);
          }
      }



    // vector of residuals
    // *******************

    r.reset(b.dim());
    for (Index i=1; i<=M; i++)
      {
        r(i) = -b(i);
        for (Index jj=1; jj<=N0; jj++)
          {
            const Index j = perm(jj);
            r(i) += A(i,j)*x0(j);
          }
      }

    // inverse matrix (cofactors)
    // **************

    // for A = LDL', and Z = inv(A), we can compute Z recursively from
    // Z = inv(D)inv(L) - (I - L')Z

    Q0.reset(N);
    Q0.set_zero();
    for (Index column=N0; column>=1; column--)
      {
        const Index j = perm(column);

        // the diagonal element 

        Float zii = Float(1)/mat(j,j);
        for (Index kk=column+1; kk<=N0; kk++)
          {
            const Index k = perm(kk);
            zii -= mat(j,k)*Q0(k,j);        // z_ii -= u_jk*z_kj
          }
        Q0(j,j) = zii;

        // elements above the diagonal

        for (Index row=column-1; row>=1; row--)
          {
            const Index i=perm(row);
            Float zij = Float();
            for (Index kk=row+1; kk<=N0; kk++)
              {
                const Index k = perm(kk);
                zij -= mat(i,k)*Q0(k,j);    // z_ij -= u_ik*z_kj
              }
            Q0(i,j) = zij;
          }
      }

    // vector of unknown parameters
    // ****************************

    if (nullity == 0)
      {
        x = x0;
      }
    else
      {
        const Index N1 = nullity + 1;
        Mat<Float, Exc> G(N, N1);
        
        // matrix of linear combinations

        for (Index i=1; i<=N0; i++)
          for (Index j=1; j<=nullity; j++)
            { 
              G(perm(i),j) = mat(perm(i),perm(N0+j));
            }
                
        // backward substitution for each column

        for (Index column=1; column<=nullity; column++)  
          for (Index ii=N0-1; ii>=1; ii--)  
            {
              const Index i=perm(ii);
              for (Index jj=ii+1; jj<=N0; jj++)
                {
                  const Index j=perm(jj);
                  G(i,column) -= mat(i,j)*G(j,column);
                }
            }
        
        // negative identity matrix corresponding to fixed paramaters in x0

        for (Index i=1; i<=nullity; i++)
          for (Index j=1; j<=nullity; j++)
            {
              G(perm(N0+i), j) = (i==j ? -1 : 0);
            }

        // particular solution x0

        for (Index i=1; i<=N; i++)
          {
            G(i,N1) = x0(i);
          }

        cout << x0;

        // the particular solution minimizing subvector defined in min_x()
        // ***************************************************************

        Vec<Index, Exc> g_perm(N1);
        for (Index i=1; i<=N1; i++) g_perm(i) = i;


        if (minx_t == ALL && minx_n != N)
          {
            delete[] minx_i;
            minx_n = N;
            minx_i = new Index[N];
            for (Index i=1; i<=N; i++) minx_i[i-1] = i;
          }


        // Gramm-Schmidt orthogonalization

        for (Index column=1; column<=nullity; column++)
          {
            Float pivot = dot(G,g_perm(column),g_perm(column));
            if (pivot < s_tol) 
              throw Exception::adjustment("AdjCholDec::solve_me() --- "
                                          "bad regularization"); 
            Index ipvt  = 0;
            for (Index i=column+1; i<=nullity; i++)
              {
                const Float t = dot(G,g_perm(i),g_perm(i));
                
                if (t > pivot)
                  {
                    pivot = t;
                    ipvt  = i;
                  }
              }
            if (ipvt) std::swap(g_perm(column),g_perm(ipvt));

            const Index pc = g_perm(column);
            pivot = std::sqrt(pivot);
            for (Index i=1; i<=N; i++) G(i,pc) /= pivot;            

            for (Index col=column+1; col<=N1; col++)
              {
                const Index c = g_perm(col);
                Float dp = dot(G, pc, c);
                for (Index i=1; i<=N; i++) G(i,c) -= dp*G(i,pc);
              }
          }


        x.reset(N);
        for (Index i=1; i<=N; i++) x(i) = G(i,N1);
      }
    
    is_solved = true;
  }


}  // namespace GNU_gama

#endif
