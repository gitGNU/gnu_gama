/*  
    C++ Matrix/Vector templates (GNU Gama / gMatVec 0.9.22)
    Copyright (C) 2002  Ales Cepek <cepek@fsv.cvut.cz>

    This file is part of the gMatVec C++ Matrix/Vector template library.
    
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
 *  $Id: bandmat2.h,v 1.5 2004/06/12 19:21:59 cepek Exp $
 *  http://www.gnu.org/software/gama/
 */

#ifndef gMatVec_Symmetric_Band_Matrix_2__H_
#define gMatVec_Symmetric_Band_Matrix_2__H_

#include <gmatvec/gmatvec.h>
#include <gmatvec/choldec.h>
#include <algorithm>

/*
 * Symmetric Band Matrix
 * =====================
 *
 * Bandwidth is defined as max{ |i-j| | a_ij != 0 }
 * 
 * Upper triangular part of the matrix is stored by rows, ie
 *
 *        d*(b+1) - b*(b+1)/2   of unzero elements, 
 *
 * where `d' is the matrix dimension and `b' the bandwidth
 *
 */


namespace gMatVec {


template <class Float=double, class Exc=Exception>
class BandMat2 : public MatBase<Float, Exc>, public CholDec<Float, Exc> {
public:

  BandMat2() : band_(0), band_1(0), dim_b(0) {}
  BandMat2(Index d, Index b) 
    : MatBase<Float, Exc>(d,d,d*(b+1) - b*(b+1)/2), band_(b), 
      band_1(b+1), dim_b(d-b)
  {
  }
  
  void   reset() { row_=col_=band_=band_1=dim_b= 0; resize(0); }
  void   reset(Index d, Index b);
  Index  dim() const { return row_; }
  Index  bandWidth() const { return band_; }
  Float  operator()(Index, Index) const; 
  Float& operator()(Index, Index);
  void   cholDec();
  void   solve(Vec<Float, Exc>&) const;
  Vec<Float, Exc> operator*(const Vec<Float, Exc>&) const;
  
  const Float* operator[](Index row) const
  {
    const Float* a_ = begin() + --row*band_1;
    if (row > dim_b) {
	const Index i_  = row - dim_b;
	a_ -= i_*(i_+1)/2;
    }  
    return a_; 
  }
  Float* operator[](Index row)
  {
    Float* a_ = begin() + --row*band_1;
    if (row > dim_b) {
      const Index i_  = row - dim_b;
      a_ -= i_*(i_+1)/2;
    }  
    return a_; 
  }
  
  std::istream&  read (std::istream&);
  std::ostream&  write(std::ostream&) const;
  
private:
  
  Index   band_, band_1, dim_b;
    
};      /* class BandMat2 */
  

template <class Float, class Exc>
void BandMat2<Float, Exc>::reset(Index d, Index b)
{
  if (dim() != d || band_ != b) 
    {
      row_   = col_ = d;
      band_  = b; 
      band_1 = b+1;
      dim_b  = d-b;	
      resize(d*(b+1) - b*(b+1)/2);
    }
}


template <class Float, class Exc>
Float  BandMat2<Float, Exc>::operator()(Index r, Index s) const
{
   if (r > s) {
      Index t = r;
      r = s;
      s = t;
   }

   if (s > r+band_)
      return 0;

   s -= r;
   return *(operator[](r) + s);
}

template <class Float, class Exc>
Float& BandMat2<Float, Exc>::operator()(Index r, Index s)
{
   if (r > s) {
      Index t = r;
      r = s;
      s = t;
   }

   if (s > r+band_)
      throw Exc(BadIndex, "Float& BandMat2::operator()(Index r, Index s)");

   s -= r;
   return *(operator[](r) + s);
}

template <class Float, class Exc>
void BandMat2<Float, Exc>::cholDec()
{
  /*
   * Cholesky factorization of positive definite matrix A = L*D*trans(L)
   *
   * L is lower triangular matrix with unity diagonal; D is diagonal matrix.
   * Matrices L and D replace factored band symmetric matrix `in situ'.
   */
   using namespace std;

   Float *B = begin();
   Index  N = dim();
   Index  W = bandWidth();

   const  Float  Tol = Abs(*B*cholTol());
   Float *p;
   Index  row, k, l, n;
   Float  pivot, q;

   if (N == 0)
     throw Exc(BadRank, "BandMat2::cholDec(Float  tol) - zero dim matrix");

   for (row=1; row<=N; row++)
     {
        if((pivot = *B) < Tol)
          throw Exc(NonPositiveDefinite, "BandMat2::cholDec(Float  tol) - "
                                         "Matrix is not positive definite");

        k = min(W, N-row);             // number of of-diagonal elements
        p = B+k;                       // next row address - 1
        for (n=1; n<=k; n++)
	{
            q = B[n]/pivot;
	    for (l=n; l<=k; l++) p[l] -= q*B[l];
            p += min(W, N-row-n);
	}  
 
        B++;                           // *B++ = pivot = sqrt(pivot); 
        for (; k; k--) *B++ /= pivot;
     }
}

template <class Float, class Exc>
void BandMat2<Float, Exc>::solve(Vec<Float, Exc>& rhs) const
{
  using namespace std;
  Index i, j, k;
  Float s;
  const Float *m;
  
  // forward substitution
  for (i=2; i<=dim(); i++)
    {
      s = 0;
      for (j = i>band_ ? i-band_ : 1; j<i; j++) s += operator()(i,j)*rhs(j);
      rhs(i) -= s;
    }
  
  // inverse of diagonal
  for (i=1; i<=dim(); i++) rhs(i) /= *operator[](i);
  
  // backward substituiton
  for (i=dim()-1; i>0; i--)
    {
      s = 0;
      m = operator[](i) + 1;
      for (k=i+1; k<=min(i+band_,dim()); k++) s += *m++ * rhs(k);
      rhs(i) -= s;
    }
}


template <class Float, class Exc>
Vec<Float, Exc> 
BandMat2<Float, Exc>::operator*(const Vec<Float, Exc>& v) const
{
   Vec<Float, Exc> T(dim());
   Index i, j, k;
   Float s;

   for (i=1; i<=dim(); i++)
   {
      s = 0;

      k = i+band_;
      if (k > dim()) k = dim();
      for (j = i>band_ ? i-band_ : 1; j<=k; j++) 
        {
          s += operator()(i,j)*v(j);
        }

      T(i) = s;
   }

   return T;
}


template <class Float, class Exc>
std::istream& BandMat2<Float, Exc>::read(std::istream& inp) 
{
   int inpd, inpb;
   inp >> inpd >> inpb;
   reset(inpd, inpb);
 
   Float  *b = begin();
   for (Index i=1; i<=dim(); i++)
      for (Index j=i; j<=i+band_; j++)
         if (j <= dim())
	    inp >> *b++;

   return inp;
}


template <class Float, class Exc>
std::ostream& BandMat2<Float, Exc>::write(std::ostream& out) const
{
   int w = out.width();
   out.width(w);
   out << dim() << ' ';
   out.width(w);
   out << band_ << "\n\n";

   const Float  *b = begin();
   for (Index i=1; i<=dim(); i++, out << '\n')
     for (Index j=i; j<=i+band_; j++)
       if (j <= dim()) 
         {
           out.width(w);
           out << *b++ << ' ';
         }
   
   return out;
}


}      //  namespace gMatVec

#endif
