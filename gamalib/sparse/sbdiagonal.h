/*  
    Geodesy and Mapping C++ Library (GNU GaMa / GaMaLib)
    Copyright (C) 2002  Ales Cepek <cepek@fsv.cvut.cz>

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
 *  $Id: sbdiagonal.h,v 1.2 2002/09/13 16:21:45 cepek Exp $
 */

#ifndef GaMaLib_Symmetric_Block_Diagonal____GaMaLib_Symmetric_Block_Diagonal__
#define GaMaLib_Symmetric_Block_Diagonal____GaMaLib_Symmetric_Block_Diagonal__

#include <cstddef>
#include <cstring>
#include <algorithm>
#include <cmath>


namespace GaMaLib {

template <class Float=double, class Index=std::size_t> 

  class BlockDiagonal {    // symmetric block diagonal matrix
    
    Index   blocks_;       // number of diagonal blocks
    Float*  nonz_;         // all nonzero elements
    Index   ncnt_;         // number of all nonzeroes
    Index   size_;
    Index*  dim_;
    Index*  width_;
    Float** begin_;
    
    BlockDiagonal (const BlockDiagonal&); 
    void operator=(const BlockDiagonal&);

    public:
    
    BlockDiagonal(Index blcks, Index floats) 
    {
      dim_    = new Index [blcks+1];  // 1 based indexes
      width_  = new Index [blcks+1];
      begin_  = new Float*[blcks+2];
      nonz_   = new Float [floats];

      ncnt_ = size_ = blocks_ = 0;    // matrix is initially empty
      begin_[1] = nonz_;
    }
    
    ~BlockDiagonal() 
    {
      delete[] dim_;
      delete[] width_;
      delete[] begin_;
      delete[] nonz_;
    }
    
    Index  blocks()        const { return blocks_;   }
    Index  dim()           const { return size_;     }
    Index  nonzeroes()     const { return ncnt_;     }
    Index  dim   (Index i) const { return dim_  [i]; }
    Index  width (Index i) const { return width_[i]; }

    const Float* begin (Index i) const { return begin_[i];   } 
    const Float* end   (Index i) const { return begin_[i+1]; }
    Float*       begin (Index i)       { return begin_[i];   }
    Float*       end   (Index i)       { return begin_[i+1]; }

    BlockDiagonal* replicate() const 
    { 
      return replicate(blocks_, ncnt_);
    }

    BlockDiagonal* replicate(Index new_floats, Index new_blocks) const 
    {
      BlockDiagonal* r = new BlockDiagonal(new_blocks, new_floats);

      for (Index i=1; i<=blocks_; i++)  
      {
        r->add_block(dim(i), width(i), begin(i));
      }

      return r;
    }

    void add_block(Index bdim, Index bwidth, const Float* mem)
    {
      Index N = bdim*(bwidth+1) - (bwidth)*(bwidth+1)/2;

      blocks_++;
      size_ += bdim;
      ncnt_ += N;
      Float* b = begin_[blocks_];
      begin_[blocks_+1] = b + N; 
      memcpy(b, mem, N*sizeof(Float));

      dim_  [blocks_] = bdim;
      width_[blocks_] = bwidth;
    }


    int cholDec(Float tol = 1e-14)
    {
      Float* B;
      Float* p;
      Index  N, W, row, k, l, n;
      Float  q, pivot;

      for (Index block=1; block<=blocks_; block++)
      {
        B = begin(block);
        N = dim  (block);
        W = width(block);

        for (row=1; row<=N; row++)
          {
            if ((pivot = *B) < tol) 
              return block;                // not positive-definite

            k = min(W, N-row);             // number of of-diagonal elements
            p = B+k+1;                     // first element of next submatrix
            for (n=1; n<=k; n++)
              {
                q = B[n]/pivot;
                for (l=n; l<=k; l++) *p++ -= q*B[l];
              }
            *B++ = pivot = sqrt(pivot);    // scaling pivot row 
            for (; k; k--) *B++ /= pivot;
          }

      }   // block ...

      return 0;
    }

  };
 
}   // namespace GaMaLib

#endif



#ifdef GaMaLib_sparse_demo

#include <iostream>

using namespace std;
using namespace GaMaLib;


void write(ostream& cout, BlockDiagonal<>* bd)
{
  // cout.precision(7); 
  cout << endl;
  for (unsigned long i=1; i<=bd->blocks(); i++)
    {
      cout << i << " : [" << bd->dim(i) << " | " << bd->width(i) << "]\n";
      double* b = bd->begin(i);
      double* e = bd->end(i);
      cout << "    ";
      while (b != e) cout << *b++ << ' ';
      cout << endl;
    }
}

int main()
{
  cout << "\n---  Symmetric Block Diagonal Matrix demo  -----------------\n";
      
  double b1[] = {1.1, 1.2, 1.3};
  double b2[] = {44.2, 5.2, 66.2, 7.2, 88.2};
  double b3[] = {111.3, 12.3, 13.3, 114.3, 15.3, 116.3};
  double b4[] = {19, 1, 2, 3, 18, 1, 2, 17, 1, 16};
  double b5[] = {19, 1, 2, 18, 1, 2, 17, 1, 2, 16, 1, 5};

  BlockDiagonal<>* m1 = new BlockDiagonal<>(20, 5000);
  m1->add_block(3, 0, b1);
  m1->add_block(3, 1, b2);
  m1->add_block(3, 2, b3);
  m1->add_block(4, 3, b4);
  m1->add_block(5, 2, b5);
  write(cout, m1);
  
  BlockDiagonal<>* m2 = m1->replicate();
  m2->cholDec();
  write(cout, m2);

  delete m1;
  delete m2;
}


#endif
