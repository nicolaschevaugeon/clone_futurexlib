/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/

#ifndef _xFixedCCSLowerMatrix_
#define _xFixedCCSLowerMatrix_
#include <algorithm>
#include <cassert>
#include <iostream>
#include <iterator>
#include <set>
#include <vector>

// DOFINDEX concept : must contain two value of type size_t : globalid and localid.
// globalid is the dof number as numbered in the CCSMatrix. Per historical convention
//  (from xCSRMatrix) global indices start at 1. Localid is the local id of the dof
//  in a local matrix Per historical convention (comming from xFemMatrix), local index start at 0.
//  A comparaison (comp) function object in order to sort DOFINDEX according to globalid must be
//  provided .

template <class DOFINDEX>
bool comp(const DOFINDEX &dof1, const DOFINDEX &dof2)
{
   return (dof1.globalid < dof2.globalid);
}

//! An implementation of the compressed row storage format.
//!  Only the lower triangular part of a matrix is stored. (usefull in particular for symmetric problems)
template <class T>
class xFixedCCSLowerMatrix
{
   // Member variables
  public:
   // number of matrix row.
   int nrows;
   // number of matrix columns
   int ncols;
   // number of matrix non-zero terms.
   int nnz;
   // Pointer to an array containing the non zero value of the Matrix. values are stored columns by columns. for each columns, the
   // nonzero value are sorted by ascending row index.
   T *values;
   // rowptr and colindex permit to get the index of a value in values. All index start at 0 (c convention)
   // Pointer to an array of integer. The array must be of size ncols. rowptr[k] return the line number corresponding to values[k]
   int *rowptr;
   // Pointer to an array of integer. colindex[j] return the index of the start of col j in rowptr as well as in values. rowptr
   // size must be ncolumns+1, rowptr[ncolumns] =  nnz
   int *colindex;
   // owner is a bool telling if the object "own" the memory for the array values, rowpt, and colindex.
   // If owner is true, the array will be deleted upon matrix destruction or resizing.
   bool owner;
   double tmp;

  public:
   //! When calling this constructor, memory is alocated for the user.
   xFixedCCSLowerMatrix(const int &_nrows, const int &_ncolumns, const int &_nnz)
       : nrows(_nrows),
         ncols(_ncolumns),
         nnz(_nnz),
         values(new T[nnz]),
         rowptr(new int[nnz]),
         colindex(new int[ncols + 1]),
         owner(true)
   {
      colindex[0] = 0;
      colindex[ncols] = nnz;
      for (int i = 0; i < nnz; ++i) values[i] = 0.;
   };

   // For this constructor, the data are supposed to be allocated by the user, and directly given, to the constructor
   /*xFixedCRSMatrix(const int &_nrows, const int &_ncolumns, const int &_nnz, T* _values, int * _columns, int *
     _rowindex):nrows(_nrows), ncolumns(_ncolumns), nnz(_nnz), values(_values), columns(_columns), rowindex(_rowindex),
     owner(false){ rowindex[0] = 0; rowindex[nrows] = nnz ;
     };*/

   ~xFixedCCSLowerMatrix()
   {
      if (owner)
      {
         delete[] values;
         delete[] colindex;
         delete[] rowptr;
      }
   }

   //! indexing the matrix, const version
   const T &operator()(const int &i, const int &j) const
   {
      return (const_cast<xFixedCCSLowerMatrix<T> *>(this))->operator()(i, j);
   }

   //! indexing of the matrix
   /*! Some safety is provided via assert in debug mode if index (i,j) is zero or in the upper triangular part
     out of debug mod, be aware that the operator return a wrong adress !
   */
   T &operator()(const int &i, const int &j)
   {
      assert(i >= j);
      assert(i < nrows);
      assert(j < ncols);
      const int jstart = colindex[j];
      const int jend = colindex[j + 1];
      int *pos = lower_bound(&rowptr[jstart], &rowptr[jend], i);
      assert(pos != &rowptr[jend] && (*pos == i));
      int k = distance(rowptr, pos);
      return values[k];
   }

   //! xfem style AddMatrix operation (i and j are given starting at one (fortran convention))
   /*! The term will be added only if i>=j and if the index (i,j) does not correspond to a null term */
   void AddMatrix(int i, int j, T val)
   {
      if (i >= j) this->operator()(i - 1, j - 1) += val;
   }

   /*
   //!  Blocked AddMatrix operation (i and j are given starting at one (fortran convention))
   void AddMatrix(const int * indexi, const int *localindexi, const size_t sizei, const int * indexj, const int *localindexj,
   const size_t sizej, const T *vallocal, const int localsize){ for (int j =0; j < sizej; ++j) AddMatrix(indexi, localindexi,
   sizei, indexj[j], localindexj[j], vallocal, localsize);
   }

   void AddMatrix( const int * indexi, const int *localindexi, const size_t sizei, const int &gj, const int &lj, const T
   *vallocal, const int localsize){ const int jstart  = colindex[gj-1]; const int jend    = colindex[gj]; int *istart =
   &rowptr[jstart]; int *iend = &rowptr[jend];
     // const double *vallocalj = &vallocal[ lj*localsize];
     for (int i = 0; i < sizei; ++i){
       const int gi = indexi[i] -1;
       if (gi >= (gj-1)){
         const int li = localindexi[i];
         int *pos =  lower_bound(istart, iend ,gi);
         assert (pos !=  &rowptr[jend]&& (*pos == gi)  );
         int k = distance(rowptr, pos );
         values[k]+= vallocal[ li*localsize + lj] ;
       }
     }
   }
   */

   //! Blocked Version of Add Matrix
   /*! On entry, dofsi and dofsj are vector of DOFINDEX.
     valloc is a pointer to an array of value, corresponding to a square matrix of size localsize, stored row by row.
     The values of the local matrix are added in the CCSMatrix column by columns. A gain in performance can be acchieved since the
     indices in dofsi can be sorted for faster retrievial of the global value. dofsi contain the globals and locals row index.
     dofsj contain the globals and locals columns index.
   */
   template <class DOFINDEX>
   void AddMatrix(const std::vector<DOFINDEX> &dofsi, const std::vector<DOFINDEX> &dofsj, const T *vallocal, const int localsize)
   {
      std::vector<DOFINDEX> dofsis = dofsi;
      sort(dofsis.begin(), dofsis.end(), comp<DOFINDEX>);
      const size_t sizej = dofsj.size();
      for (size_t j = 0; j < sizej; ++j) AddMatrixSorted(dofsis, dofsj[j], vallocal, localsize);
   }

   //! Blocked version of AddMatrix ... Follow up
   /*!  dofsi contain the globals and locals row index. dofj contain the columns number in global and local matrix.
     dofsi is supposed to be sorted according to global index. */
   template <class datadof>
   void AddMatrixSorted(const std::vector<datadof> &dofsi, const datadof &dofj, const T *vallocal, const int localsize)
   {
      const size_t lj = dofj.localid;
      const size_t gj = dofj.globalid - 1;
      const size_t jstart = colindex[gj];
      const size_t jend = colindex[gj + 1];
      const size_t sizei = dofsi.size();
      int *istart = &rowptr[jstart];
      int *iend = &rowptr[jend];
      typename std::vector<datadof>::const_iterator iti = dofsi.begin();
      typename std::vector<datadof>::const_iterator itiend = dofsi.end();
      for (; iti != itiend; ++iti)
      {
         size_t gi = (*iti).globalid - 1;
         if (gi >= gj) break;
      }
      for (; iti != itiend; ++iti)
      {
         const size_t gi = (*iti).globalid - 1;
         const size_t li = (*iti).localid;
         int *pos = lower_bound(istart, iend, gi);
         assert(pos != &rowptr[jend] && (*pos == gi));
         int k = distance(rowptr, pos);
         values[k] += vallocal[lj * localsize + li];
         istart = pos + 1;
      }
   }

   void print()
   {
      for (int j = 0; j < ncols; ++j)
      {
         cout << j << endl;
         int startj = colindex[j];
         int endj = colindex[j + 1];
         for (int k = startj; k < endj; ++k)
         {
            cout << "A(" << rowptr[k] << ", " << j << ")= " << values[k] << endl;
         }
      }
   }
};

/*!
  Starting from a GRAPH, build rowptr and colindex.
  The graph on entry is supposed to be a random access container of container. Each random entry of the graph represent the   a
  column of the matrix. For each column of the matrix the corresponding container is supposed to contain the non-zero line of the
  corresponding column, in ascending order. the array columns is supposed to have at least size  nnz where nnz is the total number
  of non zero term in that Matrix. the array rowindex is supposed to have a size n+1, where n is the number of column of the
  matrix.
  !*/
template <class GRAPH>
void graphToCCSLow(const GRAPH &graph, int *rowptr, int *colindex)
{
   int ncols = graph.size();
   colindex[0] = 0;
   for (size_t j = 0; j < ncols; ++j)
   {
      colindex[j + 1] = colindex[j] + graph[j].size();
      std::copy(graph[j].begin(), graph[j].end(), &rowptr[colindex[j]]);
   }
}
#endif
