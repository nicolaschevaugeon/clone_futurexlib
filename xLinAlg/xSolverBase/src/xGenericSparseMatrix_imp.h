/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
 */

#ifndef _XGENERICSPARSEMATRIXIMP_H
#define _XGENERICSPARSEMATRIXIMP_H

namespace xlinalg
{
//
template <typename T, typename DEFINED, typename PATTERN, typename STORAGE_TYPE, typename INDEXING>
template <typename G>
xGenericSparseMatrix<T, DEFINED, PATTERN, STORAGE_TYPE, INDEXING>::xGenericSparseMatrix(G &graph)
    : n(0), m(0), nnz(0), index1(nullptr), index2(nullptr), data(nullptr), connected_status(false)
{
   // check if graph is compatible with PATTERN
   if (xPolicyGenericSparseMatrix<PATTERN, STORAGE_TYPE, INDEXING>::checkGraph(graph.isSym()))
   {
      std::ostringstream oss;
      oss << " Graph storage and Pattern of this container are incompatible !\n";
      throw xGenericSparseMatrixException(oss.str(), __FILE__, __LINE__, __DATE__, __TIME__);
   }

   // initialize matrix structure whith graph
   initMemoryWithGraph(graph);

   return;
}

template <typename T, typename DEFINED, typename PATTERN, typename STORAGE_TYPE, typename INDEXING>
template <typename TO, typename DEFINEDO, typename PATTERNO, typename STORAGE_TYPEO, typename INDEXINGO>
xGenericSparseMatrix<T, DEFINED, PATTERN, STORAGE_TYPE, INDEXING>::xGenericSparseMatrix(
    const xGenericSparseMatrix<TO, DEFINEDO, PATTERNO, STORAGE_TYPEO, INDEXINGO> &other, int *sel_n, int *sel_m, bool pack)
    : n(0), m(0), nnz(0), index1(nullptr), index2(nullptr), data(nullptr), connected_status(false)
{
   // local
   int i, j, k, l, ke, le, nn, nm;
   xPolicyGenericSparseMatrixCopyBase *pp;
   T val;

   // use a non const reference as we know what we are doing here. Says we know that we are not going to
   // modify other
   xGenericSparseMatrix<TO, DEFINEDO, PATTERNO, STORAGE_TYPEO, INDEXINGO> &other_nc =
       const_cast<xGenericSparseMatrix<TO, DEFINEDO, PATTERNO, STORAGE_TYPEO, INDEXINGO> &>(other);

   // Connecting to "other"
   int on = other_nc.n;
   int om = other_nc.m;
   int onnz = other_nc.nnz;
   int *oindex1 = other_nc.index1;
   int *oindex2 = other_nc.index2;
   T *odata = other_nc.data;

   // find dimension of new matrix and generate packed index
   nn = 0;
   nm = 0;
   std::vector<int> sel_n_loc;
   std::vector<int> id_pack_n;
   if (sel_n)
   {
      sel_n_loc.resize(on, 0);
      id_pack_n.resize(on + 1, -1);
      int *p = &id_pack_n[1];
      if (pack)
      {
         for (i = 0; i < on; ++i)
         {
            if (sel_n[i] > 0)
            {
               sel_n_loc[i] = i + 1;
               p[i] = nn;
               ++nn;
            }
            else if (sel_n[i] < 0)
            {
               assert(sel_n[-sel_n[i] - 1] > 0);
               sel_n_loc[i] = -sel_n[i];
            }
         }
      }
      else
      {
         for (i = 0; i < on; ++i)
         {
            p[i] = i;
            if (sel_n[i] > 0)
            {
               sel_n_loc[i] = i + 1;
            }
            else if (sel_n[i] < 0)
            {
               assert(sel_n[-sel_n[i] - 1] > 0);
               sel_n_loc[i] = -sel_n[i];
            }
         }
         nn = on;
      }
   }
   else
      nn = on;
   std::vector<int> sel_m_loc;
   std::vector<int> id_pack_m;
   if (sel_m)
   {
      sel_m_loc.resize(om, 0);
      id_pack_m.resize(om + 1, -1);
      int *p = &id_pack_m[1];
      if (pack)
      {
         for (j = 0; j < om; ++j)
         {
            if (sel_m[j] > 0)
            {
               sel_m_loc[j] = j + 1;
               p[j] = nm;
               ++nm;
            }
            else if (sel_m[j] < 0)
            {
               assert(sel_m[-sel_m[j] - 1] > 0);
               sel_m_loc[j] = -sel_m[j];
            }
         }
      }
      else
      {
         for (j = 0; j < om; ++j)
         {
            p[j] = j;
            if (sel_m[j] > 0)
            {
               sel_m_loc[j] = j + 1;
            }
            else if (sel_m[j] < 0)
            {
               assert(sel_m[-sel_m[j] - 1] > 0);
               sel_m_loc[j] = -sel_m[j];
            }
         }
         nm = om;
      }
   }
   else
      nm = om;

   // generate apropriate policy
   if (sel_n && sel_m)
      pp = new xPolicyGenericSparseMatrixCopy<PATTERN, PATTERNO, true, true>(sel_n_loc.data(), sel_m_loc.data(), &id_pack_n[0],
                                                                             &id_pack_m[0]);
   else if (sel_n)
      pp = new xPolicyGenericSparseMatrixCopy<PATTERN, PATTERNO, true, false>(sel_n_loc.data(), sel_m, &id_pack_n[0], nullptr);
   else if (sel_m)
      pp = new xPolicyGenericSparseMatrixCopy<PATTERN, PATTERNO, false, true>(sel_n, sel_m_loc.data(), nullptr, &id_pack_m[0]);
   else
      pp = new xPolicyGenericSparseMatrixCopy<PATTERN, PATTERNO, false, false>(sel_n, sel_m, nullptr, nullptr);

   // create a graph to store new non null term of matrix
   xGraphMatrix g(nn, nm, xTraitsGenericSparseMatrixPattern<PATTERN>::symGraph());

   // fill graph
   ke = xPolicyGenericSparseMatrix<PATTERNO, STORAGE_TYPEO, INDEXINGO>::firstIndexEnd(on, om, onnz, oindex1, oindex2);
   for (k = 0; k < ke; ++k)
   {
      l = xPolicyGenericSparseMatrix<PATTERNO, STORAGE_TYPEO, INDEXINGO>::secondIndexStart(k);
      le = xPolicyGenericSparseMatrix<PATTERNO, STORAGE_TYPEO, INDEXINGO>::secondIndexEnd(k);
      for (; l < le; ++l)
      {
         i = xPolicyGenericSparseMatrix<PATTERNO, STORAGE_TYPEO, INDEXINGO>::indexI(k, l) - 1;
         j = xPolicyGenericSparseMatrix<PATTERNO, STORAGE_TYPEO, INDEXINGO>::indexJ(k, l) - 1;
         pp->addGraph(g, i, j);
      }
   }
   g.countNNZ();

   // initialize matrix structure whith graph
   initMemoryWithGraph(g);

   // fill matrix terms
   ke = xPolicyGenericSparseMatrix<PATTERNO, STORAGE_TYPEO, INDEXINGO>::firstIndexEnd(on, om, onnz, oindex1, oindex2);
   for (k = 0; k < ke; ++k)
   {
      l = xPolicyGenericSparseMatrix<PATTERNO, STORAGE_TYPEO, INDEXINGO>::secondIndexStart(k);
      le = xPolicyGenericSparseMatrix<PATTERNO, STORAGE_TYPEO, INDEXINGO>::secondIndexEnd(k);
      for (; l < le; ++l)
      {
         i = xPolicyGenericSparseMatrix<PATTERNO, STORAGE_TYPEO, INDEXINGO>::indexI(k, l) - 1;
         j = xPolicyGenericSparseMatrix<PATTERNO, STORAGE_TYPEO, INDEXINGO>::indexJ(k, l) - 1;
         val = *(
             static_cast<T *>(&(((TO *)odata)[xPolicyGenericSparseMatrix<PATTERNO, STORAGE_TYPEO, INDEXINGO>::indexVal(k, l)])));
         pp->addMatrix(i, j,
                       std::bind(&xGenericSparseMatrix<T, DEFINED, PATTERN, STORAGE_TYPE, INDEXING>::AddMatrix, this,
                                 std::placeholders::_1, std::placeholders::_2, val));
      }
   }

   if (pp) delete pp;

   return;
}

template <typename T, typename DEFINED, typename PATTERN, typename STORAGE_TYPE, typename INDEXING>
xGenericSparseMatrix<T, DEFINED, PATTERN, STORAGE_TYPE, INDEXING>::xGenericSparseMatrix(
    const xGenericSparseMatrix<T, DEFINED, PATTERN, STORAGE_TYPE, INDEXING> &other, bool trans)
    : n(0), m(0), nnz(0), index1(nullptr), index2(nullptr), data(nullptr), connected_status(false)
{
   // local
   int i, j, k, l, ke, le, nn, nm;

   // use a non const reference as we know what we are doing here. Says we know that we are not going to
   // modify other
   xGenericSparseMatrix<T, DEFINED, PATTERN, STORAGE_TYPE, INDEXING> &other_nc =
       const_cast<xGenericSparseMatrix<T, DEFINED, PATTERN, STORAGE_TYPE, INDEXING> &>(other);

   // Connecting to "other"
   int on = other_nc.n;
   int om = other_nc.m;
   int onnz = other_nc.nnz;
   int *oindex1 = other_nc.index1;
   int *oindex2 = other_nc.index2;
   T *odata = other_nc.data;

   // set dimension of new matrix
   if (trans)
   {
      nn = om;
      nm = on;
   }
   else
   {
      nn = on;
      nm = om;
   }

   // create a graph to store new non null term of matrix
   xGraphMatrix g(nn, nm, xTraitsGenericSparseMatrixPattern<PATTERN>::symGraph());

   // fill graph
   ke = xPolicyGenericSparseMatrix<PATTERN, STORAGE_TYPE, INDEXING>::firstIndexEnd(on, om, onnz, oindex1, oindex2);
   if (trans)
   {
      for (k = 0; k < ke; ++k)
      {
         l = xPolicyGenericSparseMatrix<PATTERN, STORAGE_TYPE, INDEXING>::secondIndexStart(k);
         le = xPolicyGenericSparseMatrix<PATTERN, STORAGE_TYPE, INDEXING>::secondIndexEnd(k);
         for (; l < le; ++l)
         {
            j = xPolicyGenericSparseMatrix<PATTERN, STORAGE_TYPE, INDEXING>::indexI(k, l) - 1;
            i = xPolicyGenericSparseMatrix<PATTERN, STORAGE_TYPE, INDEXING>::indexJ(k, l) - 1;
            g.add(i, j);
         }
      }
   }
   else
   {
      for (k = 0; k < ke; ++k)
      {
         l = xPolicyGenericSparseMatrix<PATTERN, STORAGE_TYPE, INDEXING>::secondIndexStart(k);
         le = xPolicyGenericSparseMatrix<PATTERN, STORAGE_TYPE, INDEXING>::secondIndexEnd(k);
         for (; l < le; ++l)
         {
            i = xPolicyGenericSparseMatrix<PATTERN, STORAGE_TYPE, INDEXING>::indexI(k, l) - 1;
            j = xPolicyGenericSparseMatrix<PATTERN, STORAGE_TYPE, INDEXING>::indexJ(k, l) - 1;
            g.add(i, j);
         }
      }
   }
   g.countNNZ();

   // initialize matrix structure whith graph
   initMemoryWithGraph(g);

   // fill matrix terms
   ke = xPolicyGenericSparseMatrix<PATTERN, STORAGE_TYPE, INDEXING>::firstIndexEnd(on, om, onnz, oindex1, oindex2);
   if (trans)
   {
      for (k = 0; k < ke; ++k)
      {
         l = xPolicyGenericSparseMatrix<PATTERN, STORAGE_TYPE, INDEXING>::secondIndexStart(k);
         le = xPolicyGenericSparseMatrix<PATTERN, STORAGE_TYPE, INDEXING>::secondIndexEnd(k);
         for (; l < le; ++l)
         {
            j = xPolicyGenericSparseMatrix<PATTERN, STORAGE_TYPE, INDEXING>::indexI(k, l);
            i = xPolicyGenericSparseMatrix<PATTERN, STORAGE_TYPE, INDEXING>::indexJ(k, l);
            AddMatrix(i, j, odata[xPolicyGenericSparseMatrix<PATTERN, STORAGE_TYPE, INDEXING>::indexVal(k, l)]);
         }
      }
   }
   else
   {
      for (k = 0; k < ke; ++k)
      {
         l = xPolicyGenericSparseMatrix<PATTERN, STORAGE_TYPE, INDEXING>::secondIndexStart(k);
         le = xPolicyGenericSparseMatrix<PATTERN, STORAGE_TYPE, INDEXING>::secondIndexEnd(k);
         for (; l < le; ++l)
         {
            i = xPolicyGenericSparseMatrix<PATTERN, STORAGE_TYPE, INDEXING>::indexI(k, l);
            j = xPolicyGenericSparseMatrix<PATTERN, STORAGE_TYPE, INDEXING>::indexJ(k, l);
            AddMatrix(i, j, odata[xPolicyGenericSparseMatrix<PATTERN, STORAGE_TYPE, INDEXING>::indexVal(k, l)]);
         }
      }
   }

   return;
}

template <typename T, typename DEFINED, typename PATTERN, typename STORAGE_TYPE, typename INDEXING>
xGenericSparseMatrix<T, DEFINED, PATTERN, STORAGE_TYPE, INDEXING>::~xGenericSparseMatrix()
{
   // in case of probleme  this destructor may be called with unintialyzed memory => checking
   if (index1 != nullptr) delete[] index1;
   if (index2 != nullptr) delete[] index2;
   if (data != nullptr) delete[] data;

   /*
      if (connected_status)
      std::cout << "Warning : This containner have been connected with a solver or something else at least one time. Check that
      you don't use this other means now ! "  << std::endl;
    */
}

template <typename T, typename DEFINED, typename PATTERN, typename STORAGE_TYPE, typename INDEXING>
template <typename G>
void xGenericSparseMatrix<T, DEFINED, PATTERN, STORAGE_TYPE, INDEXING>::initMemoryWithGraph(G &g)
{
   // local
   int k, l, err = 0, offset = 0;
   int *idx;

   // set dimension
   n = g.getN();
   m = g.getM();
   nnz = g.getNNZ();

   // allocating memory
   try
   {
      xPolicyGenericSparseMatrix<PATTERN, STORAGE_TYPE, INDEXING>::allocMemory(n, m, nnz, &index1, &index2, extra1);
      data = new T[nnz];
   }
   catch (std::bad_alloc &e)
   {
      std::cout << "Standard exception: " << e.what() << std::endl;
      std::ostringstream oss;
      oss << " Space for matrix  was not allocated normaly !\n";
      throw xGenericSparseMatrixException(oss.str(), __FILE__, __LINE__, __DATE__, __TIME__);
   }

   // initialise value to zero as allocation doesn't
   resetMatrixToZero();

   // set matrix structure by following graph structure
   // loop on colonnes
   for (k = 0; k < m; ++k)
   {
      idx = g.getCol(k);
      l = g.getSizeCol(k);
      if ((err = xPolicyGenericSparseMatrix<PATTERN, STORAGE_TYPE, INDEXING>::setStorage(k, l, idx, index1, index2, &extra1[0],
                                                                                         offset)))
         break;
   }

   // check
   if (err)
   {
      std::ostringstream oss;
      oss << " Defining matrix structure from computed matrix graph object failled.\n";
      throw xGenericSparseMatrixException(oss.str(), __FILE__, __LINE__, __DATE__, __TIME__);
   }

   // second pass if needed
   if (xPolicyGenericSparseMatrix<PATTERN, STORAGE_TYPE, INDEXING>::twoPass)
   {
      offset = n;
      for (k = 0; k < m; ++k)
      {
         idx = g.getCol(k);
         l = g.getSizeCol(k);
         if ((err = xPolicyGenericSparseMatrix<PATTERN, STORAGE_TYPE, INDEXING>::setStoragePassTwo(k, l, idx, index1, index2,
                                                                                                   &extra1[0], offset)))
            break;
      }

      // check
      if (err)
      {
         std::ostringstream oss;
         oss << " Defining matrix structure from computed matrix graph object failled.\n";
         throw xGenericSparseMatrixException(oss.str(), __FILE__, __LINE__, __DATE__, __TIME__);
      }
   }

   // finalyse storage creation
   xPolicyGenericSparseMatrix<PATTERN, STORAGE_TYPE, INDEXING>::endSetStorage(n, m, index1, index2, &extra1[0], offset);

   return;
}
template <typename T, typename DEFINED, typename PATTERN, typename STORAGE_TYPE, typename INDEXING>
void xGenericSparseMatrix<T, DEFINED, PATTERN, STORAGE_TYPE, INDEXING>::getMemoryAccess(int &ns, int &ms, int &nnzs, int **idx1,
                                                                                        int **idx2, T **dat,
                                                                                        STORAGE_TYPE &storage, INDEXING &indexing)
{
   // check that it is not already connected
   if (connected_status)
      std::cout << "Warning : This containner is already connected with a solver or something else ! " << std::endl;

   // set dimension
   ns = n;
   ms = m;
   nnzs = nnz;

   // set pointeur
   *idx1 = index1;
   *idx2 = index2;
   *dat = data;

   // set status
   connected_status = true;

   return;
}

template <typename T, typename DEFINED, typename PATTERN, typename STORAGE_TYPE, typename INDEXING>
void xGenericSparseMatrix<T, DEFINED, PATTERN, STORAGE_TYPE, INDEXING>::AddMatrix(const int &i, const int &j, const T &val)
{
   // fortran numbering by convention
   //#pragma omp critical
   ((T *)data)[xPolicyGenericSparseMatrix<PATTERN, STORAGE_TYPE, INDEXING>::getIndexFromFortran(i, j, index1, index2,
                                                                                                &extra1[0])] += val;
}

template <typename T, typename DEFINED, typename PATTERN, typename STORAGE_TYPE, typename INDEXING>
T &xGenericSparseMatrix<T, DEFINED, PATTERN, STORAGE_TYPE, INDEXING>::operator()(const int &i, const int &j)
{
   return ((T *)data)[xPolicyGenericSparseMatrix<PATTERN, STORAGE_TYPE, INDEXING>::getIndex(i, j, index1, index2, &extra1[0])];
}

template <typename T, typename DEFINED, typename PATTERN, typename STORAGE_TYPE, typename INDEXING>
T xGenericSparseMatrix<T, DEFINED, PATTERN, STORAGE_TYPE, INDEXING>::operator()(const int &i, const int &j) const
{
   return ((T *)data)[xPolicyGenericSparseMatrix<PATTERN, STORAGE_TYPE, INDEXING>::getIndex(i, j, index1, index2, &extra1[0])];
}

template <typename T, typename DEFINED, typename PATTERN, typename STORAGE_TYPE, typename INDEXING>
void xGenericSparseMatrix<T, DEFINED, PATTERN, STORAGE_TYPE, INDEXING>::resetMatrixToZero()
{
   T v_zero = xtool::xDataType<T>::zero();
   for (int j = 0; j < nnz; ++j)
   {
      data[j] = v_zero;
   }
}
template <typename T, typename DEFINED, typename PATTERN, typename STORAGE_TYPE, typename INDEXING>
void xGenericSparseMatrix<T, DEFINED, PATTERN, STORAGE_TYPE, INDEXING>::saveState()
{
   if (stored_data.size() == 0)
   {
      stored_data.resize(nnz);
   }
   memcpy(&stored_data[0], data, nnz * sizeof(T));
}
template <typename T, typename DEFINED, typename PATTERN, typename STORAGE_TYPE, typename INDEXING>
void xGenericSparseMatrix<T, DEFINED, PATTERN, STORAGE_TYPE, INDEXING>::restorState()
{
   if (stored_data.size() == 0)
   {
      std::ostringstream oss;
      oss << "Storage of matrix value was not done as vector is of null size ! You must call at least once saveState before "
             "calling this methode";
      throw xGenericSparseMatrixException(oss.str(), __FILE__, __LINE__, __DATE__, __TIME__);
   }
   memcpy(data, &stored_data[0], nnz * sizeof(T));
}
template <typename T, typename DEFINED, typename PATTERN, typename STORAGE_TYPE, typename INDEXING>
template <typename... VARG>
int xGenericSparseMatrix<T, DEFINED, PATTERN, STORAGE_TYPE, INDEXING>::gemvDispatch(const T &alpha, const T *X, T *Y,
                                                                                    bool sparse_X, bool transpose, VARG &... args)
{
   //===============PLEASE READ THIS BEFORE ANY CHANGE IN THIS method====================
   // As noted bellow std:function has been banished from this method implementation
   //
   // By itself std:function is perfect to hold lambda function set by a if : std::function.... bar; if something bar=a lambda
   // otherwise bar=an other lambda But it have an intrinsic cost when used. If the lambda computation cost a lot it is invisible.
   // But when lambda computation do not cost anything like it is the case with prod and test_v_zero, the std:function  intrinsic
   // cost comes into play.
   //
   // As prod and test_v_zero are used rather intensively in algo keeping direct use of lambda alleviate any unnecessary extra
   // cost. In this way, compiler optimization may even remove completely prod() and test_v_zero() when they corresponds to
   // nothing...
   //
   //====================================================================================

   // error code
   int err;

   const T v_zero = xtool::xDataType<T>::zero();
   const T v_un = xtool::xDataType<T>::one();

   // if matrix product have to be scaled
   if (alpha != v_un)
   {
      // Note: please do not use std::function to simplify code. Here direct use of lambda prevent performance loss. That may be
      // huge ...
      auto prod = [&alpha](T &val) -> void { val *= alpha; };

      // if X is considered as sparse (add X value testing)
      if (sparse_X)
      {
         // Note: please do not use std::function to simplify code. Here direct use of lambda prevent performance loss. That may
         // be huge ...
         auto test_v_zero = [&v_zero](T &val) -> bool { return val == v_zero; };

         // if matrix is transposed
         if (transpose)
         {
            err = xPolicyGenericSparseMatrixGemv<T, PATTERN, STORAGE_TYPE, INDEXING>::gemv1Transposed(
                n, m, data, X, Y, index1, index2, &extra1[0], prod, test_v_zero, args...);
         }
         // if matrix is not transposed
         else
         {
            err = xPolicyGenericSparseMatrixGemv<T, PATTERN, STORAGE_TYPE, INDEXING>::gemv1(
                n, m, data, X, Y, index1, index2, &extra1[0], prod, test_v_zero, args...);
         }
      }
      // if X is considered as dense (no  X value testing)
      else
      {
         // Note: please do not use std::function to simplify code. Here direct use of lambda prevent performance loss. That may
         // be huge ...
         auto test_v_zero = [](T &val) -> bool { return false; };

         // if matrix is transposed
         if (transpose)
         {
            err = xPolicyGenericSparseMatrixGemv<T, PATTERN, STORAGE_TYPE, INDEXING>::gemv1Transposed(
                n, m, data, X, Y, index1, index2, &extra1[0], prod, test_v_zero, args...);
         }
         // if matrix is not transposed
         else
         {
            err = xPolicyGenericSparseMatrixGemv<T, PATTERN, STORAGE_TYPE, INDEXING>::gemv1(
                n, m, data, X, Y, index1, index2, &extra1[0], prod, test_v_zero, args...);
         }
      }
   }
   // if matrix product is not scaled
   else
   {
      // Note: please do not use std::function to simplify code. Here direct use of lambda prevent performance loss. That may be
      // huge ...
      auto prod = [](T &val) -> void { ; };

      // if X is considered as sparse (add X value testing)
      if (sparse_X)
      {
         // Note: please do not use std::function to simplify code. Here direct use of lambda prevent performance loss. That may
         // be huge ...
         auto test_v_zero = [&v_zero](T &val) -> bool { return val == v_zero; };

         // if matrix is transposed
         if (transpose)
         {
            err = xPolicyGenericSparseMatrixGemv<T, PATTERN, STORAGE_TYPE, INDEXING>::gemv1Transposed(
                n, m, data, X, Y, index1, index2, &extra1[0], prod, test_v_zero, args...);
         }
         // if matrix is not transposed
         else
         {
            err = xPolicyGenericSparseMatrixGemv<T, PATTERN, STORAGE_TYPE, INDEXING>::gemv1(
                n, m, data, X, Y, index1, index2, &extra1[0], prod, test_v_zero, args...);
         }
      }
      // if X is considered as dense (no  X value testing)
      else
      {
         // Note: please do not use std::function to simplify code. Here direct use of lambda prevent performance loss. That may
         // be huge ...
         auto test_v_zero = [](T &val) -> bool { return false; };

         // if matrix is transposed
         if (transpose)
         {
            err = xPolicyGenericSparseMatrixGemv<T, PATTERN, STORAGE_TYPE, INDEXING>::gemv1Transposed(
                n, m, data, X, Y, index1, index2, &extra1[0], prod, test_v_zero, args...);
         }
         // if matrix is not transposed
         else
         {
            err = xPolicyGenericSparseMatrixGemv<T, PATTERN, STORAGE_TYPE, INDEXING>::gemv1(
                n, m, data, X, Y, index1, index2, &extra1[0], prod, test_v_zero, args...);
         }
      }
   }

   return err;
}
template <typename T, typename DEFINED, typename PATTERN, typename STORAGE_TYPE, typename INDEXING>
void xGenericSparseMatrix<T, DEFINED, PATTERN, STORAGE_TYPE, INDEXING>::gemv(int size_X, int size_Y, T alpha, T beta, T *X, T *Y,
                                                                             bool sparse_X, bool transpose)
{
   // dimension checking
   if (transpose)
   {
      if (size_X < n || size_Y < m)
      {
         std::ostringstream oss;
         oss << "Dimension are not correct\nX is of size " << size_X << " and should be >= " << n << "\n";
         oss << "Y is of size " << size_Y << " and should be >=" << m << "\n";
         oss << "Matrix is considered transposed in gemv\n";
         throw xGenericSparseMatrixException(oss.str(), __FILE__, __LINE__, __DATE__, __TIME__);
      }
   }
   else
   {
      if (size_X < m || size_Y < n)
      {
         std::ostringstream oss;
         oss << "Dimension are not correct\nX is of size " << size_X << " and should be >= " << m << "\n";
         oss << "Y is of size " << size_Y << " and should be >=" << n << "\n";
         throw xGenericSparseMatrixException(oss.str(), __FILE__, __LINE__, __DATE__, __TIME__);
      }
   }

   // error code
   int err;

   // zero and one of type T
   const T v_zero = xtool::xDataType<T>::zero();
   const T v_un = xtool::xDataType<T>::one();

   // if beta is not one Y have to be scaled => blas
   if (beta != v_un)
   {
      xlinalg::xCPPBlasDef<T>::scal(size_Y, beta, Y);
   }

   // if alpha is null no matrix product to do
   if (alpha == v_zero) return;

   // do the computation via the dispatcher
   err = gemvDispatch(alpha, X, Y, sparse_X, transpose);

   if (err < 0)
   {
      std::ostringstream oss;
      oss << "Problem while doing multiplication of sparse matrix";
      throw xGenericSparseMatrixException(oss.str(), __FILE__, __LINE__, __DATE__, __TIME__);
   }
}

template <typename T, typename DEFINED, typename PATTERN, typename STORAGE_TYPE, typename INDEXING>
void xGenericSparseMatrix<T, DEFINED, PATTERN, STORAGE_TYPE, INDEXING>::gemv(T alpha, T beta, const xDistVector<T> &X,
                                                                             xDistVector<T> &Y, bool sparse_X, bool transpose)
{
   const xDistIndex &idx_X = X.getDistIndex();
   const xDistIndex &idx_Y = Y.getDistIndex();

   const int size_X = idx_X.getGlobalIndexSize();
   const int size_Y = idx_Y.getGlobalIndexSize();

   // dimension checking
   if (transpose)
   {
      if (size_X < n || size_Y < m)
      {
         std::ostringstream oss;
         oss << "Dimension are not correct\nX is of size " << size_X << " and should be >= " << n << "\n";
         oss << "Y is of size " << size_Y << " and should be >=" << m << "\n";
         oss << "Matrix is considered transposed in gemv\n";
         throw xGenericSparseMatrixException(oss.str(), __FILE__, __LINE__, __DATE__, __TIME__);
      }
   }
   else
   {
      if (size_X < m || size_Y < n)
      {
         std::ostringstream oss;
         oss << "Dimension are not correct\nX is of size " << size_X << " and should be >= " << m << "\n";
         oss << "Y is of size " << size_Y << " and should be >=" << n << "\n";
         throw xGenericSparseMatrixException(oss.str(), __FILE__, __LINE__, __DATE__, __TIME__);
      }
   }

   // error code
   int err;

   // zero and one of type T
   T v_zero = xtool::xDataType<T>::zero();
   T v_un = xtool::xDataType<T>::one();

   // Y distributed vector status
   bool old_Y_status = Y.isInInsertModeOn();
   if (!old_Y_status) Y.switchInsertModeOn();

   // if beta is not one Y have to be scaled
   if (beta != v_un)
   {
      Y.switchInsertModeOff();
      Y.scal(beta);
      Y.switchInsertModeOn();
   }

   // if alpha is null no matrix product to do
   if (alpha == v_zero) return;

   // X distributed vector status
   xDistVector<T> &X_ = const_cast<xDistVector<T> &>(X);
   bool old_X_status = X.isInInsertModeOn();
   if (old_X_status) X_.switchInsertModeOff();
   // insure correct computation: Matrix term may be partially defined on proc but associated X term must be the sum of all
   // contribution
   X_.switchToGlobalValue();

   // do the computation via the dispatcher
   err = gemvDispatch(alpha, X.data(), Y.data(), sparse_X, transpose, idx_X, idx_Y);

   if (err < 0)
   {
      std::ostringstream oss;
      oss << "Problem while doing multiplication of sparce matrix";
      throw xGenericSparseMatrixException(oss.str(), __FILE__, __LINE__, __DATE__, __TIME__);
   }

   // reset state
   if (!old_Y_status) Y.switchInsertModeOff();
   if (old_X_status) X_.switchInsertModeOn();
}

template <typename T, typename DEFINED, typename PATTERN, typename STORAGE_TYPE, typename INDEXING>
void xGenericSparseMatrix<T, DEFINED, PATTERN, STORAGE_TYPE, INDEXING>::gemm(char *trans_B, int nb_row_OPB, int nb_col_OPB,
                                                                             int nb_row_C, int ldb, int ldc, T alpha, T beta,
                                                                             T *B, T *C)
{
   // dimension checking
   if (nb_row_OPB < m || nb_row_C != n || ldc < n)
   {
      std::ostringstream oss;
      oss << "Dimension are not correct\nop(B) is of size " << nb_row_OPB << "x" << nb_col_OPB << " and should be " << m << "x"
          << nb_col_OPB << "\n";
      oss << "C is of size " << nb_row_C << "x" << nb_col_OPB << " and should be " << n << "x" << nb_col_OPB << "\n";
      oss << "ldc is " << ldc << " and should be >= to " << n << "\n";
      throw xGenericSparseMatrixException(oss.str(), __FILE__, __LINE__, __DATE__, __TIME__);
   }

   if (*trans_B == 'N')
   {
      // dimension checking
      if (ldb < m)
      {
         std::ostringstream oss;
         oss << "Dimension are not correct\nop(B) is of size " << nb_row_OPB << "x" << nb_col_OPB << " and should be " << m << "x"
             << nb_col_OPB << "\n";
         oss << "C is of size " << nb_row_C << "x" << nb_col_OPB << " and should be " << n << "x" << nb_col_OPB << "\n";
         oss << "ldc is " << ldc << " and should be >= to " << n << "\n";
         oss << "ldb is " << ldb << " and should be >= to " << m << "\n";
         throw xGenericSparseMatrixException(oss.str(), __FILE__, __LINE__, __DATE__, __TIME__);
      }

      // loop on column of B,C to do the sparce product using gemv
      // nota : far from optimal but give a quick access to matrix x matrix operation
      for (int i = 0; i < nb_col_OPB; ++i)
      {
         gemv(nb_row_OPB, nb_row_C, alpha, beta, &B[i * ldb], &C[i * ldc]);
      }
   }
   else
   {
      throw xGenericSparseMatrixException("to be done op(B)=Bt", __FILE__, __LINE__, __DATE__, __TIME__);
   }
}

template <typename T, typename DEFINED, typename PATTERN, typename STORAGE_TYPE, typename INDEXING>
template <typename DEFINEDB, typename PATTERNB, typename STORAGE_TYPEB, typename INDEXINGB, typename DEFINEDC, typename PATTERNC,
          typename STORAGE_TYPEC, typename INDEXINGC>
void xGenericSparseMatrix<T, DEFINED, PATTERN, STORAGE_TYPE, INDEXING>::gemm(
    T alpha, T beta, const xGenericSparseMatrix<T, DEFINEDB, PATTERNB, STORAGE_TYPEB, INDEXINGB> &B,
    xGenericSparseMatrix<T, DEFINEDC, PATTERNC, STORAGE_TYPEC, INDEXINGC> &C)
{
   // local
   int i, j, k, l, ke, le;
   int idx, nb_color;

   // zero and one of type T
   T v_zero = xtool::xDataType<T>::zero();
   T v_un = xtool::xDataType<T>::one();

   // if alpha is null no product
   if (alpha == v_zero)
   {
      if (beta != v_un) C.scaleMatrix(beta);
      return;
   }
   // warpe rhs nature,alpha,sparseX and transpose parameter into a lambda to symplify A.X calling sequence
   std::function<int(T * X, T * Y)> prod = [this, &alpha](T *X, T *Y) -> int {
      return this->gemvDispatch(alpha, X, Y, true, false);
   };

   // Connecting to "B"
   int Bn = B.n;
   int Bm = B.m;
   int Bnnz = B.nnz;
   int *Bindex1 = B.index1;
   int *Bindex2 = B.index2;

   // Connecting to "C"
   int Cn = C.n;
   int Cm = C.m;
   int Cnnz = C.nnz;
   int *Cindex1 = C.index1;
   int *Cindex2 = C.index2;
   T *Cdata = C.data;

   // check dimension
   if (n != Cn)
      throw xGenericSparseMatrixException("gemm: C matrix must have the same number of row as A (instance calling method)",
                                          __FILE__, __LINE__, __DATE__, __TIME__);
   if (Bm != Cm)
      throw xGenericSparseMatrixException("gemm: C matrix must have the same number of column as B", __FILE__, __LINE__, __DATE__,
                                          __TIME__);
   if (m != Bn)
      throw xGenericSparseMatrixException(
          "gemm: A (instance calling method) matrix must have a  number of column equal to the number of row of B", __FILE__,
          __LINE__, __DATE__, __TIME__);

   // create A graph corresponding to A pattern
   xGraphMatrix Ag(n, m, xTraitsGenericSparseMatrixPattern<PATTERN>::symGraph());
   // create B graph corresponding to B pattern
   xGraphMatrix Bg(Bn, Bm, xTraitsGenericSparseMatrixPattern<PATTERNB>::symGraph());
   // create C graph corresponding to pattern corresponding to product and addition (if beta not null)
   xGraphMatrix Cg(Cn, Cm, xTraitsGenericSparseMatrixPattern<PATTERNC>::symGraph());

   // fill A graph
   ke = xPolicyGenericSparseMatrix<PATTERN, STORAGE_TYPE, INDEXING>::firstIndexEnd(n, m, nnz, index1, index2);
   for (k = 0; k < ke; ++k)
   {
      l = xPolicyGenericSparseMatrix<PATTERN, STORAGE_TYPE, INDEXING>::secondIndexStart(k);
      le = xPolicyGenericSparseMatrix<PATTERN, STORAGE_TYPE, INDEXING>::secondIndexEnd(k);
      for (; l < le; ++l)
      {
         i = xPolicyGenericSparseMatrix<PATTERN, STORAGE_TYPE, INDEXING>::indexI(k, l) - 1;
         j = xPolicyGenericSparseMatrix<PATTERN, STORAGE_TYPE, INDEXING>::indexJ(k, l) - 1;
         Ag.add(i, j);
      }
   }
   Ag.countNNZ();
   // fill B graph
   ke = xPolicyGenericSparseMatrix<PATTERNB, STORAGE_TYPEB, INDEXINGB>::firstIndexEnd(Bn, Bm, Bnnz, Bindex1, Bindex2);
   for (k = 0; k < ke; ++k)
   {
      l = xPolicyGenericSparseMatrix<PATTERNB, STORAGE_TYPEB, INDEXINGB>::secondIndexStart(k);
      le = xPolicyGenericSparseMatrix<PATTERNB, STORAGE_TYPEB, INDEXINGB>::secondIndexEnd(k);
      for (; l < le; ++l)
      {
         i = xPolicyGenericSparseMatrix<PATTERNB, STORAGE_TYPEB, INDEXINGB>::indexI(k, l) - 1;
         j = xPolicyGenericSparseMatrix<PATTERNB, STORAGE_TYPEB, INDEXINGB>::indexJ(k, l) - 1;
         Bg.add(i, j);
      }
   }
   Bg.countNNZ();

   // color container
   xGraphMatrix color_g(Cn, Cm);
   std::vector<std::vector<int>> color(Bm);
   nb_color = 0;

   // do symbolic phase
   // =================

   if (xTraitsGenericSparseMatrixPattern<PATTERN>::symGraph() && xTraitsGenericSparseMatrixPattern<PATTERNB>::symGraph())
   {
      throw xGenericSparseMatrixException("gemm: A and B matrix are symmetric, not coded yet", __FILE__, __LINE__, __DATE__,
                                          __TIME__);
   }
   else if (xTraitsGenericSparseMatrixPattern<PATTERN>::symGraph())
   {
      throw xGenericSparseMatrixException("gemm: A matrix is symmetric, not coded yet", __FILE__, __LINE__, __DATE__, __TIME__);
   }
   else if (xTraitsGenericSparseMatrixPattern<PATTERNB>::symGraph())
   {
      throw xGenericSparseMatrixException("gemm: B matrix is symmetric, not coded yet", __FILE__, __LINE__, __DATE__, __TIME__);
   }
   else
   {
      // loop on column of B
      for (j = 0; j < Bm; ++j)
      {
         int *idxJB = Bg.getCol(j);
         int nbJB = Bg.getSizeCol(j);
         // loop on row of jth B column : correspond to column of A
         for (k = 0; k < nbJB; ++k)
         {
            idx = idxJB[k];
            int *idxJA = Ag.getCol(idx);
            int nbJA = Ag.getSizeCol(idx);
            // add all terms of A column to C as corresponding B term is not null
            Cg.addLinesUnsymBlock(nbJA, j, idxJA);
            // feed current potential color
            color_g.addLinesUnsymBlock(nbJA, nb_color, idxJA);
         }
         int *idxJcc = color_g.getCol(nb_color);
         int nbJcc = color_g.getSizeCol(nb_color);
         bool new_color = true;
         if (nbJcc)
         {
            int min_cc = idxJcc[0];
            int max_cc = idxJcc[nbJcc - 1];

            // loop on previous color to see if current column is a new color or not
            for (l = 0; l < nb_color; ++l)
            {
               int *idxJc = color_g.getCol(l);
               int nbJc = color_g.getSizeCol(l);
               assert(nbJc > 0);
               if (idxJc[0] > max_cc || idxJc[nbJc - 1] < min_cc)
                  new_color = false;
               else
               {
                  new_color = false;
                  for (i = 0, k = 0; i < nbJcc && k < nbJc;)
                  {
                     if (idxJcc[i] < idxJc[k])
                        ++i;
                     else if (idxJcc[i] > idxJc[k])
                        ++k;
                     else
                     {
                        new_color = true;
                        break;
                     }
                  }
               }
               if (!new_color)
               {
                  color_g.addLinesUnsymBlock(nbJcc, l, idxJcc);
                  color_g.clearCol(nb_color);
                  color[l].push_back(j);
                  break;
               }
            }
         }
         else
         {
            new_color = false;
         }
         // if new color keep current has new and increment number of color
         if (new_color)
         {
            color[nb_color].push_back(j);
            ++nb_color;
         }
      }
   }

   // add to C graph original C pattern if beta not null
   if (beta != v_zero)
   {
      ke = xPolicyGenericSparseMatrix<PATTERNC, STORAGE_TYPEC, INDEXINGC>::firstIndexEnd(Cn, Cm, Cnnz, Cindex1, Cindex2);
      for (k = 0; k < ke; ++k)
      {
         l = xPolicyGenericSparseMatrix<PATTERNC, STORAGE_TYPEC, INDEXINGC>::secondIndexStart(k);
         le = xPolicyGenericSparseMatrix<PATTERNC, STORAGE_TYPEC, INDEXINGC>::secondIndexEnd(k);
         for (; l < le; ++l)
         {
            i = xPolicyGenericSparseMatrix<PATTERNC, STORAGE_TYPEC, INDEXINGC>::indexI(k, l) - 1;
            j = xPolicyGenericSparseMatrix<PATTERNC, STORAGE_TYPEC, INDEXINGC>::indexJ(k, l) - 1;
            Cg.add(i, j);
         }
      }
   }

   Cg.countNNZ();

   // initialize C matrix structure with graph
   xGenericSparseMatrix<T, DEFINEDC, PATTERNC, STORAGE_TYPEC, INDEXINGC> newC(Cg);

   // do numeric phase
   // =================
   //
   // loop on color to pack B column, compute and unpack result into C
   // NOTE : to simplify here we use gemv1 existing policy. We treat column by column.
   // With advanced policies we would have a matrix vector operation => L2 blas but all packing would
   // have to be done before product
   std::vector<T> X(Bn);
   std::vector<T> Y(n, v_zero);
   int err = 0;
   for (l = 0; l < nb_color; ++l)
   {
      std::fill(X.begin(), X.end(), v_zero);
      // packing B column corresponding to  current color l
      for (auto j : color[l])
      {
         int *idxJB = Bg.getCol(j);
         int nbJB = Bg.getSizeCol(j);
         int ji = xPolicyGenericSparseMatrixIndexing<INDEXINGB>::CtoI(j);
         for (k = 0; k < nbJB; ++k)
         {
            i = idxJB[k];
            X[i] = B(xPolicyGenericSparseMatrixIndexing<INDEXINGB>::CtoI(i), ji);
         }
      }
      // do sparse product
      err = prod(X.data(), Y.data());
      if (err < 0)
      {
         std::ostringstream oss;
         oss << "Problem while doing multiplication of sparse matrix : color " << k << " ";
         throw xGenericSparseMatrixException(oss.str(), __FILE__, __LINE__, __DATE__, __TIME__);
      }
      // unpack into C
      for (auto j : color[l])
      {
         int *idxJC = Cg.getCol(j);
         int nbJC = Cg.getSizeCol(j);
         int ji = xPolicyGenericSparseMatrixIndexing<INDEXINGC>::CtoI(j);
         for (k = 0; k < nbJC; ++k)
         {
            i = idxJC[k];
            newC(xPolicyGenericSparseMatrixIndexing<INDEXINGC>::CtoI(i), ji) = Y[i];
            Y[i] = v_zero;
         }
      }
   }

   // add to newC original C value
   if (beta != v_zero)
   {
      if (beta == v_un)
      {
         ke = xPolicyGenericSparseMatrix<PATTERNC, STORAGE_TYPEC, INDEXINGC>::firstIndexEnd(Cn, Cm, Cnnz, Cindex1, Cindex2);
         for (k = 0; k < ke; ++k)
         {
            l = xPolicyGenericSparseMatrix<PATTERNC, STORAGE_TYPEC, INDEXINGC>::secondIndexStart(k);
            le = xPolicyGenericSparseMatrix<PATTERNC, STORAGE_TYPEC, INDEXINGC>::secondIndexEnd(k);
            for (; l < le; ++l)
            {
               i = xPolicyGenericSparseMatrix<PATTERNC, STORAGE_TYPEC, INDEXINGC>::indexI(k, l);
               j = xPolicyGenericSparseMatrix<PATTERNC, STORAGE_TYPEC, INDEXINGC>::indexJ(k, l);
               newC.AddMatrix(i, j, Cdata[xPolicyGenericSparseMatrix<PATTERNC, STORAGE_TYPEC, INDEXINGC>::indexVal(k, l)]);
            }
         }
      }
      else
      {
         ke = xPolicyGenericSparseMatrix<PATTERNC, STORAGE_TYPEC, INDEXINGC>::firstIndexEnd(Cn, Cm, Cnnz, Cindex1, Cindex2);
         for (k = 0; k < ke; ++k)
         {
            l = xPolicyGenericSparseMatrix<PATTERNC, STORAGE_TYPEC, INDEXINGC>::secondIndexStart(k);
            le = xPolicyGenericSparseMatrix<PATTERNC, STORAGE_TYPEC, INDEXINGC>::secondIndexEnd(k);
            for (; l < le; ++l)
            {
               i = xPolicyGenericSparseMatrix<PATTERNC, STORAGE_TYPEC, INDEXINGC>::indexI(k, l);
               j = xPolicyGenericSparseMatrix<PATTERNC, STORAGE_TYPEC, INDEXINGC>::indexJ(k, l);
               newC.AddMatrix(i, j, Cdata[xPolicyGenericSparseMatrix<PATTERNC, STORAGE_TYPEC, INDEXINGC>::indexVal(k, l)] * beta);
            }
         }
      }
   }

   // scrach old C  by new one
   if (C.index1 != nullptr) delete[] C.index1;
   C.index1 = newC.index1;
   newC.index1 = nullptr;

   if (C.index2 != nullptr) delete[] C.index2;
   C.index2 = newC.index2;
   newC.index2 = nullptr;

   if (C.data != nullptr) delete[] C.data;
   C.data = newC.data;
   newC.data = nullptr;

   C.extra1 = newC.extra1;
   C.nnz = newC.nnz;

   return;
}

template <typename T, typename DEFINED, typename PATTERN, typename STORAGE_TYPE, typename INDEXING>
template <typename DEFINEDB, typename PATTERNB, typename STORAGE_TYPEB, typename INDEXINGB, typename DEFINEDC, typename PATTERNC,
          typename STORAGE_TYPEC, typename INDEXINGC>
void xGenericSparseMatrix<T, DEFINED, PATTERN, STORAGE_TYPE, INDEXING>::AoperB(
    int offset_r, int offset_c, const xGenericSparseMatrix<T, DEFINEDB, PATTERNB, STORAGE_TYPEB, INDEXINGB> &B,
    xGenericSparseMatrix<T, DEFINEDC, PATTERNC, STORAGE_TYPEC, INDEXINGC> &C, std::function<T(const T &, const T &)> oper)
{
   // local
   int i, j, jo, k, l, ke, le;

   // zero and one of type T
   T v_zero = xtool::xDataType<T>::zero();

   // Connecting to "B"
   int Bn = B.n;
   int Bm = B.m;
   int Bnnz = B.nnz;
   int *Bindex1 = B.index1;
   int *Bindex2 = B.index2;

   // Connecting to "C"
   int Cn = C.n;
   int Cm = C.m;

   // check dimension
   if (n != Cn)
      throw xGenericSparseMatrixException("AoperB: C matrix must have the same number of row as A (instance calling method)",
                                          __FILE__, __LINE__, __DATE__, __TIME__);
   if (m != Cm)
      throw xGenericSparseMatrixException("AoperB: C matrix must have the same number of column as A (instance calling method)",
                                          __FILE__, __LINE__, __DATE__, __TIME__);

   if (offset_r + Bn > Cn)
      throw xGenericSparseMatrixException(
          "AoperB: B matrix must be placed with offset_r so that it is not merging with row outside  A (instance calling method)",
          __FILE__, __LINE__, __DATE__, __TIME__);
   if (offset_c + Bm > Cm)
      throw xGenericSparseMatrixException(
          "AoperB: B matrix must be placed with offset_c so that it is not merging with column outside  A (instance calling "
          "method)",
          __FILE__, __LINE__, __DATE__, __TIME__);

   // create A graph corresponding to A pattern
   xGraphMatrix Ag(n, m, xTraitsGenericSparseMatrixPattern<PATTERN>::symGraph());
   // create B graph corresponding to B pattern
   xGraphMatrix Bg(Bn, Bm, xTraitsGenericSparseMatrixPattern<PATTERNB>::symGraph());
   // create C graph corresponding to pattern corresponding to formal addition
   xGraphMatrix Cg(Cn, Cm, xTraitsGenericSparseMatrixPattern<PATTERNC>::symGraph());

   // fill A graph
   ke = xPolicyGenericSparseMatrix<PATTERN, STORAGE_TYPE, INDEXING>::firstIndexEnd(n, m, nnz, index1, index2);
   for (k = 0; k < ke; ++k)
   {
      l = xPolicyGenericSparseMatrix<PATTERN, STORAGE_TYPE, INDEXING>::secondIndexStart(k);
      le = xPolicyGenericSparseMatrix<PATTERN, STORAGE_TYPE, INDEXING>::secondIndexEnd(k);
      for (; l < le; ++l)
      {
         i = xPolicyGenericSparseMatrix<PATTERN, STORAGE_TYPE, INDEXING>::indexI(k, l) - 1;
         j = xPolicyGenericSparseMatrix<PATTERN, STORAGE_TYPE, INDEXING>::indexJ(k, l) - 1;
         Ag.add(i, j);
      }
   }
   Ag.countNNZ();
   // fill B graph
   ke = xPolicyGenericSparseMatrix<PATTERNB, STORAGE_TYPEB, INDEXINGB>::firstIndexEnd(Bn, Bm, Bnnz, Bindex1, Bindex2);
   for (k = 0; k < ke; ++k)
   {
      l = xPolicyGenericSparseMatrix<PATTERNB, STORAGE_TYPEB, INDEXINGB>::secondIndexStart(k);
      le = xPolicyGenericSparseMatrix<PATTERNB, STORAGE_TYPEB, INDEXINGB>::secondIndexEnd(k);
      for (; l < le; ++l)
      {
         i = xPolicyGenericSparseMatrix<PATTERNB, STORAGE_TYPEB, INDEXINGB>::indexI(k, l) - 1;
         j = xPolicyGenericSparseMatrix<PATTERNB, STORAGE_TYPEB, INDEXINGB>::indexJ(k, l) - 1;
         Bg.add(i, j);
      }
   }
   Bg.countNNZ();

   // do symbolic phase
   // =================
   std::vector<int> idxJB_offset(Bn);
   if (xTraitsGenericSparseMatrixPattern<PATTERN>::symGraph() && xTraitsGenericSparseMatrixPattern<PATTERNB>::symGraph())
   {
      throw xGenericSparseMatrixException("AoperB: A and B matrix are symmetric, not coded yet", __FILE__, __LINE__, __DATE__,
                                          __TIME__);
   }
   else if (xTraitsGenericSparseMatrixPattern<PATTERN>::symGraph())
   {
      throw xGenericSparseMatrixException("AoperB: A matrix is symmetric, not coded yet", __FILE__, __LINE__, __DATE__, __TIME__);
   }
   else if (xTraitsGenericSparseMatrixPattern<PATTERNB>::symGraph())
   {
      throw xGenericSparseMatrixException("AoperB: B matrix is symmetric, not coded yet", __FILE__, __LINE__, __DATE__, __TIME__);
   }
   else
   {
      // loop on column of A
      for (j = 0; j < m; ++j)
      {
         int *idxJA = Ag.getCol(j);
         int nbJA = Ag.getSizeCol(j);
         // add all terms of A column to C
         if (nbJA) Cg.addLinesUnsymBlock(nbJA, j, idxJA);
      }
      // loop on column of B
      for (j = 0; j < Bm; ++j)
      {
         jo = j + offset_c;
         int *idxJB = Bg.getCol(j);
         int nbJB = Bg.getSizeCol(j);
         if (nbJB)
         {
            // offset row index
            std::transform(idxJB, idxJB + nbJB, idxJB_offset.data(), [&offset_r](const int i) { return i + offset_r; });
            // add all terms of B column to C
            Cg.addLinesUnsymBlock(nbJB, jo, idxJB_offset.data());
         }
      }
   }

   Cg.countNNZ();

   // initialize C matrix structure with graph
   xGenericSparseMatrix<T, DEFINEDC, PATTERNC, STORAGE_TYPEC, INDEXINGC> newC(Cg);

   // do numeric phase
   // =================
   //
   // loop on column of A/C
   xGenericSparseMatrix<T, DEFINED, PATTERN, STORAGE_TYPE, INDEXING> &A = *this;
   for (j = 0; j < m; ++j)
   {
      int *idxJA = Ag.getCol(j);
      int nbJA = Ag.getSizeCol(j);
      jo = j - offset_c;
      int jA = xPolicyGenericSparseMatrixIndexing<INDEXING>::CtoI(j);
      int jB = xPolicyGenericSparseMatrixIndexing<INDEXINGB>::CtoI(jo);
      int jC = xPolicyGenericSparseMatrixIndexing<INDEXINGC>::CtoI(j);

      // if B is not present B is considered as zero
      if (!(jo < Bm) || (jo < 0))
      {
         for (i = 0; i < nbJA; ++i)
         {
            k = idxJA[i];
            newC(xPolicyGenericSparseMatrixIndexing<INDEXINGC>::CtoI(k), jC) =
                oper(A(xPolicyGenericSparseMatrixIndexing<INDEXING>::CtoI(k), jA), v_zero);
         }
      }
      // if a mix of A and B is present
      else
      {
         int ia, ib;
         int *idxJB = Bg.getCol(jo);
         int nbJB = Bg.getSizeCol(jo);
         if (nbJB) std::transform(idxJB, idxJB + nbJB, idxJB_offset.data(), [&offset_r](const int i) { return i + offset_r; });
         for (i = 0, k = 0; i < nbJA && k < nbJB;)
         {
            ia = idxJA[i];
            ib = idxJB_offset[k];

            // one non null term in A and one null term in B at (ia,j) location
            if (ia < ib)
            {
               newC(xPolicyGenericSparseMatrixIndexing<INDEXINGC>::CtoI(ia), jC) =
                   oper(A(xPolicyGenericSparseMatrixIndexing<INDEXING>::CtoI(ia), jA), v_zero);
               ++i;
            }
            // one null term in A and one non null term in B at (ib,j) location
            else if (ia > ib)
            {
               newC(xPolicyGenericSparseMatrixIndexing<INDEXINGC>::CtoI(ib), jC) =
                   oper(v_zero, B(xPolicyGenericSparseMatrixIndexing<INDEXINGB>::CtoI(idxJB[k]), jB));
               ++k;
            }
            // one non null term in A and one non null term in B at (ia/ib,j) location
            else
            {
               newC(xPolicyGenericSparseMatrixIndexing<INDEXINGC>::CtoI(ib), jC) =
                   oper(A(xPolicyGenericSparseMatrixIndexing<INDEXING>::CtoI(ia), jA),
                        B(xPolicyGenericSparseMatrixIndexing<INDEXINGB>::CtoI(idxJB[k]), jB));
               ++i;
               ++k;
            }
         }
         for (; i < nbJA;)
         {
            ia = idxJA[i];

            // one non null term in A and one null term in B at (ia,j) location
            newC(xPolicyGenericSparseMatrixIndexing<INDEXINGC>::CtoI(ia), jC) =
                oper(A(xPolicyGenericSparseMatrixIndexing<INDEXING>::CtoI(ia), jA), v_zero);
            ++i;
         }
         for (; k < nbJB;)
         {
            ib = idxJB_offset[k];

            // one null term in A and one non null term in B at (ib,j) location
            newC(xPolicyGenericSparseMatrixIndexing<INDEXINGC>::CtoI(ib), jC) =
                oper(v_zero, B(xPolicyGenericSparseMatrixIndexing<INDEXINGB>::CtoI(idxJB[k]), jB));
            ++k;
         }
      }
   }

   // scrach old C  by new one
   if (C.index1 != nullptr) delete[] C.index1;
   C.index1 = newC.index1;
   newC.index1 = nullptr;

   if (C.index2 != nullptr) delete[] C.index2;
   C.index2 = newC.index2;
   newC.index2 = nullptr;

   if (C.data != nullptr) delete[] C.data;
   C.data = newC.data;
   newC.data = nullptr;

   C.extra1 = newC.extra1;
   C.nnz = newC.nnz;

   return;
}

template <typename T, typename DEFINED, typename PATTERN, typename STORAGE_TYPE, typename INDEXING>
template <typename S>
void xGenericSparseMatrix<T, DEFINED, PATTERN, STORAGE_TYPE, INDEXING>::gemvfac(S &solver, int ncol, double alpha, double beta,
                                                                                xCSRVector &X, xCSRVector &Y, xCSRVector &Z)
{
   if (ncol > 1)
   {
      throw xGenericSparseMatrixException("to be done ncol>1", __FILE__, __LINE__, __DATE__, __TIME__);
   }
   else
   {
      // doing Y = alpha.K-1.X + beta.Y is equivalent to solve K.Z=X and add alpha.Z to beta.Y

      // zero and one of type T
      // nota : T is use insteade of double to break compilation process if gemvfac is use with a other type as double : not
      // implemented nor thinked
      T v_zero = xtool::xDataType<T>::zero();
      T v_un = xtool::xDataType<T>::one();

      // Y = alpha.K-1.X + beta.Y
      // if beta = 0 it reduced to  Y = alpha.K-1.X
      if (beta == v_zero)
      {
         // calculate directely in Y
         solver.solve(X, Y);

         // scale by alpha if not unit
         if (alpha != v_un) xlinalg::xCPPBlasDef<T>::scal(n, alpha, Y.GetArray());
      }
      else
      {
         T *y = Y.GetArray();

         // first find Z
         solver.solve(X, Z);

         // second scale y if beta not unit
         if (beta != v_un) xlinalg::xCPPBlasDef<T>::scal(n, beta, y);

         // third use axpy blas
         xlinalg::xCPPBlasDef<T>::axpy(n, alpha, Z.GetArray(), y);
      }
   }
}
template <typename T, typename DEFINED, typename PATTERN, typename STORAGE_TYPE, typename INDEXING>
void xGenericSparseMatrix<T, DEFINED, PATTERN, STORAGE_TYPE, INDEXING>::printMatrixMarket(std::ostream &os)
{
   // see for example http://math.nist.gov/MatrixMarket/formats.html to understand loop below
   std::map<std::string, std::string> xDataTypeToMatrixMarketType{
       {"real_single", "real"}, {"real_double", "real"}, {"complex_single", "complex"}, {"complex_double", "complex"}};

   auto itMMsType = xDataTypeToMatrixMarketType.find(xtool::xDataType<T>::stype());
   if (itMMsType == xDataTypeToMatrixMarketType.end())
   {
      std::cout << xtool::xDataType<T>::stype() << " equivalent type for MatrixMarket format not defined.\n";
      std::cout << "Update printMatrixMarket method !\n" << std::flush;
      throw;
   }

   // write banner
   os << "%%MatrixMarket matrix coordinate " << itMMsType->second << " " << xTraitsGenericSparseMatrixPattern<PATTERN>::spatern()
      << std::endl;
   os << "% generated by xFem xGenericSparseMatrix container" << std::endl;
   os << "% Sparse matrix" << std::endl;
   // write data line
   os << n << " " << m << " " << nnz << std::endl;
   // following data line
   // loop on n,m dimension folowing storage and pattern policy
   forEach([&os](const int i, const int j, const T &val) { os << i << " " << j << " " << val << std::endl; });
}

template <typename T, typename DEFINED, typename PATTERN, typename STORAGE_TYPE, typename INDEXING>
void xGenericSparseMatrix<T, DEFINED, PATTERN, STORAGE_TYPE, INDEXING>::fillZeroDiag(const T &val, const int norm)
{
   const bool debug = false;

   T *p = (T *)data;
   T val_norm;
   T v_zero = xtool::xDataType<T>::zero();

   switch (norm)
   {
      case 3:
      {
         val_norm = normMax();
         break;
      }
      default:
      {
         std::ostringstream oss;
         oss << " Only norm 3 implemented so far !\n";
         throw xGenericSparseMatrixException(oss.str(), __FILE__, __LINE__, __DATE__, __TIME__);
      }
   }

   // automatic thresold calculation
   // nota : comes from MUMPS with nom=2 (i.e. infinit norme)
   T treshold = ((T)(1.e-5)) * val_norm;
   treshold *= xtool::xDataType<T>::epsilonNumeric();
   if (debug) std::cout << std::scientific << " ||A||= " << val_norm << " thresold=" << treshold << std::endl;

   // corective value
   val_norm *= val;

   // loop on diag termes
   int k;
   xPolicyGenericSparseMatrix<PATTERN, STORAGE_TYPE, INDEXING>::initDiagSearch(n, m, nnz, index1, index2, &extra1[0]);
   for (k = 0; k < m; ++k)
   {
      const int l = xPolicyGenericSparseMatrix<PATTERN, STORAGE_TYPE, INDEXING>::diagIndex(k);
      if (l < 0)
      {
         std::ostringstream oss;
         oss << " No diagonal term (" << k << "," << k << ") !\n";
         throw xGenericSparseMatrixException(oss.str(), __FILE__, __LINE__, __DATE__, __TIME__);
      }
      const T diag_val = p[l];
      if ((diag_val < v_zero && diag_val > -treshold) || (diag_val < treshold))
      {
         p[l] = val_norm;
         if (debug) std::cout << "Diagonal terme " << k << " was modified from " << diag_val << " to " << val_norm << std::endl;
      }
   }
}

template <typename T, typename DEFINED, typename PATTERN, typename STORAGE_TYPE, typename INDEXING>
T xGenericSparseMatrix<T, DEFINED, PATTERN, STORAGE_TYPE, INDEXING>::normMax()
{
   T v_zero = xtool::xDataType<T>::zero();

   int nmax = xlinalg::xCPPBlasDef<T>::iamax(n, data);

   if (!nmax)
   {
      std::ostringstream oss;
      oss << " It is imposible to find the maximum absolute value of A !\n";
      throw xGenericSparseMatrixException(oss.str(), __FILE__, __LINE__, __DATE__, __TIME__);
   }

   T val_norm = ((T *)data)[nmax - 1];
   if (val_norm < v_zero) val_norm = -val_norm;

   return val_norm;
}

template <typename T, typename DEFINED, typename PATTERN, typename STORAGE_TYPE, typename INDEXING>
void xGenericSparseMatrix<T, DEFINED, PATTERN, STORAGE_TYPE, INDEXING>::scaleMatrix(T scale)
{
   if (scale != xtool::xDataType<T>::one()) xlinalg::xCPPBlasDef<T>::scal(nnz, scale, data);
}

template <typename T, typename DEFINED, typename PATTERN, typename STORAGE_TYPE, typename INDEXING>
void xGenericSparseMatrix<T, DEFINED, PATTERN, STORAGE_TYPE, INDEXING>::operator+=(
    const xGenericSparseMatrix<T, DEFINED, PATTERN, STORAGE_TYPE, INDEXING> &rhs)
{
   // test graph ! (index1, 2)
   xlinalg::xCPPBlasDef<T>::axpy(nnz, xtool::xDataType<T>::one(), rhs.data, data);
}

template <typename T, typename DEFINED, typename PATTERN, typename STORAGE_TYPE, typename INDEXING>
void xGenericSparseMatrix<T, DEFINED, PATTERN, STORAGE_TYPE, INDEXING>::toDense(int nbr, int nbc, T *dense_matrix, int *ids_r,
                                                                                int *ids_c,
                                                                                std::function<int(int, int)> denseIndex)
{
   int i, j, k, l, ke, le, idxd;
   T val;

   // reset to zero
   std::fill(dense_matrix, dense_matrix + nbr * nbc, xtool::xDataType<T>::zero());

   // generate apropriate policy
   xPolicyDenseMatrixCopyBase *pp{nullptr};
   if (ids_r && ids_c)
      pp = new xPolicyDenseMatrixCopy<PATTERN, true, true>(nbr, nbc, ids_r, ids_c, denseIndex);
   else if (ids_r)
      pp = new xPolicyDenseMatrixCopy<PATTERN, true, false>(nbr, nbc, ids_r, nullptr, denseIndex);
   else if (ids_c)
      pp = new xPolicyDenseMatrixCopy<PATTERN, false, true>(nbr, nbc, nullptr, ids_c, denseIndex);
   else
      pp = new xPolicyDenseMatrixCopy<PATTERN, false, false>(nbr, nbc, nullptr, nullptr, denseIndex);

   // fill matrix terms
   ke = xPolicyGenericSparseMatrix<PATTERN, STORAGE_TYPE, INDEXING>::firstIndexEnd(n, m, nnz, index1, index2);
   for (k = 0; k < ke; ++k)
   {
      l = xPolicyGenericSparseMatrix<PATTERN, STORAGE_TYPE, INDEXING>::secondIndexStart(k);
      le = xPolicyGenericSparseMatrix<PATTERN, STORAGE_TYPE, INDEXING>::secondIndexEnd(k);
      for (; l < le; ++l)
      {
         i = xPolicyGenericSparseMatrix<PATTERN, STORAGE_TYPE, INDEXING>::indexI(k, l);
         j = xPolicyGenericSparseMatrix<PATTERN, STORAGE_TYPE, INDEXING>::indexJ(k, l);
         val = *(static_cast<T *>(&(((T *)data)[xPolicyGenericSparseMatrix<PATTERN, STORAGE_TYPE, INDEXING>::indexVal(k, l)])));
         idxd = pp->denseSelectedIndex(i, j);
         if (idxd > -1) dense_matrix[idxd] = val;
      }
   }

   delete pp;
}
template <typename T, typename DEFINED, typename PATTERN, typename STORAGE_TYPE, typename INDEXING>
template <typename F>
void xGenericSparseMatrix<T, DEFINED, PATTERN, STORAGE_TYPE, INDEXING>::forEach(const F &func)
{
   int k, l, ke, le;
   ke = xPolicyGenericSparseMatrix<PATTERN, STORAGE_TYPE, INDEXING>::firstIndexEnd(n, m, nnz, index1, index2);
   for (k = 0; k < ke; ++k)
   {
      l = xPolicyGenericSparseMatrix<PATTERN, STORAGE_TYPE, INDEXING>::secondIndexStart(k);
      le = xPolicyGenericSparseMatrix<PATTERN, STORAGE_TYPE, INDEXING>::secondIndexEnd(k);
      for (; l < le; ++l)
      {
         func(xPolicyGenericSparseMatrix<PATTERN, STORAGE_TYPE, INDEXING>::indexI(k, l),
              xPolicyGenericSparseMatrix<PATTERN, STORAGE_TYPE, INDEXING>::indexJ(k, l),
              data[xPolicyGenericSparseMatrix<PATTERN, STORAGE_TYPE, INDEXING>::indexVal(k, l)]);
      }
   }
}

template <typename T, typename DEFINED, typename PATTERN, typename STORAGE_TYPE, typename INDEXING>
int xGenericSparseMatrix<T, DEFINED, PATTERN, STORAGE_TYPE, INDEXING>::getN() const
{
   return n;
}
template <typename T, typename DEFINED, typename PATTERN, typename STORAGE_TYPE, typename INDEXING>
int xGenericSparseMatrix<T, DEFINED, PATTERN, STORAGE_TYPE, INDEXING>::getM() const
{
   return m;
}
template <typename T, typename DEFINED, typename PATTERN, typename STORAGE_TYPE, typename INDEXING>
int xGenericSparseMatrix<T, DEFINED, PATTERN, STORAGE_TYPE, INDEXING>::getNNZ() const
{
   return nnz;
}

}  // namespace xlinalg
#endif
