/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>

namespace xlinalg
{
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// xLinearSystemSolverSuperLUBase class implementation
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// xLinearSystemSolverSuperLU class implementation
///////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template <typename T>
xLinearSystemSolverSuperLU<T>::xLinearSystemSolverSuperLU()
    : status(0),
      keep_data(nullptr),
      keep_orig(nullptr),
      lwork(0),
      work(nullptr),
      etree(nullptr),
      R(nullptr),
      C(nullptr),
      ferr(nullptr),
      berr(nullptr),
      rpg(0.0),
      rcond(0.0)
{
   equed[0] = ' ';
   // set default parameter
   setParam(Equil, 0, YES);
   setParam(ColPerm, 0, COLAMD);
   setParam(DiagPivotThresh, 0, 1.0);
   setParam(Trans, 0, NOTRANS);
   setParam(IterRefine, 0, NOREFINE);
   setParam(SymmetricMode, 0, NO);
   setParam(PivotGrowth, 0, NO);
   setParam(ConditionNumber, 0, NO);
   setParam(PrintStat, 0, NO);
}

template <typename T>
xLinearSystemSolverSuperLU<T>::~xLinearSystemSolverSuperLU()
{
   // memory have to be cleaned in interface solver space
   deallocateMemoryInterface();
   cleanFactor();
}
template <typename T>
void xLinearSystemSolverSuperLU<T>::allocateMemoryInterface(const int &n, const int &m, const int &nnz)
{
   // local
   int i;

   try
   {
      etree = new int[m];
      R = new data_sub_type[n];
      C = new data_sub_type[m];
      perm_r = new int[n];
      if (ext_perm_c)
      {
         perm_c = ext_perm_c;
         ext_perm_c = nullptr;
         previous_ext_perm_c = true;
      }
      else
      {
         perm_c = new int[m];
         for (i = 0; i < m; ++i) perm_c[i] = 0;
      }
      if (keep) keep_data = new T[nnz];
   }
   catch (std::bad_alloc &e)
   {
      std::cout << "Standard exception: " << e.what() << std::endl;
      std::ostringstream oss;
      oss << " SuperLU interface data space was not allocated normaly ?!\n";
      throw xLinearSystemSolverSuperLUException(oss.str(), __FILE__, __LINE__, __DATE__, __TIME__);
   }

   // set to zero
   if (n < m)
   {
      for (i = 0; i < n; ++i)
      {
         etree[i] = perm_r[i] = 0;
         R[i] = C[i] = (data_sub_type)(0.);
      }
      for (; i < m; ++i)
      {
         etree[i] = 0;
         C[i] = (data_sub_type)(0.);
      }
   }
   else
   {
      for (i = 0; i < m; ++i)
      {
         etree[i] = perm_r[i] = 0;
         R[i] = C[i] = (data_sub_type)(0.);
      }
      for (; i < n; ++i)
      {
         perm_r[i] = 0;
         R[i] = (data_sub_type)(0.);
      }
   }

   // to fullfil generic driver interface RHS and SOL have to be given with consistant data => init here to empty
   // with zero colone (no solve asked)
   setPointerSuperMatrixDense<T>(RHS, n, 0, n, nullptr);
   setPointerSuperMatrixDense<T>(SOL, m, 0, m, nullptr);
}

template <typename T>
void xLinearSystemSolverSuperLU<T>::deallocateMemoryInterface()
{
   try
   {
      if (etree) delete[] etree;
      if (R) delete[] R;
      if (C) delete[] C;
      if (ferr) delete[] ferr;
      if (berr) delete[] berr;
      if (perm_r) delete[] perm_r;
      if (perm_c)
      {
         if (!previous_ext_perm_c)
         {
            delete[] perm_c;
         }
         perm_c = nullptr;
      }
      etree = nullptr;
      R = nullptr;
      C = nullptr;
      ferr = nullptr;
      berr = nullptr;
      perm_r = nullptr;
      if (keep_data)
      {
         delete[] keep_data;
         keep_data = nullptr;
      }
   }
   catch (std::bad_alloc &e)
   {
      std::cout << "Standard exception: " << e.what() << std::endl;
      std::ostringstream oss;
      oss << " SuperLU interface data space was not deallocated normaly ?!\n";
      throw xLinearSystemSolverSuperLUException(oss.str(), __FILE__, __LINE__, __DATE__, __TIME__);
   }

   previous_ext_perm_c = false;
}

template <typename T>
template <typename M>
void xLinearSystemSolverSuperLU<T>::init()
{
   // test sym from traits of matrix
   // this is done to lock what can be done with SuperLU interface
   // it's not used, as for now, SuperLU is only dedicated to nonsingulare non symetric storage  matrix
   // Matric may be either nonsingular of definite positive, symetrique or not, but it must have a unsymetrique storage or
   // being diagonal
   // interface implementation effort may be done in the future to see if SLU_SYL may be used ....
   // also diagonal dominant traits may be added in the future to use symetric usage  for symbolique phase (SymmetricMode). tmp
   // will have some sense then
   int tmp = xlinalg::xTraitsLinearSystemSolverSuperLUMatrixType<typename M::matrix_pattern, typename M::matrix_defined>::sym;
   tmp++;  // to avoid -Wunused-variable warning
   // xlinalg::xTraitsLinearSystemSolverSuperLUMatrixType < typename M::matrix_pattern, typename M::matrix_defined > :: sym;

   // if SuperLU instance have already been set we are in the case of a re - init
   if (status & XLINEARSOLVERSUPERLU_INIT)
   {
      // if user memory have already been  connected => connection is lost => status lose XLINEARSOLVERSUPERLU_ALLOC
      // done naturaly below
      // memory have also to be cleaned in interface solver space
      //
      if (verbose) std::cout << "Begin reinit SuperLU instance" << std::endl;

      deallocateMemoryInterface();
      cleanFactor();
      cleanStore(A);

      // set default parameter from package
      setDefaultParameter();
      // reinit from previous settings
      reInitParameter();
   }
   else
   {
      if (verbose) std::cout << "Begin init SuperLU instance" << std::endl;
   }

   // allocating interface space
   allocateMemoryInterface(nrows, ncol, nzdata);

   // set status whatever it was
   status = XLINEARSOLVERSUPERLU_INIT;
}

template <typename T>
template <typename M>
void xLinearSystemSolverSuperLU<T>::connectMatrix(M &matrix, int *external_col_perm)
{
   ext_perm_c = external_col_perm;
   connectMatrix(matrix);
}
template <typename T>
template <typename M>
void xLinearSystemSolverSuperLU<T>::connectMatrix(M &matrix)
{
   // temporary local
   int n, m, nnz;
   // temporary local pointeur
   int *index1;
   int *index2;
   T *orig_data, *data;
   matrix_storage storage;
   matrix_indexing indexing;

   // attache "matrix space" to the solver
   matrix.getMemoryAccess(n, m, nnz, &index1, &index2, &orig_data, storage, indexing);

   // save dimentions
   nrows = n;
   ncol = m;
   nzdata = nnz;

   // initialise or re-initialise SuperLu instance
   init<M>();

   // use intermediate container to compute A-1
   if (keep)
   {
      // set computational container to keep_data
      data = keep_data;
      // save original container adress
      keep_orig = orig_data;
   }
   // use original container (attached one) to compute A-1
   else
   {
      data = orig_data;
   }

   // link information in A container
   setPointerSuperMatrixSparse<T, matrix_storage>(A, n, m, nnz, index1, index2, data);

   // update status
   status |= XLINEARSOLVERSUPERLU_ALLOC;

   // verbode stuff
   if (verbose)
   {
      std::cout << "Matrix size " << n << " x " << m << std::endl;
      std::cout << "Number of non zeros " << nnz << std::endl;
   }

   // do by default the symbolic phase
   symb();
}

template <typename T>
void xLinearSystemSolverSuperLU<T>::symb()
{
   if ((status & XLINEARSOLVERSUPERLU_INIT) && (status & XLINEARSOLVERSUPERLU_ALLOC))
   {
      // check
      if (status != (XLINEARSOLVERSUPERLU_SYMB - 1))
      {
         std::ostringstream oss;
         oss << " SuperLU symbolic phase is not intended, for now, to be used without re-initialisation of the SuperLU instance "
                "!\n";
         oss << " Reconnect a matrix will do that ! Otherwise you have to add your own implementation\n";
         throw xLinearSystemSolverSuperLUException(oss.str(), __FILE__, __LINE__, __DATE__, __TIME__);
      }

      // update status
      // nota : we may argue on setting here the status but for now it's the most simple way to folow
      // other interface pattern
      status |= XLINEARSOLVERSUPERLU_SYMB;
   }
   else
   {
      std::ostringstream oss;
      oss << " Symbolic phase is not possible without a previous SuperLU instance initialized and a matrix structure connected "
             "!\n";
      throw xLinearSystemSolverSuperLUException(oss.str(), __FILE__, __LINE__, __DATE__, __TIME__);
   }
}

template <typename T>
void xLinearSystemSolverSuperLU<T>::fact()
{
   // check
   if ((status & XLINEARSOLVERSUPERLU_INIT) && (status & XLINEARSOLVERSUPERLU_ALLOC) && (status & XLINEARSOLVERSUPERLU_SYMB))
   {
      // numeric and others may have already be done >> clean
      // only factor,RHS,SOL and equed have to be cleaned
      // etree, ... are the same has connect/symb is the only way to change matrix structure
      if (status != (XLINEARSOLVERSUPERLU_FACT - 1))
      {
         cleanFactor();
         equed[0] = ' ';

         // to fullfil generic driver interface RHS and SOL have to be given with consistant data => init here to empty
         // with zero colone (no solve asked)
         setPointerSuperMatrixDense<T>(RHS, nrows, 0, nrows, nullptr);
         setPointerSuperMatrixDense<T>(SOL, ncol, 0, ncol, nullptr);

#ifdef DEBUG
         if (!keep)
         {
            std::cout
                << "Warning : Initial matrix attached to this SuperLU instance may have been modified by a previous factorisation"
                << std::endl;
            std::cout << "          Be shur that this new factorisation is based on correct value of A" << std::endl;
         }
#endif
      }

      // if keep option enable :
      //    if in context of first factorisation after conecting matrix : keep_data have to be filled by values of attached matrix
      //    if in context of a new factorisation with same matrix dimension : keep_data have to be filled by values of attached
      //    matrix as previous computation may
      //                                                                      have change keep_data values and user may have
      //                                                                      change value of attached matrix
      //
      if (keep)
      {
         memcpy(keep_data, keep_orig, nzdata * sizeof(T));
      }

      // set option to factorise
      setParam(Fact, 0, DOFACT);

      // Initialize the statistics variables.
      initStat();

      // doing the work
      if (verbose) std::cout << "Begin factorisation with SuperLU instance" << std::endl;

      // do the job
      int info = expertDriver();

      // checking exit status
      if (info)
      {
         std::ostringstream oss;
         oss << " SuperLU factorisation phase failed !\n\tinfo=" << info << std::endl;
         throw xLinearSystemSolverSuperLUException(oss.str(), __FILE__, __LINE__, __DATE__, __TIME__);
      }

      // verbode stuff
      if (verbose)
      {
         printPerf();
         getInternalData(PivotGrowth, 0, info);
         if (info == YES) std::cout << "Reciprocal pivot growth = " << rpg << std::endl;
         getInternalData(ConditionNumber, 0, info);
         if (info == YES) std::cout << "Reciprocal condition number  = " << rcond << std::endl;
      }

#ifdef DEBUG
      // checking reordoring of the matrix
      if ((equed[0] != 'N') && !keep)
      {
         std::cout << "Warning : Initial matrix attached to this SuperLU instance have been modified" << std::endl;
         std::cout << "          To get back to the original matrix use R and/or C according to equed=" << equed[0] << std::endl;
      }
#endif

      // print the statistics variables. (if asked by the user)
      printStat();

      // clear the statistics variables.
      clearStat();

      status |= XLINEARSOLVERSUPERLU_FACT;
   }
   else
   {
      std::ostringstream oss;
      oss << "factorisation phase is not possible without a previous SuperLU matrix structure connected !\n";
      throw xLinearSystemSolverSuperLUException(oss.str(), __FILE__, __LINE__, __DATE__, __TIME__);
   }
}

template <typename T>
template <typename V>
void xLinearSystemSolverSuperLU<T>::solve(const V &B, V &X)
{
   // check
   if ((status & XLINEARSOLVERSUPERLU_INIT) && (status & XLINEARSOLVERSUPERLU_ALLOC) && (status & XLINEARSOLVERSUPERLU_SYMB))
   {
      // local
      int i;

      // set rhs depending on number of colonne
      int nrhs = xPolicyLinearSystemSolverSuperLURHS<T, V>::nbColRhs(B);

      // use a intermediate container Y to keep B unchanged as SuperLU may scale B
      V Y;

      // duplicate B in Y and attache to rhs in SuperLU structure
      T *rhs = (T *)xPolicyLinearSystemSolverSuperLURHS<T, V>::copyBInX(B, Y);
      setPointerSuperMatrixDense<T>(RHS, nrows, nrhs, nrows, rhs);

      // attach X in SuperLu structure
      setPointerSuperMatrixDense<T>(SOL, ncol, nrhs, ncol, xPolicyLinearSystemSolverSuperLURHS<T, V>::getPointer(X));

      // allocate extra vector ferr,berr which depend on number of rhs
      // nota :
      //  after using this methode ferr and berr exist. They are not cleanned to give the user the oportunity to
      //  have a look at them via getData
      try
      {
         // cleaning extra vector
         if (ferr) delete[] ferr;
         if (berr) delete[] berr;
         // allocating vectors
         ferr = new data_sub_type[nrhs];
         berr = new data_sub_type[nrhs];
      }
      catch (std::bad_alloc &e)
      {
         std::cout << "Standard exception: " << e.what() << std::endl;
         std::ostringstream oss;
         oss << " SuperLU extra vector for solve phase were not allocated normaly ?!\n";
         throw xLinearSystemSolverSuperLUException(oss.str(), __FILE__, __LINE__, __DATE__, __TIME__);
      }

      // initialise to zero ferr and berr
      {
         data_sub_type zero = 0.;
         for (i = 0; i < nrhs; ++i) ferr[i] = berr[i] = zero;
      }

      // if already factorized
      if (status & XLINEARSOLVERSUPERLU_FACT)
      {
         setParam(Fact, 0, FACTORED);
         if (verbose) std::cout << "Begin solving with SuperLU , factor already calculated" << std::endl;
      }
      else
      {
         // if keep option enable :
         //    if in context of first factorisation after conecting matrix : keep_data have to be filled by values of attached
         //    matrix no other context possible here
         //
         if (keep)
         {
            memcpy(keep_data, keep_orig, nzdata * sizeof(T));
         }
         setParam(Fact, 0, DOFACT);
         if (verbose) std::cout << "Begin factorisation and solving with SuperLU " << std::endl;
      }

      // Initialize the statistics variables.
      initStat();

      // doing the work
      int info = expertDriver();

      // checking exit status
      if (info)
      {
         std::ostringstream oss;
         if (status & XLINEARSOLVERSUPERLU_FACT)
            oss << " SuperLU solving phase failed !\n";
         else
            oss << "SuperLU factorisation and solving phase failed !\n";
         oss << "  info=" << info;
         throw xLinearSystemSolverSuperLUException(oss.str(), __FILE__, __LINE__, __DATE__, __TIME__);
      }

      // verbode stuff
      if (verbose)
      {
         if (!(status & XLINEARSOLVERSUPERLU_FACT))
         {
            printPerf();
            getInternalData(PivotGrowth, 0, info);
            if (info == YES) std::cout << "Reciprocal pivot growth = " << rpg << std::endl;
            getInternalData(ConditionNumber, 0, info);
            if (info == YES) std::cout << "Reciprocal condition number  = " << rcond << std::endl;
         }
         getInternalData(IterRefine, 0, info);
         if (info != NOREFINE)
         {
            getInternalData(RefineSteps, 0, info);
            printf("Iterative Refinement:\n");
            printf("%8s%8s%16s%16s\n", "rhs", "Steps", "FERR", "BERR");
            for (i = 0; i < nrhs; ++i) printf("%8d%8d%16e%16e\n", i + 1, info, ferr[i], berr[i]);
         }
      }

#ifdef DEBUG
      // checking reordoring of the matrix
      if (!(status & XLINEARSOLVERSUPERLU_FACT) && (equed[0] != 'N') && !keep)
      {
         std::cout << "Warning : Initial matrix attached to this SuperLU instance have been modified" << std::endl;
         std::cout << "          To get back to the original matrix use R and/or C according to equed=" << equed[0] << std::endl;
      }
#endif

      // print the statistics variables. (if asked by the user)
      printStat();

      // clear the statistics variables.
      clearStat();

      // update status
      if (status & XLINEARSOLVERSUPERLU_FACT)
         status |= XLINEARSOLVERSUPERLU_SOLV;
      else
         status |= XLINEARSOLVERSUPERLU_SOLV | XLINEARSOLVERSUPERLU_FACT;
   }
   else
   {
      std::ostringstream oss;
      oss << " solving phase is not possible without a previous matrix connected !\n";
      throw xLinearSystemSolverSuperLUException(oss.str(), __FILE__, __LINE__, __DATE__, __TIME__);
   }
}

template <typename T>
int xLinearSystemSolverSuperLU<T>::expertDriver()
{
   int info = 0;
   gssvx<T>(options, A, perm_c, perm_r, etree, &equed[0], (void *)R, (void *)C, L, U, work, lwork, RHS, SOL, (void *)(&rpg),
            (void *)&rcond, (void *)ferr, (void *)berr,
#if SUPERLU_VERSION == 5
            Glu,
#endif
            mem_usage, stat, &info);
   return info;
}

template <typename T>
template <typename P>
void xLinearSystemSolverSuperLU<T>::getData(int param_tab, int param_index, P &param_ref)
{
   // depending on param_tab
   switch (param_tab)
   {
      case GET_VECTC:
      {
         param_ref = C[param_index];
         break;
      }
      case GET_VECTR:
      {
         param_ref = R[param_index];
         break;
      }
      case GET_FERR:
      {
         param_ref = ferr[param_index];
         break;
      }
      case GET_BERR:
      {
         param_ref = berr[param_index];
         break;
      }
      case GET_RCOND:
      {
         param_ref = rcond;
         break;
      }
      case GET_RPG:
      {
         param_ref = rpg;
         break;
      }
      default:
      {
         std::ostringstream oss;
         oss << "Invocation of getData use wrong param_tab ; see xLinearSolverSuperLU.h enum_param_tab to use a good one\n "
                "Remeber that param_ref is typed\n";
         throw xLinearSystemSolverSuperLUException(oss.str(), __FILE__, __LINE__, __DATE__, __TIME__);
      }
   }
}

}  // namespace xlinalg
