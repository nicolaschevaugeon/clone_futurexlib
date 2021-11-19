/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/

#ifndef _ASSEMBLER_TEST_H
#define _ASSEMBLER_TEST_H

#include <vector>

#include "xAssemblerDispatcher.h"
#include "xAssemblerTraitPolicy.h"
#include "xDataType.h"
#include "xValue.h"

namespace xlinalg
{
class xCSRMatrix;
class xCSRVector;
}  // namespace xlinalg

namespace xfem
{
/*!
the following define the Assembler Concept. The function needed by the Assemble algorithm from xAlgorithm.h
Note that most member function are template, which prevent any of them to be virtual : No dynamic polymorphism allowed !
This way we can assemble any type of local data.
The only alternative to be able to do that was to use some kind of type earasure, for example by writting the interface in term of
std::function <double (int, int) > to have a uniform signature whatever the type of local Matrix/ Vector ...  The problem is that
then the access to a given term was taking too long, and no latter optimisation like block assembly for instance would be possible


class xAssemblerConcept {
public:
 template<typename LOCALMATRIX, typename ITERDOF>
 void assemble_matrix(ITERDOF firsti, ITERDOF lasti, ITERDOF firstj, ITERDOF lastj, const LOCALMATRIX& matr);
 template<typename LOCALMATRIX, typename ITERDOF>
 void assemble_matrix(ITERDOF first, ITERDOF last, const LOCALMATRIX& matr);
 template <typename LOCALVECTOR, typename ITERDOF>
 void assemble_vector(ITERDOF first, ITERDOF last, const LOCALVECTOR& vect);
 template <typename LOCALSCALAR>
 void assemble_scalar(const LOCALSCALAR& scal);
 void setCoeff(const double& c) { coeff = c; }
 };

 Here LOCALMATRIX concept must provides at least :
 A parentheses operator that takes two integers and return the term associated to these arguments in the matrix
 In one xAssembler (xAssemblerLumpedBlockEqu) two extra methods getNbRow() and getNbCol() are also required. To clarified this
 concept we may consider these methods as mandatory ones, for all xAssembler. Or xAssemblerLumpedBlockEqu may be changed to avoid
 use of getNbCol()

 Here LOCALVECTOR concept must provides :
 A [] operator that takes one integer and return the term associated to this argument in the vector

!*/

// Defining here the assemble methode give the possibility to optimize the loop on elementary matrix matr by just dispaching
// needed part of the matr (for sym, only half of the matix,...) This concerne only the "potential" symetric matrix as invocated
// by Assemble with only trial field pointeur value iterator. In case of a receiving global symetric matrix this implicitely act
// as if matr is symetric and it's upper triangle is defined. In case of a receiving global unsymetric matrix this implicitely act
// as if matr is unsymetric and it's fully defined. In case of a receiving global diagonal matrix this implicitely act as if matr
// is at least fully defined on it' diagonal. It's up to xFemMatrix matrix class to fullfill these requirements with for now a
// indexation of a full matrix. If this last point is changed policies and immplementation shoul be reconsidered
//

/*!This template class permits to define the data for the different implementation of the xAssembler concept defined in this file
 * : xAssemblerBasic, xAssembleBasicAndTranspose, xAssemblerLumped, xAssemblerLumpedEqu  !*/
template <class M = xlinalg::xCSRMatrix, class V = xlinalg::xCSRVector, class S = double>
class xAssemblerBaseData
{
  public:
   typedef M matrix_type;
   typedef V vector_type;
   typedef S scalar_type;
   xAssemblerBaseData(M* A, V* b, S* s) : mat(A), vec(b), sca(s), coeff(1.) {}
   void setCoeff(const S& c) { coeff = c; }
   void setTarget(M& A, V& b, S& s)
   {
      mat = &A;
      vec = &b;
      sca = &s;
   }
   void setTarget(M& A, V& b)
   {
      mat = &A;
      vec = &b;
      sca = nullptr;
   }
   void setTarget(M& A)
   {
      mat = &A;
      vec = nullptr;
      sca = nullptr;
   }
   void setTarget(V& b)
   {
      mat = 0;
      vec = &b;
      sca = 0;
   }
   void setTarget(S& s)
   {
      mat = 0;
      vec = 0;
      sca = &s;
   }

  protected:
   M* mat;
   V* vec;
   S* sca;
   S coeff;
};

/*! This is a basic implementation of the xAssembler concept. !*/
template <class M = xlinalg::xCSRMatrix, class V = xlinalg::xCSRVector, class S = double>
class xAssemblerBasic : public xAssemblerBaseData<M, V, S>
{
  public:
   xAssemblerBasic(M& A, V& b, S& s) : xAssemblerBaseData<M, V, S>(&A, &b, &s) {}
   xAssemblerBasic(M& A, V& b) : xAssemblerBaseData<M, V, S>(&A, &b, nullptr) {}
   xAssemblerBasic(M& A) : xAssemblerBaseData<M, V, S>(&A, nullptr, nullptr) {}
   xAssemblerBasic(V& b) : xAssemblerBaseData<M, V, S>(nullptr, &b, nullptr) {}
   xAssemblerBasic(S& s) : xAssemblerBaseData<M, V, S>(0, 0, &s) {}
   xAssemblerBasic() : xAssemblerBaseData<M, V, S>(nullptr, nullptr, nullptr) {}

   /// Functions from the xAssembleConcept
   template <typename LOCALMATRIX, typename ITERDOF>
   void assemble_matrix(ITERDOF firsti, ITERDOF lasti, ITERDOF firstj, ITERDOF lastj, const LOCALMATRIX& matr)
   {
      ITERDOF firstj_bak = firstj;
      for (int i = 0; firsti != lasti; ++firsti, ++i)
      {
         firstj = firstj_bak;
         for (int j = 0; firstj != lastj; ++firstj, ++j)
         {
            auto val = this->coeff * matr(i, j);
            if (!xPolicyAssemblerOnZero<typename M::matrix_assembly_on_zero>::isNull(val))
               xMatVecScaAssembler<M, V, S>::assemble((*firsti)->getState(), (*firstj)->getState(), val, this->mat, this->vec,
                                                      this->sca);
         }
      }
   }

   template <typename LOCALMATRIX, typename ITERDOF>
   void assemble_matrix(ITERDOF first, ITERDOF last, const LOCALMATRIX& matr)
   {
      ITERDOF firsti = first, firstj, lastj;
      int j;
      for (int i = 0; firsti != last; ++firsti, ++i)
      {
         xPolicyAssemblerLoop<typename M::matrix_pattern>::firstForLoop(first, firsti, firstj, i, j);
         lastj = xPolicyAssemblerLoop<typename M::matrix_pattern>::lastForLoop(last, firsti);
         for (; firstj != lastj; ++firstj, ++j)
         {
            auto val = this->coeff * matr(i, j);
            if (!xPolicyAssemblerOnZero<typename M::matrix_assembly_on_zero>::isNull(val))
               xMatVecScaAssembler<M, V, S>::assemble((*firsti)->getState(), (*firstj)->getState(), val, this->mat, this->vec,
                                                      this->sca);
         }
      }
   }

   template <typename LOCALVECTOR, typename ITERDOF>
   void assemble_vector(ITERDOF first, ITERDOF last, const LOCALVECTOR& vect)
   {
      for (int i = 0; first != last; ++first, ++i)
      {
         auto val = vect[i] * this->coeff;
         if (val != xtool::xDataType<S>::zero())
            xMatVecScaAssembler<M, V, S>::assemble((*first)->getState(), val, this->vec, this->sca);
      }
   }

   template <typename LOCALSCALAR>
   void assemble_scalar(const LOCALSCALAR& scal)
   {
      auto val = scal * this->coeff;
      if (val != xtool::xDataType<S>::zero()) xMatVecScaAssembler<M, V, S>::assemble(val, this->sca);
   }
   /// End  Function for the xAssembleConcept
};

template <class M = xlinalg::xCSRMatrix, class V = xlinalg::xCSRVector, class S = double>
class xAssemblerBasicAndTranspose : public xAssemblerBaseData<M, V, S>
{
  public:
   xAssemblerBasicAndTranspose(M& A, V& b, S& s) : xAssemblerBaseData<M, V, S>(&A, &b, &s) {}
   xAssemblerBasicAndTranspose(M& A, V& b) : xAssemblerBaseData<M, V, S>(&A, &b, nullptr) {}
   xAssemblerBasicAndTranspose(M& A) : xAssemblerBaseData<M, V, S>(&A, nullptr, nullptr) {}
   xAssemblerBasicAndTranspose(V& b) : xAssemblerBaseData<M, V, S>(0, &b, 0) {}
   xAssemblerBasicAndTranspose(S& s) : xAssemblerBaseData<M, V, S>(0, 0, &s) {}
   xAssemblerBasicAndTranspose() : xAssemblerBaseData<M, V, S>(0, 0, 0) {}

   /// Functions from the xAssembleConcept
   template <typename LOCALMATRIX, typename ITERDOF>
   void assemble_matrix(ITERDOF firsti, ITERDOF lasti, ITERDOF firstj, ITERDOF lastj, const LOCALMATRIX& matr)
   {
      ITERDOF firstj_bak = firstj;
      for (int i = 0; firsti != lasti; ++firsti, ++i)
      {
         firstj = firstj_bak;
         for (int j = 0; firstj != lastj; ++firstj, ++j)
         {
            S val = this->coeff * matr(i, j);
            if (!xPolicyAssemblerOnZero<typename M::matrix_assembly_on_zero>::isNull(val))
            {
               AssembleTermMatrix(*firsti, *firstj, val);
            }
         }
      }
   }

   template <typename LOCALMATRIX, typename ITERDOF>
   void assemble_matrix(ITERDOF first, ITERDOF last, const LOCALMATRIX& matr)
   {
      ITERDOF firsti = first, firstj, lastj;
      int j;
      for (int i = 0; firsti != last; ++firsti, ++i)
      {
         xPolicyAssemblerLoop<typename M::matrix_pattern>::firstForLoop(first, firsti, firstj, i, j);
         lastj = xPolicyAssemblerLoop<typename M::matrix_pattern>::lastForLoop(last, firsti);
         for (; firstj != lastj; ++firstj, ++j)
         {
            S val = this->coeff * matr(i, j);
            if (!xPolicyAssemblerOnZero<typename M::matrix_assembly_on_zero>::isNull(val))
               AssembleTermMatrix(*firsti, *firstj, val);
         }
      }
   }

   template <typename LOCALVECTOR, typename ITERDOF>
   void assemble_vector(ITERDOF first, ITERDOF last, const LOCALVECTOR& vect)
   {
      for (int i = 0; first != last; ++first, ++i)
      {
         S val = vect[i] * this->coeff;
         if (val != xtool::xDataType<S>::zero())
            xMatVecScaAssembler<M, V, S>::assemble((*first)->getState(), val, this->vec, this->sca);
      }
   }
   template <typename LOCALSCALAR>
   void assemble_scalar(const LOCALSCALAR& scal)
   {
      S val = scal * this->coeff;
      if (val != xtool::xDataType<S>::zero()) xMatVecScaAssembler<M, V, S>::assemble(val, this->sca);
   }

  private:
   void AssembleTermMatrix(const xValue<S>* dofi, const xValue<S>* dofj, S val)
   {
      if (xTraitsAssembler<typename M::matrix_pattern>::do_transpose)
      {
         xMatVecScaAssembler<M, V, S>::assemble(dofi->getState(), dofj->getState(), val, this->mat, this->vec, this->sca);
         xMatVecScaAssembler<M, V, S>::assemble(dofj->getState(), dofi->getState(), val, this->mat, this->vec, this->sca);
      }
      else
      {
         xMatVecScaAssembler<M, V, S>::assemble(dofi->getState(), dofj->getState(), val, this->mat, this->vec, this->sca);
      }
   }
   /// End  Function for the xAssembleConcept
};

/// classical lumped assembler, all entries are summed on the diagonal
template <class M = xlinalg::xCSRMatrix, class V = xlinalg::xCSRVector, class S = double>
class xAssemblerLumped : public xAssemblerBaseData<M, V, S>
{
  public:
   xAssemblerLumped(M& A) : xAssemblerBaseData<M, V, S>(&A, nullptr, nullptr) {}
   xAssemblerLumped(V& b) : xAssemblerBaseData<M, V, S>(0, &b, 0) {}
   xAssemblerLumped() : xAssemblerBaseData<M, V, S>(0, 0, 0) {}

   template <typename LOCALMATRIX, typename ITERDOF>
   void assemble_matrix(ITERDOF firsti, ITERDOF lasti, ITERDOF firstj, ITERDOF lastj, const LOCALMATRIX& matr)
   {
      ITERDOF firstj_bak = firstj;
      for (int i = 0; firsti != lasti; ++firsti, ++i)
      {
         firstj = firstj_bak;
         for (int j = 0; firstj != lastj; ++firstj, ++j)
         {
            S val = this->coeff * matr(i, j);
            if (!xPolicyAssemblerOnZero<typename M::matrix_assembly_on_zero>::isNull(val))
               xMatVecScaAssembler<M, V, S>::assemble((*firsti)->getState(), (*firstj)->getState(), val, this->mat, this->vec,
                                                      this->sca);
         }
      }
   }
   template <typename LOCALMATRIX, typename ITERDOF>
   void assemble_matrix(ITERDOF first, ITERDOF last, const LOCALMATRIX& matr)
   {
      ITERDOF firsti = first, firstj, lastj;
      int j;
      for (int i = 0; firsti != last; ++firsti, ++i)
      {
         xPolicyAssemblerLoop<typename M::matrix_pattern>::firstForLoop(first, firsti, firstj, i, j);
         lastj = xPolicyAssemblerLoop<typename M::matrix_pattern>::lastForLoop(last, firsti);
         for (; firstj != lastj; ++firstj, ++j)
         {
            S val = this->coeff * matr(i, j);
            if (!xPolicyAssemblerOnZero<typename M::matrix_assembly_on_zero>::isNull(val))
               xMatVecScaAssembler<M, V, S>::assemble((*firsti)->getState(), (*firsti)->getState(), val, this->mat, this->vec,
                                                      this->sca);
         }
      }
   }
   template <typename LOCALVECTOR, typename ITERDOF>
   void assemble_vector(ITERDOF first, ITERDOF last, const LOCALVECTOR& vect)
   {
      for (int i = 0; first != last; ++first, ++i)
      {
         S val = vect[i] * this->coeff;
         if (val != xtool::xDataType<S>::zero())
            xMatVecScaAssembler<M, V, S>::assemble((*first)->getState(), val, this->vec, this->sca);
      }
   }
   template <typename LOCALSCALAR>
   void assemble_scalar(const LOCALSCALAR& scal)
   {
      S val = scal * this->coeff;
      if (val != xtool::xDataType<S>::zero()) xMatVecScaAssembler<M, V, S>::assemble(val, this->sca);
   }
};

/// same as xAssemblerLumped, but coeffcients of lumped matrix are assembled in vector
template <class M = xlinalg::xCSRMatrix, class V = xlinalg::xCSRVector, class S = double>
class xAssemblerLumpedInVector : public xAssemblerBaseData<M, V, S>
{
  public:
   xAssemblerLumpedInVector(M& A) : xAssemblerBaseData<M, V, S>(&A, nullptr, nullptr) {}
   xAssemblerLumpedInVector(V& b) : xAssemblerBaseData<M, V, S>(0, &b, 0) {}
   xAssemblerLumpedInVector() : xAssemblerBaseData<M, V, S>(0, 0, 0) {}
   template <typename LOCALMATRIX, typename ITERDOF>
   void assemble_matrix(ITERDOF firsti, ITERDOF lasti, ITERDOF firstj, ITERDOF lastj, const LOCALMATRIX& matr)
   {
      ITERDOF firstj_bak = firstj;
      for (int i = 0; firsti != lasti; ++firsti, ++i)
      {
         firstj = firstj_bak;
         for (int j = 0; firstj != lastj; ++firstj, ++j)
         {
            S val = this->coeff * matr(i, j);
            if (!xPolicyAssemblerOnZero<typename M::matrix_assembly_on_zero>::isNull(val))
               xMatVecScaAssembler<M, V, S>::assemble((*firsti)->getState(), (*firstj)->getState(), val, this->mat, this->vec,
                                                      this->sca);
         }
      }
   }
   template <typename LOCALMATRIX, typename ITERDOF>
   void assemble_matrix(ITERDOF first, ITERDOF last, const LOCALMATRIX& matr)
   {
      ITERDOF firsti = first, firstj, lastj;
      int j;
      for (int i = 0; firsti != last; ++firsti, ++i)
      {
         xPolicyAssemblerLoop<typename M::matrix_pattern>::firstForLoop(first, firsti, firstj, i, j);
         lastj = xPolicyAssemblerLoop<typename M::matrix_pattern>::lastForLoop(last, firsti);
         for (; firstj != lastj; ++firstj, ++j)
         {
            S val = this->coeff * matr(i, j);
            if (!xPolicyAssemblerOnZero<typename M::matrix_assembly_on_zero>::isNull(val))
               xMatVecScaAssembler<M, V, S>::assemble((*firsti)->getState(), val, this->vec, this->sca);
         }
      }
   }
   template <typename LOCALVECTOR, typename ITERDOF>
   void assemble_vector(ITERDOF first, ITERDOF last, const LOCALVECTOR& vect)
   {
      std::cout << "xAssemblerLumpedInVector is not supposed to be used to assemble vectors !!! " << std::endl;
      throw -1;
   }
   template <typename LOCALSCALAR>
   void assemble_scalar(const LOCALSCALAR& scal)
   {
      S val = scal * this->coeff;
      if (val != xtool::xDataType<S>::zero()) xMatVecScaAssembler<M, V, S>::assemble(val, this->sca);
   }
};

/// special lumped assembler, it assembles the same thing at each node.
/// It does the same job as  xAssemblerLumped except when the element
/// has a void part
template <class M = xlinalg::xCSRMatrix, class V = xlinalg::xCSRVector, class S = double>
class xAssemblerLumpedEqu : public xAssemblerLumped<M, V, S>
{
  public:
   xAssemblerLumpedEqu(M& A) : xAssemblerLumped<M, V, S>(A) {}
   xAssemblerLumpedEqu(V& b) : xAssemblerLumped<M, V, S>(b) {}
   xAssemblerLumpedEqu() : xAssemblerLumped<M, V, S>() {}

   template <typename LOCALMATRIX, typename ITERDOF>
   void assemble_matrix(ITERDOF firsti, ITERDOF lasti, ITERDOF firstj, ITERDOF lastj, const LOCALMATRIX& matr)
   {
      ITERDOF firstj_bak = firstj;
      for (int i = 0; firsti != lasti; ++firsti, ++i)
      {
         firstj = firstj_bak;
         for (int j = 0; firstj != lastj; ++firstj, ++j)
         {
            S val = this->coeff * matr(i, j);
            if (!xPolicyAssemblerOnZero<typename M::matrix_assembly_on_zero>::isNull(val))
               xMatVecScaAssembler<M, V, S>::assemble((*firsti)->getState(), (*firstj)->getState(), val, this->mat, this->vec,
                                                      this->sca);
         }
      }
   }

   template <typename LOCALMATRIX, typename ITERDOF>
   void assemble_matrix(ITERDOF first, ITERDOF last, const LOCALMATRIX& matr)
   {
      S sum = xtool::xDataType<S>::zero();
      int i;
      ITERDOF firsti = first, firstj;
      for (i = 0; firsti != last; ++firsti, ++i)
      {
         firstj = first;
         for (int j = 0; firstj != last; ++firstj, ++j)
         {
            sum += this->coeff * matr(i, j);
         }
      }
      sum /= (S)i;
      if (!xPolicyAssemblerOnZero<typename M::matrix_assembly_on_zero>::isNull(sum)) return;
      for (firsti = first; firsti != last; ++firsti)
         xMatVecScaAssembler<M, V, S>::assemble((*firsti)->getState(), (*firsti)->getState(), sum, this->mat, this->vec,
                                                this->sca);
   }

   template <typename LOCALVECTOR, typename ITERDOF>
   void assemble_vector(ITERDOF first, ITERDOF last, const LOCALVECTOR& vect)
   {
      S sum = xtool::xDataType<S>::zero();
      int i;
      ITERDOF firsti = first, firstj;
      for (i = 0; firsti != last; ++firsti, ++i)
      {
         sum += this->coeff * vect[i];
      }
      sum /= (S)i;
      for (firsti = first; firsti != last; ++firsti)
         if (sum != xtool::xDataType<S>::zero())
            xMatVecScaAssembler<M, V, S>::assemble((*firsti)->getState(), sum, this->vec, this->sca);
      return;
   }

   template <typename LOCALSCALAR>
   void assemble_scalar(const LOCALSCALAR& scal)
   {
      S val = scal * this->coeff;
      if (val != xtool::xDataType<S>::zero()) xMatVecScaAssembler<M, V, S>::assemble(val, this->sca);
   }
};

/// special lumped assembler, for block matrices
/// takes each block of the assembled matrix, compute the average and assemble it on the corresponding diagonal coefficients
/// not supposed to be used to assemble vectors
template <class M = xlinalg::xCSRMatrix, class V = xlinalg::xCSRVector, class S = double>
class xAssemblerLumpedBlockEqu : public xAssemblerLumped<M, V, S>
{
  public:
   xAssemblerLumpedBlockEqu(M& A) : xAssemblerLumped<M, V, S>(A) {}
   xAssemblerLumpedBlockEqu(V& b) : xAssemblerLumped<M, V, S>(b) {}
   xAssemblerLumpedBlockEqu() : xAssemblerLumped<M, V, S>() {}

   template <typename LOCALMATRIX, typename ITERDOF>
   void assemble_matrix(ITERDOF firsti, ITERDOF lasti, ITERDOF firstj, ITERDOF lastj, const LOCALMATRIX& matr)
   {
      ITERDOF firstj_bak = firstj;
      for (int i = 0; firsti != lasti; ++firsti, ++i)
      {
         firstj = firstj_bak;
         for (int j = 0; firstj != lastj; ++firstj, ++j)
         {
            S val = this->coeff * matr(i, j);
            if (!xPolicyAssemblerOnZero<typename M::matrix_assembly_on_zero>::isNull(val))
               xMatVecScaAssembler<M, V, S>::assemble((*firsti)->getState(), (*firstj)->getState(), val, this->mat, this->vec,
                                                      this->sca);
         }
      }
   }

   template <typename LOCALMATRIX, typename ITERDOF>
   void assemble_matrix(ITERDOF first, ITERDOF last, const LOCALMATRIX& matr)
   {
      if (matr.getNbRow() != matr.getNbCol())
      {
         std::cout << "Error : xAssemblerLumpedBlockEqu works only for square matrices ! " << std::endl;
         throw -1;
      }
      std::vector<S> sum(matr.getNbRow());

      ITERDOF firsti = first, firstj;
      for (int i = 0; firsti != last; ++firsti, ++i)
      {
         S val = xtool::xDataType<S>::zero();
         std::vector<int> indexes;

         firstj = first;
         for (int j = 0; firstj != last; ++firstj, ++j)
         {
            val += this->coeff * matr(i, j);
            indexes.push_back(j);
         }
         val /= indexes.size();
         for (int k = 0; k < indexes.size(); k++)
         {
            sum[indexes[k]] += val;
         }
      }
      int l = 0;
      for (firsti = first; firsti != last; ++firsti, ++l)
         xMatVecScaAssembler<M, V, S>::assemble((*firsti)->getState(), (*firsti)->getState(), sum[l], this->mat, this->vec,
                                                this->sca);
   }

   template <typename LOCALVECTOR, typename ITERDOF>
   void assemble_vector(ITERDOF first, ITERDOF last, const LOCALVECTOR& vect)
   {
      std::cout << "xAssemblerLumpedBlockEqu is not supposed to be used to assemble vectors !!! " << std::endl;
      throw -1;
   }

   template <typename LOCALSCALAR>
   void assemble_scalar(const LOCALSCALAR& scal)
   {
      S val = scal * this->coeff;
      if (val != xtool::xDataType<S>::zero()) xMatVecScaAssembler<M, V, S>::assemble(val, this->sca);
   }
};
/// does the same thing as xAssemblerLumpedBlockEqu, but the coefficients of the lumped matrix are put in a
/// vector
template <class M = xlinalg::xCSRMatrix, class V = xlinalg::xCSRVector, class S = double>
class xAssemblerLumpedBlockEquInVector : public xAssemblerLumped<M, V, S>
{
  public:
   xAssemblerLumpedBlockEquInVector(M& A) : xAssemblerLumped<M, V, S>(A) {}
   xAssemblerLumpedBlockEquInVector(V& b) : xAssemblerLumped<M, V, S>(b) {}
   xAssemblerLumpedBlockEquInVector() : xAssemblerLumped<M, V, S>() {}

   template <typename LOCALMATRIX, typename ITERDOF>
   void assemble_matrix(ITERDOF firsti, ITERDOF lasti, ITERDOF firstj, ITERDOF lastj, const LOCALMATRIX& matr)
   {
      ITERDOF firstj_bak = firstj;
      for (int i = 0; firsti != lasti; ++firsti, ++i)
      {
         firstj = firstj_bak;
         for (int j = 0; firstj != lastj; ++firstj, ++j)
         {
            S val = this->coeff * matr(i, j);
            if (!xPolicyAssemblerOnZero<typename M::matrix_assembly_on_zero>::isNull(val))
               xMatVecScaAssembler<M, V, S>::assemble((*firsti)->getState(), (*firstj)->getState(), val, this->mat, this->vec,
                                                      this->sca);
         }
      }
   }

   template <typename LOCALMATRIX, typename ITERDOF>
   void assemble_matrix(ITERDOF first, ITERDOF last, const LOCALMATRIX& matr)
   {
      if (matr.getNbRow() != matr.getNbCol())
      {
         std::cout << "Error : xAssemblerLumpedBlockEquInVector works only for square matrices ! " << std::endl;
         throw -1;
      }
      std::vector<S> sum(matr.getNbRow());

      ITERDOF firsti = first, firstj;
      for (int i = 0; firsti != last; ++firsti, ++i)
      {
         S val = xtool::xDataType<S>::zero();
         std::vector<int> indexes;

         firstj = first;
         for (int j = 0; firstj != last; ++firstj, ++j)
         {
            val += this->coeff * matr(i, j);
            indexes.push_back(j);
         }
         val /= indexes.size();
         for (int k = 0; k < indexes.size(); k++)
         {
            sum[indexes[k]] += val;
         }
      }
      int l = 0;
      for (firsti = first; firsti != last; ++firsti, ++l)
         xMatVecScaAssembler<M, V, S>::assemble((*firsti)->getState(), sum[l], this->vec, this->sca);
   }

   template <typename LOCALVECTOR, typename ITERDOF>
   void assemble_vector(ITERDOF first, ITERDOF last, const LOCALVECTOR& vect)
   {
      std::cout << "xAssemblerLumpedBlockEquInVector is not supposed to be used to assemble vectors !!! " << std::endl;
      throw -1;
   }

   template <typename LOCALSCALAR>
   void assemble_scalar(const LOCALSCALAR& scal)
   {
      S val = scal * this->coeff;
      if (val != xtool::xDataType<S>::zero()) xMatVecScaAssembler<M, V, S>::assemble(val, this->sca);
   }
};

}  // namespace xfem

#endif
