/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/
#ifdef _XDISTVECTOR_H
#ifndef _XDISTVECTOR_IMP_H
#define _XDISTVECTOR_IMP_H

#include "xBlasDef.h"
#include "xDataType.h"
#include "xMPIDataType.h"

namespace xlinalg
{
template <typename VT>
inline xDistVector<VT> &xDistVector<VT>::operator=(const xDistVector<VT> &rhs)
{
   if (this != &rhs)
   {
      if (&dist_index != &rhs.dist_index)
      {
         // assignement of vector with two different indexing is not proposed for now
         throw -1;
      }
      assert(status & INSERTOFF);
      // do copy to maintaint dat status => exchanger remain correct
      std::copy(rhs.dat.begin(), rhs.dat.end(), dat.begin());
      status = rhs.status;
      // switch to INSERTOFF to respect initial state. Done only if rhs status is !INSERTOFF
      if (!(status & INSERTOFF)) switchInsertModeOff();
   }
   return (*this);
}

template <typename VT>
inline void xDistVector<VT>::AddVal(const int &i, const VT &val)
{
   assert(!(status & INSERTOFF));
   // i : FORTRAN NUMBERING
   dat[dist_index.getPackedIndex(i)] += val;
}

template <typename VT>
inline VT &xDistVector<VT>::getVal(xDistIndex::idx_t i)
{
   assert(!(status & INSERTOFF));
   return dat.at(dist_index.getPackedIndex(i));
}

template <typename VT>
inline const VT &xDistVector<VT>::getVal(xDistIndex::idx_t i) const
{
   return dat.at(dist_index.getPackedIndex(i));
}

template <typename VT>
inline VT &xDistVector<VT>::operator[](xDistIndex::idx_t i)
{
   assert(!(status & INSERTOFF));
   return dat[i];
}

template <typename VT>
inline const VT &xDistVector<VT>::operator[](xDistIndex::idx_t i) const
{
   return dat[i];
}

template <typename VT>
inline VT &xDistVector<VT>::at(xDistIndex::idx_t i)
{
   assert(!(status & INSERTOFF));
   return dat.at(i);
}

template <typename VT>
inline const VT &xDistVector<VT>::at(xDistIndex::idx_t i) const
{
   return dat.at(i);
}

template <typename VT>
inline const VT *xDistVector<VT>::data() const
{
   return dat.data();
}

template <typename VT>
inline VT *xDistVector<VT>::data()
{
   assert(!(status & INSERTOFF));
   return dat.data();
}

template <typename VT>
inline typename xDistVector<VT>::iterator xDistVector<VT>::begin()
{
   assert(!(status & INSERTOFF));
   return dat.begin();
}

template <typename VT>
inline typename xDistVector<VT>::iterator xDistVector<VT>::end()
{
   assert(!(status & INSERTOFF));
   return dat.end();
}

template <typename VT>
inline typename xDistVector<VT>::const_iterator xDistVector<VT>::begin() const
{
   return dat.begin();
}

template <typename VT>
inline typename xDistVector<VT>::const_iterator xDistVector<VT>::end() const
{
   return dat.end();
}

template <typename VT>
inline void xDistVector<VT>::Zerodata()
{
   std::fill(dat.begin(), dat.end(), xtool::xDataType<VT>::zero());
   MPI_Barrier(dist_index.getComm());
}

template <typename VT>
inline void xDistVector<VT>::Printdata(std::ostream &out) const
{
   out << "###########################\n";
   out << "Vector of dimension " << dat.size() << " is:\n";
   std::copy(dat.begin(), dat.end(), std::ostream_iterator<VT>(out, " "));
   out << "\n";
   out << "###########################\n";
}

template <typename VT>
inline xDistVector<VT> &xDistVector<VT>::operator+=(const xDistVector<VT> &rhs)
{
   axpy(xtool::xDataType<VT>::one(), rhs);
   return *this;
}

template <typename VT>
inline xDistVector<VT> xDistVector<VT>::operator+(const xDistVector<VT> &rhs) const
{
   if (&dist_index != &rhs.dist_index)
   {
      // addition of vector with two different indexing is not proposed for now
      throw -1;
   }
   xDistVector<VT> Y(*this);
   Y.axpy(xtool::xDataType<VT>::one(), rhs);
   return Y;
}

template <typename VT>
inline xDistVector<VT> &xDistVector<VT>::operator-=(const xDistVector<VT> &rhs)
{
   axpy(-xtool::xDataType<VT>::one(), rhs);
   return *this;
}

template <typename VT>
inline xDistVector<VT> xDistVector<VT>::operator-(const xDistVector<VT> &rhs) const
{
   if (&dist_index != &rhs.dist_index)
   {
      // addition of vector with two different indexing is not proposed for now
      throw -1;
   }
   xDistVector<VT> Y(*this);
   Y.axpy(-xtool::xDataType<VT>::one(), rhs);
   return Y;
}

template <typename VT>
inline xDistVector<VT> xDistVector<VT>::operator*(const VT &alpha) const
{
   xDistVector<VT> res(*this);
   res.scal(alpha);
   return res;
}

template <typename VT>
inline VT xDistVector<VT>::operator*(const xDistVector<VT> &rhs) const
{
   return dot(rhs);
}

template <typename VT>
inline double xDistVector<VT>::nrm2() const
{
   VT res = dot(*this);
   return sqrt(res);
}

template <>
inline double xDistVector<std::complex<double>>::nrm2() const
{
   std::complex<double> res = dot(*this);
   assert(fabs(res.imag()) < std::numeric_limits<double>::epsilon());
   return sqrt(res.real());
}

//===== xKeyManagerDistVector =====================================================================================
template <typename VT>
xDistVector<VT>::xKeyManagerDistVector::xKeyManagerDistVector(const xDistIndex &dist_index_, xDistVector<VT>::Vector &data_)
    : dist_index(dist_index_), dat(data_)
{
}

template <typename VT>
auto xDistVector<VT>::xKeyManagerDistVector::localObjectKey(const xDistIndex::idx_t &lo) -> information_key_t
{
   return &dat[dist_index.getPackedIndex(lo)];
}

template <typename VT>
xtool::xConstPartitionObject<xDistIndex::idx_t> xDistVector<VT>::xKeyManagerDistVector::getConstPartitionObject(
    const xDistIndex::idx_t &lo)
{
   return dist_index.getConstPartitionObject(lo);
}

template <typename VT>
auto xDistVector<VT>::xKeyManagerDistVector::remoteObjectKey(const xtool::xRemoteObject<xDistIndex::idx_t> &ro,
                                                             const xDistIndex::idx_t &lo) -> information_key_t
{
   // Here we cheat : normally we should have given remote address of the vector cell corresponding to lo like in localObjectKey
   // we give a proper pointer address. But this would have required to create a partition manager just for this instance vector
   // and exchange address. Not pleasant at all.
   // In fact all what we need here is something giving the right order in remote proc of impacted cell in dat .....
   // For that we use address of remote index of dist_index instance corresponding to lo given by ro. Has those address represent
   // pack_id address these are well ordered. Remote will use for keys it's local vector address constructed by use of
   // getPackedIndex, says the pack_id address ....
   //
   // To make this remote address a key we need to cast it.
   return reinterpret_cast<information_key_t>(const_cast<xDistIndex::idx_t *>(ro.getObjectAddress()));
}
//===== xInfoManagerDistVectorReduce  =====================================================================================
template <typename VT>
auto xDistVector<VT>::xInfoManagerDistVectorReduce::getInfo(information_key_t key, int sendto) -> information_t
{
   VT v = *key;
   *key = xtool::xDataType<VT>::zero();
   return v;
}

template <typename VT>
void xDistVector<VT>::xInfoManagerDistVectorReduce::setInfo(information_key_t key, const information_t &info, int receivedfrom)
{
   *key += info;
}
//===== xInfoManagerDistVectorSet  ========================================================================================
template <typename VT>
auto xDistVector<VT>::xInfoManagerDistVectorSet::getInfo(information_key_t key, int sendto) -> information_t
{
   return *key;
}

template <typename VT>
void xDistVector<VT>::xInfoManagerDistVectorSet::setInfo(information_key_t key, const information_t &info, int receivedfrom)
{
   *key = info;
}

//===== xInfoManagerDistVectorDoubleReduce =====================================================================================
template <typename VT>
xDistVector<VT>::xInfoManagerDistVectorDoubleReduce::xInfoManagerDistVectorDoubleReduce(xDistVector<VT>::Vector &data_,
                                                                                        xDistVector<VT>::Vector &other_)
    : dat(data_), other(other_)
{
}

template <typename VT>
auto xDistVector<VT>::xInfoManagerDistVectorDoubleReduce::getInfo(information_key_t key, int sendto) -> information_t
{
   const size_t i = key - &dat[0];
   information_t v = {{*key, other[i]}};
   *key = xtool::xDataType<VT>::zero();
   other[i] = xtool::xDataType<VT>::zero();
   return v;
}

template <typename VT>
void xDistVector<VT>::xInfoManagerDistVectorDoubleReduce::setInfo(information_key_t key, const information_t &info,
                                                                  int receivedfrom)
{
   *key += info[0];
   other[key - &dat[0]] += info[1];
}

//===== xDistVector<VT> ===========================================================================================
/*
 * Coherence table
 *   x: set                         ISon : insert mode "on"
 *   I: impossible                  ISoff: insert mode "off"
 *   P: possible                    G    : global state             |xxxxxxxxxx|xxxxxx|
 *                                  L    : local state              |x0tx0xxywv|u0x0g0|
 *                                  R    : reduced state            |xxxxxxxxxx|000000|
 *                                                                   owned       remote
 *
 *         | ISon | ISoff |  G  |  L  |  R  |
 *   ISon  |   x  |   I   |  I  |  P  | I/P*|   P* : possible due to scatter
 *   ISoff |      |   x   |  P  |  P  |  P  |
 *    G    |      |       |  x  |  I  |  I  |
 *    L    |      |       |  I  |  x  |  P  |
 *    R    |      |       |  I  |  x  |  x  |
 *
 *
 */
template <typename VT>
xDistVector<VT>::xDistVector(const xDistIndex &dist_index_)
    : dat(dist_index_.getPackedIndexSize(), xtool::xDataType<VT>::zero()),
      dist_index(dist_index_),
      reduce_keys(dist_index.getComm()),
      set_keys(dist_index.getComm()),
      status(RESET),
      proc_id(-1)
{
   // create key manager
   xKeyManagerDistVector key_manager(dist_index, dat);
   // create keys
   reduce_keys.accumulateKeysOwnerGather(dist_index.begin(), dist_index.end(), key_manager);
   set_keys.accumulateKeysOwnerScatter(dist_index.begin(), dist_index.end(), key_manager);
   // getting rank in world
   MPI_Comm_rank(dist_index.getComm(), &proc_id);
}

template <typename VT>
xDistVector<VT>::xDistVector(const xDistVector<VT> &other)
    : dat(other.dat),
      dist_index(other.dist_index),
      reduce_keys(dist_index.getComm()),
      set_keys(dist_index.getComm()),
      status(other.status),
      proc_id(other.proc_id)
{
   // create key manager
   xKeyManagerDistVector key_manager(dist_index, dat);
   // create keys : here keys have changed from other has dat is not the same
   reduce_keys.accumulateKeysOwnerGather(dist_index.begin(), dist_index.end(), key_manager);
   set_keys.accumulateKeysOwnerScatter(dist_index.begin(), dist_index.end(), key_manager);
}

template <typename VT>
xDistVector<VT> *xDistVector<VT>::clone() const
{
   xDistVector<VT> *p = new xDistVector<VT>(dist_index);
   p->switchInsertModeOff();
   return p;
}

template <typename VT>
VT xDistVector<VT>::dot(const xDistVector<VT> &y) const
{
   xDistVector<VT> &x = const_cast<xDistVector<VT> &>(*this);
   xDistVector<VT> &y_ = const_cast<xDistVector<VT> &>(y);
   if (&x.dist_index != &y_.dist_index)
   {
      // dot product of vector with two different indexing is not proposed for now
      throw -1;
   }
   assert(status & INSERTOFF);
   assert(y_.status & INSERTOFF);
   int N = dist_index.getPackedLocalIndexSize();

   // if none are Reduced/Global optimize transfer if not the same vector
   if (!(status & FLAGRG) && !(y_.status & FLAGRG))
   {
      // if same vector one transfer is enough
      if (this == &y)
      {
         y_.reduceOnOwner();
      }
      // not identical
      else
      {
         x.doubleReduceOnOwner(y_);
         /*old
          *
            // should used new exchanger which exchange both vector value in on shot
            // for now dirty double exchange
            x.reduceOnOwner();
            y_.reduceOnOwner();
          */
      }
   }
   // only y reduced/global, change x
   else if (!(status & FLAGRG))
      x.reduceOnOwner();
   // only x reduced/global, change y
   else if (!(y_.status & FLAGRG))
      y_.reduceOnOwner();

   assert(status & FLAGRG);
   assert(y_.status & FLAGRG);

   // Now we can do local product
   // ----- WARNING ! -----
   // For complex dot is the hermitian scalar product
   // v = x^H . y (i.e., dot is in fact zdotc)
   VT v = xCPPBlasDef<VT>::dot(N, x.dat.data(), y_.dat.data());

   // sum local product to get result on all proc
   MPI_Allreduce(MPI_IN_PLACE, &v, 1, xtool::xMPIDataType<VT>(), MPI_SUM, dist_index.getComm());

   return v;
}

template <typename VT>
void xDistVector<VT>::axpy(const VT &a, const xDistVector<VT> &x)
{
   xDistVector<VT> &x_ = const_cast<xDistVector<VT> &>(x);
   if (&dist_index != &x_.dist_index)
   {
      // addition of vector with two different indexing is not proposed for now
      throw -1;
   }
   assert(status & INSERTOFF);
   assert(x_.status & INSERTOFF);
   // by default sum on all values
   int N = dat.size();

   bool propagate = false;

   // if not identical
   if (this != &x)
   {
      // x_ is in global status
      if (x_.status & GLOBAL)
      {
         // this is in local status, reduced or not
         // Only sum on local/owned term is needed.
         if (!(status & GLOBAL)) N = dist_index.getPackedLocalIndexSize();
      }
      // x_ is in reduced status
      else if (x_.status & REDUCED)
      {
         // this is in global status
         // Only sum on local/owned term is needed but Y have to be reset in global mode to propagate result
         // this must be temporally flagged as reduced
         if (status & GLOBAL)
         {
            propagate = true;
            N = dist_index.getPackedLocalIndexSize();
            status ^= GLOBAL;
            status |= REDUCED;
         }
      }
      // x_ is not in global nor reduced status
      else
      {
         // this is in reduce status
         // Only this state have to be changed
         if (status & REDUCED)
         {
            status ^= REDUCED;
         }
         // this is in GLOBAL status
         // Here we have not much choice, two communication is involved :
         // A) Changing x in global (2 comm) compute and back in reduce
         // or
         // B) Changing this in local  compute and back in global (2 comm)
         // B solution preserve x => choose it
         // this must be temporally flagged as not reduced as computation may break reduction.
         // Then when passing this to global state a reduce will be done
         else if (status & GLOBAL)
         {
            switchToLocalValue();
            status ^= REDUCED;
            propagate = true;
         }
      }
   }

   xCPPBlasDef<VT>::axpy(N, a, &x_.dat[0], &(dat[0]));

   if (propagate) switchToGlobalValue();
}

template <typename VT>
void xDistVector<VT>::componentProduct(const xDistVector<VT> &x, const xDistVector<VT> &y)
{
   xDistVector<VT> &x_ = const_cast<xDistVector<VT> &>(x);
   xDistVector<VT> &y_ = const_cast<xDistVector<VT> &>(y);
   if ((&dist_index != &x_.dist_index) || (&dist_index != &y_.dist_index))
   {
      // vectors with different indexing is not proposed for now
      throw -1;
   }
   assert(status & INSERTOFF);
   assert(x_.status & INSERTOFF);
   assert(y_.status & INSERTOFF);

   // x_ and y_ have to be in compatible state so that sum of product contribution correspond to correct product

   // reseting to insert mode off "this" => local status
   // No cleaning as it is scratched by operation below
   status = INSERTOFF;

   char oper = 0;

   // Remove impossible combination -----------------------
   // if x not in global nor reduced mode
   if (!(x_.status & FLAGRG))
   {
      // if y is also not in global nor reduced mode
      // => x is turned in global state (avoid reducing x and y, only reduce and scatter x)
      if (!(y_.status & FLAGRG)) x_.switchToGlobalValue();
      // if y is in reduced mode
      // => y is turned in global state
      else if (y_.status & REDUCED)
         y_.switchToGlobalValue();
   }
   // if y is not in global nor reduced mode but x is in reduced or global state
   // => x is turned in global state what ever it was
   else if (!(y_.status & FLAGRG))
      x_.switchToGlobalValue();

   // analyze combination ---------------------------------
   if ((x_.status & REDUCED) || (y_.status & REDUCED))
      oper = 1;
   else if ((x_.status & GLOBAL) && (y_.status & GLOBAL))
      oper = 2;
   else
      oper = 3;

   // do operation
   switch (oper)
   {
      case 1:
      {
         status |= REDUCED;
         xDistIndex::idx_t n = dist_index.getPackedLocalIndexSize();
         std::transform(x_.dat.begin(), x_.dat.begin() + n, y_.dat.begin(), dat.begin(), std::multiplies<VT>());
         std::fill(dat.begin() + n, dat.end(), xtool::xDataType<VT>::zero());
         break;
      }
      case 2:
      {
         status |= GLOBAL;
      }
      case 3:
      {
         std::transform(x_.dat.begin(), x_.dat.end(), y_.dat.begin(), dat.begin(), std::multiplies<VT>());
         break;
      }
   }
}

template <typename VT>
void xDistVector<VT>::scal(const VT &alpha)
{
   assert(status & INSERTOFF);
   xCPPBlasDef<VT>::scal(dat.size(), alpha, &dat[0]);
}

template <typename VT>
xDistVector<VT> &xDistVector<VT>::operator=(const VT &alpha)
{
   assert(status & INSERTOFF);
   if (status & GLOBAL)
   {
      for_each(dat.begin(), dat.end(), [&alpha](VT &v) { v = alpha; });
   }
   else
   {
      xDistIndex::idx_t n = dist_index.getPackedLocalIndexSize();
      for_each(dat.begin(), dat.begin() + n, [&alpha](VT &v) { v = alpha; });
      if (!(status & REDUCED)) std::fill(dat.begin() + n, dat.end(), xtool::xDataType<VT>::zero());
   }
   return *this;
}

template <typename VT>
void xDistVector<VT>::doubleReduceOnOwner(xDistVector<VT> &other)
{
   assert(status & INSERTOFF);
   assert(other.status & INSERTOFF);
   assert(!(status & FLAGRG) && !(other.status & FLAGRG));
   xInfoManagerDistVectorDoubleReduce double_reduce_info(dat, other.dat);
   xtool::exchangeInformation(reduce_keys, double_reduce_info);
   status |= REDUCED;
   other.status |= REDUCED;
}

template <typename VT>
void xDistVector<VT>::reduceOnOwner()
{
   assert(status & INSERTOFF);
   switchToLocalValue();
   if (!(status & REDUCED))
   {
      xtool::exchangeInformation(reduce_keys, reduce_info);
   }
   status |= REDUCED;
}

template <typename VT>
void xDistVector<VT>::switchToGlobalValue()
{
   assert(status & INSERTOFF);
   if (!(status & GLOBAL))
   {
      if (!(status & REDUCED))
      {
         xtool::exchangeInformation(reduce_keys, reduce_info);
      }
      else
         status ^= REDUCED;
      xtool::exchangeInformation(set_keys, set_info);
      status |= GLOBAL;
   }
}

template <typename VT>
void xDistVector<VT>::switchToLocalValue()
{
   assert(status & INSERTOFF);
   if (status & GLOBAL)
   {
      status ^= GLOBAL;
      int n = dist_index.getPackedLocalIndexSize();
      for (auto it = dat.begin() + n, ite = dat.end(); ite != it; ++it) *it = xtool::xDataType<VT>::zero();
      status |= REDUCED;
   }
}

template <typename VT>
void xDistVector<VT>::gather(std::vector<VT> &V, int target_proc) const
{
   // to by pass const
   xDistVector<VT> &x = const_cast<xDistVector<VT> &>(*this);

   // auto commute insert mode status
   char old_status = status;
   if (!(status & INSERTOFF)) x.switchInsertModeOff();

   // if not reduced and not global reduce on owner so that all owned value be correct and the next gather
   // grab correct values
   if (!(status & FLAGRG)) x.reduceOnOwner();  // if not GLOBAL and not REDUCED turn it to REDUCED

   if (dist_index.isInIncreasingIndexing())
   {
      if (proc_id == target_proc && static_cast<xDistIndex::idx_t>(V.size()) != dist_index.getGlobalIndexSize())
         V.resize(dist_index.getGlobalIndexSize());
      MPI_Gatherv(reinterpret_cast<void *>(x.dat.data()), dist_index.getPackedLocalIndexSize(), xtool::xMPIDataType<VT>(),
                  V.data(), reinterpret_cast<int *>(const_cast<xDistIndex::idx_t *>(dist_index.getAllPackedLocalIndexSize())),
                  reinterpret_cast<int *>(const_cast<xDistIndex::idx_t *>(dist_index.getAllPackedLocalIndexDisp())),
                  xtool::xMPIDataType<VT>(), target_proc, dist_index.getComm());
   }
   else
   {
      throw -1;
   }
   // reset vector to its initial state if it was insert mode "On" (loose REDUCED state but this mode impose it)
   if (!(old_status & INSERTOFF)) x.switchInsertModeOn();
}

template <typename VT>
void xDistVector<VT>::scatter(std::vector<VT> &V, int source_proc)
{
   assert((proc_id == source_proc) ? V.size() == static_cast<size_t>(dist_index.getGlobalIndexSize()) : true);

   if (dist_index.isInIncreasingIndexing())
   {
      MPI_Scatterv(reinterpret_cast<void *>(V.data()),
                   reinterpret_cast<int *>(const_cast<xDistIndex::idx_t *>(dist_index.getAllPackedLocalIndexSize())),
                   reinterpret_cast<int *>(const_cast<xDistIndex::idx_t *>(dist_index.getAllPackedLocalIndexDisp())),
                   xtool::xMPIDataType<VT>(), reinterpret_cast<void *>(dat.data()), dist_index.getPackedLocalIndexSize(),
                   xtool::xMPIDataType<VT>(), source_proc, dist_index.getComm());
   }
   else
   {
      throw -1;
   }

   // Adapt to status
   //
   // if was in insert mode "On"
   if (!(status & INSERTOFF))
   {
      // normally GLOBALE mode impossible
      assert(!(status & GLOBAL));

      // It may be in REDUCE mode (previous scatter)
      // in this case nothing to do. Scatterv have
      // smashed reduced value and not owned stay null.
      // Otherwise we must clean not owned and pass in REDUCED
      if (!(status & REDUCED))
      {
         // fake to be in global state and clean not owned information => pass in reduced state
         status |= GLOBAL;
         status |= INSERTOFF;  // to bypass assert in switchToLocalValue
         switchToLocalValue();
         status ^= INSERTOFF;
      }
   }
   // if was in insert mode "Off"
   else
   {
      // if was in global state
      if (status & GLOBAL)
      {
         // continue scattering information everywhere => conform to global state
         xtool::exchangeInformation(set_keys, set_info);
      }
      // if was not in reduce state
      else if (!(status & REDUCED))
      {
         // fake to be in global state an clean not owned information => pass in reduced state
         status |= GLOBAL;
         switchToLocalValue();
      }
   }
}

template <typename VT>
void xDistVector<VT>::switchInsertModeOff()
{
   assert(!(status & INSERTOFF));
   assert(!(status & GLOBAL));
   // assert(!( status&REDUCED )); in fact a scatter leave things in REDUCED
   status |= INSERTOFF;
}
template <typename VT>
void xDistVector<VT>::switchInsertModeOffFromUserGlobalSetting()
{
   assert(!(status & INSERTOFF));
   assert(!(status & GLOBAL));
   status |= INSERTOFF;
   status |= GLOBAL;
}

template <typename VT>
void xDistVector<VT>::switchInsertModeOn()
{
   switchToLocalValue();
   status = RESET;
}

template <typename VT>
bool xDistVector<VT>::isInInsertModeOn() const
{
   return (!(status & INSERTOFF));
}

template <typename VT>
const xDistIndex &xDistVector<VT>::getDistIndex() const
{
   return dist_index;
}

// =================== friend function ===================================
// result = v1 + c2*v2
template <typename VT>
void add(const xDistVector<VT> &v1, const VT &c2, const xDistVector<VT> &v2, xDistVector<VT> &result)
{
   // TODO: to be optimized : remove intermediate step
   // be carefull:
   //  result=v1;
   //  result.axpy(c2,v2);
   //  may not work if result is in fact a reference to v2 !!!
   result = v1 + v2 * c2;
   return;
}

// v1 = v1 + c2*v2
template <typename VT>
void add(xDistVector<VT> &v1, VT c2, const xDistVector<VT> &v2)
{
   v1.axpy(c2, v2);
   return;
}

//! result = c1*v1 + c2*v2
template <typename VT>
void add(const VT &c1, const xDistVector<VT> &v1, const VT &c2, const xDistVector<VT> &v2, xDistVector<VT> &result)
{
   // TODO: to be optimized : remove intermediate step
   result = v1 * c1 + v2 * c2;
   return;
}

//! result = alpha(v1 + v2)
template <typename VT>
void add(const VT &alpha, const xDistVector<VT> &v1, const xDistVector<VT> &v2, xDistVector<VT> &result)
{
   // TODO: to be optimized : remove intermediate step
   result = v1 + v2;
   result.scal(alpha);
   return;
}

//! result = v1 + v2 + v3
template <typename VT>
void add(const xDistVector<VT> &v1, const xDistVector<VT> &v2, const xDistVector<VT> &v3, xDistVector<VT> &result)
{
   // TODO: to be optimized : remove intermediate step
   result = v1 + v2 + v3;
   return;
}

//! result = v1 - v2
template <typename VT>
void subtract(const xDistVector<VT> &v1, const xDistVector<VT> &v2, xDistVector<VT> &result)
{
   result = v1 - v2;
   return;
}

//! return the inner product of v1 and v2
template <typename VT>
VT InnerProduct(const xDistVector<VT> &v1, const xDistVector<VT> &v2)
{
   return v1.dot(v2);
}

}  // namespace xlinalg

#endif
#endif
