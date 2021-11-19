/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.
*/
#ifndef XSPARSEVECTOR_H
#define XSPARSEVECTOR_H

//std
#include <vector>
#include <cassert>

namespace xlinalg
{
/// This class is a rather crud sparse vector interface. It fulfill FastMarching transport concept requirements.
/*! It have been created from TLS modesValuesVector and has a replacement, it must continue to respect initial
 * interface. Fill free to add methods but remember that internally index must be in ascending order. This
 * is why add method have to be kept in private part to avoid user mistakes for now.
 * From original implementation this class is now template on the type of data treated.
 * 
 */
template <typename T>
class xSparseVector
{
    public:
        typedef typename std::vector < int >::const_iterator iter_idx_t;
        typedef typename std::vector < T >::const_iterator iter_val_t;

        // constructor
        xSparseVector();
        xSparseVector(int index, T val);
        xSparseVector(const xSparseVector &in);

        // arithmetic
        xSparseVector & operator*=(const T &scal);
        xSparseVector & operator/=(const T &scal);
        xSparseVector & operator=(const xSparseVector &in);
        xSparseVector operator+(const xSparseVector & rhs) const;
        xSparseVector operator*(const T &rhs) const;
        xSparseVector operator/(const T &rhs) const;

        // basic iterator
        iter_idx_t beginIdx() const;
        iter_idx_t endIdx() const;
        iter_val_t beginVal() const;
        iter_val_t endVal() const;

        // size
        size_t nnz() const;

        // controled reset
        template <typename ITI,typename ITV>
        void resetBy(ITI bi,ITI bie,ITV bv, ITV bve);
        void clear();

    private:
        void add(int index, T val);
        std::vector < int > idx; 
        std::vector < T > values; 
};

} //end namespace xlinalg

#include  "xSparseVector_imp.h"
#endif
