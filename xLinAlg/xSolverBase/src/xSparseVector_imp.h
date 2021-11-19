/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
 */

#ifdef XSPARSEVECTOR_H
#ifndef XSPARSEVECTOR_IMP__H
#define XSPARSEVECTOR__IMP_H

namespace xlinalg
{

template < typename T >
xSparseVector < T >::xSparseVector(){};

template < typename T >
xSparseVector < T >::xSparseVector(int index, T val) : idx(1,index),values(1,val){};

template < typename T >
xSparseVector < T >::xSparseVector(const xSparseVector < T > &in) = default;

template < typename T >
void xSparseVector < T >::add(int index, T val)
{
    idx.push_back(index);
    values.push_back(val);
}
template < typename T >
xSparseVector < T > & xSparseVector < T >::operator*=(const T &scal)
{
    for_each(values.begin(), values.end(), [&scal](T & val){ val *= scal; } );
    return *this;
}
template < typename T >
xSparseVector < T > & xSparseVector < T >::operator/=(const T &scal)
{
    for_each(values.begin(), values.end(), [&scal](T & val){ val /= scal; } );
    return *this;
}
template < typename T >
xSparseVector < T > & xSparseVector < T >::operator=(const xSparseVector < T > &in)
    = default;
template < typename T >
xSparseVector < T > xSparseVector < T >::operator+(const xSparseVector < T > & rhs) const
{
    auto it_idx = idx.begin();
    auto it_idx_end = idx.end();
    auto it_idx_rhs = rhs.idx.begin();
    auto it_idx_rhs_end = rhs.idx.end();
    auto it_val = values.begin();
    auto it_val_rhs = rhs.values.begin();
    xSparseVector < T > tmp;
    bool have_it = ( it_idx != it_idx_end );
    bool have_it_rhs = ( it_idx_rhs != it_idx_rhs_end );
    while (have_it || have_it_rhs)
    {

        if (have_it && have_it_rhs)
        {
            if (*it_idx < *it_idx_rhs)
            {
                tmp.add(*it_idx,*it_val);
                ++it_idx;
                ++it_val;
                have_it = ( it_idx != it_idx_end );
            }
            else if (*it_idx > *it_idx_rhs)
            {
                tmp.add(*it_idx_rhs,*it_val_rhs);
                ++it_idx_rhs;
                ++it_val_rhs;
                have_it_rhs = ( it_idx_rhs != it_idx_rhs_end );
            }
            else
            {
                tmp.add(*it_idx,( *it_val )+( *it_val_rhs ));
                ++it_idx;
                ++it_val;
                ++it_idx_rhs;
                ++it_val_rhs;
                have_it = ( it_idx != it_idx_end );
                have_it_rhs = ( it_idx_rhs != it_idx_rhs_end );
            }
        }
        else if (have_it)
        {
            tmp.idx.insert(tmp.idx.end(),it_idx,it_idx_end);
            tmp.values.insert(tmp.values.end(),it_val,values.end());
            break;
        }
        else
        {
            assert(have_it_rhs);
            tmp.idx.insert(tmp.idx.end(),it_idx_rhs,it_idx_rhs_end);
            tmp.values.insert(tmp.values.end(),it_val_rhs,rhs.values.end());
            break;
        }
    }

    return tmp;
}
template < typename T >
xSparseVector < T > xSparseVector < T >::operator*(const T &rhs) const
{
    xSparseVector < T > lhs(*this);
    return lhs *= rhs;
}
template < typename T >
xSparseVector < T > xSparseVector < T >::operator/(const T &rhs) const
{
    xSparseVector < T > lhs(*this);
    return lhs /= rhs;
}
template < typename T >
xSparseVector < T > operator*(const T &lhs,  xSparseVector < T > rhs)
{
    return rhs *= lhs;
}

template < typename T >
auto xSparseVector < T >::beginIdx() const->iter_idx_t
{ return idx.cbegin(); }

template < typename T >
auto xSparseVector < T >::endIdx() const->iter_idx_t
{ return idx.cend(); }

template < typename T >
auto xSparseVector < T >::beginVal() const->iter_val_t
{ return values.cbegin(); }

template < typename T >
auto xSparseVector < T >::endVal() const->iter_val_t
{ return values.cend(); }

template < typename T >
size_t xSparseVector < T >::nnz() const
{
    assert(idx.size() == values.size());
    return idx.size();
}
template < typename T >
template < typename ITI,typename ITV >
void xSparseVector < T >::resetBy(ITI bi,ITI bie,ITV bv, ITV bve)
{
    size_t nb = std::distance(bi,bie);
    assert(std::distance(bv,bve) == nb);
    if (nb)
    {
        idx.resize(nb);
        values.resize(nb);
        int *pi=idx.data();
        *pi=*bi;
        for (auto id : xtool::make_range(bi+1,bie))
        {
            if (id > *pi)
            {
                *(++pi)=id;
            }
            else
            {
                pi=nullptr;
                break;
            }
        }
        assert(pi);
        std::copy(bv,bve,values.begin());
    }
    else
    {
        idx.clear();
        values.clear();
    }
}
template < typename T >
void xSparseVector < T >::clear()
{
    idx.clear();
    values.clear();
    return;
}

}         //end namespace
#endif
#endif
