#ifndef _XVECTOR_H_
#error Do NOT include xvector_imp.h alone
#endif

#include "xDataType.h"

/// left multiply by a scalar
template <typename VT>
inline xVector<VT> operator*(const VT &l ,const  xVector<VT> &r);

/// output stream operator (write ...)
template <typename VT>
inline std::ostream& operator<<(std::ostream& s, const xVector<VT>& v);

/// input stream operator (read ...)
template <class charT, class traits, typename VT>
inline
std::basic_istream<charT, traits> & operator >> (std::basic_istream<charT, traits> & strm, xVector<VT>& v);

template <typename VT>
inline xVector<VT> operator*(const VT &l ,const  xVector<VT> &r){
    return r*l;
}


template <typename VT>
inline VT xVector<VT>::mag() const
{
    return sqrt(pos[0]*pos[0]+pos[1]*pos[1]+pos[2]*pos[2]);
}

template <typename VT>
inline VT  xVector<VT>::normValue()
{
    const VT n = mag();
    if(n == xtool::xDataType<VT>::zero())return n;
    (*this)/=n;
    return n;
}

template <typename VT>
inline xVector<VT>&  xVector<VT>::norm()
{
    const VT n = mag();
    if(n == xtool::xDataType<VT>::zero())return *this;
    (*this)/=n;
    return *this;
}

template <typename VT>
inline VT & xVector<VT>::operator () (int i)
{
    if(i>=3)throw -1;
    return pos[i];
}

template <typename VT>
inline const VT& xVector<VT>::operator () (int i) const
{
    return pos[i];
}

template <typename VT>
inline VT & xVector<VT>::operator [] (int i)
{

    if(i>=3)throw -1;

    return pos[i];
}

template <typename VT>
inline const VT& xVector<VT>::operator [] (int i) const
{

    if(i>=3)throw -1;

    return pos[i];
}

template <typename VT>
inline xVector<VT> xVector<VT>::operator - () const
{
    return xVector(-pos[0], -pos[1], -pos[2]);
}

template <typename VT>
inline xVector<VT> xVector<VT>::operator - (const xVector<VT> &other) const
{
    return (xVector(*this) -= other);
}

// cross product
template <typename VT>
inline xVector<VT> xVector<VT>::operator % (const xVector<VT> &other) const
{
    return xVector(	pos[1]*other.pos[2]-pos[2]*other.pos[1],
            pos[2]*other.pos[0]-pos[0]*other.pos[2],
            pos[0]*other.pos[1]-pos[1]*other.pos[0]);
}

template <typename VT>
inline xVector<VT> xVector<VT>::operator + (const xVector<VT> &other) const
{
    return (xVector<VT>(*this) += other);
}

template <typename VT>
inline xVector<VT> xVector<VT>::operator * (const VT &other) const
{
    return (xVector<VT>(*this) *= other);
}

template <typename VT>
inline xVector<VT> xVector<VT>::operator / (const VT &other) const
{
    return (xVector<VT>(*this) /= other);
}

template <typename VT>
inline xVector<VT> & xVector<VT>::operator += (const xVector<VT> &other)
{
    pos[0]+=other.pos[0];pos[1]+=other.pos[1];pos[2]+=other.pos[2];
    return *this;
}

template <typename VT>
inline xVector<VT> & xVector<VT>::operator /= (const VT &other)
{
    const VT d = 1./other;
    (*this) *= d;
    return *this;
}

template <typename VT>
inline xVector<VT> & xVector<VT>::operator -= (const xVector<VT> &other)
{
    pos[0]-=other.pos[0];pos[1]-=other.pos[1];pos[2]-=other.pos[2];
    return *this;
}

template <typename VT>
inline xVector<VT> & xVector<VT>::operator *= (const VT &other)
{
    pos[0]*=other;pos[1]*=other;pos[2]*=other;
    return *this;
}

template <typename VT>
inline double xVector<VT>::angleRad(const xVector<VT> &v) const
{
    // have to copy not to modify original vectors
    xVector<VT> y(v);
    xVector<VT> x(*this);
    x.norm();
    y.norm();
    const double cosA = x * y;
    const double sinA = ( x % y).mag();
    return atan2(sinA,cosA);
}

template <typename VT>
inline double xVector<VT>::angleDeg(const xVector<VT> &v) const
{
    return angleRad(v) * (180./M_PI);
}

template <typename VT>
inline std::ostream& operator<<(std::ostream& s, const xVector<VT>& v)
{
    s << v(0) << " " << v(1) << " " << v(2);
    return s;
}


template <class charT, class traits, typename VT>
inline
std::basic_istream<charT, traits> & operator >> (std::basic_istream<charT, traits> & strm, xVector<VT>& v){
    VT x{xtool::xDataType<VT>::zero()},y{xtool::xDataType<VT>::zero()},z{xtool::xDataType<VT>::zero()};
    //char a;
    int nopen=0;
    // std::cout << "reading vector" << std::endl;
    // strm>> std::ws;
    //std::cout << "ppeke " << strm.peek()<< std::endl;
    while ((strm.peek()=='{') ||(strm.peek()==' ' )||(strm.peek()=='(' )){
        if ((strm.peek()=='{')||(strm.peek()=='(' )) nopen++;
        strm.ignore();
    }
    strm>> x;
    // std::cout << " x " <<  x << std::endl;
    while ((strm.peek()==',') ||(strm.peek()==' ' ) ||(strm.peek()==';' )){
        strm.ignore();
    }
    strm>> y;
    //std::cout << "y " <<  y << std::endl;
    while ((strm.peek()==',') ||(strm.peek()==' ' ) ||(strm.peek()==';' )){
        strm.ignore();
    }
    strm>>z;
    //std::cout << " z " <<  z << std::endl;
    while (nopen!=0&&strm.good()){
        if ((strm.peek()=='}')||(strm.peek()==')' )) nopen--;
        strm.ignore();
    }
    if (strm.fail()) {
        std::cout << "a parentezes was open but never closed while reading a xVector" << std::endl;
        throw;
    }
    v= xVector<VT>(x,y,z);
    return strm;
}



template<typename VT>
xVector<VT>& xVector<VT>::operator *= ( const xTensor2<VT> &other )
{
    // RIGHT multiply of a matrix to a vector: xTensor2 * xVector
    xVector<VT> m(*this);
    //flatten loops...
    pos[0] = other(0,0) * m.pos[0] + other(0,1) * m.pos[1] + other(0,2) * m.pos[2];
    pos[1] = other(1,0) * m.pos[0] + other(1,1) * m.pos[1] + other(1,2) * m.pos[2];
    pos[2] = other(2,0) * m.pos[0] + other(2,1) * m.pos[1] + other(2,2) * m.pos[2];
    return *this;
}

template<typename VT>
xVector<VT> xVector<VT>::operator * ( const xTensor2<VT>& M ) const
{
    // LEFT multiply of a vector to a matrix: xVector * xTensor2
    return !M * (*this);
}




