/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/

#ifndef H_Products
#define H_Products

#include <iostream>
#include "xGenericOperations.h"// DANS XTOOL !!!
#include "xVector.h"
#include "xTensor2.h"
#include "xTensor2Sym.h"
#include "xTensor4.h"


namespace xtensor{


template<typename VT = double>
class xTensorProduct : public std::binary_function<xVector<VT>, xVector<VT>, xTensor2<VT> >
{
public:
    typedef typename std::binary_function<xVector<VT>, xVector<VT>, xTensor2<VT> >::result_type result_type;
    typedef typename std::binary_function<xVector<VT>, xVector<VT>, xTensor2<VT> >::first_argument_type first_argument_type;
    typedef typename std::binary_function<xVector<VT>, xVector<VT>, xTensor2<VT> >::second_argument_type second_argument_type;
    result_type operator()(const first_argument_type& f, const second_argument_type& s) const {return tensor_product(f,s);  }
};

template<typename VT = double>
class xTensorContraction : public std::binary_function<xTensor2<VT>, xTensor2<VT>, double>
{
public:
    typedef typename std::binary_function<xTensor2<VT>, xTensor2<VT>, double>::result_type result_type;
    typedef typename std::binary_function<xTensor2<VT>, xTensor2<VT>, double>::first_argument_type first_argument_type;
    typedef typename std::binary_function<xTensor2<VT>, xTensor2<VT>, double>::second_argument_type second_argument_type;
    result_type operator()(const first_argument_type& f, const second_argument_type& s) const {return f.contract(s);  }
};

template <class T, typename VT = double>
class xContraction : public std::binary_function<T, T, VT>
{
public:
    typedef typename std::binary_function<T, T, VT>::result_type result_type;
    typedef typename std::binary_function<T, T, VT>::first_argument_type first_argument_type;
    typedef typename std::binary_function<T, T, VT>::second_argument_type second_argument_type;
    
    result_type operator()(const first_argument_type& f, const second_argument_type& s) const {return f*s;  }
};

template <typename VT>
class xContraction<xTensor2<VT>, VT > : public std::binary_function<xTensor2<VT>, xTensor2<VT>, VT>
{
public:
    typedef typename std::binary_function<xTensor2<VT>, xTensor2<VT>, VT>::result_type result_type;
    typedef typename std::binary_function<xTensor2<VT>, xTensor2<VT>, VT>::first_argument_type first_argument_type;
    typedef typename std::binary_function<xTensor2<VT>, xTensor2<VT>, VT>::second_argument_type second_argument_type;
    
    result_type operator()(const first_argument_type& f, const second_argument_type& s) const {return f.contract(s);  }
};

template <typename VT>
class xContraction<xTensor2Sym<VT>, VT > : public std::binary_function<xTensor2Sym<VT>, xTensor2Sym<VT>, double>
{
public:
    typedef typename std::binary_function<xTensor2Sym<VT>, xTensor2Sym<VT>, double>::result_type result_type;
    typedef typename std::binary_function<xTensor2Sym<VT>, xTensor2Sym<VT>, double>::first_argument_type first_argument_type;
    typedef typename std::binary_function<xTensor2Sym<VT>, xTensor2Sym<VT>, double>::second_argument_type second_argument_type;
    
    result_type operator()(const first_argument_type& f, const second_argument_type& s) const {return f.contract(s);  }
};

template< typename T >
struct xCastToComplex : public std::unary_function< T, std::complex<T> >
{
public:
    typedef std::unary_function< T, std::complex<T> >               base_type;
    typedef typename base_type::argument_type         argument_type;
    typedef typename base_type::result_type           result_type;

    result_type operator()( argument_type& t ) const { return result_type(t); }
    const result_type operator()( const argument_type& t ) const { return result_type(t); }
};

template< typename T >
struct xCastToVectorComplex : public std::unary_function< xtensor::xVector<T>, xtensor::xVector<std::complex<T> > >
{
public:
    typedef std::unary_function< xtensor::xVector<T>, xtensor::xVector<std::complex<T> > >               base_type;
    typedef typename base_type::argument_type         argument_type;
    typedef typename base_type::result_type           result_type;

    result_type operator()( argument_type& t ) const { return result_type(t(0), t(1), t(2)); }
    const result_type operator()( const argument_type& t ) const { return result_type(t(0), t(1), t(2)); }
};

template< typename T >
struct xCastToTensor2Complex : public std::unary_function< xtensor::xTensor2<T>, xtensor::xTensor2<std::complex<T> > >
{
public:
    typedef std::unary_function< xtensor::xTensor2<T>, xtensor::xTensor2<std::complex<T> > >               base_type;
    typedef typename base_type::argument_type         argument_type;
    typedef typename base_type::result_type           result_type;

    result_type operator()( argument_type& t ) const { return result_type(t(0,0), t(0,1), t(0,2),
                                                                          t(1,0), t(1,1), t(1,2),
                                                                          t(2,0), t(2,1), t(2,2)); }
    const result_type operator()( const argument_type& t ) const { return result_type(t(0,0), t(0,1), t(0,2),
                                                                                      t(1,0), t(1,1), t(1,2),
                                                                                      t(2,0), t(2,1), t(2,2)); }
};

template <typename VT = double>
class xExtractCompVector : public std::unary_function<xVector<VT>, VT>
{
public:
    xExtractCompVector(const int l) : indice(l) {}
    VT operator()(const xVector<VT>& vect) const { return vect(indice); }

private :
    const int indice;
};

template <typename VT = double>
class xAddCompVector : public std::unary_function<VT, xVector<VT> >
{
public:
    xAddCompVector(const int l) : indice(l) {}
    xVector<VT> operator()(const VT& s) const
    {
        xVector<VT> vect;
        vect(indice)=s;
        return vect;
    }
private:
    const int indice;
};

template <typename VT = double>
class xExtractNormVector :  public std::unary_function<xVector<VT>, VT>{
public :
    //xExtractNormVector() {}
    VT operator()(const xVector<VT>& vect) const {return vect.mag();}
};

template <typename VT = double>
class xNormalizeVector : public std::unary_function<xVector<VT>, xVector<VT> >{
public :
    xVector<VT> operator()(const xVector<VT>& vect) const {
        xVector<VT> normed(vect);
        normed.norm();
        return normed;
    }
}; 


//template <class T1, class T2,class T3>
template <typename VT = double>
class xExtractCompTensor : public std::unary_function<xTensor2<VT>, VT>
{
public:
    xExtractCompTensor(const int l, const int c) : ligne(l), col(c) {}
    VT operator()(const xTensor2<VT>& tens) const { return tens(ligne,col); }

private :
    const int ligne;
    const int col;

};

template <typename VT = double>
struct xNorm : std::unary_function<xVector<VT>, VT>
{
    VT operator () ( xVector<VT> const &arg ) const
    {
        return arg.mag();
    }
};

template <typename VT = double>
struct xNormalize : std::unary_function<xVector<VT>, xVector<VT> >
{
    xVector<VT> operator () ( xVector<VT> const &arg ) const
    {
        xVector<VT> tmp(arg);
        return tmp.norm();
    }
};

template <typename VT = double>
class xSymmetrize : public std::unary_function<xTensor2<VT>, xTensor2<VT> >
{
public: 
    xTensor2<VT>& operator()(xTensor2<VT>& in) const
    {
        in(0,1) = 0.5 * (in(0,1) + in(1,0)); in(1,0) = in(0,1);
        in(0,2) = 0.5 * (in(0,2) + in(2,0)); in(2,0) = in(0,2);
        in(1,2) = 0.5 * (in(1,2) + in(2,1)); in(2,1) = in(1,2);
        return in;
    }
};

template <typename VT = double>
class xStrainPlaneStress : public std::unary_function<xTensor2<VT>, xTensor2<VT> >
{
public: 
    xStrainPlaneStress(VT nu)
    {
        k = -nu/(1-nu);
    }
    xTensor2<VT>& operator()(xTensor2<VT>& in) const
    {
        in(0,1) = 0.5 * (in(0,1) + in(1,0)); in(1,0) = in(0,1);
        in(0,2) = 0.5 * (in(0,2) + in(2,0)); in(2,0) = in(0,2);
        in(1,2) = 0.5 * (in(1,2) + in(2,1)); in(2,1) = in(1,2);
        in(2,2) = k * (in(0,0)+in(1,1));
        return in;
    }
private :
    double k;
};

template <typename VT = double>
class xSymmetrizeStrong : public std::unary_function<xTensor2<VT>, xTensor2Sym<VT>>
{
public: 
    xTensor2Sym<VT> operator()(const xTensor2<VT>& in) const
    {
        return xTensor2Sym<VT>(in);
    }
};

template <typename VT = double>
struct xAntiSymmetrize : public std::unary_function<xTensor2<VT>, xTensor2<VT> >
{
    xTensor2<VT>& operator()(xTensor2<VT>& in) const
    {
        in(0,1) = 0.5 * (in(0,1) - in(1,0)); in(1,0) = - in(0,1);
        in(0,2) = 0.5 * (in(0,2) - in(2,0)); in(2,0) = - in(0,2);
        in(1,2) = 0.5 * (in(1,2) - in(2,1)); in(2,1) = - in(1,2);
        in(0,0) = 0.;
        in(1,1) = 0.;
        in(2,2) = 0.;
        return in;
    }
};

template <typename VT = double>
class xDeviatoric : public std::unary_function<xTensor2<VT>,xTensor2<VT> >
{
public: 
    xTensor2<VT>&  operator()(xTensor2<VT>& in) const
    {
        VT volumetric= in(0,0)+in(1,1)+in(2,2);
        in(0,0) -= volumetric/3.;
        in(1,1) -= volumetric/3.;
        in(2,2) -= volumetric/3.;
        return in;
    }
};

template <typename VT = double>
class xTrace : public std::unary_function<xTensor2<VT>, VT>
{
public: 
    VT   operator()(xTensor2<VT>& in) const
    {
        return in(0,0)+in(1,1)+in(2,2);
    }
};

/// Use to compute the second invariant of a tensor c  (IIc =  0.5* tr (c^T*c) )
template <typename VT = double>
class xSecondInvariant : public std::unary_function<xTensor2<VT>, VT>
{
public:
    VT operator()(xTensor2<VT>& in) const
    {
        return 0.5*(pow(in(0,0),2)+pow(in(1,1),2)+pow(in(2,2),2)
                    +pow(in(0,1),2)+pow(in(0,2),2)+pow(in(1,2),2)
                    +pow(in(1,0),2)+pow(in(2,0),2)+pow(in(2,1),2));
    }
};

template <typename VT = double>
class xTranspose : public std::unary_function<xTensor2<VT>, xTensor2<VT> >
{
public: 
    xTensor2<VT>&  operator()(xTensor2<VT>& in) const
    {
        std::swap(in(0,1),in(1,0));
        std::swap(in(0,2),in(2,0));
        std::swap(in(1,2),in(2,1));
        return in;
    }
};


/// Use to compute Von Mises Norm of a xTensor2<VT>
template <typename VT = double>
class xVonMisesNorm : public std::unary_function<xTensor2<VT>, VT>
{
public:
    VT operator()(xTensor2<VT>& in) const
    {
        return in.vonMisesNorm();
    }
};

/// Used to compute  one Principal component
/// Only for double
class xPrincipalComponent : public std::unary_function<xTensor2<double>, double>
{
public:
    xPrincipalComponent(const int component): c(component) {}
    double operator()(xTensor2<double>& in) const
    {
        double val[3];
        try
        {
            in.getAnalyticalEigenvaluesAndDerivative(val);
            std::sort(val,val+3);
            if (fabs(val[c]) < std::numeric_limits<double>::min()) val[c] = 0.0;
        }
        catch(int e)
        {
            val[0]=val[1]=val[2]=0.;
            if (e!=-1) throw -1;

        }
        return val[c];
    }
private :
    const unsigned  int c;
};


} // end of namespace


#endif


