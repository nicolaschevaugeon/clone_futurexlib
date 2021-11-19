/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/

//////////////////////////////////////////////////////////////////////
#ifndef _XTensor2Sym_H_
#define _XTensor2Sym_H_

#include <iostream>
#include "xVector.h"
//GILLES: #include "mTensor2.h"
#include "xTensor4.h"

namespace xtensor {

template<typename VT>
class xTensor4;

class xTensor4Isotropic;

/**
     Simple class for a 2nd order tensor in 3d which
     has a 3x3 matrix representation.
  */

template<typename VT = double>
class xTensor2Sym
{
    friend class xTensor4<VT>;
    friend class xTensor4Isotropic;
    template<typename UT>
    friend void axpy_tensor2Sym (const UT& a, const xTensor2Sym<UT>& x, xTensor2Sym<UT>& y);
    VT pos[6];
    //xVector pos[3];
public:
    //inline xTensor2Sym(const xVector &c1, const xVector &c2, const xVector &c3);
    inline xTensor2Sym(const VT& init = 0.0)
    {
        VT *p = &pos[0];
        for (unsigned int i=0; i<6; ++i)
            *p = init;
    }


    /*Greg Add*/

    inline VT delta(int i, int j) const {return ((i == j)?1.0:0.0);}
    


    // New Implementation
    //
    //
    inline xTensor2Sym(const xTensor2Sym<VT>& in)
    {
        const VT *r = &in.pos[0];
        VT *l = &pos[0];
        for(unsigned int i=0;i<6;++i)
        {
            *l = *r;
            ++l; ++r;
        }
    }
    
    inline xTensor2Sym(const xTensor2<VT>& in)
    {
        for (unsigned int i=0; i< 3; ++i)
            pos[i] = in(i,i);
        pos[3] = 0.5* (in(0,1) + in(1,0));
        pos[4] = 0.5* (in(0,2) + in(2,0));
        pos[5] = 0.5* (in(1,2) + in(2,1));
    }

    
    
    inline xTensor2Sym(const xTensor2Sym<VT>& in, VT scal)
    {
        const VT *r = &in.pos[0];
        VT *l = &pos[0];
        for(unsigned int i=0;i<6;++i)
        {
            *l = *r * scal;
            ++l; ++r;
        }
    }
    
    inline xTensor2Sym &operator *= (const VT &scalar)
    {
        VT * l = &pos[0];
        for (unsigned int i=0; i< 6;  ++i){
            (*l) *= scalar;
            ++l;
        }
        return *this;
    }
    
    inline xTensor2Sym operator * (const VT &scalar) const
    {
        xTensor2Sym<VT> other(*this);
        other *=scalar;
        return other;
    }
    
    inline xTensor2Sym &operator /= (const VT &scalar)
    {
        VT * l = &pos[0];
        for (unsigned int i=0; i< 6;  ++i){
            *l /= scalar;
            ++l;
        }
        return *this;
    }
    
    inline xTensor2Sym &operator += (const xTensor2Sym<VT> &other)
    {
        VT *l = &pos[0];
        const VT *r  = &other.pos[0];
        for(unsigned int i=0;i<6;++i)
        {
            *l +=  *r;
            ++l; ++r;
        }
        return *this;
    }

    inline xTensor2Sym &operator -= (const xTensor2Sym<VT> &other)
    {
        VT *l = &pos[0];
        const VT *r  = &other.pos[0];
        for(unsigned int i=0;i<6;++i)
        {
            *l -=  *r;
            ++l; ++r;
        }
        return *this;
    }
    
    inline xTensor2Sym operator +  (const xTensor2Sym<VT> &other) const
    {
        xTensor2Sym<VT> m(*this) ;
        m += other;
        return m;
    }
    
    inline xTensor2Sym operator -  (const xTensor2Sym<VT> &other) const
    {
        xTensor2Sym<VT> m(*this);
        m -= other;
        return m;
    }
    
    inline double contract (const xTensor2Sym<VT> &other) const
    {
        VT x= 0.;
        const VT * l = &pos[0];
        const VT * r = &other.pos[0];

        for(unsigned int i=0;i<3;++i)
        {
            x += (*l++) * (*r++);
        }
        for(unsigned int i=3;i<6;++i)
        {
            x += 2.* (*l++) * (*r++);
        }

        return x;
    }
    
    inline xTensor2Sym operator * (const xTensor4Isotropic &other) const
    {
        xTensor2Sym<VT> m( *this, 2*other.mu );
        VT lam =  other.lam;
        VT t = trace()*lam;
        for(unsigned int i=0;i<3;i++){
            m.pos[i] +=  t;
        }
        return m;
    }
    
    inline xTensor2Sym & operator *= (const xTensor4Isotropic &other)
    {
        *this *= 2*other.mu;
        VT lam =  other.lam;
        VT t = trace()*lam;
        for(unsigned int i=0;i<3;i++){
            pos[i] +=  t;
        }
        return *this;
    }
    

    /*    inline xTensor2(const Trellis_Util::mTensor2& in)
    {
      for(unsigned int i=0;i<3;i++)
    {
      for(unsigned int j=0;j<3;j++)
        {
          pos[i][j] = in(i,j);
        }
    }
    }*/


    double & operator() ( int, int );
    const double& operator() ( int, int ) const;
    inline VT trace () const
    {
        return pos[0] + pos[1] + pos[2];
    }

    //double trace2 () const;
    //long eigen   (xVector e[3],double v[3]) const;
    long eigen2d (xVector<VT> e[3],VT v[3]) const;//NE FONCTIONNE PAS ??? IL FAUDRAIT DES REFERENCES SUR LES ARGUMENTS ???
};



template<typename VT>
inline double &xTensor2Sym<VT>::operator () (int i, int j)
{
    if (i>j)
    {
        unsigned int tmp = j;
        j=i;
        i=tmp;
    }
    if (i==j) return pos[i];
    else return pos[3+i+j];
}

template<typename VT>
inline const double& xTensor2Sym<VT>::operator () (int i, int j) const
{
    if (i>j)
    {
        unsigned int tmp = j;
        j=i;
        i=tmp;
    }
    if (i==j) return pos[i];
    else return pos[3+i+j];
}



// Fonctions de sortie :

template<typename VT>
inline std::ostream& operator<<(std::ostream& s, const xTensor2Sym<VT>& t)
{
    s << t(0,0) << " " << t(0,1) << " " << t(0,2) << std::endl;
    s << t(1,0) << " " << t(1,1) << " " << t(1,2) << std::endl;
    s << t(2,0) << " " << t(2,1) << " " << t(2,2) << std::endl;
    return s;
}

template<typename UT>
inline void axpy_tensor2Sym (const UT& a, const xTensor2Sym<UT>& x, xTensor2Sym<UT>& y)
{
    UT *py = &y.pos[0];
    const UT *px = &x.pos[0];
    for (unsigned int i =0; i<6 ;++i){
        *py += a * (*px);
        ++py; ++px;
    }
}


//Only specialized for double
//Note: inline is used in order to avoid multiple definitions linking error (cf https://stackoverflow.com/a/35017825)
template<>
inline long xTensor2Sym<double>::eigen2d(xVector<double> e[3], double v[3]) const
{
    e[2] = xVector<double>(0,0,1);

    double b = -pos[0]-pos[1];
    double c = (pos[0]*pos[1]-pos[3]*pos[3]);
    double delta = b*b-4*c;

    if (delta < 0){
        const double tolerance = std::numeric_limits<double>::epsilon()*100.;
        if (delta>-tolerance)
            delta=0.;
        else
            return 0;
    }

    v[0] = (-b+sqrt(delta))/2.;
    v[1] = (-b-sqrt(delta))/2.;
    v[2] = 1.0;

    long nbEigen = 2;

    if (fabs(v[1]) > fabs(v[0]))
    {
        double temp = v[0];
        v[0] = v[1];
        v[1] = temp;
    }

    double result[4];
    unsigned int nb_vec=0;
    while(1)
    {
        double a[4] = {pos[0]-v[nb_vec],pos[3],
                   pos[3],pos[1]-v[nb_vec]};
        double eps = 1.e-8;
        unsigned int nb = 0;
        while (1)
        {
            nb = NullSpace (a,result,eps,2);
            if (nb != 0)break;
            eps *= 2.0;
        }
        unsigned int kk=0;
        for (unsigned int i=nb_vec;i<nb+nb_vec;i++)
        {
            e[i] = xVector<double> (result[0+kk*2],result[1+kk*2],0.0);
            e[i].norm();
            kk++;
        }
        nb_vec += nb;
        if (nb_vec == 2)return nbEigen;
        if (nb > 2)throw -1;
    }
    throw -1;
}


} // end of namespace

#endif
