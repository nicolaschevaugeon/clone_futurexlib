/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/

#ifndef _XTENSOR3_H_
#define _XTENSOR3_H_

#include "xTensor2.h"
#include "xVector.h"

namespace xtensor {


template<typename VT = double>
class xTensor3
{
    VT pos[3][3][3];
public:
    xTensor3 ()= default;
    xTensor3 (const VT c[3][3][3]);

    inline VT operator () (int i,int j, int k) const;
    inline VT & operator () (int i,int j, int k);

    inline xVector<VT> operator * (const xTensor2<VT> &other) const
    {
        xVector<VT> v(0.0);
        for(int i=0;i<3;i++)
        {
            // 	v(i) = 0.0;//Inutile, deja initialise a zero
            for(int j=0;j<3;j++)
            {
                for(int k=0;k<3;k++)
                {
                    v(i) += pos[i][j][k] * other(j,k) ;
                }
            }
        }
        return v;
    }

    inline xTensor2<VT> operator * (const xVector<VT> &other) const //tens3*vec
    {
        xTensor2<VT> m(0.0);
        for(int i=0;i<3;i++)
        {
            for(int j=0;j<3;j++)
            {
                //           m(i,j) = 0.0;//Inutile, deja initialise a zero
                for(int k=0;k<3;k++)
                {
                    m(i,j) += pos[k][i][j] * other(k) ;
                }
            }
        }
        return m;
    }


    inline xTensor2<VT> multtransp(const xVector<VT> &other) const // tens3^t*vect
    {
        xTensor2<VT> m(0.0);
        for(int i=0;i<3;i++)
        {
            for(int j=0;j<3;j++)
            {
                //           m(i,j) = 0.0;//Inutile, deja initialise a zero
                for(int k=0;k<3;k++)
                {
                    m(i,j) += pos[k][i][j] * other(k) ;
                }
            }
        }
        return m;
    }
};

template<typename VT>
inline xTensor3<VT>::xTensor3 (const VT c[3][3][3])
{
    int i, j, k;
    for (i = 0; i < 3; ++i){
        for (j = 0; j < 3; ++j){
            for (k = 0; k < 3; ++k){
                pos[i][j][k] = c[i][j][k];
            }
        }
    }
}

template<typename VT>
inline VT xTensor3<VT>::operator () (int i,int j, int k) const
{
    return pos[i][j][k];
}

template<typename VT>
inline VT & xTensor3<VT>::operator () (int i,int j, int k)
{
    return pos[i][j][k];
}


template<typename VT>
inline std::ostream& operator<<(std::ostream& s, const xtensor::xTensor3<VT>& t)
{
    s << t(0,0,0) << " " << t(0,0,1) << " " << t(0,0,2) << std::endl;
    s << t(0,1,0) << " " << t(0,1,1) << " " << t(0,1,2) << std::endl;
    s << t(0,2,0) << " " << t(0,2,1) << " " << t(0,2,2) << std::endl;
    s << std::endl;
    s << t(1,0,0) << " " << t(1,0,1) << " " << t(1,0,2) << std::endl;
    s << t(1,1,0) << " " << t(1,1,1) << " " << t(1,1,2) << std::endl;
    s << t(1,2,0) << " " << t(1,2,1) << " " << t(1,2,2) << std::endl;
    s << std::endl;
    s << t(2,0,0) << " " << t(2,0,1) << " " << t(2,0,2) << std::endl;
    s << t(2,1,0) << " " << t(2,1,1) << " " << t(2,1,2) << std::endl;
    s << t(2,2,0) << " " << t(2,2,1) << " " << t(2,2,2) << std::endl;
    return s;
}

template<typename VT>
inline  xVector<VT> operator *(xTensor2<VT> &t1,xTensor3<VT> &t2)
{
    xVector<VT> m;
    //   cout << t2 << endl;
    for(int i=0;i<3;i++){
        //       m(i) = 0.0;//Inutile, deja initialise a zero
        for(int k=0;k<3;k++){
            for(int l=0;l<3;l++){
                m(i) += t1(k,l) * t2(i,k,l);
            }
        }
    }
    //       cout << endl << m << endl << endl;
    return m;
}

template<typename VT>
inline  xTensor2<VT> operator *(xVector<VT> &t1,xTensor3<VT> &t2)
{
    xTensor2<VT> m;
    //   cout << t2 << endl;
    for(int i=0;i<3;i++){
        for(int j=0;j<3;j++){
            //         m(i,j) = 0.0;//Inutile, deja initialise a zero
            for(int k=0;k<3;k++){
                m(i,j) += t1(k) * t2(k,i,j);
            }
        }
    }
    //       cout << endl << m << endl << endl;
    return m;
}


} // end of namespace






#endif
