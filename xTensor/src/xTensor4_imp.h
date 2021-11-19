#ifndef _XTensor4_H_
#error Do NOT include xTensor4_imp.h alone
#endif

#include "xDataType.h"


// #########################################################################################
// IMPLEMENTATION of inlined functions
// #########################################################################################

inline xTensor4Isotropic::xTensor4Isotropic () : lam(0.), mu(0.) {}

inline xTensor4Isotropic::xTensor4Isotropic (const double& E, const double& nu)
{
    lam = nu*E/((1.+nu)*(1.-2.*nu));
    mu  = E/(2.*(1.+nu));
}

inline xTensor2<double>  xTensor4Isotropic::operator * (const xTensor2<double> &other) const
{
    xTensor2<double> m(0.);
    double t = other.trace();
    for(unsigned int i=0;i<3;i++){
        for(unsigned int j=0;j<3;j++){
            m(i,j) = lam * t * this->delta(i,j) + mu * (other.pos[i][j] + other.pos[j][i]);
        }
    }
    return m;
}

inline xTensor4Isotropic xTensor4Isotropic::operator * (const double &other) const
{
    xTensor4Isotropic m(0., 0.);
    m.lam = lam*other;
    m.mu = mu*other;
    return m;
}

inline xTensor4Isotropic &xTensor4Isotropic::operator *= (const double &other)
{
    this->lam *= other;
    this->mu *= other;
    return *this;
}


inline bool choleskyL(const xTensor4Isotropic &A4, xTensor4<double> &L, const double prec ){
    const double lambda = A4.getLambda();
    const double mu  = A4.getMu();
    const double Asquare = lambda+2*mu;
    if (Asquare < 0) return false;
    const double A = sqrt(Asquare);
    const double B = lambda/A ;
    if (mu < 0) return false;
    const double C = sqrt(mu);
    const double Dsquare = Asquare - B*B;
    if ( Dsquare < 0) return false;
    const double D = sqrt(Dsquare);
    const double E = (lambda - B*B)/D;
    const double Fsquare = Dsquare - E*E;
    if ( Fsquare < 0) return false;
    const double F = sqrt(Fsquare);
    L = xTensor4<double>(0.);
    L(0,0,0,0) = A;
    L(1,1,0,0) = B;
    L(2,2,0,0) = B;
    L(0,1,0,1) = C;
    L(1,0,0,1) = C;
    L(0,2,0,2) = C;
    L(2,0,0,2) = C;
    L(1,2,1,2) = C;
    L(2,1,1,2) = C;
    L(1,1,1,1) = D;
    L(2,2,1,1) = E;
    L(2,2,2,2) = F;
    return true;
}

inline bool choleskyU(const xTensor4Isotropic &A4, xTensor4<double> &U, const double prec){
    const double lambda = A4.getLambda();
    const double mu  = A4.getMu();
    const double Asquare = lambda+2*mu;
    if (Asquare < 0) return false;
    const double A = sqrt(Asquare);
    const double B = lambda/A ;
    if (mu < 0) return false;
    const double C = sqrt(mu);
    const double Dsquare = Asquare - B*B;
    if ( Dsquare < 0) return false;
    const double D = sqrt(Dsquare);
    const double E = (lambda - B*B)/D;
    const double Fsquare = Dsquare - E*E;
    if ( Fsquare < 0) return false;
    const double F = sqrt(Fsquare);
    U = xTensor4<double>(0.);
    U(0,0,0,0) = A;
    U(0,0,1,1) = B;
    U(0,0,2,2) = B;
    U(0,1,0,1) = C;
    U(0,1,1,0) = C;
    U(0,2,0,2) = C;
    U(0,2,2,0) = C;
    U(1,2,1,2) = C;
    U(1,2,2,1) = C;
    U(1,1,1,1) = D;
    U(1,1,2,2) = E;
    U(2,2,2,2) = F;
    return true;
}

template<typename VT>
inline xTensor4<VT>::xTensor4 () {}

template<typename VT>
inline xTensor4<VT>::xTensor4(const VT& init){
    VT *p = data();
    for (unsigned int i= 0; i< 81; ++i) *p++ = init;
}

template<typename VT>
inline xTensor4<VT>::xTensor4 (
        const VT &a0000, const VT &a0001,  const VT &a0002,
        const VT &a0010, const VT &a0011,  const VT &a0012,
        const VT &a0020, const VT &a0021,  const VT &a0022,

        const VT &a0100, const VT &a0101,  const VT &a0102,
        const VT &a0110, const VT &a0111,  const VT &a0112,
        const VT &a0120, const VT &a0121,  const VT &a0122,

        const VT &a0200, const VT &a0201,  const VT &a0202,
        const VT &a0210, const VT &a0211,  const VT &a0212,
        const VT &a0220, const VT &a0221,  const VT &a0222,

        const VT &a1000, const VT &a1001,  const VT &a1002,
        const VT &a1010, const VT &a1011,  const VT &a1012,
        const VT &a1020, const VT &a1021,  const VT &a1022,

        const VT &a1100, const VT &a1101,  const VT &a1102,
        const VT &a1110, const VT &a1111,  const VT &a1112,
        const VT &a1120, const VT &a1121,  const VT &a1122,

        const VT &a1200, const VT &a1201,  const VT &a1202,
        const VT &a1210, const VT &a1211,  const VT &a1212,
        const VT &a1220, const VT &a1221,  const VT &a1222,

        const VT &a2000, const VT &a2001,  const VT &a2002,
        const VT &a2010, const VT &a2011,  const VT &a2012,
        const VT &a2020, const VT &a2021,  const VT &a2022,

        const VT &a2100, const VT &a2101,  const VT &a2102,
        const VT &a2110, const VT &a2111,  const VT &a2112,
        const VT &a2120, const VT &a2121,  const VT &a2122,

        const VT &a2200, const VT &a2201,  const VT &a2202,
        const VT &a2210, const VT &a2211,  const VT &a2212,
        const VT &a2220, const VT &a2221,  const VT &a2222

        ):pos
{
{
{{a0000, a0001, a0002}, {a0010, a0011, a0012}, {a0020, a0021, a0022}},
{{a0100, a0101, a0102}, {a0110, a0111, a0112}, {a0120, a0121, a0122}},
{{a0200, a0201, a0202}, {a0210, a0211, a0212}, {a0220, a0221, a0222}}
              },

{
{{a1000, a1001, a1002}, {a1010, a1011, a1012}, {a1020, a1021, a1022}},
{{a1100, a1101, a1102}, {a1110, a1111, a1112}, {a1120, a1121, a1122}},
{{a1200, a1201, a1202}, {a1210, a1211, a1212}, {a1220, a1221, a1222}}
              },

{
{{a2000, a2001, a2002}, {a2010, a2011, a2012}, {a2020, a2021, a2022}},
{{a2100, a2101, a2102}, {a2110, a2111, a2112}, {a2120, a2121, a2122}},
{{a2200, a2201, a2202}, {a2210, a2211, a2212}, {a2220, a2221, a2222}}
              }

              }{}

template<typename VT>
inline xTensor4<VT>::xTensor4(const xTensor4& in){
    std::copy(reinterpret_cast<const VT *>(in.data()), reinterpret_cast<const VT *>(in.data())+81, reinterpret_cast<VT *>(data()));
}

template<typename VT>
inline xTensor4<VT> &xTensor4<VT>::operator = (const xTensor4 &other){
    std::copy(reinterpret_cast<const VT *>(other.data()), reinterpret_cast<const VT *>(other.data())+81, reinterpret_cast<VT *>(data()));
    return *this;
}

template<typename VT>
inline  xTensor4<VT>  xTensor4<VT>::operator -  (const xTensor4 &other) const
{
    xTensor4 m(*this);
    m -= other;
    return m;
}

template<typename VT>
inline xTensor4<VT> xTensor4<VT>::operator +  (const xTensor4 &other) const
{
    xTensor4 m(*this);
    m += other;
    return m;
}

template<typename VT>
inline xTensor4<VT>  xTensor4<VT>::operator *  (const VT &scalar) const
{
    xTensor4 m(*this);
    m *= scalar;
    return m;
}

template<typename VT>
inline xTensor4<VT> xTensor4<VT>::operator /  (const VT &scalar) const
{
    xTensor4 m(*this);
    m *= (1./scalar);
    return m;
}

template<typename VT>
inline xTensor2<VT> xTensor4<VT>::operator * (const xTensor2<VT> &other) const
{
    xTensor2<VT> m;
    //	Version With Only 1 loop.
    const VT *T4data = data();
    const VT *T2data = other.pos[0];
    VT * T2res = m.pos[0];
    for(unsigned int i=0;i<9;i++){
        *T2res =
                T4data[0]*T2data[0] + T4data[1]*T2data[1] + T4data[2]*T2data[2]
                +T4data[3]*T2data[3] + T4data[4]*T2data[4] + T4data[5]*T2data[5]
                +T4data[6]*T2data[6] + T4data[7]*T2data[7] + T4data[8]*T2data[8];
        T4data += 9;
        ++T2res;
    }
    return m;
}

template<typename VT>
inline xTensor4<VT> &xTensor4<VT>::operator /= (const VT &scalar)
{
    (*this) *= (1./scalar);
    return *this;
}

template<typename VT>
inline xTensor4<VT> &xTensor4<VT>::operator *= (const VT &scalar)
{

    for ( VT *dpos = data(); dpos < data()+81; ++dpos ) *dpos *= scalar;
    return *this;
}


/// return the transpose of this.
template<typename VT>
inline xTensor4<VT> xTensor4<VT>::operator !() const{
    xTensor4<VT> B;
    VT *pB = B.data();
    const VT *pA = data();
    for (unsigned int i=0; i < 9; ++i){
        for (unsigned int j=0; j < 9; ++j){
            pB[9*i+j] = pA[9*j+i];
        }
    }
    return B;
};



template<typename VT>
inline xTensor4<VT>::xTensor4 (const VT E, const VT nu)
{
    const VT lam = nu*E/((1.+nu)*(1.-2.*nu));
    const VT mu  = E/(2.*(1.+nu));
    unsigned int i, j, k, l;
    for (i = 0; i < 3; ++i){
        for (j = 0; j < 3; ++j){
            for (k = 0; k < 3; ++k){
                for (l = 0; l < 3; ++l){
                    pos[i][j][k][l] = lam * this->delta(i,j) * this->delta(k,l)
                            + mu * (this->delta(i,k) * this->delta(j,l) + this->delta(i,l) * this->delta(j,k));
                }
            }
        }
    }
    return;
}

template<typename VT>
inline xTensor4<VT>::xTensor4 (const VT E, const VT nu, const VT beta )
{
    // generate some kind of cubic or isotropic transverse material by nullifying
    // shear termes of the standard isotropic hook tensor
    //
    // nota : ugly implementation. Should be more somthing like setting
    // directely the right termes.
    // ok for now
    const VT lam = nu*E/(( 1.+nu )*( 1.-2.*nu ));
    const VT mu = E/( 2.*( 1.+nu ));
    const VT mu_null = beta*mu;
    unsigned int i, j, k, l;
    // generate standard isotropic hook law
    for (i = 0; i < 3; ++i)
    {
        for (j = 0; j < 3; ++j)
        {
            for (k = 0; k < 3; ++k)
            {
                for (l = 0; l < 3; ++l)
                {
                    pos[i][j][k][l] = lam * this->delta(i,j) * this->delta(k,l)
                            + mu * ( this->delta(i,k) * this->delta(j,l) + this->delta(i,l) * this->delta(j,k));
                }
            }
        }
    }
    // nullify (beta) shear terms
    //pos[2][1][2][1] = mu_null;
    //pos[1][2][1][2] = mu_null;
    //pos[2][1][1][2] = mu_null;
    //pos[1][2][2][1] = mu_null;
    pos[2][0][2][0] = mu_null;
    pos[0][2][0][2] = mu_null;
    pos[2][0][0][2] = mu_null;
    pos[0][2][2][0] = mu_null;
    pos[0][1][0][1] = mu_null;
    pos[1][0][1][0] = mu_null;
    pos[0][1][1][0] = mu_null;
    pos[1][0][0][1] = mu_null;

    return;
}

template<typename VT>
inline xTensor4<VT>::xTensor4 (const VT c[3][3][3][3])
{
    std::copy(reinterpret_cast<const VT *>(c), reinterpret_cast<const VT *>(c)+81, reinterpret_cast<VT *>(pos));
    return;
}

//Specialize only for double
template<>
inline xTensor4<double> xTensor4<double>::rotate(const xTensor2<double> &rot)
{
    // nota : ugly implementation. Should use the fact that only 21 termes
    // are independante.
    // A try to view tensor as a 81 termes vector and use a blas operation like Dpqrs=R*Cijkl
    // with R a 81x81 matrice filled by rot matrice was not succesfull : the blas was
    // giving better performance then this implemention but creation of R was to expansive.
    // Maybe by taking into acompte the minor and major symetri .....
    // ok for now
    xTensor4<double> res(xtool::xDataType<double>::zero());
    unsigned int i, j, k, l, p, q, r, s;
    double tmp,vi,vj,vk;
    for (p = 0; p < 3; ++p)
    {
        for (q = 0; q < 3; ++q)
        {
            for (r = 0; r < 3; ++r)
            {
                for (s = 0; s < 3; ++s)
                {
                    tmp = xtool::xDataType<double>::zero();
                    for (i = 0; i < 3; ++i)
                    {
                        vi = rot(p,i);
                        for (j = 0; j < 3; ++j)
                        {
                            vj = vi*rot(q,j);
                            for (k = 0; k < 3; ++k)
                            {
                                vk = vj*rot(r,k);
                                for (l = 0; l < 3; ++l)
                                {
                                    tmp += vk*rot(s,l)*pos[i][j][k][l];
                                }
                            }
                        }
                    }
                    res(p,q,r,s) = tmp;
                }
            }
        }
    }

    return res;
}

template<typename VT>
inline xTensor4<VT> & xTensor4<VT>::operator += (const xTensor4 &other)
{
    VT *l = data();
    const VT *r  = other.data();
    for(unsigned int i=0;i<81;++i)
    {
        *l++ +=  *r++;
    }
    return *this;
}

template<typename VT>
inline xTensor4<VT> & xTensor4<VT>::operator -= (const xTensor4 &other){
    VT *l = data();
    const VT *r  = other.data();
    for(unsigned int i=0;i<81;++i)
    {
        *l++ -=  *r++;
    }
    return *this;
}

template<typename VT>
inline VT xTensor4<VT>::normValue() const{
    const VT *pA = data();
    VT tmp = xtool::xDataType<VT>::zero();
    for (unsigned int i= 0; i < 81; ++i) {
        const VT ai = *(pA+i);
        tmp += ai*ai;
    }
    return sqrt(tmp);
}

template<typename VT>
inline VT * xTensor4<VT>::data () {return pos[0][0][0];}

template<typename VT>
inline const VT * xTensor4<VT>::data () const {return pos[0][0][0];}

template<typename VT>
inline const VT &xTensor4<VT>::operator () (int i, int j, int k, int l) const
{
    return pos[i][j][k][l];
}

template<typename VT>
inline VT & xTensor4<VT>::operator () (int i, int j, int k, int l)
{
    return pos[i][j][k][l];
}

template<typename VT>
inline xTensor4<VT> operator*(const xTensor4<VT> &A, const xTensor4<VT> &B){
    xTensor4<VT> tmp(xtool::xDataType<VT>::zero());
    VT *ptmp = &tmp(0,0,0,0);
    const VT *pA = &A(0,0,0,0);
    const VT *pB = &B(0,0,0,0);
    for (unsigned int i=0; i < 9; ++i){
        for (unsigned int k=0; k < 9; ++k){
            for (unsigned int j=0; j < 9; ++j){
                ptmp[9*i+j] += pA[9*i+k]*pB[9*k+j];
            }
        }
    }
    return tmp;
};

template<typename VT>
inline std::ostream& operator<<(std::ostream& s, const xTensor4<VT>& t){
    for (unsigned int i= 0; i< 3; ++i){
        for (unsigned int j= 0; j< 3; ++j){
            for (unsigned int k= 0; k< 3; ++k){
                for (unsigned int l= 0; l< 3; ++l){
                    s << t(i,j,k,l) << " ";
                }
            }
            s << std::endl;
        }
    }
    return s;
}




inline bool choleskyL(const xTensor4<double> &A, xTensor4<double> &L, const double prec){
    const double *pA = &A(0,0,0,0);
    const double norm = A.normValue();
    const double eps= prec*norm;
    double  d[] = {0.,0.,0.,0.,0.,0.,0.,0.,0.};

    L ={
        1.,0.,0.,0.,0.,0.,0.,0.,0.,
        0.,1.,0.,0.,0.,0.,0.,0.,0.,
        0.,0.,1.,0.,0.,0.,0.,0.,0.,
        0.,0.,0.,1.,0.,0.,0.,0.,0.,
        0.,0.,0.,0.,1.,0.,0.,0.,0.,
        0.,0.,0.,0.,0.,1.,0.,0.,0.,
        0.,0.,0.,0.,0.,0.,1.,0.,0.,
        0.,0.,0.,0.,0.,0.,0.,1.,0.,
        0.,0.,0.,0.,0.,0.,0.,0.,1.
    };

    double *pL = &L(0.,0.,0.,0.);
    for (unsigned int j=0; j < 9; ++j){
        d[j] = pA[9*j+ j];
        for (unsigned int k = 0; k < j; ++k ){
            d[j] -= pL[9*j+k]*pL[9*j+k]*d[k];
        }
        if (fabs(d[j]) >  eps)
            for(unsigned int i=j+1; i <9; ++i){
                pL[9*i+j] = pA[9*j +i];
                for (unsigned int k =0; k< j ; ++k){
                    pL[9*i+j] -= pL[9*i+k]*pL[9*j+k]*d[k];
                }
                pL[9*i+j] /= d[j] ;
            }
        else{
            d[j] = 0.;
            for(unsigned int i=j+1; i <9; ++i){
                pL[9*j+i] = 0.;
            }
        }
    }
    for (unsigned int j =0; j < 9; ++j){
        if (d[j] < 0.) return false;
        else {
            const double piv =sqrt(d[j]);
            for (unsigned int i = j; i<9; ++i){
                pL[9*i+j] *= piv;
            }
        }
    }
    return true;
}


inline bool choleskyU(const xTensor4<double> &A, xTensor4<double> &U, const double prec){
    const double *pA = &A(0,0,0,0);
    const double norm = A.normValue();
    const double eps= prec*norm;
    double  d[] = {0.,0.,0.,0.,0.,0.,0.,0.,0.};
    U ={
        1.,0.,0.,0.,0.,0.,0.,0.,0.,
        0.,1.,0.,0.,0.,0.,0.,0.,0.,
        0.,0.,1.,0.,0.,0.,0.,0.,0.,
        0.,0.,0.,1.,0.,0.,0.,0.,0.,
        0.,0.,0.,0.,1.,0.,0.,0.,0.,
        0.,0.,0.,0.,0.,1.,0.,0.,0.,
        0.,0.,0.,0.,0.,0.,1.,0.,0.,
        0.,0.,0.,0.,0.,0.,0.,1.,0.,
        0.,0.,0.,0.,0.,0.,0.,0.,1.
    };

    double *pU = &U(0.,0.,0.,0.);
    for (unsigned int i=0; i < 9; ++i){
        d[i] = pA[9*i+ i];
        for (unsigned int k = 0; k < i; ++k ){
            d[i] -= pU[9*k+i]*pU[9*k+i]*d[k];
        }
        if (fabs(d[i]) >  eps)
            for(unsigned int j=i+1; j <9; ++j){
                pU[9*i+j] = pA[9*i +j];
                for (unsigned int k =0; k< i ; ++k){
                    pU[9*i+j] -= pU[9*k+i]*pU[9*k+j]*d[k];
                }
                pU[9*i+j] /= d[i] ;
            }
        else{
            d[i] = 0.;
            for(unsigned int j=i+1; j <9; ++j){
                pU[9*i+j] = 0.;
            }
        }
    }
    for (unsigned int i =0; i < 9; ++i){
        if (d[i] < 0.) return false;
        else {
            const double piv =sqrt(d[i]);
            for (unsigned int j = i; j<9; ++j){
                pU[9*i+j] *= piv;
            }
        }
    }
    return true;
}

inline xTensor4AnisoPlaneStrain::xTensor4AnisoPlaneStrain(): nu(0.), E(0.){};

inline void xTensor4AnisoPlaneStrain::updateCoeffs(const double beta){
    const double lam = nu*E/((1.+nu)*(1.-2.*nu));
    const double mu  = E/(2.*(1.+nu));
    coefs(0,0) = lam + 2*mu;
    coefs(0,1) = lam;
    coefs(1,0) = lam;
    coefs(1,1) = lam + 2*mu;
    coefs(2,2) = beta * mu;
    return;
}

inline xTensor4AnisoPlaneStrain::xTensor4AnisoPlaneStrain(const xTensor2<double> & _coefs, double _nu):coefs(_coefs), nu(_nu) {};

inline xTensor4AnisoPlaneStrain::xTensor4AnisoPlaneStrain(const double _E, const double _nu): nu(_nu), E(_E){
    updateCoeffs();
    return;
}

inline xTensor4AnisoPlaneStrain::xTensor4AnisoPlaneStrain(const double _E, const double _nu, const double beta ): nu(_nu), E(_E){
    updateCoeffs(beta);
    return;
}

inline xTensor2<double>  xTensor4AnisoPlaneStrain::operator * (const xTensor2<double> &other) const{
    // the 3 unknown components of a "2D" (plane strain) xTensor2<double>
    xVector<double> vec,res;
    vec(0) = other(0,0);
    vec(1) = other(1,1);
    vec(2) = other(0,1)+other(1,0);
    res = coefs*vec;
    xTensor2<double> m(0.);
    m(0,0) = res(0);
    m(1,1) = res(1);
    m(0,1) = res(2);
    m(1,0) = res(2);
    m(2,2) = nu*(res(0)+res(1));// plane strain assumption !!!
    return m;
};

inline xTensor4AnisoPlaneStrain  xTensor4AnisoPlaneStrain::operator * (const double &other) const{
    xTensor4AnisoPlaneStrain m(coefs,nu);
    m.coefs*= other;
    return m;
};

inline xTensor4AnisoPlaneStrain  xTensor4AnisoPlaneStrain::rotate(const xTensor2<double>& rot) const{
    xTensor2<double> Q;
    Q(0,0) = rot(0,0)*rot(0,0);
    Q(0,1) = rot(1,0)*rot(1,0);
    Q(0,2) = rot(0,0)*rot(1,0);
    Q(1,0) = rot(0,1)*rot(0,1);
    Q(1,1) = rot(1,1)*rot(1,1);
    Q(1,2) = rot(0,1)*rot(1,1);
    Q(2,0) = 2.*rot(0,0)*rot(0,1);
    Q(2,1) = 2.*rot(1,0)*rot(1,1);
    Q(2,2) = rot(0,0)*rot(1,1)+rot(0,1)*rot(1,0);
    xTensor2<double> Pinv(Q);
    Pinv.transpose();
    return xTensor4AnisoPlaneStrain((Pinv*coefs)*Q, nu);
}

inline void xTensor4AnisoPlaneStress::updateCoeffs(){
    const double a1 = E/(1.-(nu*nu));
    coefs(0,0) = a1;
    coefs(0,1) = a1*nu;
    coefs(1,0) = a1*nu;
    coefs(1,1) = a1;
    coefs(2,2) = a1*0.5*(1.-nu);
    return;
}

inline xTensor4AnisoPlaneStress::xTensor4AnisoPlaneStress(): nu(0.), E(0.){};

inline xTensor4AnisoPlaneStress::xTensor4AnisoPlaneStress(const xTensor2<double> & _coefs, double _nu):coefs(_coefs), nu(_nu) {};

inline xTensor4AnisoPlaneStress::xTensor4AnisoPlaneStress(const double _E, const double _nu): nu(_nu), E(_E){
    updateCoeffs();
    return;
}

inline xTensor2<double>  xTensor4AnisoPlaneStress::operator * (const xTensor2<double> &other) const{
    // the 3 unknown components of a "2D" (plane stress) xTensor2<double>
    xVector<double> vec,res;
    vec(0) = other(0,0);
    vec(1) = other(1,1);
    vec(2) = other(0,1)+other(1,0);
    res = coefs*vec;
    xTensor2<double> m(0.);
    m(0,0) = res(0);
    m(1,1) = res(1);
    m(0,1) = res(2);
    m(1,0) = res(2);
    m(2,2) = 0.;// plane stress assumption !!!
    return m;
};

inline xTensor4AnisoPlaneStress  xTensor4AnisoPlaneStress::operator * (const double &other) const{
    xTensor4AnisoPlaneStress m(coefs,nu);
    m.coefs*= other;
    return m;
};

