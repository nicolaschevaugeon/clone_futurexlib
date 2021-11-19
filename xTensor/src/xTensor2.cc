/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/
///////////////////////////////////////////////////////////////

#include <algorithm>
#include <functional>
#include <cmath>
#include <cstdio>
#include <limits>
#include <cassert>
#include "xTensor2.h"
#include "xTensor4.h"
#include "xVector.h"


namespace xtensor {

//Implementation only valid for double, so define it here
//Why ? to avoid later circular dependencies...
template <>
xTensor2<double> xTensor2<double>::operator * (const xTensor4Isotropic &other) const
{
    xTensor2<double> m( *this, 2*other.mu );
    m.symmetrize();
    double lam =  other.lam;
    double t = trace()*lam;
    for(unsigned int i=0;i<3;i++){
        m.pos[i][i] +=  t;
    }
    return m;
}



//Implementation only valid for double, so define it here
//Why ? to avoid later circular dependencies...
template <>
xTensor2<double> xTensor2<double>::operator * (const xTensor4AnisoPlaneStrain &other) const
{
    // the 3 unknown components of a "2D" (plane strain) xTensor2
    xVector<double> vec,res;
    vec(0) = pos[0][0];
    vec(1) = pos[1][1];
    vec(2) = (pos[0][1]+pos[1][0]);
    xTensor2 coefs_tr(other.coefs);
    coefs_tr.transpose();
    res = coefs_tr*vec;
    xTensor2 m(0.0);
    m(0,0) = res(0);
    m(1,1) = res(1);
    m(0,1) = res(2);
    m(1,0) = res(2);
    m(2,2) = other.nu*(res(0)+res(1));// plane strain assumption !!!
    return m;
}



//Implementation only valid for double, so define it here
//Why ? to avoid later circular dependencies...
template <>
xTensor2<double> xTensor2<double>::operator * (const xTensor4AnisoPlaneStress &other) const
{
    // the 3 unknown components of a "2D" (plane stress) xTensor2
    xVector<double> vec,res;
    vec(0) = pos[0][0];
    vec(1) = pos[1][1];
    vec(2) = (pos[0][1]+pos[1][0]);
    xTensor2 coefs_tr(other.coefs);
    coefs_tr.transpose();
    res = coefs_tr*vec;
    xTensor2 m(0.0);
    m(0,0) = res(0);
    m(1,1) = res(1);
    m(0,1) = res(2);
    m(1,0) = res(2);
    m(2,2) = 0.;// plane stress assumption !!!
    return m;
}



long FindCubicRoots(const double coeff[4], double x[3])
{
    double a1 = coeff[2] / coeff[3];
    double a2 = coeff[1] / coeff[3];
    double a3 = coeff[0] / coeff[3];

    double Q = (a1 * a1 - 3 * a2) / 9.;
    double R = (2. * a1 * a1 * a1 - 9. * a1 * a2 + 27. * a3) / 54.;
    double Qcubed = Q * Q * Q;
    double d = Qcubed - R * R;

    //    printf ("d = %22.15e Q = %12.5E R = %12.5E Qcubed %12.5E\n",d,Q,R,Qcubed);

    /// three roots, 2 equal
    if(Qcubed == 0.0 || fabs ( Qcubed - R * R ) < 1.e-8 * (fabs ( Qcubed) + fabs( R * R)) )
    {
        double theta;
        if (Qcubed <= 0.0)theta = acos(1.0);
        else if (R / sqrt(Qcubed) > 1.0)theta = acos(1.0);
        else if (R / sqrt(Qcubed) < -1.0)theta = acos(-1.0);
        else theta = acos(R / sqrt(Qcubed));
        double sqrtQ = sqrt(Q);
        //	printf("sqrtQ = %12.5E teta=%12.5E a1=%12.5E\n",sqrt(Q),theta,a1);
        x[0] = -2 * sqrtQ * cos( theta           / 3) - a1 / 3;
        x[1] = -2 * sqrtQ * cos((theta + 2 * M_PI) / 3) - a1 / 3;
        x[2] = -2 * sqrtQ * cos((theta + 4 * M_PI) / 3) - a1 / 3;
        return (3);
    }

    /* Three real roots */
    if (d >= 0.0) {
        double theta = acos(R / sqrt(Qcubed));
        double sqrtQ = sqrt(Q);
        x[0] = -2 * sqrtQ * cos( theta           / 3) - a1 / 3;
        x[1] = -2 * sqrtQ * cos((theta + 2 * M_PI) / 3) - a1 / 3;
        x[2] = -2 * sqrtQ * cos((theta + 4 * M_PI) / 3) - a1 / 3;
        return (3);
    }

    /* One real root */
    else {
        printf("IMPOSSIBLE !!!\n");

        double e = pow(sqrt(-d) + fabs(R), 1. / 3.);
        if (R > 0)
            e = -e;
        x[0] = (e + Q / e) - a1 / 3.;
        return (1);
    }
}

#define MAXN 32
#define R(i,j)	result[n*(i)+(j)]

long NullSpace(const double *a, double *result, double eps, long n)
{
    int r[MAXN], c[MAXN];
    long i, j, k;
    int jj=0;
    int kk=0;
    int t;
    double max, temp;
    int ec;

    for (i = 0; i < n; i++)
        r[i] = c[i] = -1;			/* Reset row and column pivot indices */

    // copy the input matrix if not in place
    if (result != a)
        for (i = 0; i < n*n; i++)
            result[i] = a[i];
    // rest of algorithm is in place wrt result[]

    for (i = 0; i < n; i++) {
        /* Find the biggest element in the remaining submatrix
       * for the next full pivot.
       */
        max = 0.0;
        for (k = 0; k < n; k++) {
            if (r[k] < 0) {
                for (j = 0; j < n; j++) {
                    if ((c[j] < 0) && ((temp = fabs(R(k, j))) > max)) {
                        kk = k;
                        jj = j;
                        max = temp;
                    }
                }
            }
        }
        if (max < eps)
            break;		/* Consider this and all subsequent pivots to be zero */

        c[jj] = kk;					/* The row */
        r[kk] = jj;					/*	      and column of the next pivot */

        temp = 1.0 / R(kk, jj);
        R(kk, jj) = 1.0;
        for (j = 0; j < n; j++)		/* Should this be for j != jj ? */
            R(kk, j) *= temp;		/* Row equilibration */

        for (k = 0; k < n; k++) {	/* Row elimination */
            if (k == kk)
                continue;			/* Don't do a thing to the pivot row */
            temp = R(k, jj);
            R(k, jj) = 0.0;
            for (j = 0; j < n; j++) {
                R(k, j) -= temp * R(kk, j);	/* Subtract row kk from row k */
                if (fabs(R(k, j)) < eps)
                    R(k, j) = 0.0;	/* Flush to zero if too small */
            }
        }
    }

    /* Sort into a truncated triangular matrix */
    for (j = 0; j < n; j++) {		/* For all columns... */
        while ((c[j] >= 0) && (j != c[j])) {
            for (k = 0; k < n; k++) {
                if (r[k] < 0) {
                    /* Aha! a null column vector */
                    temp = R(k, j);	/* Get it on top */
                    R(k, j) = R(k, c[j]);
                    R(k, c[j]) = temp;
                }
            }
            t = c[j];				/* Twiddle until pivots are on the diagonal */
            c[j] = c[t];
            c[t] = t;
        }
    }

    /* Copy the null space vectors into the top of the A matrix */
    ec = 0;
    for (k = 0; k < n; k++) {
        if (r[k] < 0) {
            R(k, k) = 1.0;			/* Set the pivot equal to 1 */
            if (ec != k) {
                for (j = 0; j < n; j++) {
                    R(ec, j) = R(k, j);
                }
            }
            ec++;
        }
    }
    /* The first  ec  rows of the matrix  a  are the vectors which are
     * orthogonal to the columns of the matrix  a.
     */
    return (ec);
}




//Specialize only for double
template <>
xVector<double> xTensor2<double>::operator / (const xVector<double> & b) const
{
    xVector<double> res;

    double detm = det();

    if (detm == 0.0){
        throw -1;
    }

    double ud = 1. / (detm);

    res[0] = b[0] * (pos[1][1] * pos[2][2] - pos[1][2] * pos[2][1]) -
            pos[0][1] * (b[1] * pos[2][2] - pos[1][2] * b[2]) +
            pos[0][2] * (b[1] * pos[2][1] - pos[1][1] * b[2]);

    res[1] = pos[0][0] * (b[1] * pos[2][2] - pos[1][2] * b[2]) -
            b[0] * (pos[1][0] * pos[2][2] - pos[1][2] * pos[2][0]) +
            pos[0][2] * (pos[1][0] * b[2] - b[1] * pos[2][0]);

    res[2] = pos[0][0] * (pos[1][1] * b[2] - b[1] * pos[2][1]) -
            pos[0][1] * (pos[1][0] * b[2] - b[1] * pos[2][0]) +
            b[0] * (pos[1][0] * pos[2][1] - pos[1][1] * pos[2][0]);

    for (unsigned int i = 0; i < 3; i++)
        res[i] *= ud;

    return res;
}

//Specialize only for double
template <>
xTensor2<double> xTensor2<double>::invert() const
{
    xTensor2<double> inv;

    double detm = det();

    if (detm == 0.0){
        throw -1;
    }

    // Modifie le 03/10/2003 : invert donne l'inverse de la matrice
    //                         et non plus la transposee de l'inverse
    //                         erreur dans AOMD ?

    double ud = 1. / (detm);
    inv.pos[0][0] = ud * (pos[1][1] * pos[2][2] - pos[1][2] * pos[2][1]);
    inv.pos[1][0] = -ud * (pos[1][0] * pos[2][2] - pos[1][2] * pos[2][0]);
    inv.pos[2][0] = ud * (pos[1][0] * pos[2][1] - pos[1][1] * pos[2][0]);
    inv.pos[0][1] = -ud * (pos[0][1] * pos[2][2] - pos[0][2] * pos[2][1]);
    inv.pos[1][1] = ud * (pos[0][0] * pos[2][2] - pos[0][2] * pos[2][0]);
    inv.pos[2][1] = -ud * (pos[0][0] * pos[2][1] - pos[0][1] * pos[2][0]);
    inv.pos[0][2] = ud * (pos[0][1] * pos[1][2] - pos[0][2] * pos[1][1]);
    inv.pos[1][2] = -ud * (pos[0][0] * pos[1][2] - pos[0][2] * pos[1][0]);
    inv.pos[2][2] = ud * (pos[0][0] * pos[1][1] - pos[0][1] * pos[1][0]);
    return inv;
}


struct greater_abs
{
    bool operator () (const double &a, const double &b)
    {
        return fabs(a) > fabs(b);
    }
};

//Specialize only for double
template <>
long xTensor2<double>::eigen (xVector<double> e[3], double v[3]) const
{
    /// characteristic polynomial of T : find v root of
    /// v^3 - I1 v^2 + I2 T + I3 = 0
    /// I1 : first invariant , trace(T)
    /// I2 : second invariant , 1/2 (I1^2 -trace(T^2))
    /// I3 : third invariant , det T
    double I[4];
    I[3] = 1.0;
    I[2] = - trace();
    I[1] = 0.5 * (I[2]*I[2] - trace2());
    I[0] = - det();

    //    printf (" %lf x^3 +  %lf x^2 + %lf x + %lf = 0\n",
    //    	  I[3],I[2],I[1],I[0]);

    long nbEigen = FindCubicRoots (I,v);

    std::sort(v,v+3, greater_abs() );

    //    printf ("nbEigen = %d %12.5E %12.5E %12.5E\n",nbEigen,v[0],v[1],v[2]);

    double result[12];
    unsigned int nb_vec=0;
    while(1)
    {
        double a[9] = {pos[0][0]-v[nb_vec],pos[0][1],pos[0][2],
                       pos[1][0],pos[1][1]-v[nb_vec],pos[1][2],
                       pos[2][0],pos[2][1],pos[2][2]-v[nb_vec]};

        double eps = 1.e-3;
        unsigned int nb = 0;
        while (1)
        {
            nb = NullSpace (a,result,eps,3);
            if (nb != 0)break;
            eps *= 2.0;
        }
        unsigned int kk=0;
        for (unsigned int i=nb_vec;i<nb+nb_vec;i++)
        {
            e[i] = xVector<double> (result[0+kk*3],result[1+kk*3],result[2+kk*3]);
            e[i].norm();
            kk++;
            if (i == 2)return nbEigen;
        }
        nb_vec += nb;
        if (nb_vec == 3)return nbEigen;
        if (nb > 3)throw -1;
    }
    //    printf (" %lf x^3 +  %lf x^2 + %lf x + %22.15E = 0\n",
    //	  I[3],I[2],I[1],I[0]);
    throw -1;
}

//Specialize only for double
template <>
long xTensor2<double>::eigen2d (xVector<double> e[3], double v[3]) const
{

    double a = 1.0;
    double b = - pos[0][0] - pos[1][1];
    double c = (pos[0][0] * pos[1][1] - pos[1][0] * pos[0][1]);

    e[2] = xVector<double>(0,0,1);

    double delta = b*b - 4 * a * c;

    if (delta < 0){
        const double tolerance = std::numeric_limits<double>::epsilon()*100.;
        if (delta>-tolerance)
            delta=0.;
        else
            return 0;
    }


    v[0] = (-b+sqrt(delta))/(2.*a);
    v[1] = (-b-sqrt(delta))/(2.*a);
    v[2] = 1.0;

    long nbEigen = 2;

    if (fabs(v[1]) > fabs(v[0]))
    {
        double temp = v[0];
        v[0] = v[1];
        v[1] = temp;
    }

    //    printf ("nbEigen = %d %12.5E %12.5E %12.5E\n",nbEigen,v[0],v[1],v[2]);

    double result[4];
    unsigned int nb_vec=0;
    while(1)
    {
        double a[4] = {pos[0][0]-v[nb_vec],pos[0][1],
                       pos[1][0],pos[1][1]-v[nb_vec]};
        double eps = 1.e-8;
        unsigned int nb = 0;
        while (1)
        {
            nb = NullSpace (a,result,eps,2);
            if (nb != 0)break;
            eps *= 2.0;
            //	    printf ("esp = %12.5E\n",eps);
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

//Specialize only for double
template <>
void xTensor2<double>::getAnalyticalEigenvaluesAndDerivative( double lambda[3], double der[][3][3],  double dersec[][3][3][3][3])
{


    // index
    unsigned int i,j,k,l,m,n,ki;

    // numerical constante
    const double zero = 0.;
    const double one = 1.;
    const double two = 2.;
    const double three = 3.;
    const double six = 6.;
    const double oneontwo = 1./2.;
    const double oneonthree = 1./3.;
    // const double oneonfor = 1./4.;
    //const double oneonsix = 1./6.;
    //const double oneonnine = 1./9.;
    //const double twoonnine = 2./9.;
    //const double oneontwentyseven = 1./27.;
    const double twoontwentyseven = 2./27.;
    const double twoonthree = 2./3.;
    //const double threeontwo = 3./2.;
    const double fouronthree = 4./3.;
    const double epsilonn = std::numeric_limits < double >::epsilon();

    // check as implementation is limited : to obtaine dersec der have to be done !
    if (dersec && !der)
    {
        std::cout<<"FILE "<<__FILE__<<" LINE "<<__LINE__<<std::endl;
        std::cout<<"Sorry but you ask for second derivative without derivative\n this is not treated yet\n for now stop\n";
        throw -1;
    }

    // calculation of normalisation factor
    const double nf = std::max(fabs(pos[0][0]),std::max(fabs(pos[0][1]),std::max(fabs(pos[0][2]),std::max(fabs(pos[1][1]),std::max(fabs(pos[1][2]),fabs(pos[2][2]))))));

    // if nf is zero we shouldn't have call this function !
    assert(nf != zero);

    // normalized strain
    xTensor2 strainn(*this,1/nf);

    // normalized invariant calculation
    const double i1 = strainn.trace();
    const double i2 = oneontwo*( i1*i1 - strainn.trace2());
    const double i3 = strainn.det();


    // equation of order 3 witch gives eigenvalues is :
    //
    //    lambda^3 - invariant1.lambda^2 + invariant2.lambda -invariant3 = 0
    //
    // as we normalize tensor we solve in fact :
    //
    //    lambdan^3 - i1.lambdan^2 + i2.lambdan -i3 = 0
    //
    //    with lambdan=lambda/nf
    //
    // it is rewrited with bi as :
    //
    // lambdan^3 + b2.lambdan^2 + b1.lambdan + b0 = 0  (1)
    //
    // with :
    //    b2 = -i1
    //    b1 =  i2
    //    b0 = -i3
    //


    // intermediate value
    const double i1i1 = i1*i1;
    const double b1 = i2;
    const double b2 = -i1;
    const double de = b2*oneonthree;

    //  we set x as lambdan=x-b2/3
    //
    //  this lead (1) to become
    //
    //  x^3 + p.x = q    (2)
    //
    // with :
    //   p = b1 -b2^2 / 3
    //   q = (9.b1.b2 - 2.b2^3 - 27.b0 )/27
    //
    const double p = b1-i1i1*oneonthree;
    const double q = i3+de*b1+i1i1*i1*twoontwentyseven;

    //
    //   setting :
    //        Q = p/3
    //        R = q/2
    //
    //   we can obtaine the polynomial discriminant D :
    //     D = Q^3 + R^2
    //
    //     depending on signe of D we have for (2) :
    //        D > 0 => one real root and two complex conjugates
    //        D = 0 => all roots are real with at least two equale
    //        D < 0 => all roots are real and not equal
    //
    //   nota : somme particular case are treated befor (p and/or q null)
    //   nota : D is effectively calculated below for the general case
    //
    const double Q = p*oneonthree;
    const double R = q*oneontwo;
    const double Q3 = Q*Q*Q;
    const double R2 = R*R;

    // set test on p,q null
    // nota : here we use Q3 for p test as general case divide a quantity by Q^3 so
    // this is a more strict test as it validate p numerical zero in a larger zone
    // and Q^3 numerical zero
    // same for q and R we use R^2
    const bool p_is_null = ( fabs(Q3) < epsilonn );
    const bool q_is_null = ( R2 < epsilonn );

    // in eigen values calculation we gone keep if D is null or not
    // for derivatives calculation
    bool D_is_null = false;

    //===========================================================================================
    //     eigen values calculation
    //===========================================================================================
    //
    // --------------------------------------------
    // particular case :
    //  p=0 q#0 => D > 0
    //          => (2) become x^3=q
    //          => on real root  lambdan=q^(1/3)-b2/3 an 2  complex conjugates
    if ( p_is_null && !q_is_null)
    {
        std::cout<<"FILE "<<__FILE__<<" LINE "<<__LINE__<<std::endl;
        std::cout<<"It appends !!! invariants lead to complex eigenvalues !!! context p=0 and q#0\n";
        std::cout<<this<<std::endl<<strainn;
        std::cout<<"i1 ="<<std::scientific<<i1<<std::endl;
        std::cout<<"i2 ="<<std::scientific<<i2<<std::endl;
        std::cout<<"i3 ="<<std::scientific<<i3<<std::endl;
        throw -1;
    }
    // --------------------------------------------
    // particular case :
    //  p=0 q=0 => D = 0
    //          => (2) become x^3=0
    //          => three identical real root  lambdan=-b2/3
    else if ( p_is_null && q_is_null)
    {
        // eigenvalue
        const double r = -de;
        lambda[0] = lambda[1] = lambda[2] = r;

        D_is_null = true;
    }
    // --------------------------------------------
    // particular case :
    //  p#0 q=0 => D = Q^3
    //          => (2) become x(x^2+p)=0
    //   if p<0 ,then D<0, three differents real root :
    //                       lambdan=-b2/3
    //                       lambdan=(-p)^(1/2)-b2/3
    //                       lambdan=-(-p)^(1/2)-b2/3
    //   if p>0, then D>0,one real root (-b2/3) and two complex conjugates : imposible in this mecanical context
    //
    else if ( !p_is_null && q_is_null)
    {
        // Stopping if complexe
        //  nota : here we test against  zero as p_is_null allready eliminate close to zero case
        if (p > zero)
        {
            std::cout<<"FILE "<<__FILE__<<" LINE "<<__LINE__<<std::endl;
            std::cout<<"It appends !!! invariants lead to  complex eigenvalues !!! context q=0 and p>0\n";
            std::cout<<this<<std::endl<<strainn;
            std::cout<<"i1 ="<<std::scientific<<i1<<std::endl;
            std::cout<<"i2 ="<<std::scientific<<i2<<std::endl;
            std::cout<<"i3 ="<<std::scientific<<i3<<std::endl;
            throw -1;
        }
        // eigenvalue
        const double r = sqrt(-p);
        lambda[0] = -de;
        lambda[1] = r-de;
        lambda[2] = -r-de;

    }
    // --------------------------------------------
    // general case with p#0 and q#0 :
    //
    //  D may still be every thing, evean if particular test have been eliminated,
    //  so we have to test it
    //
    else
    {
        const double D = Q3+R2;

        //  +++++++++++++++++++++++
        //  D > 0
        //  We still have to eliminate de positive case
        if ( D > epsilonn)
        {
            std::cout<<"FILE "<<__FILE__<<" LINE "<<__LINE__<<std::endl;
            std::cout<<"It appends !!! invariants lead to  complex eigenvalues !!! context q#0 and p#0\n";
            std::cout<<this<<std::endl<<strainn;
            std::cout<<"i1 ="<<std::scientific<<i1<<std::endl;
            std::cout<<"i2 ="<<std::scientific<<i2<<std::endl;
            std::cout<<"i3 ="<<std::scientific<<i3<<std::endl;
            throw -1;
        }
        //
        //  +++++++++++++++++++++++
        //  D = 0
        //
        //   (2) have 3 reals solutions and 2 are equale :
        //
        //       x1 = x2 = -sign(q)(-p/3)^(1/2) = -R^(1/3)
        //       x3= -3q/p = 2*sign(q)*(-p/3)^(1/2) = 2.R^(1/3) = -2*x1
        //
        //       where sign(q)=1 if q>0 and -1 if q<0
        //
        //       witch gives ei=xi-b2/3
        //
        // nota :
        //   This can be viewed as a particular case of equation below for D<0 with teta = 0 or pi
        //   but to avoid numerical probleme with "acos" defined in [-1,1]  we treat this case appart
        //
        // nota :
        //   We  gone a use R expretion as they comme from the general expretion or roots of a cubic
        //   polynome :
        //       x1 = -1/2.(S+T)+1/2.i.sqrt(3).(S-T)
        //       x2 = -1/2.(S+T)-1/2.i.sqrt(3).(S-T)
        //       x3 = S+T
        //
        //     with
        //        i^2=-1
        //        S=(R+sqrt(D))^(1/3)
        //        T=(R-sqrt(D))^(1/3)
        //        here D=0 => S=T=R^(1/3) and we get above expretion X1=X2=-R^(1/3), X3=2.R^(1/3)
        //
        //    using R^(1/3) impose anyway to check signe of R as pow(R,1/3) is nan if R<0
        //
        //
        else if ( D > -epsilonn )
        {

            // intermediate expretion
            const double CR = ( q < zero ) ? -pow(-R,oneonthree) : pow(R,oneonthree);
            const double r = -CR-de;

            // eigenvalue
            lambda[0] = two*CR-de;
            lambda[1] = lambda[2] = r;

            D_is_null = true;
        }
        //  +++++++++++++++++++++++
        //  D < 0
        //
        //    using R and Q the general solution of (2) is then :
        //
        //     x1 = 2*(-Q)^(1/2)*cos (teta/3)
        //     x2 = 2*(-Q)^(1/2)*cos ((teta+2*pi)/3)
        //     x3 = 2*(-Q)^(1/2)*cos ((teta+4*pi)/3)
        //
        //     with teta=acos(R.(-Q)^(-3/2))
        //
        //     witch gives lambdan=xi-b2/3
        //
        //     D<0  imply
        //          Q^3<0 => Q<0 => -Q>0 => sqrt(-Q) ok
        //          Q^3<-R^2 => R^2/(-Q)^3<1  => -1<R/(-Q)^(3/2)<1 => acos and acos' ok
        else
        {


            // intermediate expretion
            const double SMQ = sqrt(-Q);
            const double SMQ2 = two*SMQ;
            const double K = R/sqrt(-Q3);
            // teta calulation
            // epression R/-Q^(3/2) is from test above in ]-1,1[ => no check, only a assert
            assert (( K < one ) && ( K > -one ));
            // here it's alpha=teta/3 wich is used
            const double alpha = oneonthree*acos(K);
            const double beta = M_PI*( twoonthree )+alpha;
            const double gama = M_PI*( fouronthree )+alpha;
            const double calpha = cos(alpha);
            const double cbeta = cos(beta);
            const double cgama = cos(gama);

            // eigenvalue
            lambda[0] = calpha*SMQ2-de;
            lambda[1] = cbeta*SMQ2-de;
            lambda[2] = cgama*SMQ2-de;


        } // end D<0

    } // end general case p#0 q#0

    if (der)
    {
        //===========================================================================================
        //     first derivatives calculation
        //===========================================================================================
        //
        // the general form of the first derivative of eigen value by invariants is obtained from (1) by it's derivation
        //
        // d(lambdan)/d(ik).(3.lambdan^2+2.b2.lambdan+b1) = -d(b2)/d(ik).lambdan^2 - d(b1)/d(ik).lambdan - d(b0)/d(ik) (3)
        //
        // with  k=1,2,3 (3 invariants)
        //
        // this gives d(lambdan)/d(ik) only if (3.lambdan^2+2.b2.lambdan+b1)#0
        //
        // the question is then : is it possible to have eign values solution of (1) wich satisfy 3.lambdan^2+2.b2.lambdan+b1=0  ?
        //
        //      setting equation (L as unknow)
        //      3.L^2+2.b2.L+b1=0 (3b)
        //
        // solution of (3b) is given with determinant calculation of (3b) : delta=-4p/3
        //
        //  * if p>0 => delta<0 and (3b) have no solution which is good but isn't usefull as solutions of (1) are complexe
        //
        //  * if p<0 => delta>0 and (3b) have 2 solutions :
        //
        //           L1= sqrt(-Q) - b2/3
        //           L2= -sqrt(-Q) - b2/3
        //
        //         p<0 imply D<=0
        //
        //         + D<0 => lambdanj=2.cos(aj).sqrt(-Q)-b2/3 j=1,2,3 (eigen values)
        //
        //                  with aj= alpha if j=1
        //                          beta  if j=2
        //                          gama  if j=3
        //
        //                 lambdanj=L1 => cos(aj)=1/2
        //                 lambdanj=L2 => cos(aj)=-1/2
        //
        //                 this is possible only if teta=0 or pi which correspond to D=0
        //                 => derivative is then obtain from (3) by dividing by (3.lambdan^2+2.b2.lambdan+b1)#0
        //
        //         + D=0 => one of the root is lambdan=sign(q).sqrt(-Q) - b2/3 wich corespond to L1 or L2 depending on signe of q
        //           here all the derivatives are imposible to obtaine with (3) ?!?!?
        //
        //  * if p=0 => delta=0 and (3b) have 1 solutions :
        //
        //         L= - b2/3
        //
        //         when p=0 if q#0 (1) have complexe solution wich is not treated here
        //         when p=0 and q=0, D=0 and (1) have a unique real triple solution wich is -b2/3 !!!!!
        //           here de derivatives is imposible to obtaine with (3) ?!?!?
        //
        //
        //


        // intermediate derivative table
        double der_l_i[3][3];
        double der_i_e[3][3][3];
        double *u,*v;//*w;
        //

        // for now no solution
        if (D_is_null)
        {
            std::cout<<"FILE "<<__FILE__<<" LINE "<<__LINE__<<std::endl;
            std::cout<<"It appends !!! we are asking for derivatives and D=0 !\n";
            std::cout<<" pertubate to get derivatives !!! \n";
            std::cout<<this<<std::endl<<strainn;
            std::cout<<"i1 ="<<std::scientific<<i1<<std::endl;
            std::cout<<"i2 ="<<i2<<std::endl;
            std::cout<<"i3 ="<<i3<<std::endl;
            throw -1;
        }

        // initialize intermediate derivative table
        //
        u = &der_i_e[0][0][0];
        v = &der[0][0][0];
        for (i = 0; i < 27; ++i)
        {
            u[i] = zero;
            v[i] = zero;
        }

        // compute derivatives against invariants
        //
        // for L1 :
        const double L1 = lambda[0];
        const double r1 = one/(( three*L1+two*b2 )*L1+b1 );

        der_l_i[0][1] = -L1*r1;
        der_l_i[0][0] = -der_l_i[0][1]*L1;
        der_l_i[0][2] = r1;

        // for L2 :
        const double L2 = lambda[1];
        const double r2 = one/(( three*L2+two*b2 )*L2+b1 );

        der_l_i[1][1] = -L2*r2;
        der_l_i[1][0] = -der_l_i[1][1]*L2;
        der_l_i[1][2] = r2;

        // for L3 :
        const double L3 = lambda[2];
        const double r3 = one/(( three*L3+two*b2 )*L3+b1 );

        der_l_i[2][1] = -L3*r3;
        der_l_i[2][0] = -der_l_i[2][1]*L3;
        der_l_i[2][2] = r3;


        // normalisezd strain componant
        const double e11 = strainn(0,0);
        const double e12 = strainn(0,1);
        const double e13 = strainn(0,2);
        const double e22 = strainn(1,1);
        const double e23 = strainn(1,2);
        const double e33 = strainn(2,2);


        der_i_e[0][0][0] = 1.0;
        der_i_e[0][1][1] = 1.0;
        der_i_e[0][2][2] = 1.0;

        der_i_e[1][0][0] = e22+e33;
        der_i_e[1][0][1] = -e12;
        der_i_e[1][0][2] = -e13;
        der_i_e[1][1][0] = -e12;
        der_i_e[1][1][1] = e11+e33;
        der_i_e[1][1][2] = -e23;
        der_i_e[1][2][0] = -e13;
        der_i_e[1][2][1] = -e23;
        der_i_e[1][2][2] = e11+e22;

        der_i_e[2][0][0] = e22*e33-e23*e23;
        der_i_e[2][0][1] = e13*e23-e12*e33;
        der_i_e[2][0][2] = e12*e23-e13*e22;
        der_i_e[2][1][0] = e13*e23-e12*e33;
        der_i_e[2][1][1] = e11*e33-e13*e13;
        der_i_e[2][1][2] = e12*e13-e11*e23;
        der_i_e[2][2][0] = e12*e23-e13*e22;
        der_i_e[2][2][1] = e12*e13-e11*e23;
        der_i_e[2][2][2] = e11*e22-e12*e12;


        // nota : normalisation is  automatiquely  taken into acount for first derivatives :
        //
        //   d(lambdan)/d(epsilon_normalized_ij) = d(lambda)/d(epsilon_ij)
        //

        // computing derivatives product
        double r;
        for (n = 0; n < 3; ++n)
        {
            u = &der[n][0][0];
            for (k = 0; k < 3; ++k)
            {
                r = der_l_i[n][k];
                v = &der_i_e[k][0][0];

                // daxpy ?
                for ( i = 0; i < 9; ++i)
                {
                    u[i] += r*v[i];

                }
            }
        }





        if (dersec)
        {
            //===========================================================================================
            //     second derivative computation
            //===========================================================================================
            //
            // by derivation of (3) we get the general form of the second derivative
            //
            //
            // d2(lambdan)/d(ik)d(il).(3.lambdan^2+2.b2.lambdan+b1) = -d(lambdan)/d(il) . ( 2.lambdan.d(b2)/d(ik) + d(b1)/d(ik) )
            //                                                        -d(lambdan)/d(ik) . ( 2.lambdan.d(b2)/d(il) + d(b1)/d(il) )
            //                                                        -d(lambdan)/d(ik) .  d(lambdan)/d(il) .( 6.lambdan + 2.b2 )  (4)
            //
            //
            // with  k=1,2,3 (3 invariants)
            //       l=1,2,3 (3 invariants)
            //
            // again this is available only if (3b) is not null
            //
            //
            // intermediate derivative table
            double der2_l_i[3][3][3];
            double der2_i_e[3][3][3][3][3];
            //
            // compute second derivatives against invariants
            //
            // for L1 :
            double expr0 = -two*L1;
            double expr3 = six*L1+two*b2;

            double expr30 = der_l_i[0][0]*expr3;
            der2_l_i[0][0][0] = -( der_l_i[0][0]*( expr0+expr0+expr30 ) )*r1;
            der2_l_i[0][0][1] = -( der_l_i[0][1]*expr0+der_l_i[0][0]+der_l_i[0][1]*expr30 )*r1;
            der2_l_i[0][0][2] = -( der_l_i[0][2]*expr0+der_l_i[0][2]*expr30 )*r1;

            double expr31 = der_l_i[0][1]*expr3;
            der2_l_i[0][1][0] = der2_l_i[0][0][1];
            der2_l_i[0][1][1] = -( der_l_i[0][1]*( two+expr31 ))*r1;
            der2_l_i[0][1][2] = -( der_l_i[0][2]+der_l_i[0][2]*expr31 )*r1;

            der2_l_i[0][2][0] = der2_l_i[0][0][2];
            der2_l_i[0][2][1] = der2_l_i[0][1][2];
            der2_l_i[0][2][2] = -( der_l_i[0][2]*der_l_i[0][2]*expr3 )*r1;

            // for L2 :
            expr0 = -two*L2;
            expr3 = six*L2+two*b2;

            expr30 = der_l_i[1][0]*expr3;
            der2_l_i[1][0][0] = -( der_l_i[1][0]*( expr0+expr0+expr30 ))*r2;
            der2_l_i[1][0][1] = -( der_l_i[1][1]*expr0+der_l_i[1][0]+der_l_i[1][1]*expr30 )*r2;
            der2_l_i[1][0][2] = -( der_l_i[1][2]*expr0+der_l_i[1][2]*expr30 )*r2;

            expr31 = der_l_i[1][1]*expr3;
            der2_l_i[1][1][0] = der2_l_i[1][0][1];
            der2_l_i[1][1][1] = -( der_l_i[1][1]*( two+expr31 ))*r2;
            der2_l_i[1][1][2] = -( der_l_i[1][2]+der_l_i[1][2]*expr31 )*r2;

            der2_l_i[1][2][0] = der2_l_i[1][0][2];
            der2_l_i[1][2][1] = der2_l_i[1][1][2];
            der2_l_i[1][2][2] = -( der_l_i[1][2]*der_l_i[1][2]*expr3 )*r2;

            // for L3 :
            expr0 = -two*L3;
            expr3 = six*L3+two*b2;

            expr30 = der_l_i[2][0]*expr3;
            der2_l_i[2][0][0] = -( der_l_i[2][0]*( expr0+expr0+expr30 ))*r3;
            der2_l_i[2][0][1] = -( der_l_i[2][1]*expr0+der_l_i[2][0]+der_l_i[2][1]*expr30 )*r3;
            der2_l_i[2][0][2] = -( der_l_i[2][2]*expr0+der_l_i[2][2]*expr30 )*r3;

            expr31 = der_l_i[2][1]*expr3;
            der2_l_i[2][1][0] = der2_l_i[2][0][1];
            der2_l_i[2][1][1] = -( der_l_i[2][1]*( two+expr31 ))*r3;
            der2_l_i[2][1][2] = -( der_l_i[2][2]+der_l_i[2][2]*expr31 )*r3;

            der2_l_i[2][2][0] = der2_l_i[2][0][2];
            der2_l_i[2][2][1] = der2_l_i[2][1][2];
            der2_l_i[2][2][2] = -( der_l_i[2][2]*der_l_i[2][2]*expr3 )*r3;

            // initialize intermediate derivative tables
            n = 243;
            u=&der2_i_e[0][0][0][0][0];
            for (i=0;i<n;++i) u[i]=zero;

            der2_i_e[1][0][0][1][1] = 1.0;
            der2_i_e[1][0][0][2][2] = 1.0;
            der2_i_e[1][1][0][0][1] = -1.0;
            der2_i_e[1][2][0][0][2] = -1.0;
            der2_i_e[1][0][1][1][0] = -1.0;
            der2_i_e[1][1][1][0][0] = 1.0;
            der2_i_e[1][1][1][2][2] = 1.0;
            der2_i_e[1][2][1][1][2] = -1.0;
            der2_i_e[1][0][2][2][0] = -1.0;
            der2_i_e[1][1][2][2][1] = -1.0;
            der2_i_e[1][2][2][0][0] = 1.0;
            der2_i_e[1][2][2][1][1] = 1.0;

            der2_i_e[2][0][0][1][1] = e33;
            der2_i_e[2][0][0][2][1] = -e23;
            der2_i_e[2][0][0][1][2] = -e23;
            der2_i_e[2][0][0][2][2] = e22;
            der2_i_e[2][1][0][0][1] = -e33;
            der2_i_e[2][1][0][2][1] = e13;
            der2_i_e[2][1][0][0][2] = e23;
            der2_i_e[2][1][0][2][2] = -e12;
            der2_i_e[2][2][0][0][1] = e23;
            der2_i_e[2][2][0][1][1] = -e13;
            der2_i_e[2][2][0][0][2] = -e22;
            der2_i_e[2][2][0][1][2] = e12;
            der2_i_e[2][0][1][1][0] = -e33;
            der2_i_e[2][0][1][2][0] = e23;
            der2_i_e[2][0][1][1][2] = e13;
            der2_i_e[2][0][1][2][2] = -e12;
            der2_i_e[2][1][1][0][0] = e33;
            der2_i_e[2][1][1][2][0] = -e13;
            der2_i_e[2][1][1][0][2] = -e13;
            der2_i_e[2][1][1][2][2] = e11;
            der2_i_e[2][2][1][0][0] = -e23;
            der2_i_e[2][2][1][1][0] = e13;
            der2_i_e[2][2][1][0][2] = e12;
            der2_i_e[2][2][1][1][2] = -e11;
            der2_i_e[2][0][2][1][0] = e23;
            der2_i_e[2][0][2][2][0] = -e22;
            der2_i_e[2][0][2][1][1] = -e13;
            der2_i_e[2][0][2][2][1] = e12;
            der2_i_e[2][1][2][0][0] = -e23;
            der2_i_e[2][1][2][2][0] = e12;
            der2_i_e[2][1][2][0][1] = e13;
            der2_i_e[2][1][2][2][1] = -e11;
            der2_i_e[2][2][2][0][0] = e22;
            der2_i_e[2][2][2][1][0] = -e12;
            der2_i_e[2][2][2][0][1] = -e12;
            der2_i_e[2][2][2][1][1] = e11;

            // nota : normalisation is have to be taken into acount for second derivatives :
            //
            //   d2(lambdan)/d(epsilon_normalized_ij)d(epsilon_normalized_kl) = nf . d2(lambda)/d(epsilon_ij)d(epsilon_kl)
            //

            // computing second derivative products
            //dscal_(&i,  &zero, der2_i_e, &one);
            // ddot_(int* N,  double * x, int *incx, double *y, int *incy)

            /*
       * test with blas : little slower and with  bugge
       double h;
       for (m = 0; m < 3; ++m)
       {
       for (n = 0; n < 3; ++n)
       {

       i = 81;
       j = 1;
       daxpy_(&i, &der_l_i[n][m], &der2_i_e[m][0][0][0][0], &j, &dersec[n][0][0][0][0], &j);
       for ( ki = 0; ki < 3; ++ki)
       {
       h = der2_l_i[n][m][ki];
       v = &der_i_e[ki][0][0];
       for (i = 0; i < 3; ++i)
       {

       for (j = 0; j < 3; ++j)
       {
       r = der_i_e[m][i][j]*h;
       u = &dersec[n][i][j][0][0];
       for ( k = 0; k < 9; k++)
       {
       u[k] += v[k]*r;

       }
       }
       }
       }
       }
       }
       r = one/nf;
       i = 243;
       j = 1;
       dscal_(&i,  &r,&dersec[0][0][0][0][0], &j);
      */
            for (n = 0; n < 3; ++n)
            {
                for (i = 0; i < 3; ++i)
                {

                    for (j = 0; j < 3; ++j)
                    {
                        for ( k = 0; k < 3; ++k)
                        {
                            for ( l = 0; l < 3; ++l)
                            {
                                r = zero;
                                for (m = 0; m < 3; ++m)
                                {
                                    r += der_l_i[n][m]*der2_i_e[m][i][j][k][l];
                                    for ( ki = 0; ki < 3; ++ki)
                                    {
                                        r += der2_l_i[n][m][ki]*der_i_e[ki][k][l]*der_i_e[m][i][j];
                                    }
                                }

                                // normalization procedure
                                dersec[n][i][j][k][l] = r/nf;
                            }
                        }
                    }
                }
            }

            for ( k = 0; k < 3; ++k) //symmetric tensors
            {
                for ( i = 0; i < 3; ++i)
                {
                    for ( j = 0; j < 3; ++j)
                    {
                        dersec[k][i][j][1][0] += dersec[k][i][j][0][1];
                        dersec[k][i][j][2][0] += dersec[k][i][j][0][2];
                        dersec[k][i][j][2][1] += dersec[k][i][j][1][2];
                        dersec[k][i][j][1][0] /= 2.;
                        dersec[k][i][j][2][0] /= 2.;
                        dersec[k][i][j][2][1] /= 2.;
                        dersec[k][i][j][0][1] = dersec[k][i][j][1][0];
                        dersec[k][i][j][0][2] = dersec[k][i][j][2][0];
                        dersec[k][i][j][1][2] = dersec[k][i][j][2][1];
                    }
                }
            }

        } // end second derivative
    }    // end first derivative

    //transform lambdan in lambda
    lambda[0] *= nf;
    lambda[1] *= nf;
    lambda[2] *= nf;

    return;

}

} // end of namespace

