/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
 */

#include <iostream>
#include <sstream>
#include <fstream>
#include <complex>

// eXlibris_tools
#include "xMPIEnv.h"
// xlinalg
#include "xDistIndex.h"
#include "xDistVector.h"

/*
 * Testing xDistVector class with following example :
 *  On 3 procs each proc have the following set of Fortran index
 *
 *       p0        p1     p2
 *                        16
 *       1                1
 *       2
 *                 3      3
 *                 5
 *                 13
 *                 15
 *       8         8      8
 *       10
 *
 *
 *  Tested : all  public method
 *
 */

using namespace std;


/* function that do the test: template on arithmetic
 */
template < typename T >
void do_test(ofstream &ofs, xlinalg::xDistIndex & idx, xlinalg::xDistIndex & idxi,const T & facti,
             const T &coef1,
             const T &coef13,
             const T &coef05,
             const T &coef2,
             const T &coef27,
             const T &coef3,
             const T &coef2340,
             const T &coef35,
             int proc_id)
{
    ofs <<"Test with "<<xtool::xDataType < T >::stype()<<" arithmetic"<<endl;
    // create 4 vectors with idx  : x is a unit vector initially sprayed on all proc
    //                              y is made of Fortran index
    //                              z is uninitialized
    // create 1 vectors with idxi : w is made of Fortran index
    xlinalg::xDistVector < T > x(idx);
    xlinalg::xDistVector < T > y(idx);
    xlinalg::xDistVector < T > z(idx);
    xlinalg::xDistVector < T > w(idxi);

    for (auto i : idx) y[idx.getPackedIndex(i)] = facti*static_cast < double >( i );
    for (auto i : idxi) w[idxi.getPackedIndex(i)] = facti*static_cast < double >( i );
    switch (proc_id)
    {
        case 0 :
         {
             x.AddVal(1,coef05);
             x.getVal(2) = coef1;
             x[2] = coef13;
             x.at(3) = coef1;
             break;
         }
        case 1 :
         {
             T *p = x.data();
             for (int i = 0; i < 5; ++i) p[i] = coef1;
             x.getVal(3) = coef05;
             x.getVal(8) = coef13;
             break;
         }
        case 2 :
         {
             for (auto i : idx) x[idx.getPackedIndex(i)] = coef05;
             x.getVal(16) = coef1;
             x.getVal(8) = coef13;
             break;
         }
    }

    ofs <<"Test retriving distIndex from x : "<<( &x.getDistIndex() == &idx )<<endl;
    ofs <<"Initial x state : "<<endl;
    x.Printdata(ofs);
    // do some computation -> insert off
    x.switchInsertModeOff();
    // test dot with same vector, one transfert
    ofs <<"dot  : x.x  "<<x.dot(x)<<endl;
    ofs <<"x state after dot : "<<endl;
    x.Printdata(ofs);
    x.switchToGlobalValue();
    ofs <<"x state after passing in global state : "<<endl;
    x.Printdata(ofs);

    z.switchInsertModeOff();
    z = x;
    ofs <<"z copy of X when in global state : "<<endl;
    z.Printdata(ofs);
    ofs <<"x state after being copied : "<<endl;
    const xlinalg::xDistVector < T > & xc = x;
    std::copy(xc.begin(), xc.end(), std::ostream_iterator < T >(ofs, " "));
    ofs <<endl;

    ofs <<"Initial y state : "<<endl;
    y.Printdata(ofs);

    std::vector < T > A;
    y.switchInsertModeOff();
    y.gather(A,2);
    ofs <<"y state after gather : "<<endl;
    y.Printdata(ofs);
    ofs <<"A state : "<<endl;
    std::copy(A.begin(), A.end(), std::ostream_iterator < T >(ofs, " "));
    ofs <<endl;

    std::vector < T > B;
    w.switchInsertModeOff();
    ofs <<"Initial w state : "<<endl;
    w.Printdata(ofs);
    w.gather(B,1);
    ofs <<"w state after gather : "<<endl;
    w.Printdata(ofs);
    ofs <<"B state : "<<endl;
    std::copy(B.begin(), B.end(), std::ostream_iterator < T >(ofs, " "));
    ofs <<endl;

    int k = 20;
    for (auto & i : B ) i = facti*static_cast < double >( ++k );
    ofs <<"new B state : "<<endl;
    std::copy(B.begin(), B.end(), std::ostream_iterator < T >(ofs, " "));
    ofs <<endl;
    w.switchInsertModeOn();
    w.scatter(B,1);
    ofs <<"w state after scatter in insert mode on : "<<endl;
    w.Printdata(ofs);
    w.switchInsertModeOff();
    w.scatter(B,1);
    ofs <<"w state after scatter in insert mode off (auto switching modes) : "<<endl;
    w.Printdata(ofs);
    w.switchToGlobalValue();
    w.scatter(B,1);
    ofs <<"w state after scatter in insert mode off and global value(auto switching modes) : "<<endl;
    w.Printdata(ofs);
    w.switchInsertModeOn();
    w.switchInsertModeOff();
    w.scatter(B,1);
    ofs <<"w state after scatter in insert mode off and not in reduced state (auto switching modes) : "<<endl;
    w.Printdata(ofs);


    // test axpy with this in special internal reduced mode and arg in globale mode
    y.axpy(coef2,x);
    ofs <<"y state after y=2*x+y  : "<<endl;
    y.Printdata(ofs);

    x.switchToLocalValue();
    ofs <<"x state after passing in local state : "<<endl;
    x.Printdata(ofs);


    y.switchToGlobalValue();
    // change z
    z.switchInsertModeOn();
    for (auto i : idx) z[idx.getPackedIndex(i)] = facti*static_cast < double >( i );
    ofs <<"z state after modification  : "<<endl;
    z.Printdata(ofs);
    z.switchInsertModeOff();
    // test axpy with this in global mode and arg in local mode
    y.axpy(-xtool::xDataType < T >::one(),z);
    ofs <<"z state after y=-1*z+y  : "<<endl;
    z.Printdata(ofs);
    ofs <<"y state after y=-1*z+y  : "<<endl;
    y.Printdata(ofs);

    // test dot with two vector, one global and on local
    ofs <<"dot  : y.z  "<<y.dot(z)<<endl;
    ofs <<"z state after dot y.z  : "<<endl;
    z.Printdata(ofs);
    ofs <<"y state after dot y.z  : "<<endl;
    y.Printdata(ofs);

    // test axpy with this in global mode and arg in special internal reduced mode
    y.axpy(-xtool::xDataType < T >::one(),z);
    ofs <<"z state after y=-1*z+y  : "<<endl;
    z.Printdata(ofs);
    ofs <<"y state after y=-1*z+y  : "<<endl;
    y.Printdata(ofs);


    // test axpy with this in special internal reduced mode and arg in local
    x.switchInsertModeOn();
    x[1] += coef1;
    x.switchInsertModeOff();
    z.axpy(-xtool::xDataType < T >::one(),x);
    ofs <<"x state after z=-1*x+z  : "<<endl;
    x.Printdata(ofs);
    ofs <<"z state after z=-1*x+z  : "<<endl;
    z.Printdata(ofs);


    // check nrm2
    ofs <<"nrm2  : z  "<<z.nrm2()<<endl;
    ofs <<"z state after nrm2  : "<<endl;
    z.Printdata(ofs);

    z.switchToGlobalValue();
    // test axpy with this in global mode and arg in global mode (no special treatement)
    y.axpy(-xtool::xDataType < T >::one(),z);
    ofs <<"z state after y=-1*z+y  : "<<endl;
    z.Printdata(ofs);
    ofs <<"y state after y=-1*z+y  : "<<endl;
    y.Printdata(ofs);

    // x in local mode we can do a collective or local scale
    x.scal(coef35);
    if (proc_id) x.scal(coef2);
    ofs <<"x state after scaling by 3.5 and locally by 2.  : "<<endl;
    x.Printdata(ofs);
    // in global mode only collective scal
    x.switchToGlobalValue();
    x.scal(coef27);
    ofs <<"x state after scaling by 2./7.  : "<<endl;
    x.Printdata(ofs);


    // here we don't care about state as all these operators rely on all ready tested dot,axpy and scal
    // check operator + - *
    x.switchToLocalValue();
    y.switchToLocalValue();
    z.switchToLocalValue();
    ofs <<"x,y,z state before test on operator + - *  : "<<endl;
    x.Printdata(ofs);
    y.Printdata(ofs);
    z.Printdata(ofs);
    z = ( x+y )*(( z*y )/coef2340 )-( z*coef3 );
    ofs <<"x,y,z state after test on operator + - * : "<<endl;
    x.Printdata(ofs);
    y.Printdata(ofs);
    z.Printdata(ofs);

    // check operator += -=
    z -= y;
    y += x;
    ofs <<"y,z state after test on operator += -= : "<<endl;
    y.Printdata(ofs);
    z.Printdata(ofs);

    // test dot with two vector in local mode (double reduce testing)
    y.switchInsertModeOn();
    z.switchInsertModeOn();
    for (auto i : idx) y[idx.getPackedIndex(i)] = facti*static_cast < double >( i );
    for (auto i : idx) z[idx.getPackedIndex(i)] = facti*static_cast < double >( 16-i );
    y.switchInsertModeOff();
    z.switchInsertModeOff();
    ofs <<"z state before dot y.z  : "<<endl;
    z.Printdata(ofs);
    ofs <<"y state before dot y.z  : "<<endl;
    y.Printdata(ofs);
    ofs <<"dot  : y.z  "<<y.dot(z)<<endl;
    ofs <<"z state after dot y.z  : "<<endl;
    z.Printdata(ofs);
    ofs <<"y state after dot y.z  : "<<endl;
    y.Printdata(ofs);

    // test dot with two vector, one local and on global
    y.switchInsertModeOn();
    z.switchInsertModeOn();
    for (auto i : idx) y[idx.getPackedIndex(i)] = facti*static_cast < double >( i );
    for (auto i : idx) z[idx.getPackedIndex(i)] = facti*static_cast < double >( 16-i );
    y.switchInsertModeOff();
    z.switchInsertModeOff();
    z.switchToGlobalValue();
    ofs <<"z globale state before dot y.z  : "<<endl;
    z.Printdata(ofs);
    ofs <<"y local state before dot y.z  : "<<endl;
    y.Printdata(ofs);
    ofs <<"dot  : y.z  "<<y.dot(z)<<endl;
    ofs <<"z state after dot y.z  : "<<endl;
    z.Printdata(ofs);
    ofs <<"y state after dot y.z  : "<<endl;
    y.Printdata(ofs);


    // test componentProduct ====
    x.switchToGlobalValue();
    y.switchToGlobalValue();
    z.switchToGlobalValue();
    ofs <<"global z x and y state before componentProduct z = x .*. y: "<<endl;
    z.Printdata(ofs);
    x.Printdata(ofs);
    y.Printdata(ofs);
    z.componentProduct(x,y);
    ofs <<"z state after componentProduct z = x .*. y: "<<endl;
    z.Printdata(ofs);
    // put y in local mode having all value in reduced mode
    y.switchInsertModeOn();
    y.switchInsertModeOff();
    ofs <<" x global  y local before componentProduct z = x .*. y: "<<endl;
    x.Printdata(ofs);
    y.Printdata(ofs);
    z.componentProduct(x,y);
    ofs <<"z state after componentProduct z = x .*. y: "<<endl;
    z.Printdata(ofs);
    // put x in local mode having all value in reduced mode
    x.switchInsertModeOn();
    x.switchInsertModeOff();
    ofs <<" x local  y local before componentProduct z = x .*. y: "<<endl;
    x.Printdata(ofs);
    y.Printdata(ofs);
    z.componentProduct(x,y);
    ofs <<"z state after componentProduct z = x .*. y: "<<endl;
    z.Printdata(ofs);
    ofs <<"x state after componentProduct z = x .*. y: "<<endl;
    x.Printdata(ofs);
    // put x in global state and y in reduced
    x.switchToGlobalValue();
    y.switchToGlobalValue();
    y.switchToLocalValue();
    ofs <<" x global  y reduced before componentProduct z = x .*. y: "<<endl;
    x.Printdata(ofs);
    y.Printdata(ofs);
    z.componentProduct(x,y);
    ofs <<"z state after componentProduct z = x .*. y: "<<endl;
    z.Printdata(ofs);
    // put x in local mode and y stay in reduced
    x.switchInsertModeOn();
    x.switchInsertModeOff();
    ofs <<" x local not in reduced mode but with reduced values,  y reduced before componentProduct z = x .*. y: "<<endl;
    x.Printdata(ofs);
    y.Printdata(ofs);
    z.componentProduct(x,y);
    ofs <<"z state after componentProduct z = x .*. y: "<<endl;
    z.Printdata(ofs);

    // test switchInsertModeOffFromUserGlobalSetting
    std::vector < T > C(idx.getPackedIndexSize());
    x.switchToGlobalValue();
    std::copy(xc.begin(), xc.end(), C.begin());
    ofs <<"x state in global mode: "<<endl;
    x.Printdata(ofs);
    ofs <<"C state (must be the same as x) : "<<endl;
    std::copy(C.begin(), C.end(), std::ostream_iterator < T >(ofs, " "));
    ofs <<endl;
    y.switchInsertModeOn();
    std::copy(C.begin(), C.end(), y.begin());
    ofs <<"y state in insert mode on with global value from x via C: "<<endl;
    y.Printdata(ofs);
    y.switchInsertModeOffFromUserGlobalSetting();
    ofs <<"y state in insert mode off and global value state (should be the same as above): "<<endl;
    y.Printdata(ofs);
    y.switchToLocalValue();
    ofs <<"y state in insert mode off and local value (check global was set switchInsertModeOffFromUserGlobalSetting): "<<endl;
    y.Printdata(ofs);

    // test = with rhs in insert mode on
    // coverage test: no way to see that z state being turn to insert mode on by x state and imediately
    // switch back to insert mode off
    x.switchInsertModeOn();
    for (auto i : idx) x[idx.getPackedIndex(i)] = facti*static_cast < double >( i );
    z = x;
    ofs <<"z copy of X when in insert mode on : "<<endl;
    z.Printdata(ofs);
    ofs <<"x state after being copied : "<<endl;
    std::copy(xc.begin(), xc.end(), std::ostream_iterator < T >(ofs, " "));
    ofs <<endl;
    ofs <<"z state after being copied from x : "<<endl;
    z.Printdata(ofs);

    // force same value in instance
    z=facti;
    ofs <<"z state after being set to arbitrary value : "<<endl;
    z.Printdata(ofs);
    z.switchToGlobalValue();
    z=facti+facti;
    ofs <<"z in globale state after being set to arbitrary value : "<<endl;
    z.Printdata(ofs);


    // throws
    try
    {
        w = x;
    }
    catch (int e)
    {
        ofs <<"throw with = with # xDistIndex ok"<<endl;
    }
    try
    {
        w += z;
    }
    catch (int e)
    {
        ofs <<"throw with axpy and # xDistIndex ok"<<endl;
    }
    try
    {
        z = w + z;
    }
    catch (int e)
    {
        ofs <<"throw with + and # xDistIndex ok"<<endl;
    }
    try
    {
        z = w - z;
    }
    catch (int e)
    {
        ofs <<"throw with - and # xDistIndex ok"<<endl;
    }
    try
    {
        z = w * z;
    }
    catch (int e)
    {
        ofs <<"throw with dot and # xDistIndex ok"<<endl;
    }
    try
    {
        z.componentProduct(w,x);
    }
    catch (int e)
    {
        ofs <<"throw with componentProduct and # xDistIndex ok"<<endl;
    }
}



int main(int argc, char *argv[])
{
    // initialize mpi universe
    xtool::xMPIEnv::init(argc,argv);

    // local variable
    MPI_Comm world = MPI_COMM_WORLD;
    int proc_id,nb_proc;

    // get rank for world
    MPI_Comm_rank(world, &proc_id);
    MPI_Comm_size(world, &nb_proc);


    if (nb_proc != 3)
    {
        cout << argv[0] << " Need 3 mpi process to run "<< endl;
        MPI_Abort(world,-1);
    }

    // open reference
    ofstream ofs("proc_"+std::to_string(proc_id)+"_ref.txt",ios::out | ios::trunc);

    // output
    ofs<<"In proc "<<proc_id<<" Starting"<<endl;

    // create index to be used with vector instances
    // note : this is the same input as xDistIndex test, so this part is already checked.
    xlinalg::xDistIndex idx(world);
    xlinalg::xDistIndex idxi(world);
    switch (proc_id)
    {
        case 0 :
         {
             idx.insertIndex(1);
             idx.insertIndex(2);
             idx.insertIndex(8);
             idx.insertToFrom(1,2);
             idx.insertIndex(10);
             idx.insertToFrom(8,1);
             idx.insertToFrom(8,2);
             idx.finalize(4,10,true); // here we say that it is in increasing order but it is false. Just to see the difference with other index when using gather
             idxi.insertIndex(1);
             idxi.insertIndex(2);
             idxi.insertIndex(3);
             idxi.insertIndex(4);
             idxi.insertToFrom(1,2);
             idxi.insertToFrom(3,1);
             idxi.insertToFrom(3,2);
             idxi.finalize(4,4,true);
             break;
         }
        case 1 :
         {
             idx.insertIndex(3);
             idx.insertIndex(5);
             idx.insertIndex(13);
             idx.insertIndex(15);
             idx.insertIndex(8);
             idx.insertToFrom(3,2);
             idx.insertToFrom(8,0);
             idx.insertToFrom(8,2);
             idx.finalize(4,15,true); // here we say that it is in increasing order but it is false. Just to see the difference with other index when using gather
             idxi.insertIndex(5);
             idxi.insertIndex(6);
             idxi.insertIndex(7);
             idxi.insertIndex(8);
             idxi.insertIndex(3);
             idxi.insertToFrom(5,2);
             idxi.insertToFrom(3,0);
             idxi.insertToFrom(3,2);
             idxi.finalize(4,8,true);
             break;
         }
        case 2 :
         {
             idx.insertIndex(16);
             idx.insertIndex(1);
             idx.insertIndex(8);
             idx.insertIndex(3);
             idx.insertToFrom(1,0);
             idx.insertToFrom(3,1);
             idx.insertToFrom(8,0);
             idx.insertToFrom(8,1);
             idx.finalize(1,16,true); // here we say that it is in increasing order but it is false. Just to see the difference with other index when using gather
             idxi.insertIndex(9);
             idxi.insertIndex(1);
             idxi.insertIndex(5);
             idxi.insertIndex(3);
             idxi.insertToFrom(1,0);
             idxi.insertToFrom(3,0);
             idxi.insertToFrom(3,1);
             idxi.insertToFrom(5,1);
             idxi.finalize(1,9,true);
             break;
         }
    }

    do_test < float >(ofs,idx,idxi,1.f,1.f,1.f/3.f,0.5f,2.f,2.f/7.f,3.f,2340.f, 3.5f, proc_id);
    do_test < double >(ofs,idx,idxi,1.,1.,1./3.,0.5,2.,2./7.,3.,2340.,3.5, proc_id);
    do_test < complex < double > >(ofs,idx,idxi,complex < double >(1.,-1.),complex < double >(1.,-1.),complex < double >(1./3.,-1./3.),complex < double >(0.5,-0.5),
                                   complex < double >(2.,0.),complex < double >(2./7.,0.),complex < double >(3.,0.),complex < double >(2340.,0.), complex < double >(3.5,-3.5),proc_id);

    ofs.close();
    return xtool::xMPIEnv::finalize();
}
