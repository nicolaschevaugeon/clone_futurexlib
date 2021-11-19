/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.
*/

#include <iostream>
#include <sstream>
#include <fstream>

// eXlibris_tools
#include "xMPIEnv.h"
// xlinalg
#include "xDistIndex.h"
#include "xDistVector.h"
#include "xDistBlockDenseMatrix.h"

/*
 * Testing xDistBlockDenseMatrix class with following example :
 *  The 14x14 matrix is despatched on 4 procs, each proc having the following
 *  set of dense block matrix at start
 *
 *      P0 1    2   3   4     5    6                         P0 owne: 1-6 connect to 1: 3-4 connect to 2: 3,5-6
 *      1  1    2   3   4     5    6
 *      2  2    4   6   8    10   12
 *      3  3    6  9/3 12/2 15/2 18/2
 *      4  4    8   6  16/2  20   24
 *      5  5   10  7.5  20  25/2 30/2
 *      6  6   12   9   24  15   36/2
 *
 *       P1  7   8     9     3    4                          P1 owne: 7-9 connect to 0: 3-4 connect to 2: 3,8-9
 *       7  49   56   63    21   28
 *       8  56  64/2 72/2  24/2  32
 *       9  63   36  81/2  27/2  36
 *       3  21  24/2 27/2   9/3 12/2
 *       4  28   32   36     6  16/2
 *
 *       P2  10   11    12   3    5    6     8     9           P2 owne: 10-12 connect to 0: 3,5-6 connect to 1: 3,8-9 connect to 3: 11-12
 *       10 100  110   120  30   50   60    80    90
 *       11 110 121/2 132/2 33   55   66    88    99
 *       12 120   66  144/2 36   60   72    96   108
 *        3  30   33    36  9/3 15/2 18/2  24/2  27/2
 *        5  50   55    60  7.5 25/2 30/2   40    45
 *        6  60   66    72   9   15  36/2   48    54
 *        8  80   88    96 24/2   40   48  64/2  72/2
 *        9  90   99   108 27/2   45   54  72/2  81/2
 *
 *
 *       P3  13  14   11     12                              P3 owne: 13-14 connect to 2: 11-12
 *       13 169 182  143    156
 *       14 182 196  154    168
 *       11 143 154 121/2  132/2
 *       12 156 168 132/2  144/2
 *
 * whitch is representing :
 *
 *       1  2  3  4  5  6  7  8   9  10  11  12  13   14
 *    1  1  2  3  4  5  6  0  0   0   0   0   0   0    0
 *    2  2  4  6  8 10 12  0  0   0   0   0   0   0    0
 *    3  3  6  9 12 15 18 21 24  27  30  33  36   0    0
 *    4  4  8 12 16 20 24 28 32  36   0   0   0   0    0
 *    5  5 10 15 20 25 30  0 40  45  50  55  60   0    0
 *    6  6 12 18 24 30 36  0 48  54  60  66  72   0    0
 *    7  0  0 21 28  0  0 49 56  63   0   0   0   0    0
 *    8  0  0 24 32 40 48 56 64  72  80  88  96   0    0
 *    9  0  0 27 36 45 54 63 72  81  90  99 108   0    0
 *   10  0  0 30  0 50 60  0 80  90 100 110 120   0    0
 *   11  0  0 33  0 55 66  0 88  99 110 121 132 143  154
 *   12  0  0 36  0 60 72  0 96 108 120 132 144 156  168
 *   13  0  0  0  0  0  0  0  0   0   0 143 156 169  182
 *   14  0  0  0  0  0  0  0  0   0   0 154 168 182  196
 *
 */

using namespace std;

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


    if (nb_proc != 4)
    {
        cout << argv[0] << " Need 4 mpi process to run "<< endl;
        MPI_Abort(world,-1);
    }

    // open reference
    ofstream ofs("proc_"+std::to_string(proc_id)+"_ref.txt",ios::out | ios::trunc);

    // output
    ofs<<"In proc "<<proc_id<<" Starting"<<endl;

    // create index to be used with matrix instances
    xlinalg::xDistIndex idx(world);
    switch (proc_id)
    {
        case 0 :
         {
             //      P0 1    2   3   4     5    6                         P0 owne: 1-6 connect to 1: 3-4 connect to 2: 3,5-6
             idx.insertIndex(1);
             idx.insertIndex(2);
             idx.insertIndex(3);
             idx.insertIndex(4);
             idx.insertIndex(5);
             idx.insertIndex(6);
             idx.insertToFrom(3,1);
             idx.insertToFrom(3,2);
             idx.insertToFrom(4,1);
             idx.insertToFrom(5,2);
             idx.insertToFrom(6,2);
             idx.finalize(6,6,true);
             break;
         }
        case 1 :
         {
             //       P1  7   8     9     3    4                          P1 owne: 7-9 connect to 0: 3-4 connect to 2: 3,8-9
             idx.insertIndex(7);
             idx.insertIndex(8);
             idx.insertIndex(9);
             idx.insertIndex(3);
             idx.insertIndex(4);
             idx.insertToFrom(3,0);
             idx.insertToFrom(4,0);
             idx.insertToFrom(3,2);
             idx.insertToFrom(8,2);
             idx.insertToFrom(9,2);
             idx.finalize(3,9,true);
             break;
         }
        case 2 :
         {
             //       P2  10   11    12   3    5    6    8    9           P2 owne: 10-12 connect to 0: 3,5-6 connect to 1: 3,8-9 connect to 3: 11-12
             idx.insertIndex(10);
             idx.insertIndex(11);
             idx.insertIndex(12);
             idx.insertIndex(3);
             idx.insertIndex(5);
             idx.insertIndex(6);
             idx.insertIndex(8);
             idx.insertIndex(9);
             idx.insertToFrom(3,0);
             idx.insertToFrom(3,1);
             idx.insertToFrom(5,0);
             idx.insertToFrom(6,0);
             idx.insertToFrom(8,1);
             idx.insertToFrom(9,1);
             idx.insertToFrom(11,3);
             idx.insertToFrom(12,3);
             idx.finalize(3,12,true);
             break;
         }
        case 3 :
         {
             //      P3  13  14   11     12                              P3 owne: 13-14 connect to 2: 11-12
             idx.insertIndex(13);
             idx.insertIndex(14);
             idx.insertIndex(11);
             idx.insertIndex(12);
             idx.insertToFrom(11,2);
             idx.insertToFrom(12,2);
             idx.finalize(2,14,true);
             break;
         }
    }

    // create empty dense block matrix based on index distribution
    xlinalg::xDistBlockDenseMatrix < double,xTraitMatrixDefinitePositive,xTraitMatrixUnSym > A(idx);

    int ng,nl;
    double *data;
    A.getMemoryAccess(ng,nl,&data);
    for (auto idi : idx)
    {
        const int i = idx.getPackedIndex(idi);
        for (auto idj : idx)
        {
            const int j = idx.getPackedIndex(idj);
            data[i+nl*j] = idi*idj;
        }
    }
    switch (proc_id)
    {
        case 0 :
         {
             /*          0    1   2   3     4    5
              *      P0  1    2   3   4     5    6
              *   0   1  1    2   3   4     5    6
              *   1   2  2    4   6   8    10   12
              *   2   3  3    6  9/3 12/2 15/2 18/2
              *   3   4  4    8   6  16/2  20   24
              *   4   5  5   10  7.5  20  25/2 30/2
              *   5   6  6   12   9   24  15   36/2
              *      */
             for (int j = 3; j < 6; ++j)
             {
                 data[2+nl*j] /= 2.;
                 data[j+nl*2] /= 2.;
             }
             data[2+nl*2] /= 3.;
             data[3+nl*3] /= 2.;
             data[4+nl*4] /= 2.;
             data[4+nl*5] /= 2.;
             data[5+nl*4] /= 2.;
             data[5+nl*5] /= 2.;
             break;
         }
        case 1 :
         {
             /*
              *            0   1     2     3    4
              *        P1  7   8     9     3    4
              *   0    7  49   56   63    21   28
              *   1    8  56  64/2 72/2  24/2  32
              *   2    9  63   36  81/2  27/2  36
              *   3    3  21  24/2 27/2   9/3 12/2
              *   4    4  28   32   36     6  16/2
              *       */
             data[1+nl*1] /= 2.;
             data[1+nl*2] /= 2.;
             data[1+nl*3] /= 2.;
             data[2+nl*1] /= 2.;
             data[2+nl*2] /= 2.;
             data[2+nl*3] /= 2.;
             data[3+nl*1] /= 2.;
             data[3+nl*2] /= 2.;
             data[3+nl*3] /= 3.;
             data[3+nl*4] /= 2.;
             data[4+nl*3] /= 2.;
             data[4+nl*4] /= 2.;
             break;
         }
        case 2 :
         {
             /*             0    1     2   3    4    5    6    7
              *        P2  10   11    12   3    5    6    8    9
              *  0     10 100  110   120  30   50   60   80   90
              *  1     11 110 121/2 132/2 33   55   66   88   99
              *  2     12 120   66  144/2 36   60   72   96  108
              *  3      3  30   33    36  9/3 15/2 18/2  24/2 27/2
              *  4      5  50   55    60  7.5 25/2 30/2  40   45
              *  5      6  60   66    72   9   15  36/2  48   54
              *  6      8  80   88    96 24/2   40   48  64/2 72/2
              *  7      9  90   99   108 27/2   45   54  72/2 81/2
              */
             data[1+nl*1] /= 2.;
             data[1+nl*2] /= 2.;
             data[2+nl*1] /= 2.;
             data[2+nl*2] /= 2.;
             data[3+nl*3] /= 3.;
             data[3+nl*4] /= 2.;
             data[3+nl*5] /= 2.;
             data[3+nl*6] /= 2.;
             data[3+nl*7] /= 2.;
             data[4+nl*3] /= 2.;
             data[4+nl*4] /= 2.;
             data[4+nl*5] /= 2.;
             data[5+nl*3] /= 2.;
             data[5+nl*4] /= 2.;
             data[5+nl*5] /= 2.;
             data[6+nl*3] /= 2.;
             data[6+nl*6] /= 2.;
             data[6+nl*7] /= 2.;
             data[7+nl*3] /= 2.;
             data[7+nl*6] /= 2.;
             data[7+nl*7] /= 2.;
             break;
         }
        case 3 :
         {
             /*             0   1    2      3
              *        P3  13  14   11     12
              *  0     13 169 182  143    156
              *  1     14 182 196  154    168
              *  2     11 143 154 121/2  132/2
              *  3     12 156 168 132/2  144/2
              */
             data[2+nl*2] /= 2.;
             data[2+nl*3] /= 2.;
             data[3+nl*2] /= 2.;
             data[3+nl*3] /= 2.;
             break;
         }
    }
    A.Printdata(ofs);

    // find reduce on owner
    int nb_owned = idx.getPackedLocalIndexSize();
    for (auto idi : idx)
    {
        const int i = idx.getPackedIndex(idi);
        if (i < nb_owned)
        {
            for (auto idj : idx)
            {
                const int j = idx.getPackedIndex(idj);
                // based on fact that we input half(third) term
                if (j < nb_owned && ( data[i+nl*j] != idi*idj ) )
                    ofs<<" term "<<idi<<" "<<idj<<" have to be reduced"<<std::endl;
            }
        }
    }

    // create vectors
    xlinalg::xDistVector< > x(idx);
    xlinalg::xDistVector< > y(idx);
    for (auto idi : idx)
    {
        const int i = idx.getPackedIndex(idi);
        if (i < nb_owned)
            x[i] = 1./idi;
    }

    x.Printdata(ofs);

    // do matrix vector operation y=A.x
    double alpha = 1.;
    double beta = 0.;
    A.gemv(alpha,beta,x,y);
    y.Printdata(ofs);

    // gather on proc 0 to check
    std::vector < double >  V;
    y.gather(V, 0);
    if (!proc_id)
    {
        ofs << "###########################"<<std::endl;
        ofs << "Global vector  is:\n";
        std::copy(V.begin(), V.end(), std::ostream_iterator < double >(ofs, " "));
        ofs << std::endl;
        ofs << "###########################"<<std::endl;

    }

    A.reduceOnOwner();

    A.Printdata(ofs);

    // check reduce on owner
    for (auto idi : idx)
    {
        const int i = idx.getPackedIndex(idi);
        if (i < nb_owned)
        {
            for (auto idj : idx)
            {
                const int j = idx.getPackedIndex(idj);
                if (j < nb_owned && ( data[i+nl*j] != idi*idj ) )
                    ofs<<"error term "<<idi<<" "<<idj<<" is not correct"<<std::endl;
            }
        }
    }

    // do matrix vector operation y=A.x wit A in reduced mode
    y.switchInsertModeOff();
    y = 0.;
    A.gemv(alpha,beta,x,y);
    y.Printdata(ofs);

    // gather on proc 0 to check
    y.gather(V, 0);
    if (!proc_id)
    {
        ofs << "###########################"<<std::endl;
        ofs << "Global vector  is:\n";
        std::copy(V.begin(), V.end(), std::ostream_iterator < double >(ofs, " "));
        ofs << std::endl;
        ofs << "###########################"<<std::endl;

    }

    ofs.close();
    return xtool::xMPIEnv::finalize();
}
