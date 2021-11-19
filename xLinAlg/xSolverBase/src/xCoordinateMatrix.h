/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/

#include <algorithm>
#include <cassert>
#include <iostream>
#include <iterator>

template <class T>
class xCoordinateMatrix
{
  public:
   xCoordinateMatrix(int _nrows, int _ncolumns, int _nnz)
       : nrows(_nrows),
         ncolumns(_ncolumns),
         nnz(_nnz),
         values(new T[nnz]),
         rows(new int[nnz]),
         columns(new int[nnz]),
         valuesadded(0)
   {
   }
   T &operator()(const int &i, const int &j)
   {
      // cout <<" values added " <<  valuesadded << endl;
      assert((i < nrows) && (j < ncolumns));
      int jcurrent = -1;
      int pos = -1;
      int *begin = rows;
      int *linestart = 0;
      int *end = &rows[valuesadded];
      while (jcurrent != j)
      {
         linestart = &rows[pos + 1];
         int *pi = find(linestart, end, i);
         // cout << begin << " " << pi << " " << end << std::endl;;
         if (pi == end)
         {
            rows[valuesadded] = i;
            columns[valuesadded] = j;
            valuesadded += 1;
            if (valuesadded > nnz)
            {
               cout << "adding more than nnz values" << endl;
               throw;
            }
            return values[valuesadded - 1];
         }

         if (pi != end)
         {
            pos = distance(begin, pi);
            jcurrent = columns[pos];
         }
      }
      return values[pos];
   };
   ~xCoordinateMatrix()
   {
      delete[] values;
      delete[] rows;
      delete[] columns;
   }
   void print()
   {
      for (int n = 0; n < nnz; ++n) cout << "A(" << rows[n] << ", " << columns[n] << ")= " << values[n] << endl;
   };

  public:
   int nrows;
   int ncolumns;
   int nnz;
   T *values;
   int *rows;
   int *columns;
   int valuesadded;
};
