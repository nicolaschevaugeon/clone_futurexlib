/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.
*/

#ifndef _GRAPHMATRIX_H
#define _GRAPHMATRIX_H


#include <iostream>
#include <vector>
#include <sstream>
#include <fstream>
#include <functional>


namespace xlinalg
{

class xGraphMatrix
{
    public:
        // typedefs
        typedef std::vector < std::vector < int > > graph_matrix_t;

        /// constructor
        //! sym member is fixed here by the user. By  default the wall graph
        //! is stored
        //! If symetrique pattern is considered (trial/test are identical in AssembleGraph) 
        //! sym may be true => we only store the lower triangular part of this symetric pattern
        //! sym may be false => we fully store this symetric pattern. This make sens as the matrix
        //!     may be unsymetric but with symetric pattern. The full graph must be create as the solver
        //!     in this case is designed for unsymetrique problem (LU decomposition)
        //! If unsymetrique pattern is considered (trial/test are different in AssembleGraph) 
        //! sym may be false => we fully store this unsymetric pattern.
        //! sym may be true => This make sens, and correspond to a global symetrique matrice for which we
        //!     add in it's lower triangular part description, a out of diagonal bloc
        //
        //! n et m are the dimenson of the matrix
        //! n and m are given here as there is no simple way for the AssembleGraph to
        //        fixe them as field given may be only part of a largeur problem
        xGraphMatrix(int n_, int m_, bool sym_ = false);


        // Public methodes
        //

        /// methode mandatory at least by AssembleGraph
        //! all these methodes use c indexing (starting at zero)
        //
        //! add, ordered, term (i,j) in the graph structure
        void add(const int i,const int j);
        //
        //! add, ordered, line terms of column j in the graph structure considering all terms
        //! nb = number of line termes in idx
        //! j = column number
        //! idx = line 
        void addLinesUnsymBlock(const int nb,const int j, int *idx);
        //
        //! add, ordered, line terms of column j in the graph structure considering only some terms depending on sym
        //! il = index of the first line term
        //! nb = number of line termes in idx
        //! j = column number
        //! idx = line 
        void addLinesSymBlock(const int il, const int nb,const int j, int *idx);
        //
        //! count non zero
        void countNNZ();
        //
        /// respond true if sym is true
        //nota : also mandatory by xGenericSparseMatrix
        bool isSym() {return sym;};
        // end mandatory at least by AssembleGraph

        // begin mandatory xValueCreatorLinkOnFrontFiltered
        //
        //! remove term from the graph
        //! use c indexing (starting at zero)
        void remove(const int i,const int j);
        //
        //! clear an entire column
        void clearCol(const int j) { matrix_struct[j].clear(); }
        // end mandatory xValueCreatorLinkOnFrontFiltered


        /// methode mandatory at least by xGenericSparseMatrix
        //! getting graph information
        int getN(){ return n; };
        int getM(){ return m; };
        int getNNZ(){ return nnz; };
        int * getCol(const int j){ return &matrix_struct[j][0]; };
        int  getSizeCol(const int j){ return matrix_struct[j].size(); };
        // end mandatory at least by xGenericSparseMatrix

        // This methode, giving  row and colone permutation reoder the whole graph generated so far
        // Size of vector containing permutation have to be given for dimension check
        void ReoderGraphFromReordering(int nr, int *perm_r, int nc, int *perm_c );

        /// compute bandwidth of the matrix
        int getBandWidth();

        /// clear graph struture
        void clear();

        // methode to output graph
        // interface style
        void print();
        // for scilab/matlab input : use  coordinate matrix storage with a arbitrary value of 1 for matrix term value
        // and fortran indexing
        void printCOO(std::ofstream & os);

        /// compute symetrique permutation for bandwith reduction using CutHillMcKee algorithm
        void getPermReverseCutHillMcKee(int nm,int *perm,int *invperm);



    private:
        typedef  std::vector<int>::size_type size_type;

        int n,m,nnz;
        // sym is given by the user to declare if the full pattern of non zero termes of the
        // matrix is stored or only the lower triangular part. 
        bool sym;
        graph_matrix_t matrix_struct;

        /// add, ordered, term (i,j) in the graph structure in lower part
        void addSym(const int i,const int j);
          
        // intermediate for sym and unsym case
        void addLinesSym(const int il,const int nb_tot,const int j, int *line_idx0);
        void addLinesUnSym(const int il,const int nb_tot,const int j, int *line_idx0);

        // function generic for addLinesUnsymBlock
        // depending on sym this is calling sym or unsym intermediate function
        std::function < void (const int,const int, const int,int *) > addLinesFunction;

};


} // end of namespace

#endif
