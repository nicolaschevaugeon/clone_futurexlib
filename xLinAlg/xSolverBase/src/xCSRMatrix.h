/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.
*/
#ifndef _CSR_MATRIX__H
#define _CSR_MATRIX__H

#include <set>
#include <vector>
#include <algorithm>
#include <iostream>
#include "xTraitsMatrix.h"



namespace xlinalg 
{

/// a struct for dynamically reallocated arrays
typedef struct {
  int nmax;
  int size; // size of element in the array
  int incr; // number of element to realloc if pool is full
  int n;   
  int isorder; 
  char *array; //array containing the data
} List_T;

void *Malloc(size_t size);
void Free(void *ptr);

List_T *List_Create(int n, int incr, int size);
void    List_Delete(List_T *liste);
void    List_Realloc(List_T *liste,int n);
void    List_Add(List_T *liste, void *data);
int     List_Nbr(List_T *liste);
int     List_Nbr0(List_T *liste);
void    List_Insert(List_T *liste, void *data, int (*fcmp)(const void *a, const void *b));
int     List_Replace(List_T *liste, void *data, int (*fcmp)(const void *a, const void *b));
void    List_Read(List_T *liste, int index, void *data);
void    List_Write(List_T *liste, int index, void *data);
void    List_Put(List_T *liste, int index, void *data);
void    List_Pop(List_T *liste);
void   *List_Pointer(List_T *liste, int index);
void    List_Tri(List_T *liste, int (*fcmp)(const void *a, const void *b));
int     List_Search(List_T *liste, void *data, int (*fcmp)(const void *a, const void *b));
int     List_ISearch(List_T *liste, void *data, int (*fcmp)(const void *a, const void *b));
int     List_Query(List_T *liste, void *data, int (*fcmp)(const void *a, const void *b));
void   *List_PQuery(List_T *liste, void *data, int (*fcmp)(const void *a, const void *b));
int     List_Suppress(List_T *liste, void *data, int (*fcmp)(const void *a, const void *b));
void    List_Reset(List_T *liste);
//Added by Nicolas to define a copy constructeur of CSR_Matrix
//Copy allocate and paste
List_T *List_Copy(List_T *in);
//this obe just paste
//List_T *List_Paste(List_T *in);


/// Sparse Matrix Class
/*!
  The xCSRMatrix class is a class to store and work with Sparse Matrix.
  note that fortan convention for indices is used ...
  first element of a matrix is 1 1, not 0 0.
  this is mainly for historical reason because the fortran sparskit package was used by default.
  it uses the struct List_T for memory managment purpose (memory will be added by "chunk" if initial size is not enough).
  By default, and as long as values are added inside the matrix, the storage is in 
  what is called here the HARWELL_BOEING_STORAGE.
  4 arrays are used for the storage scheme :

  a_  : contain the non-zeros value (structural non zero, if you add a 0.0 in the matrix, it's stored).
        no ordering are assumed here.

  ai_[k] : row index of element k in a_ (démarre à 1).

  ptr_[k]   : array, working as a linked list. ptr_[k] contain the index in a_ of the next element in the column of the matrix, 
           if no next, contain 0

  jptr_[j]  : array containing index in a_ for the first element in column j (démarre à 1).

  
  exemple : A = 1 2 5 

                0 3 0

		1 0 0

  depending on the order in which A was build, we could have :

  a_   = [1 2 5 3 1 ]  (the matrix was build line by line)

  ai_  = [1 1 1 2 3 ]

  ptr_ = [1 2 3]

  jptr_= [5 4 0 0 0 ] 
  
  or 

  a_   = [1 1 2 3 5 ]  (the matrix was build line by line )

  ai_  = [1 3 1 2 1 ]

  ptr_ = [1 3 5]

  jptr_= [2 0 4 0 0 ]


  When the matrix is fully build, this dynamic storage is not stricly needed, but could be usefull none the less.
  
  When a first call is made to PutInCompressedRow(), the array a_ is sorted to be in sorted compressed row format,
  eg : the non zero element are sorted line by line, and in each line, the elements are sorted by column.
  then jptr_ store the begining of each line in a_, and ptr_ store the column of the nth element in a_.
  ai_ is then used as a storage to go back to HARWELL_BOEING_STORAGE. Note that when going back to HB storage, the order
  of the element in a_ are not the same as it was at the beginning... the sorted CRS is keeped ...
  if no element is added, switch back to CRS cost 0.
  
  if any element is added to the matrix (with Matrixtool::xAdd for example), the matrix is switch back to HB storage, and the ordering is lost.
  PutInCompressRow() will cost some ordering.

!*/
class xCSRMatrix {


public :
  //! fixed traits for assemble/solver
  // arbitrary chosen as a unsymetric container with no zero termes (variable graph)
  // arbitrary chosen as a undefined. It's the most general case for solver and is coherent with unsymetric.
  // raver arbitrary again the Storage type is compressed row storage format= Compressed Sparce Row
  // it is arbitrary because in fact storage is a flipping pattern bettwen Harwell and compressed row 
  // from assembler point of view every way of discribing elementary matrices are/may be equivalent for Harwell
  // but for compressed row better take the optimal way
  // from solver point of view fixing Compressed Sparce Row simplify transfert of data if needed and give
  // a consistant interface as all methode like ReturnRowPointer, ... are base on CSR.
  //
  // to resume
  // change if you want matrix_pattern,matrix_defined or matrix_assembly_on_zero or matrix_storage  but it's
  // then your responsability that it works (i.e. depending on what you do implementation to do)
  // You better not change matrix_storage in all case ....
  //
  // A derived class mecanism have been use anyway to change this setting .... see below xCSRMatrixSym
  //
  typedef xTraitMatrixUnSym matrix_pattern;
  typedef xTraitMatrixSparceCSR matrix_storage;
  typedef xTraitMatrixNonSingular matrix_defined;
  typedef xTraitMatrixNoAssemblyOnZero matrix_assembly_on_zero;

  /// Define storage scheme
  /*! 
    HARWELL_BOEING_STORAGE : The format used when building a matrix
    COMPRESSED_ROW_STORAGE : Format  for use with most solver
  !*/
  enum CSR_Matrix_Storage_t {HARWELL_BOEING_STORAGE, COMPRESSED_ROW_STORAGE, UNKNOWN_STORAGE};

  xCSRMatrix ();
  xCSRMatrix (const int n);
  
  ~xCSRMatrix( );

  void Load (const xCSRMatrix& in);
  
  //load a matrix from a stream. the stream must have been created with OutputMatrix
  void Load (std::istream& is);
  
 
  void ReAllocate(const int n, const int nbz, const bool putzero=true);
  void SetSize(const int& n);


  // methode needed to connect this matrix to a solver or whatever
  template <typename STORAGE_TYPE,typename INDEXING >
  void getMemoryAccess(int &ns,int &ms,int &nnzs, int **idx1,int **idx2,double **dat, STORAGE_TYPE &storage,INDEXING &indexing);

  void AddNonZeroForGraph(std::vector< std::set< size_t > > & graph);

  void ExecuteReordering();

  void ReorderArray( double * array );
  void InverseReorderArray( double * array);

  /// add val to element i j of the matrix.
  /*!
    if the matrix is not curently in HB format, it will be put back to HB.
    if the added element is actually a new non zero, the flag Hard_CRS_Needed will
    be set to true, and going back to CRS format will take some time.
  !*/
  void   AddMatrix(int i, int j, double val);

  /// add all non zero elements in 'in', multiply by coeff to the matrix.
  /*!
    internally call AddMatrix(int i, int j, double val), cf this function for more detail;
  !*/
  void   AddMatrix(const xCSRMatrix& in, const double& coeff = 1.0);

  /// return a copy of element i, j of the matrix.
  double GetMatrix(int i, int j) const;

  /// delete every thing in the matrix
  void   ZeroMatrix  ( );

  /// keep the matrix structure, just put all values in a_ to zero.
  void   SoftZeroMatrix ( );

  ///Adds other in the current matrix provided they have the same graph
  void  SoftAddMatrix (const xCSRMatrix& other);

  /// return the number of (structural) non zero element.
  int   ReturnNumberOfNonZero ( ) const ;
  
  /// return a pointer to the array of non zero element.
  double * ReturnValPointer ( );

  /// const version of ReturnValPointer. User can't change the valuein the returned array
  const double * ReturnValPointer ( ) const;

  
  const int    * ReturnRowPointer  ( ) const;
  int    * ReturnRowPointer  ( );
  const int    * ReturnColumIndexPointer( ) const;
  int    * ReturnColumIndexPointer( );


  //to return raw and column shifted to zero, to be compatible 
  //with for instance petsc interface
  int    * ReturnRowPointerShiftedToZero  ( ) const;
  // int    * ReturnRowPointerShiftedToZero  ( void );
  int    * ReturnColumnIndexPointerShiftedToZero( ) const;
  //  int    * ReturnColumnIndexPointerShiftedToZero( void );



  /// tell tghe matrix that the assembly is ended. Does nothing but PutInCompressedRow(); 
  void EndOfAssembly();

  /// Change ordering of element in a_ and values in jptr_, and ptr_ so that the matrix is in sorted compressed row format.
  /*! 
    if the matrix structure was not change between 2 call, internally call PutInCompressedRowSoft(), which cost nothing.
  !*/
  void PutInCompressedRow( ) const;
  
  /// To use only if already put once before in CRS, and if structure did not change. should be private ...
  void PutInCompressedRowSoft( ) const;
  
  /// Change to HarwellBoeing. No cost.
  void PutInHarwellBoeing( ) const;
  void PutInHarwellBoeingSoft( ) const;
  
  /// return the number of lines= number of column.
  int GetNbUnknown() const; //nb lines of the matrix

  inline int ReturnStorageType()
  {return storage;};

  inline void SetMatrixReorderingToFalse( ) 
  { MatrixWasReordered = false; }
  
  /// output the matrix following the HB format. A matrix can then be loaded with the member Load(ostream &)
  void OutputMatrixDebugHBFormat(std::ostream &) const;

  /// output the matrix following the compressed row  format. A matrix can then be loaded with the member Load(ostream &)
  void OutputMatrixDebugCompressedRowFormat(std::ostream &) const;


  /// output the matrix in the i j val format, suitable for taucs.
  /*! the function is not const, even if the matrix is not modified, since the iterators didn't have a const version
      at the time it was written
  */
  void OutputMatrixij(std::ostream &os);
  
  
  void OutputMatrixMatlabFormat(std::ostream &os, const std::string& name="A", const unsigned precision=16) const ;
  void OutputMatrixOctaveFormat(std::ostream &out) const;
  void OutputMatrixOctaveFormat(const std::string & filename) const;

  void cmkreord(int n, double *a, int *ja, int *ia, double *a0,
		int *ja0, int *ia0, int *init, int * iperm, int * mask,
		int * maskval, int * nlev, int * riord, int * levels);
  void sort_col(int n,double * a, int * ja, int *ia, int * iw, double * rw);


  /// an iterator on the non zero element of a xCSRMatrix
  struct iterator{
    int ncurrent;
    xCSRMatrix *themat;
    double & operator *();
    iterator & operator++();
    bool operator!=(const iterator& rhs) const;
  };
  iterator begin();
  iterator end();

  /// an iterator on the non zero element of a column of a xCSRMatrix. HARWELL_BOEING Storage is assumed
  /*!
    when iterating over a column, the matrix element on which point the iterator comes in the order in which the appear !
    member variable line current containt the to which the iterator currently point.
    operator * return a reference to the current (nonzero) matrix element.
  !*/
  struct column_iterator_HB{
    int column;
    int ncurrent;
    int line;
    xCSRMatrix *themat;
    double & operator *();
    column_iterator_HB & operator++();
    bool operator!=(const column_iterator_HB& rhs) const;
  };
  
  column_iterator_HB begin_column_HB(int j);
  
  column_iterator_HB end_column_HB(int j);

  /// iterator over the column, based on the column_iterator_HB. It's the fastest and recommanded one to iterate over column
  typedef column_iterator_HB column_iterator;

  column_iterator begin_column(int j){
    return begin_column_HB(j);
  }
  
  column_iterator end_column(int j){
    return end_column_HB(j);
  }
  
  /// an iterator on the non zero element of a line of a xCSRMatrix. HARWELL_BOEING Storage is assumed
  /*!
    when iterating over a line, the matrix element on which point the iterator comes in the order in which the appear !
    there is no way to know the actual column ... it would be too expensive with the HARWELL_BOEING FORMAT
    operator * return a reference to the current (nonzero) matrix element.
  !*/
  struct line_iterator_HB{
    //int column; //column number can't be known in an esay way while in HB format.
    int ncurrent;
    int line;
    //int maxfoundcol;
    xCSRMatrix *themat;
    double & operator *(){
      double * a = (double*)(themat->a_->array);
      return a[ncurrent-1];
    }
    
    line_iterator_HB & operator++(){
      ++ncurrent;
      int nnz = themat->ReturnNumberOfNonZero();
      //      std::cout <<  ((int*)( themat->ai_->array))[ncurrent-1] << std::endl;
      while ((ncurrent <= nnz ) && (((int*)( themat->ai_->array))[ncurrent-1] != line)) {
	++ncurrent;
	//	std::cout <<  ((int*)( themat->ai_->array))[ncurrent-1] << std::endl;
      }     
      return *this;
    }

    bool operator!=(const line_iterator_HB& rhs) const{
      return rhs.ncurrent!=ncurrent;
    }
  };
  
  line_iterator_HB begin_line_HB(int i){
    this->PutInHarwellBoeing();
    //    maxfoundcol=0;
    line_iterator_HB it;
    it.line=i;
    it.themat = this;
    it.ncurrent = 0;
    ++it;
    return it;
  };
  
  line_iterator_HB end_line_HB(int i){
    //maxfoundcol=0
    this->PutInHarwellBoeing();
    line_iterator_HB it;
    it.themat = this;
    it.line=i;
    it.ncurrent = it.themat->ReturnNumberOfNonZero()+1;
    return it;
  };


  /// an iterator on the non zero element of a line of a xCSRMatrix. COMPRESSED_ROW_STORAGE is assumed
  /*!
    when iterating over a line, the matrix element on which point the iterator comes in the order in which the appear !
  !*/
  struct line_iterator_CRS{
    //int column; //column number can't be known in an esay way while in HB format.
    int ncurrent;
    int line;
    int column;
    //int maxfoundcol;
    xCSRMatrix *themat;
    double & operator *(){
      double * a = (double*)(themat->a_->array);
      return a[ncurrent-1];
    }
    
    line_iterator_CRS & operator++(){
      ++ncurrent;
      column = ((int*)(themat->ptr_->array))[ncurrent-1];
      return *this;
    }     

    bool operator!=(const line_iterator_CRS& rhs) const{
      return rhs.ncurrent!=ncurrent;
    }
  };
  
  line_iterator_CRS begin_line_CRS(int i){
    //    maxfoundcol=0;
    this->PutInCompressedRow();
    line_iterator_CRS it;
    it.line=i;
    it.themat = this;
    it.ncurrent =((int *) (it.themat->jptr_->array))[i-1];
    it.column = ((int *)(it.themat->ptr_->array))[it.ncurrent-1];
    return it;
  };
  
  line_iterator_CRS end_line_CRS(int i){
    //maxfoundcol=0
    this->PutInCompressedRow();
    line_iterator_CRS it;
    it.themat = this;
    it.line=i;
    it.ncurrent =((int *) (it.themat->jptr_->array))[i];
    it.column = ((int *)(ptr_->array))[it.ncurrent-1];
    return it;
  }; 
  
  /// iterator over non zero element of a line
  /*!
    based on the line_iterator_CRS. It's the fastest and recommanded one to iterate over line 
    !*/

  typedef line_iterator_CRS line_iterator;

  line_iterator begin_line(int i){
    return begin_line_CRS(i);
  }
  
  line_iterator end_line(int i){
    return end_line_CRS(i);
  }
  

  /// an iterator on the non zero element of a column of a xCSRMatrix. COMPRESSED_ROW_STORAGE is assumed
  /*!
    when iterating over a column, the matrix element on which point the iterator comes in the order in which the appear !
    !*/
  struct column_iterator_CRS{
    int ncurrent;
    int line;
    int col;
    //int maxfoundcol;
    xCSRMatrix *themat;
    double & operator *(){
      double * a = (double*)(themat->a_->array);
      return a[ncurrent-1];
    }
    
    column_iterator_CRS & operator++(){
      ++line;
      ++ncurrent;
      int nnz = themat->ReturnNumberOfNonZero();
      while ((ncurrent<=nnz)&&( ((int*)(themat->ptr_->array))[ncurrent-1]!= col )) { 
	++ncurrent;
      };
      return *this;
    }  
    
    column_iterator_CRS &next(){
      int nnz = themat->ReturnNumberOfNonZero();
     
      ++line;
      ncurrent = ((int*)(themat->jptr_->array))[line-1];
      while(ncurrent<=nnz){
	int bl = ((int*)(themat->jptr_->array))[line-1];
	int el = ((int*)(themat->jptr_->array))[line];
	int * itbegl = &((int*)(themat->ptr_->array))[bl-1];
	int * itendl = &((int*)(themat->ptr_->array))[el-1];
	int *it = std::lower_bound(itbegl, itendl, col);
	if ((it!=itendl)&&(*it==col)){
	  ncurrent +=std::distance(itbegl,it);
	  return *this;
	}
	ncurrent = el;
	++line;
      }
      return *this;
    }

    bool operator!=(const column_iterator_CRS& rhs) const{
      return rhs.ncurrent!=ncurrent;
    }
  };
  
  column_iterator_CRS begin_column_CRS(int j){
    this->PutInCompressedRow();
    column_iterator_CRS it;
    it.col=j;
    it.line=0;
    it.themat = this;
    it.ncurrent =1 ;
    it.next();
    return it;
  };
  
  column_iterator_CRS end_column_CRS(int j){
    this->PutInCompressedRow();
    column_iterator_CRS it;
    it.themat = this;
    it.col=j;
    int nnz = it.themat->ReturnNumberOfNonZero();
    it.ncurrent = nnz+1;
    return it;
  };


public :

  xCSRMatrix (const xCSRMatrix& in);
  xCSRMatrix & operator=(const xCSRMatrix& in);

  //to store the matrix
public :
  mutable List_T  *a_, *ai_, *ptr_, *jptr_ ;
  mutable List_T  *swap_nnz_, *swap_nbup1_ ;

  mutable std::vector<int> row_shifted, column_shifted;

private :

  int NbLines;
  bool AllocationDone;
  
  mutable bool Hard_CRS_Needed;  //track if the structure has changed since last call to PutInCompresssedRow()
  mutable bool Hard_HB_Needed;
  mutable CSR_Matrix_Storage_t storage;

  //for the reordering
  bool MatrixWasReordered;
  int     *permr_, *permp_, *rpermr_;

  //for connection
  std::vector<int> extra1;
  template <typename PATTERN, typename STORAGE_TYPE , typename STORAGE_TYPE_E ,typename INDEXING_E >
  void transfertToConnectionParameter(int &ns,int &ms,int &nnzs, int **idx1,int **idx2,double **dat, STORAGE_TYPE_E &storage,INDEXING_E &indexing);

  void sort2(unsigned long n, double arr[], int ai[] , int aj [] ) const;
  void deblign ( int nz , int *ptr , int *jptr , int *ai) const;
  int *ivector(long nl, long nh) const;
  void free_ivector(int *v, long nl, long nh) const;
  int  cmpij(int ai,int aj,int bi,int bj) const;


  int maskdeg(int *ja, int *ia, int *nod, int *mask, int *maskval);
  int rversp(int n, int *riord);
  int perphn(int n, int *ja, int *ia, int *init, int *iperm,
	     int * mask, int * maskval, int * nlev, int * riord, int * levels);
  void exchange(int n, int *iriord, int *iperm);
 
  int dperm(int nrow, double *a, int *ja, int * ia, double * ao,
	    int * jao, int * iao, int * perm, int *qperm, int * job);

  int cperm(int nrow, double * a, int *ja, int * ia, double * ao,
	    int * jao, int * iao, int * perm, int * job);

  void rperm(int nrow,double * a,int * ja,int * ia, double * ao, 
	     int * jao, int * iao, int * perm, int * job);
 
  int bfs(int n, int *ja,int *ia, int * nfirst, int * iperm, int * mask, 
	  int *maskval, int * riord,int *levels, int *nlev);

  int sort_irv(int * itmp, double *rtmp, int * n);

  int add_lvst(int * istart, int * iend, int * nlev, int * riord,
	       int * ja, int * ia, int * mask, int * maskval);

  void Allocate(const int n);
  void Clear();

};

#if 0
class xCSRMatrixSym : public  xCSRMatrix {

public :
  // this derived class is here only to change traits and transform matrix in this container into 
  //  a  lower symetrique defined matrix type  usable by taucs => instance of transfertToConnectionParameter have to be implemented
  //  not done as for now it means a big west of memory :
  //     in transfertToConnectionParameter implementation will need to duplicate the matrix in Harwel/CSR in CSC
  //     in taucs as by default implemetation is using internal fill-in reduction a other copy of the matrix is done !
  //  this lead in the worst caste to have 3 times the entire (struct+value) matrix in memory ....   
  //  This can be reduce to 2 time if external permutation is used but it imply then to redoo the assemble with this xCSRMatrixSym
  //  as it's him which give the graph !!!!!
  //  => use xGenericSparceMatrix 
  //  
  //
  //  if you insite to use this uncoment and implemente the transfertToConnectionParameter
  //  this little part commes from old taucs interface should do what transfertToConnectionParameter instance need to do (good start)
 /*
  *
  * to copy and modify in cc in transfertToConnectionParameter correct instance

    int *colptr = ((taucs_ccs_matrix*) (taucsmat))->colptr;
    int *rowind =((taucs_ccs_matrix*)( taucsmat))->rowind;
    taucs_double * d = ((taucs_ccs_matrix*)(taucsmat))->values.d;
    int k = 0;
    int n = Mat->GetNbUnknown();
      //std::cout << k << std::endl;
      for (int j =1; j<=n ; ++j){
  
        //std::cout << " j-1 " << j-1 << std::endl;
        colptr[j-1]=k;
        matrix_type::column_iterator  itc = Mat->begin_column(j);
        matrix_type::column_iterator  itce = Mat->end_column(j);
#ifdef REORDER
        std::map<int,double> tmp;
#endif
        while (itc!=itce&&itc.line<j){
  	++itc;
        }
        while (itc!=itce){
  	if ((itc.line)>=j){
  	  // std::cout << k << " " << itc.line-1 << std::endl;
#ifdef REORDER
          tmp[itc.line-1]=*itc;
#else
   	  d[k] = *itc;
  	  rowind[k]=itc.line-1;
  	  ++k;
#endif
  	}
  	++itc;
  	
        }
#ifdef REORDER
        std::map<int,double>::iterator itm = tmp.begin();
        std::map<int,double>::iterator itme = tmp.end();
        for (;itm!=itme; ++itm)
        {
  	    rowind[k]=itm->first;
  	    d[k] = itm->second;
  	    ++k;
        }
#endif

        
      }
      // std::cout << "ga" << std::endl;
      colptr[n]=k;

      */
  //  
  //
  //
  //
  //
  //
  //
  typedef xTraitMatrixLowerSym matrix_pattern;
  typedef xTraitMatrixSparceCSR matrix_storage;
  typedef xTraitMatrixDefinitePositive matrix_defined;
  typedef xTraitMatrixNoAssemblyOnZero matrix_assembly_on_zero;

  xCSRMatrixSym (void) : xCSRMatrix() {};
  xCSRMatrixSym (const int n) : xCSRMatrix(n) {};
};
#endif

#include "xCSRMatrix_imp.h"

} // end of namespace

#endif








