/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.
*/
#include <fstream>
#include <iostream>
#include <iterator>
#include <iomanip>
#include <cstdlib> 
#include <cmath>
#include <vector>
#include <set>
#include <cassert>
#include <string>
#include <cstring>
#include "xCSRMatrix.h"

using namespace std;

namespace xlinalg
{



#define SWAP(a,b)  temp=(a);(a)=(b);(b)=temp;
#define SWAPI(a,b) tempi=(a);(a)=(b);(b)=tempi;

const unsigned long M_sort2  = 7;
const int NSTACK   = 50;

xCSRMatrix::xCSRMatrix() : 
  a_(nullptr), ai_(nullptr), ptr_(nullptr), jptr_(nullptr), swap_nnz_(nullptr), swap_nbup1_(nullptr),
  NbLines(0), AllocationDone(false),
  Hard_CRS_Needed(true), Hard_HB_Needed(false),  storage(UNKNOWN_STORAGE), MatrixWasReordered(false), 
  permr_(nullptr), permp_(nullptr), rpermr_(nullptr) {}


void xCSRMatrix::Allocate(const int n) 
{
  NbLines = n;  
  storage = HARWELL_BOEING_STORAGE ;
  a_    = List_Create (NbLines, NbLines, sizeof(double));
  ai_   = List_Create (NbLines, NbLines, sizeof(int));
  ptr_  = List_Create (NbLines, NbLines, sizeof(int));
  jptr_ = List_Create (NbLines+1, NbLines, sizeof(int));
  int i, j=0;
  for(i=0; i<NbLines; i++) List_Add (jptr_, &j);
  ZeroMatrix ( );

  // Add a zero to each entry on the diagonal -- this gives better structure
  // and avoids a bug in the matrix vector multiply
  //printf("Before addmatrix in Allocate\n");

  for (int n = 1; n <= NbLines; n++) {
    AddMatrix(n,n,0.0);
  }
  AllocationDone = true;

  delete []rpermr_; rpermr_ = nullptr;
  delete []permr_;  permr_ = nullptr;
  delete []permp_;  permp_ = nullptr;
  
  Free(swap_nnz_); swap_nnz_=nullptr;
  Free(swap_nbup1_); swap_nbup1_=nullptr;
  
  MatrixWasReordered = false;
  Hard_CRS_Needed = true;
}

void xCSRMatrix::ReAllocate(const int n, const int nbz, const bool putzero) 
{
  List_Delete(a_);
  List_Delete(ai_);
  List_Delete(ptr_);
  List_Delete(jptr_);
  delete []rpermr_; rpermr_ = nullptr;
  delete []permr_;  permr_ = nullptr;
  delete []permp_;  permp_ = nullptr;
  Free(swap_nnz_); swap_nnz_=nullptr;
  Free(swap_nbup1_); swap_nbup1_=nullptr;
  
  NbLines = n;  
  storage = HARWELL_BOEING_STORAGE ;
  a_    = List_Create (nbz, NbLines, sizeof(double));
  ai_   = List_Create (nbz, NbLines, sizeof(int));
  ptr_  = List_Create (nbz, NbLines, sizeof(int));
  jptr_ = List_Create (NbLines+1, NbLines, sizeof(int));
  int i, j=0;
  for(i=0; i<NbLines; i++) List_Add (jptr_, &j);
  ZeroMatrix ( );
  if (putzero){
    for (int n = 1; n <= NbLines; n++) 
      {
	AddMatrix(n,n,0.0);
      }
  }
  AllocationDone = true;
  MatrixWasReordered = false;
  Hard_CRS_Needed = true;
}



void xCSRMatrix::AddNonZeroForGraph(std::vector< std::set< size_t > > & graph)
{

  std::vector< std::set< size_t > >::const_iterator it  = graph.begin();
  std::vector< std::set< size_t > >::const_iterator ite  = graph.end();
  int i = 0;
  for (; it!=ite; ++it,++i){
    const std::set< size_t >& s = *it;
    std::set< size_t >::const_iterator its   = s.begin();
    std::set< size_t >::const_iterator ites  = s.end();
    for (; its!=ites; ++its){
      AddMatrix(i+1,*its+1,0.0);
    }
  }
  PutInCompressedRow();
  return;
}


void xCSRMatrix::Load(std::istream &is) {
  int nnz; //number of nonzero
  int nlines; //number of lines
  is >> nlines;
  is >> nnz;
  ReAllocate(nlines, nnz, false);
  a_->n = nnz;
  ptr_->n = nnz;
  ai_->n = nnz;
  jptr_->n = nlines;
  double *a = (double *)a_->array;
  int *ai   = (int *)ai_->array;
  int *jptr = (int *)jptr_->array;
  int *ptr  = (int *)ptr_->array;
  for (int i=0; i < nnz; i++){
    is >> a[i];
  }
  for (int i=0; i < nnz; i++){
    is >> ai[i];
  }
  for (int i=0; i < nnz; i++){
    is >> ptr[i];
  }
  for (int i=0; i < nlines; i++){
    is >> jptr[i];
  }
  Hard_CRS_Needed = true;
}

void xCSRMatrix::OutputMatrixij(std::ostream &os) {
  int nl = GetNbUnknown();
  for (int i=1; i<= nl; i++){
    column_iterator it = begin_column(i);
    column_iterator ite = end_column(i);
    while (it !=ite){
      os << it.line << " " << it.column << " " << std::scientific << std::showpoint << std::setprecision(17) << *it << std::endl;
      ++it;
    }
  }  
}

void xCSRMatrix::OutputMatrixDebugHBFormat(std::ostream &os) const
{
  CSR_Matrix_Storage_t storage_initial = storage;
  PutInHarwellBoeing();
  int Nnz = ReturnNumberOfNonZero(); //number of non null elem in matrix
  os << NbLines << " " << Nnz << std::endl;
  double *a = (double *)a_->array;
  int *ai   = (int *)ai_->array;
  int *jptr = (int *)jptr_->array;
  int *ptr  = (int *)ptr_->array;

  os << "a" << endl;
  for (int i=0; i < Nnz; i++){
    os << std::scientific << std::showpoint << std::setprecision(17)  <<  a[i] << " " ;
  }
  os << std::endl;
  os << std::endl;
  os << "row index for each value in a (ai)" << endl;
  for (int i=0; i < Nnz; i++){
    os << ai[i] << " " ;
  }
  os << std::endl;
  os << std::endl;
  
  os << " index of the beginning of each column (jptr)" << endl;
  for (int i=0; i < NbLines; i++){
    os << jptr[i] << " " ;
  }
  os << std::endl;
  os << std::endl;
  os << "column indices (ptr)" << endl;
  for (int i=0; i < Nnz; i++){
    os << ptr[i] << " " ;
  }
  os << std::endl;

  if (storage_initial != HARWELL_BOEING_STORAGE)   PutInCompressedRow();

}

void xCSRMatrix::OutputMatrixDebugCompressedRowFormat(std::ostream &os) const
{
  CSR_Matrix_Storage_t storage_initial = storage;
  PutInCompressedRow();
  int Nnz = ReturnNumberOfNonZero(); //number of non null elem in matrix
  os << NbLines << " " << Nnz << std::endl;
  double *a = (double *)a_->array;
  int *ai   = (int *)ai_->array;
  int *jptr = (int *)jptr_->array;
  int *ptr  = (int *)ptr_->array;

  os << "a" << endl;
  for (int i=0; i < Nnz; i++){
    os << std::scientific << std::showpoint << std::setprecision(17)  <<  a[i] << " " ;
  }
  os << std::endl;
  os << std::endl;
  os << "row index for each value in a (ai)" << endl;
  for (int i=0; i < Nnz; i++){
    os << ai[i] << " " ;
  }
  os << std::endl;
  os << std::endl;
  
  os << " index of the beginning of each line (jptr)" << endl;
  for (int i=0; i < NbLines+1; i++){
    os << jptr[i] << " " ;
  }
  os << std::endl;
  os << std::endl;
  os << "column indices (ptr)" << endl;
  for (int i=0; i < Nnz; i++){
    os << ptr[i] << " " ;
  }
  os << std::endl;

  if (storage_initial != COMPRESSED_ROW_STORAGE)   PutInHarwellBoeing();

}
  
xCSRMatrix::xCSRMatrix(const int n) : 
a_(nullptr), ai_(nullptr), ptr_(nullptr), jptr_(nullptr), 
swap_nnz_(nullptr), swap_nbup1_(nullptr),
NbLines(0), AllocationDone(false), 
storage(UNKNOWN_STORAGE),
MatrixWasReordered(false), 
permr_(nullptr), permp_(nullptr), rpermr_(nullptr)
{
  Allocate(n);
  return ;
  Hard_CRS_Needed = true;
}

void xCSRMatrix::SetSize(const int& n) {

  if (AllocationDone) { 
    fprintf(stderr, "Error: Matrix already allocated\n");
    assert(0); 
  }
  Allocate(n); 
  Free(swap_nnz_); swap_nnz_=nullptr;
  Free(swap_nbup1_); swap_nbup1_=nullptr;
  Hard_CRS_Needed = true;
  return ;  
}

void xCSRMatrix::Clear(){

    NbLines = 0;
    storage = UNKNOWN_STORAGE;
    AllocationDone = false;
    MatrixWasReordered = false;
    Hard_CRS_Needed = true;
    List_Delete(a_);  a_ = nullptr;
    List_Delete(ai_); ai_ = nullptr;
    List_Delete(ptr_); ptr_ = nullptr;
    List_Delete(jptr_); jptr_ = nullptr;
    delete []rpermr_; rpermr_ = nullptr;
    delete []permr_; permr_ = nullptr;
    delete []permp_; permp_ = nullptr;
    Free(swap_nnz_); swap_nnz_=nullptr;
    Free(swap_nbup1_); swap_nbup1_=nullptr;
}


xCSRMatrix::~xCSRMatrix(){
    List_Delete(a_);
    List_Delete(ai_); 
    List_Delete(ptr_);
    List_Delete(jptr_); 
    if ( MatrixWasReordered == true ) {
      delete []rpermr_;
      delete []permr_;
      delete []permp_; 
      Free(swap_nnz_); swap_nnz_=nullptr;
      Free(swap_nbup1_); swap_nbup1_=nullptr;
    }
}

int xCSRMatrix :: GetNbUnknown() const { return NbLines ; }


void  xCSRMatrix :: SoftZeroMatrix ( )
{
  double* v = ReturnValPointer();
  int n = List_Nbr(a_);
  for (int i=0; i< n; i++) {*v = 0.0; v++; }
}

// Allows the summation of two matrices having the same graph
void  xCSRMatrix :: SoftAddMatrix (const xCSRMatrix& other)
{
  double* v = ReturnValPointer();
  const double* v_other = other.ReturnValPointer();
  int n = List_Nbr(a_);
  for (int i=0; i< n; i++) {*v += *v_other; v++; v_other++;}
}


void  xCSRMatrix :: ZeroMatrix ( )
{
//warning I believe there is a bug in this function
  List_Reset (a_);
  List_Reset (ai_);
  List_Reset (ptr_);
  List_Reset (jptr_);
  int i,j=0;
  for (i=0; i< NbLines; i++) List_Add ( jptr_, &j);
  Hard_CRS_Needed = true; 
  Free(swap_nnz_); swap_nnz_=nullptr;
  Free(swap_nbup1_); swap_nbup1_=nullptr;
}

void xCSRMatrix::AddMatrix ( const  xCSRMatrix& in, const double& coeff)
{
  int i,j;
  const int Asize = NbLines;
  double val;
  for(i=1;i<=Asize;i++){
    for(j=1;j<=Asize;j++){
      val = in.GetMatrix(i,j) ;
      if (val != 0.0)  AddMatrix(i, j, coeff * val);
    }
  }
}


void  xCSRMatrix::AddMatrix ( int ic, int il, double val) {

   assert( MatrixWasReordered == false );

  // attention a la transposition 

  int     *ai, *pp, n, iptr, iptr2, jptr, *ptr, zero = 0;
  double   *a;

  if ( storage == COMPRESSED_ROW_STORAGE ) PutInHarwellBoeing();

  //printf( "A(%d,%d) = %12.8e Number of Non Zero %d \n",ic,il,val,List_Nbr(a_));

   il--;
   pp  = (int*) jptr_->array;
   ptr = (int*) ptr_->array;
   ai  = (int*) ai_->array;
   a   = (double*) a_->array;
   
   iptr = pp[il];
   iptr2 = iptr-1;
   
   while(iptr>0){
     iptr2 = iptr-1;
     jptr = ai[iptr2];
     if(jptr == ic){
       a[iptr2] += val;
       return;
     }
     iptr = ptr[iptr2];
   }
   Hard_CRS_Needed = true;
   List_Add (a_, &val);
   List_Add (ai_, &ic);
   List_Add (ptr_, &zero);
   
   // The pointers amy have been modified
   //   if there has been a reallocation in List_Add  
   
   ptr = (int*) ptr_->array;
   ai  = (int*) ai_->array;
   a   = (double*) a_->array;
   
   n = List_Nbr(a_);

   if(!pp[il]) pp[il] = n;
   else ptr[iptr2] = n;
   
}
  
double xCSRMatrix::GetMatrix ( int ic, int il) const {
  if ( storage == HARWELL_BOEING_STORAGE && MatrixWasReordered == false) { 
    double val;
    int     *ai, *pp, iptr, iptr2, jptr, *ptr; 
    double   *a;
    il--;
    pp  = (int*) jptr_->array;
    ptr = (int*) ptr_->array;
    ai  = (int*) ai_->array;
    a   = (double*) a_->array;
    iptr = pp[il];
    iptr2 = iptr-1;
    
    while(iptr>0){
      iptr2 = iptr-1;
      jptr = ai[iptr2];
      if(jptr == ic){
	assert(a != nullptr);
	val = a[iptr2];
	return val;
      }
      iptr = ptr[iptr2];
    }
    return 0.0;
  }
  else if ( storage == COMPRESSED_ROW_STORAGE && MatrixWasReordered == false) { 
    int     *row_ptr, *col_ptr;
    double   *val_ptr;
    val_ptr    = (double*) a_->array;
    row_ptr    = (int*)    jptr_->array;
    col_ptr    = (int*)    ptr_->array;
    
    for (int t=row_ptr[ic-1]-1; t<row_ptr[ic]-1; t++)
      if (col_ptr[t] == il ) return val_ptr[t];
    
    if ( 1 <=ic && ic <= NbLines && 1 <=il && il <= NbLines ) return 0.0;
    else {
      printf( "Array accessing exception -- out of bounds.");
      assert( 0 );
      throw;
    }
  }
  else if ( storage == COMPRESSED_ROW_STORAGE && MatrixWasReordered == true ) { 
    int     *row_ptr, *col_ptr;
    double   *val_ptr;
    val_ptr    = (double*) a_->array;
    row_ptr    = (int*)    jptr_->array;
    col_ptr    = (int*)    ptr_->array;
    for (int t=row_ptr[ic-1]-1; t<row_ptr[ic]-1; t++)
      if (col_ptr[t] == il ) return val_ptr[t];
    if ( 1 <=ic && ic <= NbLines && 1 <=il && il <= NbLines ) return 0.0;
    else{
      printf( "Array accessing exception -- out of bounds.");
      assert( 0 );
      throw;
    }
    
  }
  else { 
    printf( "CSR_Matrix : Matrix storage is corrupted. \n"); 
    assert (0); 
    throw;
  }
}


int xCSRMatrix::ReturnNumberOfNonZero ( ) const { 
   return  List_Nbr(a_);
}

const double * xCSRMatrix::ReturnValPointer ( ) const
{   
  return  (double*) a_->array;
}
double * xCSRMatrix::ReturnValPointer ( ) 
{   
  return  (double*) a_->array;
}

const int * xCSRMatrix::ReturnRowPointer  ( ) const
{   
  return  (int*) jptr_->array;  
}

int * xCSRMatrix::ReturnRowPointer  ( ) 
{ 
  PutInCompressedRow();
  return  (int*) jptr_->array;  
}

const int * xCSRMatrix::ReturnColumIndexPointer( ) const 
{ 
  PutInCompressedRow();
  return (int*)  ptr_->array;  
}
int * xCSRMatrix::ReturnColumIndexPointer( )  { return (int*)  ptr_->array;  }

void  xCSRMatrix :: EndOfAssembly() {

  PutInCompressedRow();

 }



//return a vector of NbLines+1 value giving the index in a 
//of the start of each row.
int * xCSRMatrix::ReturnRowPointerShiftedToZero  ( )  const
{   
  const bool debug = false;
  PutInCompressedRow();
  row_shifted.resize(NbLines+1);

  int *jptr = (int *)jptr_->array;
  transform(jptr, jptr+NbLines+1, row_shifted.begin(), std::bind2nd(std::plus<int>(), -1));


  if (debug)
    { cout << " row index shifted " << endl;
      std::copy(row_shifted.begin(), row_shifted.end(), std::ostream_iterator<int>(std::cout, " ")); 
      cout << endl;
    }
  
  return  &row_shifted[0];
}


int * xCSRMatrix::ReturnColumnIndexPointerShiftedToZero  ( )  const
{   
  const bool debug = false;
  PutInCompressedRow();

  int nbz = List_Nbr(a_);
  column_shifted.resize(nbz);
  int *ptr  = (int *)ptr_->array;
  transform(ptr, ptr+nbz, column_shifted.begin(), std::bind2nd(std::plus<int>(), -1));


  if (debug) 
    {
      cout << " column shifted " << endl;
      std::copy(column_shifted.begin(), column_shifted.end(), std::ostream_iterator<int>(std::cout, " ")); 
      cout << endl;
    }
  return  &column_shifted[0];
}




void xCSRMatrix::PutInCompressedRow()  const {
//void csr_format in sort.c
  
 if (storage == HARWELL_BOEING_STORAGE) {
   if (Hard_CRS_Needed){
     storage = COMPRESSED_ROW_STORAGE;
     //     std::cout << "### CRHARD" << std::endl;
     int    i,*ptr,*jptr,*ai,iptr,iptr2;
     double *a;
     
     ptr  = (int*)ptr_->array;
     jptr = (int*)jptr_->array;
     ai   = (int*)ai_->array;
     a    = (double*)a_->array;
     
     
     
     for(i=0;i<NbLines;i++){
       iptr = jptr[i];
       while(iptr){
	 iptr2 = iptr - 1;
	 iptr = ptr[iptr2];
	 ptr[iptr2] = i+1;
       }
     }
     
     sort2(List_Nbr(a_),a,ai,ptr);
     
     deblign(List_Nbr(a_),ptr,jptr,ai);
     
     jptr[NbLines]=List_Nbr(a_)+1;



     Hard_CRS_Needed = false; 
     Hard_HB_Needed  = true;
   
   }
   else PutInCompressedRowSoft();
 }
 
 return;
}


void xCSRMatrix::PutInCompressedRowSoft() const {
  // note : this only work if the matrix was initially in HarwellBoeing ans symetry ???
 if (Hard_CRS_Needed) {
   std::cout << "trying a PutInCompressedRowSoft while the flag Hard_CRS_needed is on" <<std::endl;
   throw ;
 }
 if (storage == HARWELL_BOEING_STORAGE) {
   //   std::cout << "### CRSOFT" << std::endl;
   storage = COMPRESSED_ROW_STORAGE;
   char* temp = ptr_->array;
   ptr_->array = swap_nnz_->array;
   swap_nnz_->array = temp;
   
   temp = jptr_->array;
   jptr_->array =  swap_nbup1_->array;
   swap_nbup1_->array = temp ;
 }
 return;
}

void xCSRMatrix::PutInHarwellBoeingSoft() const {
  if (storage == COMPRESSED_ROW_STORAGE) {
    storage = HARWELL_BOEING_STORAGE; 
    //    std::cout << "### HBSOFT" << std::endl;
    char* temp = ptr_->array;
    ptr_->array = swap_nnz_->array;
    swap_nnz_->array = temp;
  
    temp = jptr_->array;
    jptr_->array =  swap_nbup1_->array;
    swap_nbup1_->array = temp ;
  }
}

void xCSRMatrix::PutInHarwellBoeing() const {
  // note : this only work if the matrix was initially in HarwellBoeing
  // And if theMatrix is symmetric !!!
  if (!Hard_HB_Needed){ 
    PutInHarwellBoeingSoft();
    return;
  }
  
  if (storage == COMPRESSED_ROW_STORAGE) {
    storage = HARWELL_BOEING_STORAGE;
    //    std::cout << "### HBHARD" << std::endl;
   //symetric case
   /*{
     char* temp  = ptr_->array;
     ptr_->array = ai_->array;
     ai_->array = temp;
   }*/
   
   int nnz = ReturnNumberOfNonZero();
   int *ja = (int*)ptr_->array;
   int *ia = (int*)jptr_->array;
   
   Free (swap_nnz_);
   Free (swap_nbup1_);
   swap_nnz_ = nullptr;
   swap_nbup1_ = nullptr;
   
   swap_nnz_   = List_Create (nnz, NbLines, sizeof(int));
   swap_nbup1_ = List_Create (NbLines+1, NbLines, sizeof(int));
   //   List_T *lai   = List_Create (nnz, NbLines, sizeof(int));
   
   int *ptr  = (int*)swap_nnz_->array;
   int *jptr = (int*)swap_nbup1_->array;
   int *ai   = (int*)ai_->array;
   
   for (int i= 0; i< NbLines+1; i++) jptr[i] =0;
   for (int k = 0; k< nnz; k++) ptr[k] = 0;
   int i = 1;
   for (int k = 0; k< nnz; k++) {
     if (k >= ia[i]-1) ++i;
     ai[k]=i;
     int j = ja[k];
     int ptr1 = jptr[j-1];
     if (!ptr1) jptr[j-1] = k+1;
     else {
       int ptr2 = ptr[ptr1-1];
       while (ptr2) {
	 ptr1 = ptr2;
	 ptr2 = ptr[ptr1-1];
       }
       ptr[ptr1-1] = k+1;
     }
   }
   
   char *tmp = ptr_->array;
   ptr_->array  = (char* )ptr;
   swap_nnz_->array  = tmp;
   
   tmp = jptr_->array; 
   jptr_->array = (char *)jptr;
   swap_nbup1_->array = tmp;
   Hard_HB_Needed  = false;
 }
 return;
}


  
/// an iterator on the non zero element of a xCSRMatrix
double & xCSRMatrix::iterator::operator *(){
  double * a = (double*)(themat->a_->array);
  return a[ncurrent];
}
    
xCSRMatrix::iterator & xCSRMatrix::iterator::operator++(){
  ++ncurrent;
  return *this;
}
    
bool xCSRMatrix::iterator::operator!=(const iterator& rhs) const{
  return rhs.ncurrent!=ncurrent;
}
  
xCSRMatrix::iterator xCSRMatrix::begin(){
    iterator it;
    it.themat = this;
    it.ncurrent = 0;
    return it;
}
  
xCSRMatrix::iterator xCSRMatrix::end(){
  iterator it;
  it.themat = this;
  it.ncurrent = it.themat->ReturnNumberOfNonZero();
  return it;
}
  

double & xCSRMatrix::column_iterator_HB::operator *(){
  double * a = (double*)(themat->a_->array);
  return a[ncurrent-1];
}

xCSRMatrix::column_iterator_HB & xCSRMatrix::column_iterator_HB::operator++(){
      ncurrent = ((int*)(themat->ptr_->array))[ncurrent-1];
      line =((int*)( themat->ai_->array))[ncurrent-1];
      return *this;
}

bool xCSRMatrix::column_iterator_HB::operator!=(const column_iterator_HB& rhs) const{
  return rhs.ncurrent!=ncurrent;
}

// xfem::xCSRMatrix::column_iterator_HB::operator!=(xfem::xCSRMatrix::column_iterator_HB const&) const'



xCSRMatrix::column_iterator_HB xCSRMatrix::begin_column_HB(int j){ 
  // the name is begin_column and not begin, since over loading on return argument only is not possible
  xCSRMatrix::column_iterator_HB it;
  it.column=j;
  it.themat = this;
  this->PutInHarwellBoeing();
  it.ncurrent = ((int*)(it.themat->jptr_->array))[j-1];
  it.line =((int*)( it.themat->ai_->array))[it.ncurrent-1];
  return it;
}
  
xCSRMatrix::column_iterator_HB xCSRMatrix::end_column_HB(int j){
  xCSRMatrix::column_iterator_HB it=this->begin_column_HB(j);
  while (it.ncurrent) ++it;
  return it;
} 




/*
void xCSRMatrix::PutInHarwellBoeingHard() const {
{ 
  if (storage == COMPRESSED_ROW_STORAGE){
    storage = HARWELL_BOEING_STORAGE;
    //a;
    //CSR index
    int nnz = ;
    int * csr_ja  = // ReturnColumIndexPointer() ; // size nnz //allocate a copy here
    int * csr_ia  = ReturnRowPointer() ;   //size n+1
    int * hb_ptr  =                        //size nnz
    int * hb_jptr = csr_ia ;               //size n+1 // need to be copied first
    int * hb_ai   = csr_ja ;               //size nnz //ok
    
    int k = 1 ;
    line =1 ;
    while (k<=nnz){
    if (csr_ia[line-1]) <k) ++line;
    int col = csr_ja[k-1];
    hb_ai[k-1] = line;
    int nextincol = hb_jptr[col-1];
    if (!nextincol) {
	hb_jptr[col-1] = k;
      }
      else{
	while (!(nextincol = hb_ptr[nextincol-1])){
	}
	hb_ptr[nextincol-1]=k;
	} 
    }
  }
}
*/


void  xCSRMatrix::sort2(unsigned long n, double arr[], int ai[] , int aj [] )  const {

  unsigned long i,ir=n,j,k,l=1;
  int *istack,jstack=0,tempi;
  double a,temp;
  int    b,c;
    
  istack=ivector(1,NSTACK);
  for (;;) {
    if (ir-l < M_sort2) {
      for (j=l+1;j<=ir;j++) {
	a=arr[j -1];
	b=ai[j -1];
	c=aj[j -1];
	for (i=j-1;i>=1;i--) {
	  if (cmpij(ai[i -1],aj[i -1],b,c) <= 0) break;
	  arr[i+1 -1]=arr[i -1];
	  ai[i+1 -1]=ai[i -1];
	  aj[i+1 -1]=aj[i -1];
	}
	arr[i+1 -1]=a;
	ai[i+1 -1]=b;
	aj[i+1 -1]=c;
      }
      if (!jstack) {
	free_ivector(istack,1,NSTACK);
	return;
      }
      ir=istack[jstack];
      l=istack[jstack-1];
      jstack -= 2;
    } 
    else {
      k=(l+ir) >> 1;
      SWAP(arr[k -1],arr[l+1 -1])
      SWAPI(ai[k -1],ai[l+1 -1])
      SWAPI(aj[k -1],aj[l+1 -1])
      if (cmpij(ai[l+1 -1],aj[l+1 -1],ai[ir -1],aj[ir -1])>0){
	SWAP(arr[l+1 -1],arr[ir -1])
	SWAPI(ai[l+1 -1],ai[ir -1])
	SWAPI(aj[l+1 -1],aj[ir -1])
      }
      if (cmpij(ai[l -1],aj[l -1],ai[ir -1],aj[ir -1])>0){
	SWAP(arr[l -1],arr[ir -1])
	SWAPI(ai[l -1],ai[ir -1])
	SWAPI(aj[l -1],aj[ir -1])
      }
      if (cmpij(ai[l+1 -1],aj[l+1 -1],ai[l -1],aj[l -1])>0){
	SWAP(arr[l+1 -1],arr[l -1])
	SWAPI(ai[l+1 -1],ai[l -1])
	SWAPI(aj[l+1 -1],aj[l -1])
      }
      i=l+1;
      j=ir;
      a=arr[l -1];
      b=ai[l -1];
      c=aj[l -1];
      for (;;) {
	do i++; while (cmpij(ai[i -1],aj[i -1],b,c) < 0);
	do j--; while (cmpij(ai[j -1],aj[j -1],b,c) > 0);
	if (j < i) break;
	SWAP(arr[i -1],arr[j -1])
	SWAPI(ai[i -1],ai[j -1])
	SWAPI(aj[i -1],aj[j -1])
	}
      arr[l -1]=arr[j -1];
      arr[j -1]=a;
      ai[l -1]=ai[j -1];
      ai[j -1]=b;
      aj[l -1]=aj[j -1];
      aj[j -1]=c;
      jstack += 2;
      if (jstack > NSTACK) {
	fprintf(stderr,"NSTACK too small in sort2.\n");
	exit(1);
      }
      if (ir-i+1 >= j-l) {
	istack[jstack]=ir;
	istack[jstack-1]=i;
	ir=j-1;
      } 
      else {
	istack[jstack]=j-1;
	istack[jstack-1]=l;
	l=i;
      }
    }
  }
}


int  xCSRMatrix::cmpij(int ai,int aj,int bi,int bj) const {
  if(ai<bi)return -1;
  if(ai>bi)return 1;
  if(aj<bj)return -1;
  if(aj>bj)return 1;
  return 0;
}

int * xCSRMatrix::ivector(long nl, long nh)  const {
  // allocate an int vector with subscript range v[nl..nh] 
  int *v;
  
  v=(int *)malloc((size_t) ((nh-nl+2)*sizeof(int)));
  if (!v) fprintf(stderr, "allocation failure in ivector()\n");
  return v-nl+1;
}

 void xCSRMatrix :: free_ivector(int *v, long nl, long nh) const {
  // free an int vector allocated with ivector() 
  free((char*)(v+nl-1));
}


void xCSRMatrix::ExecuteReordering() {

  //  const bool debug = xdebug_flag;
  const bool debug = false;
  if ( MatrixWasReordered == false ){

  PutInCompressedRow();

  double *a;
  int    *jptr,*ai;
  int    *mask, *levels;

  int i,j,k;
  a    = (double*) a_->array;
  jptr = (int*)    jptr_->array;
  ai   = (int*)    ptr_->array; 

//  permr_  = (int*) Malloc (NbLines * sizeof(int));
//  rpermr_ = (int*) Malloc (NbLines * sizeof(int));
//  permp_  = (int*) Malloc (2 * NbLines * sizeof(int));
  permr_  = new int[NbLines];
  rpermr_ = new int[NbLines];
  permp_  = new int[2 * NbLines];

  for (i=0 ; i<NbLines ; i++) {
    permr_[i] = rpermr_[i] = permp_[i+NbLines] = i+1;
  }

  //  fprintf(stdout, "# --> RCMK renumbering\n");

  double *a_rcmk;
  int    *jptr_rcmk, *ai_rcmk;

  a_rcmk    = (double*) Malloc (List_Nbr(a_) * sizeof(double));
  jptr_rcmk = (int*)    Malloc ((NbLines + 1) * sizeof(int));
  ai_rcmk   = (int*)    Malloc (List_Nbr(a_) * sizeof(int));

  mask   = (int*) Malloc (List_Nbr(a_) * sizeof(int));
  levels = (int*) Malloc ((NbLines + 1) * sizeof(int));
  
  i = j = k = 1;
  
  cmkreord  (NbLines, a, ai, jptr, a_rcmk, ai_rcmk, jptr_rcmk, &i,
	     permr_, mask, &j, &k, rpermr_, levels);

  double *vv;
  vv = (double*) Malloc (List_Nbr(a_) * sizeof(double));

  sort_col  ( NbLines, a_rcmk, ai_rcmk, jptr_rcmk, mask, vv);
  
  Free(vv);
  Free(mask);
  Free(levels);

  for ( i=0 ; i < List_Nbr(a_) ; i++ ){  
   a[i]    = a_rcmk[i];
   ai[i]   = ai_rcmk[i];   
 }

  for ( i=0 ; i < NbLines+1 ; i++ ) jptr[i] = jptr_rcmk[i];

    Free (a_rcmk);
    Free (jptr_rcmk);
    Free (ai_rcmk);

    //  printf("The number of entries in the matrix is %d\n", List_Nbr(a_) );

    
   if (debug) {
     printf("CSR::The row pointer is\n");
     for (i = 0; i < NbLines + 1; i++) {
       printf("%d ", jptr[i]);
     }
     printf("\n");
     printf("The col pointer and value are\n");
     for (i = 0; i < List_Nbr(a_); i++) {
       printf("%d %12.5e\n", ai[i], a[i]);
     }
   }
   
   MatrixWasReordered = true ;
  }
  
}

void xCSRMatrix::ReorderArray(double * array) {

 if ( MatrixWasReordered  == true ){
 
 int i;
 double * temp;
 temp = new double [ NbLines ];

  for ( i=0;i<NbLines;i++) temp[i] = array[rpermr_[i] - 1];
 
  for ( i=0;i<NbLines;i++) array[i] = temp[i];
  delete [] temp;
 }

}

void xCSRMatrix::InverseReorderArray(double * array) {


 if ( MatrixWasReordered == true ){

 double * temp;
 temp = new double [ NbLines ];
 
 int i,j,k;

   for ( i=0;i<NbLines;i++) {
    j = permr_[i] - 1;
    k = permp_[j+NbLines] - 1;        
    temp[i] = array[k];
   }  

    for ( i=0;i<NbLines;i++) array[i] = temp[i];
 
 delete [] temp;

 }
}

void xCSRMatrix::cmkreord(int n, double *a, int *ja, int *ia, double *a0,
			   int *ja0, int *ia0, int *init, int * iperm, int * mask,
			   int * maskval, int * nlev, int * riord, int * levels) {

    int c__1 = 1;
    // System generated locals 
    int i__1;

    // Local variables 
//    extern // Subroutine  int exchange_();
    int i;
//    extern // Subroutine  int dperm_(), perphn_(), rversp_();

// -----------------------------------------------------------------------
 
//      Cuthill-McKee Reordering : call to perphn 
//                                 renumber the nodes 
// -------------------------parameters------------------------------------
 
// on entry: 
// ---------- 
// n      = number of nodes in the graph 
// ja, ia = pattern of matrix in CSR format (the ja,ia arrays of csr data 

//          structure) 
// init   = initial node 
// iperm  = integer array indicating in which order to  traverse the graph
 
//          in order to generate all connected components. 
//          The nodes will be traversed in order iperm(1),....,iperm(n) 
//          Convention: 
//          if iperm(1) .eq. 0 on entry then BFS will traverse the 
//          nodes in the  order 1,2,...,n. 

// riord  = (also an ouput argument). on entry riord contains the labels 

//          of the nfirst nodes that constitute the first level. 

// mask   = array used to indicate whether or not a node should be 
//          condidered in the graph. see maskval. 
//          mask is also used as a marker of  visited nodes. 

// maskval= consider node i only when:  mask(i) .eq. maskval 
//          maskval must be .gt. 0. 
//          thus, to consider all nodes, take mask(1:n) = 1. 
//          maskval=1 (for example) 

// on return 
// --------- 
// mask   = on return mask is restored to its initial state. 
// riord  = `reverse permutation array'. Contains the labels of the nodes 

//          constituting all the levels found, from the first level to 
//          the last. 
// levels = pointer array for the level structure. If lev is a level 
//          number, and k1=levels(lev),k2=levels(lev+1)-1, then 
//          all the nodes of level number lev are: 
//          riord(k1),riord(k1+1),...,riord(k2) 
// nlev   = number of levels found 
// -----------------------------------------------------------------------
 
//    local variables 
//     print*,'   ... Searching levels' 
    // Parameter adjustments 
    --levels;
    --riord;
    --mask;
    --iperm;
    --ia0;
    --ja0;
    --a0;
    --ia;
    --ja;
    --a;

    // Function Body 
    *maskval = 1;
    i__1 = n;
    for (i = 1; i <= i__1; ++i) {
	mask[i] = *maskval;
    }
    iperm[1] = 0;
    perphn(n, &ja[1], &ia[1], init, &iperm[1], &mask[1], maskval, nlev, &
	    riord[1], &levels[1]);
    rversp(n, &riord[1]);
//     print*,'   ... Renumbering' 
    exchange(n, &riord[1], &iperm[1]);
    dperm(n, &a[1], &ja[1], &ia[1], &a0[1], &ja0[1], &ia0[1], &iperm[1], &
	    iperm[1], &c__1);

} 

int xCSRMatrix::perphn(int n, int *ja, int *ia, int *init, int *iperm,
			  int * mask, int * maskval, int * nlev, 
			  int * riord, int * levels){
    // System generated locals 
    int i__1;

    // Local variables 
    int j, nlevp, mindeg, nfirst, deg;
//    extern // Subroutine  int bfs_();
    int nod;
//    extern integer maskdeg_();

// -----------------------------------------------------------------------
 
//     finds a pseudo-peripheral node and does a BFS search from it. 
// -----------------------------------------------------------------------
 
// see routine  dblstr for description of parameters 
// input: 
// ------- 
// ja, ia  = list pointer array for the adjacency graph 
// mask    = array used for masking nodes -- see maskval 
// maskval = value to be checked against for determing whether or 
//           not a node is masked. If mask(k) .ne. maskval then 
//           node k is not considered. 
// init    = init node in the pseudo-peripheral node algorithm. 

// output: 
// ------- 
// init    = actual pseudo-peripherial node found. 
// nlev    = number of levels in the final BFS traversal. 
// riord   = 
// levels  = 
// -----------------------------------------------------------------------
 
    // Parameter adjustments 
    --levels;
    --riord;
    --mask;
    --iperm;
    --ia;
    --ja;

    // Function Body 
    nlevp = 0;
L1:
    riord[1] = *init;
    nfirst = 1;
    bfs(n, &ja[1], &ia[1], &nfirst, &iperm[1], &mask[1], maskval, &riord[1], 
	    &levels[1], &(*nlev));

    if (*nlev > nlevp) {
	mindeg = levels[*nlev + 1] - 1;
	i__1 = levels[*nlev + 1] - 1;
	for (j = levels[*nlev]; j <= i__1; ++j) {
	    nod = riord[j];
	    deg = maskdeg(&ja[1], &ia[1], &nod, &mask[1], maskval);
	    if (deg < mindeg) {
		*init = nod;
		mindeg = deg;
	    }
	}
	nlevp = *nlev;
	goto L1;
    }
    return 0;
} 

int  xCSRMatrix::rversp(int n, int *riord) {

    // System generated locals 
    int i__1;

    // Local variables 
    int j, k;

// -----------------------------------------------------------------------
//
//     this routine does an in-place reversing of the permutation array 
//     riord -- 
// -----------------------------------------------------------------------

    // Parameter adjustments //
    --riord;

    // Function Body
    i__1 = n / 2;
    for (j = 1; j <= i__1; ++j) {
	k = riord[j];
	riord[j] = riord[n - j + 1];
	riord[n - j + 1] = k;
// L26: 
    }
    return 0;
}


int xCSRMatrix::maskdeg(int *ja, int *ia, int *nod,
			  int *mask, int *maskval) {

    // System generated locals 
    int ret_val, i__1;

    // Local variables 
    int k, deg;

// ----------------------------------------------------------------------- 
    // Parameter adjustments 
    --mask;
    --ia;
    --ja;

    // Function Body 
    deg = 0;
    i__1 = ia[*nod + 1] - 1;
    for (k = ia[*nod]; k <= i__1; ++k) {
	if (mask[ja[k]] == *maskval) {
	    ++deg;
	}
    }
    ret_val = deg;
    return ret_val;
}

void xCSRMatrix::exchange(int n, int *iriord, int *iperm) {
    // System generated locals 
    int i__1;

    // Local variables 
    int i;

//-----------------------------------------------------------------------
//  Reverse a permutation vector 
//  On entry : 
//  ---------- 
//            n      : dimension 
//            iriord : initial reordering vector 
//  On return : 
//  ----------- 
//           iperm  : permutation vector to be used with SPARSKIT (dperm .
//-----------------------------------------------------------------------
//     local variable 
    // Parameter adjustments 
    --iperm;
    --iriord;

    // Function Body 
    i__1 = n;
    for (i = 1; i <= i__1; ++i) {
	iperm[iriord[i]] = i;
    }
} 

void xCSRMatrix::sort_col(int n,double * a, int * ja, int *ia,
			    int * iw, double * rw){
    // System generated locals 
    int i__1, i__2;

    // Local variables 
    int ideb, ifin;
//    extern  int sort_irv__();
    int i, j, k;

    // Parameter adjustments 
    --rw;
    --iw;
    --ia;
    --ja;
    --a;

    // Function Body 
    i__1 = n;
    for (i = 1; i <= i__1; ++i) {
	ideb = ia[i];
	ifin = ia[i + 1] - 1;
	k = 0;
	i__2 = ifin;
	for (j = ideb; j <= i__2; ++j) {
	    ++k;
	    iw[k] = ja[j];
	    rw[k] = a[j];
	}
	i__2 = ifin - ideb + 1;
	sort_irv(&iw[1], &rw[1], &i__2);
	k = 0;
	i__2 = ifin;
	for (j = ideb; j <= i__2; ++j) {
	    ++k;
	    ja[j] = iw[k];
	    a[j] = rw[k];
	}
    }
    return;
}


int xCSRMatrix::dperm(int nrow, double *a, int *ja, int * ia, double * ao,
			int * jao, int * iao, int * perm, int *qperm, int * job) {

    int locjob;

// -----------------------------------------------------------------------
 
// This routine permutes the rows and columns of a matrix stored in CSR 
// format. i.e., it computes P A Q, where P, Q are permutation matrices. 

// P maps row i into row perm(i) and Q maps column j into column qperm(j):
 
//      a(i,j)    becomes   a(perm(i),qperm(j)) in new matrix 
// In the particular case where Q is the transpose of P (symmetric 
// permutation of A) then qperm is not needed. 
// note that qperm should be of length ncol (number of columns) but this 

// is not checked. 
// -----------------------------------------------------------------------
 
// Y. Saad, Sep. 21 1989 / recoded Jan. 28 1991. 
// -----------------------------------------------------------------------
 
// on entry: 
// ---------- 
// n 	= dimension of the matrix 
// a, ja, 
//    ia = input matrix in a, ja, ia format 
// perm 	= integer array of length n containing the permutation arrays 
// 	  for the rows: perm(i) is the destination of row i in the 
//         permuted matrix -- also the destination of column i in case 
//         permutation is symmetric (job .le. 2) 

// qperm	= same thing for the columns. This should be provided only 
//         if job=3 or job=4, i.e., only in the case of a nonsymmetric 
// 	  permutation of rows and columns. Otherwise qperm is a dummy 

// job	= integer indicating the work to be done: 
// * job = 1,2 permutation is symmetric  Ao :== P * A * transp(P) 
// 		job = 1	permute a, ja, ia into ao, jao, iao 
// 		job = 2 permute matrix ignoring real values. 
// * job = 3,4 permutation is non-symmetric  Ao :== P * A * Q 
// 		job = 3	permute a, ja, ia into ao, jao, iao 
// 		job = 4 permute matrix ignoring real values. 

// on return: 
// ----------- 
// ao, jao, iao = input matrix in a, ja, ia format 

// in case job .eq. 2 or job .eq. 4, a and ao are never referred to 
// and can be dummy arguments. 
// Notes: 
// ------- 
//  1) algorithm is in place 
//  2) column indices may not be sorted on return even  though they may be
 
//     on entry. 
// ----------------------------------------------------------------------c
 
// local variables 

//     locjob indicates whether or not real values must be copied. 

    // Parameter adjustments 
    --qperm;
    --perm;
    --iao;
    --jao;
    --ao;
    --ia;
    --ja;
    --a;

    // Function Body 
    locjob = *job % 2;

// permute rows first 

    rperm(nrow, &a[1], &ja[1], &ia[1], &ao[1], &jao[1], &iao[1], &perm[1], &
	    locjob);

// then permute columns 

    locjob = 0;

    if (*job <= 2) {
	cperm(nrow, &ao[1], &jao[1], &iao[1], &ao[1], &jao[1], &iao[1], &
		perm[1], &locjob);
    } else {
	cperm(nrow, &ao[1], &jao[1], &iao[1], &ao[1], &jao[1], &iao[1], &
		qperm[1], &locjob);
    }

    return 0;
}

int xCSRMatrix::cperm(int nrow, double * a, int *ja, int * ia, double * ao,
			int * jao, int * iao, int * perm, int * job) {

    // System generated locals 
    int i__1;

    // Local variables 
    int i, k, nnz;

// -----------------------------------------------------------------------
 
// this subroutine permutes the columns of a matrix a, ja, ia. 
// the result id written in the output matrix  ao, jao, iao. 
// cperm computes B = A P, where  P is a permutation matrix 
// that maps column j into column perm(j), i.e., on return 
//      a(i,j) becomes a(i,perm(j)) in new matrix 
// Y. Saad, May 2, 1990 / modified Jan. 28, 1991. 
// -----------------------------------------------------------------------
 
// on entry: 
// ---------- 
// nrow 	= row dimension of the matrix 

// a, ja, ia = input matrix in csr format. 

// perm	= integer array of length ncol (number of columns of A 
//         containing the permutation array  the columns: 
//         a(i,j) in the original matrix becomes a(i,perm(j)) 
//         in the output matrix. 

// job	= integer indicating the work to be done: 
// 		job = 1	permute a, ja, ia into ao, jao, iao 
//                       (including the copying of real values ao and 
//                       the array iao). 
// 		job .ne. 1 :  ignore real values ao and ignore iao. 

// ------------ 
// on return: 
// ------------ 
// ao, jao, iao = input matrix in a, ja, ia format (array ao not needed) 


// Notes: 
// ------- 
// 1. if job=1 then ao, iao are not used. 
// 2. This routine is in place: ja, jao can be the same. 
// 3. If the matrix is initially sorted (by increasing column number) 
//    then ao,jao,iao  may not be on return. 

// ----------------------------------------------------------------------c
 
// local parameters: 

    // Parameter adjustments 
    --perm;
    --iao;
    --jao;
    --ao;
    --ia;
    --ja;
    --a;

    // Function Body 
    nnz = ia[nrow + 1] - 1;
    i__1 = nnz;
    for (k = 1; k <= i__1; ++k) {
	jao[k] = perm[ja[k]];
// L100: 
    }

//     done with ja array. return if no need to touch values. 

    if (*job != 1) {
	return 0;
    }

// else get new pointers -- and copy values too. 

    i__1 = nrow + 1;
    for (i = 1; i <= i__1; ++i) {
	iao[i] = ia[i];
// L1: 
    }

    i__1 = nnz;
    for (k = 1; k <= i__1; ++k) {
	ao[k] = a[k];
//  L2:      
    }

    return 0;

}

int xCSRMatrix::bfs(int n, int * ja, int * ia, int * nfirst, int * iperm,
		      int * mask, int * maskval, int * riord, int * levels, int * nlev) {

    // System generated locals 
    int i__1;

    // Local variables 
    int iend;
//    extern // Subroutine  int add_lvst__();
    int j, ii, istart;

//    static logical permut;
    int permut;

    int nod;

// -----------------------------------------------------------------------
 
// finds the level-structure (breadth-first-search or CMK) ordering for a 

// given sparse matrix. Uses add_lvst. Allows a  set of nodes to be 
// the initial level (instead of just one node). Allows masked nodes. 
// -------------------------parameters------------------------------------
 
// on entry: 
// ---------- 
// n      = number of nodes in the graph 
// ja, ia = pattern of matrix in CSR format (the ja,ia arrays of csr data 

//          structure) 
// nfirst = number of nodes in the first level that is input in riord 
// iperm  = integer array indicating in which order to  traverse the graph
 
//          in order to generate all connected components. 
//          The nodes will be traversed in order iperm(1),....,iperm(n) 
//          Convention: 
//          if iperm(1) .eq. 0 on entry then BFS will traverse the 
//          nodes in the  order 1,2,...,n. 

// riord  = (also an ouput argument). on entry riord contains the labels 

//          of the nfirst nodes that constitute the first level. 

// mask   = array used to indicate whether or not a node should be 
//          condidered in the graph. see maskval. 
//          mask is also used as a marker of  visited nodes. 

// maskval= consider node i only when:  mask(i) .eq. maskval 
//          maskval must be .gt. 0. 
//          thus, to consider all nodes, take mask(1:n) = 1. 
//          maskval=1 (for example) 

// on return 
// --------- 
// mask   = on return mask is restored to its initial state. 
// riord  = `reverse permutation array'. Contains the labels of the nodes 

//          constituting all the levels found, from the first level to 
//          the last. 
// levels = pointer array for the level structure. If lev is a level 
//          number, and k1=levels(lev),k2=levels(lev+1)-1, then 
//          all the nodes of level number lev are: 
//          riord(k1),riord(k1+1),...,riord(k2) 
// nlev   = number of levels found 
// -----------------------------------------------------------------------
 
// Notes on possible usage 
// ------------------------- 
// 1. if you want a CMK ordering from a known node, say node init then 
//    call BFS with nfirst=1,iperm(1) =0, mask(1:n) =1, maskval =1, 
//    riord(1) = init. 
// 2. if you want the RCMK ordering and you have a preferred initial node 

//     then use above call followed by reversp(n,riord) 
// 3. Similarly to 1, and 2, but you know a good LEVEL SET to start from 

//    (nfirst = number if nodes in the level, riord(1:nfirst) contains 
//    the nodes. 
// 4. If you do not know how to select a good initial node in 1 and 2, 
//    then you should use perphn instead. 

// -----------------------------------------------------------------------
 
//     local variables -- 
    // Parameter adjustments 
    --levels;
    --riord;
    --mask;
    --iperm;
    --ia;
    --ja;

    // Function Body 
    permut = iperm[1] != 0;

//     start pointer structure to levels 

    *nlev = 0;

//     previous end 

    istart = 0;
    ii = 0;

//     current end 

    iend = *nfirst;

//     intialize masks to zero -- except nodes of first level -- 

    i__1 = *nfirst;
    for (j = 1; j <= i__1; ++j) {
	mask[riord[j]] = 0;
// L12: 
    }
// -----------------------------------------------------------------------
 
// L13: 

L1:
    ++(*nlev);
    levels[*nlev] = istart + 1;
    add_lvst(&istart, &iend, nlev, &riord[1], &ja[1], &ia[1], &mask[1], 
	    maskval);
    if (istart < iend) {
	goto L1;
    }
L2:
    ++ii;
    if (ii <= n) {
	nod = ii;

//	if (permut) {
        if (permut == 1) {
	    nod = iperm[nod];
	}
	if (mask[nod] == *maskval) {

//     start a new level 

	    istart = iend;
	    ++iend;
	    riord[iend] = nod;
	    mask[nod] = 0;
	    goto L1;
	} else {
	    goto L2;
	}
    }
// -----------------------------------------------------------------------
 
// L3: 
    levels[*nlev + 1] = iend + 1;
    i__1 = iend;
    for (j = 1; j <= i__1; ++j) {
	mask[riord[j]] = *maskval;
    }
    return 0;
} 

int xCSRMatrix::sort_irv(int * itmp, double *rtmp, int * n) {

    // System generated locals 
    int i__1, i__2;

    // Local variables 
    int jmin, i, j, it;
    double rt;
    int itmpmin;

//-----------------------------------------------------------------------
//   Tri selon une methode on ne peut plus simple */
//-----------------------------------------------------------------------
    // Parameter adjustments */
    --rtmp;
    --itmp;

    //      Function Body */
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	itmpmin = itmp[i];
	jmin = i;
	i__2 = *n;
	for (j = i + 1; j <= i__2; ++j) {
	    if (itmp[j] < itmpmin) {
		jmin = j;
		itmpmin = itmp[jmin];
	    }
	}
	it = itmp[i];
	itmp[i] = itmp[jmin];
	itmp[jmin] = it;
	rt = rtmp[i];
	rtmp[i] = rtmp[jmin];
	rtmp[jmin] = rt;
    }
    return 0;
} 

void xCSRMatrix::rperm(int nrow,double * a,int * ja,int * ia, double * ao, 
			int * jao, int * iao, int * perm, int * job) {
 
   // System generated locals 
    int i__1, i__2;

    // Local variables 
    int i, j, k, ii, ko;

//    static logical values;
   
    int  values;

// -----------------------------------------------------------------------
 
// this subroutine permutes the rows of a matrix in CSR format. 
// rperm  computes B = P A  where P is a permutation matrix. 
// the permutation P is defined through the array perm: for each j, 
// perm(j) represents the destination row number of row number j. 
// Youcef Saad -- recoded Jan 28, 1991. 
// -----------------------------------------------------------------------
 
// on entry: 
// ---------- 
// n 	= dimension of the matrix 
// a, ja, ia = input matrix in csr format 
// perm 	= integer array of length nrow containing the permutation arrays 

// 	  for the rows: perm(i) is the destination of row i in the 
//         permuted matrix. 
//         ---> a(i,j) in the original matrix becomes a(perm(i),j) 
//         in the output  matrix. 

// job	= integer indicating the work to be done: 
// 		job = 1	permute a, ja, ia into ao, jao, iao 
//                       (including the copying of real values ao and 
//                       the array iao). 
// 		job .ne. 1 :  ignore real values. 
//                     (in which case arrays a and ao are not needed nor 
//                      used). 

// ------------ 
// on return: 
// ------------ 
// ao, jao, iao = input matrix in a, ja, ia format 
// note : 
//        if (job.ne.1)  then the arrays a and ao are not used. 
// ----------------------------------------------------------------------c
//           Y. Saad, May  2, 1990                                      c 
// ----------------------------------------------------------------------c
    // Parameter adjustments 
    --perm;
    --iao;
    --jao;
    --ao;
    --ia;
    --ja;
    --a;

    // Function Body 
    values = *job == 1;

//     determine pointers for output matix. 

    i__1 = nrow;
    for (j = 1; j <= i__1; ++j) {
	i = perm[j];
	iao[i + 1] = ia[j + 1] - ia[j];
// L50: 
    }

// get pointers from lengths 

    iao[1] = 1;
    i__1 = nrow;
    for (j = 1; j <= i__1; ++j) {
	iao[j + 1] += iao[j];
// L51: 
    }

// copying 

    i__1 = nrow;
    for (ii = 1; ii <= i__1; ++ii) {

// old row = ii  -- new row = iperm(ii) -- ko = new pointer 

	ko = iao[perm[ii]];
	i__2 = ia[ii + 1] - 1;
	for (k = ia[ii]; k <= i__2; ++k) {
	    jao[ko] = ja[k];
	    if (values==1) {
		ao[ko] = a[k];
	    }
	    ++ko;
// L60: 
	}
// L100: 
    }

    return;

}

int xCSRMatrix::add_lvst(int * istart, int * iend, int * nlev, int * riord,
			   int * ja, int * ia, int * mask, int * maskval) {

    // System generated locals
    int i__1, i__2;

    // Local variables
    int i, j, k, ir, nod;

// ---------------------------------------------------------------------- 
// adds one level set to the previous sets. span all nodes of previous
// set. Uses Mask to mark those already visited.
// -----------------------------------------------------------------------

    // Parameter adjustments 
    --mask;
    --ia;
    --ja;
    --riord;

    // Function Body
    nod = *iend;
    i__1 = *iend;
    for (ir = *istart + 1; ir <= i__1; ++ir) {
	i = riord[ir];
	i__2 = ia[i + 1] - 1;
	for (k = ia[i]; k <= i__2; ++k) {
	    j = ja[k];
	    if (mask[j] == *maskval) {
		++nod;
		mask[j] = 0;
		riord[nod] = j;
	    }
// L24: 
	}
// L25: 
    }
    *istart = *iend;
    *iend = nod;
    return 0;
//  -----------------------------------------------------------------------
 
}

void  xCSRMatrix::deblign(int nz , int *ptr , int *jptr , int *ai) const {

  int ilign =1;
 
  jptr[0] = 1;
  for(int i=1; i<nz; i++) {
    if (ai[i-1] < ai[i]) {
      jptr[ilign++]=i+1;
      ai[i-1] = 0;
    }    
    else{
      ai[i-1] = i+1;
    }
  }
  ai[nz-1] = 0;
}

xCSRMatrix& xCSRMatrix::operator=(const xCSRMatrix& in) {

  //assert(0); cerr << "operator of xCSRMatrix not debugged yet!" << endl; 

  if (this != &in) {
    Clear();

    MatrixWasReordered = in.MatrixWasReordered;
    AllocationDone     = in.AllocationDone;
    storage            = in.storage;
    NbLines            = in.NbLines;
    Hard_CRS_Needed    = in.Hard_CRS_Needed;
    Hard_HB_Needed     =  in.Hard_HB_Needed;;
    if (in.AllocationDone) {
      a_    = List_Copy (in.a_);
      ai_   = List_Copy (in.ai_);
      ptr_  = List_Copy (in.ptr_);
      jptr_ = List_Copy (in.jptr_);
      if (in.swap_nnz_) swap_nnz_ = List_Copy (in.swap_nnz_);
      if (in.swap_nbup1_) swap_nbup1_ = List_Copy (in.swap_nbup1_);
      
    }
    if (in.MatrixWasReordered) {
      permr_  = new int[NbLines];
      rpermr_ = new int[NbLines];
      permp_  = new int[2 * NbLines];
      std::memcpy(permr_, in.permr_, NbLines * sizeof(int));
      std::memcpy(rpermr_, in.rpermr_, NbLines * sizeof(int));
      std::memcpy(permp_, in.permp_, 2 * NbLines * sizeof(int));
    }
  }
  return *this;
}
  

xCSRMatrix::xCSRMatrix(const xCSRMatrix& in) {

    throw;
  assert(0); cerr << "copy constructor of xCSRMatrix not debugged yet!" << endl;

  NbLines = in.NbLines;
  MatrixWasReordered = in.MatrixWasReordered;
  storage = in.storage;
  AllocationDone = in.AllocationDone;

  if (in.AllocationDone) {
    a_    = List_Copy (in.a_);
    ai_   = List_Copy (in.ai_);
    ptr_  = List_Copy (in.ptr_);
    jptr_ = List_Copy (in.jptr_);
  }
  
  if (in.MatrixWasReordered) {
    permr_  = new int[NbLines];
    rpermr_ = new int[NbLines];
    permp_  = new int[2 * NbLines];
    
    std::memcpy(permr_, in.permr_, NbLines * sizeof(int));
    std::memcpy(rpermr_, in.rpermr_, NbLines * sizeof(int));
    std::memcpy(permp_, in.permp_, 2 * NbLines * sizeof(int));

  }

}
//
void xCSRMatrix::OutputMatrixMatlabFormat(std::ostream &os, const std::string& name, const unsigned precision) const {

    xCSRMatrix* Matrix(const_cast<xCSRMatrix *>(this));

    int nl = Matrix->GetNbUnknown();
    os << "indexI_temp = [ " ;
    for (int i=1; i<= nl; i++) {
        xCSRMatrix::column_iterator it = Matrix->begin_column(i);
        xCSRMatrix::column_iterator ite = Matrix->end_column(i);
        while (it !=ite) {
            os << it.line << " " ;
            ++it;
        }
    }
    os << "]; " <<endl;

    os << "indexJ_temp = [ " ;
    for (int i=1; i<= nl; i++) {
        xCSRMatrix::column_iterator it = Matrix->begin_column(i);
        xCSRMatrix::column_iterator ite = Matrix->end_column(i);
        while (it !=ite) {
            os << it.column << " " ;
            ++it;
        }
    }
    os << "]; " << endl;

    os << "data_temp = [ " << setprecision (precision) ;
    for (int i=1; i<= nl; i++) {
        xCSRMatrix::column_iterator it = Matrix->begin_column(i);
        xCSRMatrix::column_iterator ite = Matrix->end_column(i);
        while (it !=ite) {
            os << *it << " " ;
            ++it;
        }
    }
    os << "]; " << endl;
    os << name << " = sparse(indexI_temp,indexJ_temp,data_temp,"<<nl<<","<<nl<<");" << endl;
    os << "clear indexI_temp; clear indexJ_temp; clear data_temp;" << endl;

}


void xCSRMatrix::OutputMatrixOctaveFormat(const string& filename) const
{
      std::ofstream  out(filename.c_str());
      OutputMatrixOctaveFormat(out);
      out.close();
}

void xCSRMatrix::OutputMatrixOctaveFormat(std::ostream &out) const
{
    ///double tol = 1E-15;
    double val;
    int p=out.precision(16);
    out << setiosflags(std::ios_base::scientific);
    int i,j;
    const int Asize = NbLines;
    for(i=1;i<=Asize;i++){
      for(j=1;j<=Asize;j++){
	val = GetMatrix(i,j) ;
	///if ( fabs(val) < tol ) val = 0.;
	out << val << " "; 
      }
      out << endl;
      
    }
    out << resetiosflags(std::ios_base::scientific);
    out.precision(p);
    return;
}



void xCSRMatrix::Load (const xCSRMatrix& in) {
  assert(0); //it is bugged I do not know why, God help me
  assert(AllocationDone && in.AllocationDone && NbLines == in.NbLines 
         && in.MatrixWasReordered == false);
  //we copy thr matrix
  printf("a_->nmax %d, in.a_->nmax %d\n", a_->nmax, in.a_->nmax);
  printf("ai_->nmax %d, in.ai_->nmax %d\n", ai_->nmax, in.ai_->nmax);
  printf("ptr_->nmax %d, in.ptr_->nmax %d\n", ptr_->nmax, in.ptr_->nmax);
  printf("jptr_->nmax %d, in.ajptr_->nmax %d\n", jptr_->nmax, in.jptr_->nmax);
  printf("a_->size %d, in.a_->size %d\n", a_->size, in.a_->size);
  printf("ai_->size %d, in.ai_->size %d\n", ai_->size, in.ai_->size);
  printf("ptr_->size %d, in.ptr_->size %d\n", ptr_->size, in.ptr_->size);
  printf("jptr_->size %d, in.jptr_->size %d\n", jptr_->size, in.jptr_->size);
  assert(a_->nmax == in.a_->nmax);
  assert(ai_->nmax == in.ai_->nmax);
  assert(ptr_->nmax == in.ptr_->nmax);
  assert(jptr_->nmax == in.jptr_->nmax);
  assert(a_->size == in.a_->size);
  assert(ai_->size == in.ai_->size);
  assert(ptr_->size == in.ptr_->size);
  assert(jptr_->size == in.jptr_->size);

  std::memcpy(a_->array,    in.a_->array,    a_->nmax * a_->size);
  std::memcpy(ai_->array,   in.ai_->array,   ai_->nmax * ai_->size);
  std::memcpy(ptr_->array,  in.ptr_->array,  ptr_->nmax * ptr_->size);
  std::memcpy(jptr_->array, in.jptr_->array, jptr_->nmax * jptr_->size);
  //a_    = List_Paste (in.a_); 
  //ai_   = List_Paste (in.ai_);
  //ptr_  = List_Paste (in.ptr_);
  //jptr_ = List_Paste (in.jptr_);

  MatrixWasReordered = false;
  delete []rpermr_; rpermr_ = nullptr;
  //printf("Helly6\n");
  delete []permr_; permr_ = nullptr;
  //printf("Helly7\n");
  delete []permp_; permp_ = nullptr;
  //printf("Helly8\n");
}






/////////////////////////////////////////////////////////
///////////////////////////////////////////////////////:
// listman functions


void *Malloc(size_t size)
{
  void *ptr;

  if (!size) return(nullptr);
  ptr = malloc(size);
  return(ptr);
}

void *Realloc(void *ptr, size_t size)
{
  if (!size) return(nullptr);
  ptr = realloc(ptr,size);
  return(ptr);
}

void Free(void *ptr)
{
//printf("Hellt1\n");
  if (ptr == nullptr) return;
//printf("Hellt2\n");
  free(ptr);
//printf("Hellt3\n");
  ptr = nullptr;
//printf("Hellt4\n");
}

/* - List ---------------------------------------------------------- */

List_T *List_Create(int n, int incr, int size)
{
  List_T *liste;

  if (n <= 0)  n = 1 ;
  if (incr <= 0) incr = 1;

  liste = (List_T *)Malloc(sizeof(List_T));

  liste->nmax    = 0;
  liste->incr    = incr;
  liste->size    = size;
  liste->n       = 0;
  liste->isorder = 0;
  liste->array   = nullptr;

  List_Realloc(liste,n);
  return(liste);
}

void List_Delete(List_T *liste)
{
  if (liste != nullptr) {
    Free(liste->array);
    Free(liste);
  }
}

void List_Realloc(List_T *liste,int n)
{
  char* temp;
  if (n <= 0) return;
  if (liste->array == nullptr) {
    liste->nmax = ((n - 1) / liste->incr + 1) * liste->incr;
/*if (liste->nmax * liste->size == 0) printf("dans realloc probleme avec array null\n");*/
    liste->array = (char *)Malloc(liste->nmax * liste->size);
  }
  else {
    if (n > liste->nmax) {
      liste->nmax = ((n - 1) / liste->incr + 1) * liste->incr;
/*      liste->array = (char *)Realloc(liste->array, liste->nmax * liste->size); */
/*if (liste->nmax * liste->size == 0) printf("dans realloc probleme second case\n");*/
      temp = (char *)Realloc(liste->array, liste->nmax * liste->size);
      liste->array = temp;
    }
  }
}

void List_Add(List_T *liste, void *data)
{
  liste->n++;

  List_Realloc(liste,liste->n);
  liste->isorder = 0;
  memcpy(&liste->array[(liste->n - 1) * liste->size],data,liste->size);
}

int List_Nbr(List_T *liste)
{
  return(liste->n);
}

int List_Nbr0(List_T *liste)
{
  return (liste)? liste->n : 0 ;
}

void List_Insert(List_T *liste, void *data,
		 int (*fcmp)(const void *a, const void *b))
{
  if (List_Search(liste,data,fcmp) == 0)
    List_Add(liste,data);
}

int List_Replace(List_T *liste, void *data,
		 int (*fcmp)(const void *a, const void *b))
{
  void *ptr;

  if (liste->isorder != 1) List_Tri(liste,fcmp);
  liste->isorder = 1;
  ptr = (void *) bsearch(data,liste->array,liste->n,liste->size,fcmp);
  if (ptr == nullptr) {
    List_Add(liste,data);
    return(0);
  }
  else {
    memcpy(ptr,data,liste->size);
    return (1);
  }
}

void List_Read(List_T *liste, int index, void *data)
{
  memcpy(data,&liste->array[index * liste->size],liste->size);
}

void List_Write(List_T *liste, int index, void *data)
{
  liste->isorder = 0;
  memcpy(&liste->array[index * liste->size],data,liste->size);
}

void List_Put(List_T *liste, int index, void *data)
{
  if (index >= liste->n) {
    liste->n = index + 1;
    List_Realloc(liste,liste->n);
    List_Write(liste,index,data);
  } else {
    List_Write(liste,index,data);
  }
}

void List_Pop(List_T *liste)
{
  liste->n -- ;
}

void *List_Pointer(List_T *liste, int index)
{
/*  liste->isorder = 0; */
  return(&liste->array[index * liste->size]);
}

void *List_Pointer_NoChange(List_T *liste, int index)
{
  return(&liste->array[index * liste->size]);
}

void List_Tri(List_T *liste, int (*fcmp)(const void *a, const void *b))
{
  qsort(liste->array,liste->n,liste->size,fcmp);
}

int List_Search(List_T *liste, void *data,
		 int (*fcmp)(const void *a, const void *b))
{
  void *ptr;

  if (liste->isorder != 1) { List_Tri(liste,fcmp) ; liste->isorder = 1 ; }
  ptr = (void *) bsearch(data,liste->array,liste->n,liste->size,fcmp);
  if (ptr == nullptr) return(0);
  return (1);
}

int List_ISearch(List_T *liste, void *data,
		 int (*fcmp)(const void *a, const void *b))
{
  void *ptr;

  if (liste->isorder != 1) List_Tri(liste,fcmp);
  liste->isorder = 1;
  ptr = (void *) bsearch(data,liste->array,liste->n,liste->size,fcmp);
  if (ptr == nullptr) return(-1);
  return (((long)ptr - (long)liste->array) / liste->size);
}

int List_Query(List_T *liste, void *data,
		 int (*fcmp)(const void *a, const void *b))
{
  void *ptr;

  if (liste->isorder != 1) List_Tri(liste,fcmp);
  liste->isorder = 1;
  ptr = (void *) bsearch(data,liste->array,liste->n,liste->size,fcmp);
  if (ptr == nullptr) return(0);

  memcpy(data,ptr,liste->size);
  return (1);
}



void *List_PQuery(List_T *liste, void *data,
		 int (*fcmp)(const void *a, const void *b))
{
  void *ptr;

  if (liste->isorder != 1) List_Tri(liste,fcmp);
  liste->isorder = 1;
  ptr = (void *) bsearch(data,liste->array,liste->n,liste->size,fcmp);
  return(ptr);
}

int List_Suppress(List_T *liste, void *data,
		 int (*fcmp)(const void *a, const void *b))
{
  char *ptr;
  int len;
  
  ptr = (char*)List_PQuery(liste,data,fcmp) ;
  if (ptr == nullptr) return(0);
  
  liste->n--;
  len = liste->n - (((long)ptr - (long)liste->array) / liste->size);
  if (len > 0) memmove(ptr, ptr + liste->size, len * liste->size);
  return(1);
}

void List_Reset(List_T *liste)
{
  liste->n = 0;
}


//Added by Nicolas to define a copy constructeur of CSR_Matrix

List_T *List_Copy(List_T *in) {
  List_T * liste;
  if (in == nullptr) liste = nullptr;
  else {
    liste = (List_T *)Malloc(sizeof(List_T));
    
    liste->nmax    = in->nmax;
    liste->incr    = in->incr;
    liste->size    = in->size;
    liste->n       = in->n;
    liste->isorder = in->isorder;
    liste->array   = nullptr;
    
    List_Realloc(liste,liste->nmax);
    
    memcpy(liste->array, in->array, liste->nmax * liste->size);
  }

  return liste;
}

// specialisation template of transfertToConnectionParameter :
// xCSRMatrix considered as xTraitMatrixUnSym and xTraitMatrixSparceCSR
// connecting object : xTraitMatrixSparceCOO, xTraitMatrixFindex (Mumps for example)
template <>
void xCSRMatrix::transfertToConnectionParameter<xTraitMatrixUnSym,xTraitMatrixSparceCSR,xTraitMatrixSparceCOO,xTraitMatrixFindex>(int &ns,int &ms,int &nnzs, int **idx1,int **idx2,double **dat, xTraitMatrixSparceCOO &matrix_storage_requested,xTraitMatrixFindex &indexing_requested)
{
    // not fondamental
    PutInCompressedRow();

    // set dimensions
    ns = ms = GetNbUnknown();
    // set number of nonzeros
    nnzs  = ReturnNumberOfNonZero();

    // to avoid transfert use of memory space of xCSRMatrix
    ( *dat ) = (double *) (ReturnValPointer());

    // CSR structure 
    // colone index are directely used in COO and as xCSRMatrix is in fortran indexing
    // everything is ok
    ( *idx2 ) = const_cast < int * >( ReturnColumIndexPointer());

    // get row pointer to construct row index
    int *ptri = const_cast < int * >( ReturnRowPointerShiftedToZero());

    // allocate row index using extra1 as container (will be deleted when object will be deleted)
    int *idx_r;
    extra1.resize(nnzs);
    idx_r = ( *idx1 ) = &extra1[0];

    // filling idx_r
    int i,k,l;
    for (i = 0; i < ns; )
    {
        i++;
        l = ptri[i];
        for (k = ptri[i-1]; k < l; k++)
        {
            idx_r[k] = i;
        }
    }

    return;

}
// specialisation template of transfertToConnectionParameter :
// xCSRMatrix considered as xTraitMatrixUnSym and xTraitMatrixSparceCSR
// connecting object : xTraitMatrixSparceCSR, xTraitMatrixCindex (SuperLu for example)
template <>
void xCSRMatrix::transfertToConnectionParameter<xTraitMatrixUnSym,xTraitMatrixSparceCSR,xTraitMatrixSparceCSR,xTraitMatrixCindex>(int &ns,int &ms,int &nnzs, int **idx1,int **idx2,double **dat, xTraitMatrixSparceCSR &matrix_storage_requested,xTraitMatrixCindex &indexing_requested)
{
    // not fondamental
    PutInCompressedRow();

    // set dimensions
    ns = ms = GetNbUnknown();
    // set number of nonzeros
    nnzs  = ReturnNumberOfNonZero();

    // to avoid transfert use of memory space of xCSRMatrix with shifted index as C is ascked
    ( *dat ) = (double *) (ReturnValPointer());
    ( *idx1 ) = const_cast < int * >( ReturnRowPointerShiftedToZero());
    ( *idx2 ) = const_cast < int * >( ReturnColumnIndexPointerShiftedToZero());

    return;

}

} // end of namespace
