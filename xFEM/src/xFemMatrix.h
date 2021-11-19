/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.
*/

#ifndef __FemMatrix_H
#define __FemMatrix_H

#include <vector>
#include <iostream>
#include <iterator>

namespace xfem
{

template<typename VT>
class xFemScalar {
public:
  xFemScalar() : sca(0.0), coeff(1.0) {}
  void add(const VT& val)
    {
      sca += coeff*val;
    } 
  void setCoeff(const VT& coe_acc) {coeff=coe_acc;}
  const VT& getCoeff() const {return coeff;}
  const VT& getVal() const { return sca; }
  void setVal(const VT& v) { sca = v; }
private:
  VT sca;
  VT coeff;
};

template<typename VT>
class xFemMatrix {

public:
  xFemMatrix() {coeff=1.0;}
  xFemMatrix(int a, int b) 
    { coeff =1.0; 
      resize(a,b);
    }
  void resize( int a,  int b) 
    {  
      nbrow  = a;  nbcol = b;
      int siz =nbrow*nbcol;
      vals.resize(siz, 0.0);
      std::fill(vals.begin(), vals.end(), 0.);
    }

  void setCoeff(const VT& coe_acc) {coeff=coe_acc;}
  const VT& getCoeff() const {return coeff;}

  template<typename UT>
  friend std::ostream& operator<<(std::ostream& s, const xFemMatrix<UT>& m);
  inline VT& operator() (int i, int j)       {return vals[i * nbcol + j];}
  inline const VT&  operator() (int i, int j) const {return vals[i * nbcol + j];}
  inline void add(int i, int j, const VT& v)
     {
      vals[i * nbcol + j] += coeff * v;
    }
  inline void addSym(int i, int j, const VT& v)
    {
//       vals[i * nbcol + j] += coeff * v;
//       if (i!=j) vals[j * nbcol + i] += coeff * v;

      //Only store upper triangular matrix
      if(j>=i) vals[i * nbcol + j] += coeff * v;
    }
  int getNbRow() const  {return nbrow;}
  int getNbCol() const {return nbcol;}
  VT * ReturnValPointer ( ){return &vals[0];};
  const VT * ReturnValPointer ( ) const {return &vals[0];};

  void symmetrize(){
//     std::cout<<"Syme\n";
//Symmetrize the upper-triangular matrix into a full matrix
    for(int i=0;i<nbrow; ++i){
      int inbc=i * nbcol;
      for(int j=i+1; j<nbcol; ++j){
        int jnbc=j * nbcol;
        vals[jnbc + i]=vals[inbc + j];
      }
    } 
  }

  
private:
  int nbrow;
  int nbcol;
  VT coeff;
  std::vector<VT> vals;
};

/// compute B = alpha*A +B; for double
//! It use blas axpy
// note: it is note implemented at this level to avoid xBlasDef.h header to be included at
// this level. Including this header would add a extra dependency to xSolverBase library to
// all code including this header ... Not wanted for now.
// note: it is not declared inlined with implementation in .cc as some compiler seems to 
// have problems to find out these symbols at linking stage  
void axpy(const double &alpha, const xFemMatrix<double> &A, xFemMatrix<double> &B);


// compute B = transpose(A);
template<typename VT>
inline void transpose(const xFemMatrix<VT> &A,  xFemMatrix<VT> &B){
  int nl =  A.getNbRow();
  int nc =  A.getNbCol();
  B.resize(nc, nl);
  for (int i=0; i < nl; ++i){
    for (int j=0; j < nc; ++j){
      B(j,i) = A(i,j);
    }
  }
}
 
template<typename UT>
inline std::ostream& operator<<(std::ostream& s, const xFemMatrix<UT>& m)
  {
    std::cout << " FemMatrix of size " << m.nbrow << " rows times " << m.nbcol << " columns" << std::endl; 
    for ( int i = 0; i < m.nbrow; ++i)
      {
	for ( int j = 0; j < m.nbcol; ++j)
	  {
	    std::cout << m(i,j) << " ";
	  }
	std::cout << std::endl;
      }
    return s;
  }


template<typename VT>
class xFemVector {

public:
  xFemVector(): coeff(1.0) {} 
  xFemVector(int a)  
    {
      resize(a);
    }
  void resize(int a)
    {
      nbval = a;
      vals.resize(nbval, 0.0);
      for (unsigned int i = 0; i < vals.size(); ++i) vals[i] = 0.0;
    }    
  const VT&  operator[] (int i) const {return vals[i];}
  const VT&  operator() (int i) const {return vals[i];}
  //a virer
  //a virer double& operator[] (int i)       {return vals[i];} 

  template<typename UT>
  friend std::ostream& operator<<(std::ostream& s, const xFemVector<UT>& v);
  inline void add(int i, const VT& v)
     {
      vals[i] += coeff * v;
    }
  void setCoeff(const VT& coe_acc) {coeff=coe_acc;}
  const VT& getCoeff() const {return coeff;}


private:
  int nbval;
  VT coeff;
  std::vector<VT> vals;
};

template<typename UT>
inline std::ostream& operator<<(std::ostream& s, const xFemVector<UT>& v)
  {  
    std::cout << " FemVector of size " << v.nbval << std::endl;
    std::copy(v.vals.begin(), v.vals.end(), std::ostream_iterator<UT>(s, " "));
    return s;
  }

} // end of namespace


#endif








