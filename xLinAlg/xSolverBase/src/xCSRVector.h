/* 
   This file is a part of eXlibris C++ Library
   under the GNU Lesser General Public License.
   See the NOTICE.md & LICENSE.md files for terms 
   and conditions.
*/
#ifndef _CSR_VECTOR_H
#define _CSR_VECTOR_H

#include <iostream> 
#include <fstream>
#include <string>
#include <vector>
#include <iterator>
#include <algorithm>
#include <cmath>
#include <cassert>


namespace xlinalg 
{

  class xCSRVector 
  {
  private : 
    typedef std::vector<double> Vector;
    Vector Array;
  public :
    typedef Vector::iterator iterator;
    typedef Vector::const_iterator const_iterator;
    typedef Vector::reference reference;
    typedef Vector::reference const_reference;
    typedef Vector::value_type value_type;

    inline xCSRVector( ) {};
    inline xCSRVector( const int& n )
    {
      Array.insert(Array.begin(), n, 0.0);
    }
    inline xCSRVector(const xCSRVector& rhs)
      = default;

    inline xCSRVector & operator +=(const xCSRVector& rhs){
      if (size()!=rhs.size()) throw;
      std::transform(rhs.begin(), rhs.end(), begin(), begin(), std::plus<double>());
      return (*this);
    }

    inline xCSRVector & operator -=(const xCSRVector& rhs){
      if (size()!=rhs.size()) throw;
      std::transform(begin(), end(), rhs.begin(), begin(), std::minus<double>());
      return (*this);
    } 

    inline xCSRVector & operator *=(const double scale){
      std::transform(begin(), end() , begin(), bind2nd(std::multiplies<double>(), scale));
      return (*this);
    }
  
    inline void resize(const int& n) 
    {
      Array.resize(n);
    }

    inline ~xCSRVector() = default;

    inline iterator begin() { return Array.begin();}
    inline iterator end() { return Array.end();}

    inline const_iterator begin() const { return Array.begin();}
    inline const_iterator end() const { return Array.end();}

    inline xCSRVector & operator=(const xCSRVector & rhs)
      {
	if (this != &rhs)
	  Array = rhs.Array;
	return (*this);
      }  
    inline double operator()(const int & i) const{
      return Array[i];
    }
  
    inline double & operator()(const int & i){
      return Array[i];
    }
  
    inline void duplicate(const xCSRVector& rhs)   
    {
      if (this != &rhs)  Array = rhs.Array;
    }
			
    inline double scalar_product(const xCSRVector& rhs)
    {
      double res=0;
      for (int it=1; it<=getSize(); it++)
	{	res += GetVal(it)*rhs.GetVal(it);}
      return res;
    }


    inline void AddVal( const int& i, const double& val)
    {
      // FORTRAN NUMBERING ? (oui)
      //#pragma omp critical
      Array[i-1] += val;
    }

    inline double GetVal( const int i) const
    {
      // FORTRAN NUMBERING ? (oui)
      return Array[i-1];
    }


    inline const double*  GetArray() const
    {
      return &(*Array.begin());
    }

    inline double*  GetArray()
    {
      return &(*Array.begin());
    }
  
    inline void ZeroArray()
    {
      std::fill(Array.begin(), Array.end(), 0.0);
    }
  
    inline void PrintArray( std::ostream  & out) const
    {
      out << "###########################\n";
      out << "Vector of dimension " << Array.size() << " is:\n";
      std::copy(Array.begin(), Array.end(), std::ostream_iterator<double>(out, " ")); 
      //for (unsigned int i = 0; i < Array.size() ; ++i) out << Array[i] << " ";
      out << "\n";
      out << "###########################\n";
    }
  
    inline void OutputVector( std::ostream  & out) const
    {
      out <<  Array.size() << std::endl;  
      std::copy(Array.begin(), Array.end(), std::ostream_iterator<double>(out, " ")); 
      //for (unsigned int i = 0; i < Array.size() ; ++i) out << Array[i] << " ";
      out << "\n";
    }
  
    inline void LoadVector( std::istream  & in)
    {
      double n;//Only the first element is an integer, all the other one are doubles !
      in >> n;
      Array.resize((int) n);
      for (unsigned int i = 0; i < Array.size() ; ++i) in >> Array[i];
    }

  

    inline void PrintArray(const std::string& filename) const
    {
      std::ofstream out(filename.c_str());
      PrintArray(out);
      out.close();
    }
  
    inline double& operator[](unsigned int i)
    {
      return Array[i];
    }

    inline const  double& operator[](unsigned int i) const 
    {
      return Array[i];
    }

    inline int size() const {return Array.size() ;}
    inline int getSize() const {return Array.size() ;}
    void OutputOctaveFormat(std::ostream &out) const;
    void OutputOctaveFormat(const std::string  &out) const;

    inline double two_norm() const 
    {
      double res=0;
      for (const_iterator it=Array.begin(); it!=Array.end(); ++it)
	{	res+= (*it)*(*it);}

      return sqrt(res);
    }
  };

  void   axpy(const double &a, const xCSRVector & x,  xCSRVector & y);
  double nrm2( const xCSRVector & x);
  double dot( const xCSRVector & x, const xCSRVector &y);
  void   scal(const double& alpha,    xCSRVector & x);
  
  inline xCSRVector  operator +(const xCSRVector& A, const xCSRVector &B){
    if (A.size()!=B.size()) throw;
    xCSRVector C(A);
    C+=B;
    return C;
  }

  
  
  inline xCSRVector  operator -(const xCSRVector& A, const xCSRVector &B){
    if (A.size()!=B.size()) throw;
    xCSRVector C(A);
    C-=B;
    return C;
  }


  inline xCSRVector  operator *(const xCSRVector& A, const double &alpha){
    xCSRVector C(A);
    scal(alpha, C);
    return C;
  }
  
 
  
} // end of namespace

#endif







