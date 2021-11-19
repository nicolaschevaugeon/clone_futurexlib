/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.
*/

#ifndef __STRING_MANAGER__
#define __STRING_MANAGER__


#include <map>
#include <string>
#include <set>
#include <iostream>
#include <functional>
#include <limits>

namespace xtool
{
  template <class T>
    class firstinpaircontaining{
  public:
  firstinpaircontaining(const std::string &subs_):subs(subs_){}
    bool operator()(const std::pair<std::string, T>  & pair){
      return (pair.first.find(subs) != std::string::npos)  ;
    }
    std::string subs;
  };

//---------------Hash function-------------------------------------------------------------------------------------
/// Generic hash function for xStringManager
//! For general case It use std::hash<std::string> function from C++ standard library.
//! To give a conforming type to what user decide to use in xStringManager the result of this function
//! is casted by T. This unfortunately increase collision for T type that are small.
//
//! Example for T=short:
//! std::hash<std::string> function give for "INTEG_POINT_308" string : 13975685018915129296
//! std::hash<std::string> function give for "INTEG_POINT_174" string : 17478252458999826384
//! Those 2 hash values are distinct but when they are casted to short they both give 18384 !!! 
//
//! For this short type a specialization is proposed. 
//
// 
template<typename  T>
class hashString
{
    public:
        T operator () (const std::string& s)
        {
            return  static_cast<T>(str_hash(s));
        }
    private:
        std::hash<std::string> str_hash;
};
/// hashString specialization for T being a short
//! It is base on djb2 algorithm (k=33) first reported by Dan Bernstein many years ago in comp.lang.c 
//! For example  we have 
//! for "INTEG_POINT_308" string : 3464
//! for "INTEG_POINT_174" string : 1513
//
//! but is in default with "GAUSS_1010" and "GAUSS_9" !!!!
//
template<>
class hashString <short>
{
    public:
        hashString<short>():mx(std::numeric_limits<unsigned short>::max()){}
        short operator () (const std::string& s)
        {
            unsigned short h = 0;
            for (auto c : s)
            {
                h = ((( h << 5 ) + h ) + c)%mx; /* hash * 33 + c */
            }

            return static_cast<short>(h);
        }
    private:
        unsigned short mx;

};
/// No hashing traits
//! This type is the default for class xStringManager
class noHashString {}; 

//---------------Policy--------------------------------------------------------------------------------------------
/// policy to be able to treat different strategy under the same xStringManager class
template <typename RT, typename HASH>
class  xStringManagerPolicy;
template <typename RT>
class  xStringManagerPolicy<RT,noHashString>
{
    public:
        static RT getIdPolicy(const std::string& s,size_t nb) 
        {
            return  static_cast<RT>(nb);
        }
        static void checkIdPolicy(const std::map<RT, std::string>& t,RT nb,const std::string& s) 
        {
             ; // no check in this context
        }
};
template <typename RT>
class  xStringManagerPolicy<RT,hashString<RT> >
{
    public:
        static RT getIdPolicy(const std::string& s,size_t nb) 
        {
            // The c++11 hash class is expected to give no collision in regard of the set
            // of string used in xfem.
            // Unfortunately in Xfem xKeyInfo if using a T of type  smaller then the size_t type
            // returned by C++11 hash operator it introduce problem. This first casting, that appears here with static_cast call,
            // is a first concern:  does the casting will not  increase collision ?
            // To avoid this ugly rising collision probability due to casting,  a specific hash function may have to be set to 
            // produce directly the smallest type (short).
            return hash_func(s);
        }
        static void checkIdPolicy(const std::map<RT, std::string>& t,RT nb,const std::string& s) 
        {
            auto itn = t.find(nb);
            if   (itn != t.end()) 
            {
                std::cout<<"Collision occurs !!!"<<std::endl<<"String "<<s<<" is not yet in the stringmanager and its hash number is "<<nb<<std::endl;
                std::cout<<"Unfortunately string "<<itn->second<<" is having the same hash number"<<std::endl;
                std::cout<<"One solution is to  change one of those two strings  so that they do not collied"<<std::endl;
                std::cout<<"An other solution is to  change hash function so that it do not collied for those 2 strings ... and others"<<std::endl;
                throw -2345;
            }
        }
    private:
        static hashString<RT> hash_func;
};

template < typename RT >
hashString < RT > xStringManagerPolicy < RT,hashString < RT > >::hash_func=hashString < RT >();

//-----------------------------------------------------------------------------------------------------------------
/// xString Manager class
template<typename  T,  typename H=noHashString >
class xStringManager {

public:
  typedef typename std::map<std::string, T>::const_iterator const_iterator;
  typedef typename std::map<std::string, T>::iterator             iterator;
  typedef typename std::map<T, std::string>::const_iterator const_iterator2;
  typedef typename std::map<T, std::string>::iterator             iterator2;
  
  T  getId(const std::string& s) 
  {
    const_iterator it = s_.find(s);
    if   (it != s_.end()) return it->second;

    T nb = xStringManagerPolicy<T,H>::getIdPolicy(s, s_.size());
    xStringManagerPolicy<T,H>::checkIdPolicy(t_, nb,s);

    s_[s]  = nb;
    t_[nb] =  s;
    return nb; 
  }
  
  std::set<T > getIdContaining(const std::string & s){
    std::set<T> res;
    iterator itend = s_.end();
    iterator it = s_.begin();
    firstinpaircontaining<T > pred(s);
    while (it !=itend){
      it = find_if ( it, itend, pred);
      if (it !=itend) {
	res.insert(it->second) ;
	++it;
      }
    }
    return res;
  }
    
  std::string getName(T t)  {
    const_iterator2 it = t_.find(t);
    if   (it != t_.end()) return it->second;
    std::string s = "dummy";
    t_[t]  = s;
    s_[s]  = t;
    return s; 
  }

  void clear(){s_.clear();t_.clear();}
  
public:
 typedef T int_type_t;

private:
  T n_;
  std::map<std::string, T> s_;
  std::map<T, std::string> t_;

};
} // end of namespace

#endif
