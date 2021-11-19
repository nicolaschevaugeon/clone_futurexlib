/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.
*/

#ifndef _VALKEY_H
#define _VALKEY_H

#include <iostream>
#include <functional>
#include <string>

#include "mEntity.h" 
#include "xStringManager.h"
#include <boost/functional/hash.hpp>


//! The Key of a Value
namespace xfem
{

template <typename T>
T defaultIdPhys(void);
template <typename T>
T defaultIdGeom(void);

class xValKey {

public:
  typedef size_t ids_size_t;

private:
  ids_size_t      Phys;  // Physical Nature of the val
  ids_size_t      Geom;  // Geometrical Nature of the val
  AOMD::mEntity*  Enti;  // mesh entity to which the key is associated
  int             Refe;  // additional key if needed

public:
  friend std::ostream & operator << (std::ostream & ofs, const xValKey&); 

  AOMD::mEntity*     getEnti() const  { return Enti; }
  void             setEnti(AOMD::mEntity* in) { Enti = in; }

  ids_size_t   getGeom() const                  {return Geom;}
  void         setGeom(ids_size_t in)           {Geom = in;}

  ids_size_t   getPhys() const                  {return Phys;}
  void         setPhys(ids_size_t in)           {Phys = in; }

  int          getRefe() const {return Refe; }
  void         setRefe(int in)             {Refe = in; }

   xValKey(ids_size_t a = defaultIdPhys<ids_size_t>(), ids_size_t b = defaultIdGeom<ids_size_t>()) : Phys(a), Geom(b), Enti(nullptr), Refe(-1) {}
    xValKey(ids_size_t a, ids_size_t b, AOMD::mEntity* c, int d = -1) :
     Phys(a), Geom(b), Enti(c), Refe(d) {}
    xValKey(const xValKey & rhs) = default;



  friend bool operator < (const xValKey & c1, const xValKey & c2){
    return ( compkey (c1,c2) == -1 );
  }
  friend bool operator > (const xValKey & c1, const xValKey & c2){
    return ( compkey (c1,c2) ==  1 );
  }
  friend bool operator == (const xValKey & c1, const xValKey & c2){
    return ( compkey (c1,c2) ==  0 );
  }

  friend int compkey (const xValKey & c1, const xValKey & c2){
    if(c1.Enti  > c2.Enti) return 1;
    if(c1.Enti  < c2.Enti) return -1;
    if(c1.Phys  > c2.Phys) return 1;
    if(c1.Phys  < c2.Phys) return -1;
    if(c1.Geom  > c2.Geom) return 1;
    if(c1.Geom  < c2.Geom) return -1;
    if(c1.Refe  > c2.Refe) return 1;
    if(c1.Refe  < c2.Refe) return -1;
    return 0;
  };

//    static xKeyInfo info;
};


struct xHashValKey {
  int operator()(const xValKey& key) const 
    {
      AOMD::mEntity * e = key.getEnti();
      return e->getId();
    }
};


//Far more efficient for pointwise fields
struct xHashValKeyGauss {

    int operator()(const xValKey& key) const {

        AOMD::mEntity * e = key.getEnti();
        std::size_t seed = 0;
        boost::hash_combine(seed, e->getId() );
        boost::hash_combine(seed, key.getGeom());
        return static_cast<int>(seed);
    }
};

struct xEqualValKey  {
  bool operator()(const xValKey& key1, const xValKey& key2) const 
  {return compkey(key1, key2) == 0; }
};


class xKeyInfo {

public :
typedef  xtool::xStringManager<xValKey::ids_size_t,xtool::hashString<xValKey::ids_size_t> > string_manager_t;
typedef  string_manager_t::int_type_t int_type_t; 

public:
static std::string getPhysName(int_type_t in)         {return PhysString.getName(in);} 
static int_type_t    getPhysId(const std::string& in) {return PhysString.getId(in);  }
static std::set<int_type_t> getPhysIdContaining(const std::string & in){return PhysString.getIdContaining(in); }
static std::string getGeomName(int_type_t in)         {return GeomString.getName(in);} 
static int_type_t    getGeomId(const std::string& in) {return GeomString.getId(in);  } 
static std::set<int_type_t> getGeomIdContaining(const std::string & in){return GeomString.getIdContaining(in); }
static void clearPhysName() {PhysString.clear();}
static void clearGeomName() {GeomString.clear();}

private:
static string_manager_t PhysString;
static string_manager_t GeomString;

};

template <typename T>
T defaultIdPhys(void)
{
    static_assert( std::is_same < typename xKeyInfo::int_type_t, T >::value, "T must be of  xKeyInfo::int_type_t type") ;
    return xKeyInfo::getPhysId("no_phys");
}
template <typename T>
T defaultIdGeom(void)
{
    static_assert( std::is_same < typename xKeyInfo::int_type_t, T >::value, "T must be of  xKeyInfo::int_type_t type") ;
    return xKeyInfo::getGeomId("no_geom");
}


typedef std::unary_function<xValKey, bool> key_filter_t;

class xBoth : public  std::unary_function<xValKey, bool> {
public:
  template <class T1, class T2>
  xBoth(const T1& f1, const T2& f2) :
   filter1(f1), filter2(f2) {}
  bool operator()(const xValKey& v) const  {return (filter1(v) && filter2(v));}
private:
  mutable std::function<bool (xValKey)> filter1;
  mutable std::function<bool (xValKey)> filter2;
};
struct xIsOnEntity : public key_filter_t {
  xIsOnEntity(AOMD::mEntity* e) : intern(e) {}
  bool operator()(const xValKey& v) const {return (v.getEnti() == intern);}
  AOMD::mEntity* intern;
};
struct xIsNotOnEntity : public key_filter_t {
  xIsNotOnEntity(AOMD::mEntity* e) : intern(e) {}
  bool operator()(const xValKey& v) const {return (v.getEnti() != intern);}
  AOMD::mEntity* intern;
};
struct xIsPhys : public key_filter_t {
  xIsPhys(const std::string& s)  { phys_id = xKeyInfo::getPhysId(s); }
  bool operator()(const xValKey& v) const {return (v.getPhys() == phys_id);}
  xValKey::ids_size_t phys_id;
};
struct xIsNotPhys : public key_filter_t {
  xIsNotPhys(const std::string& s)  { phys_id = xKeyInfo::getPhysId(s); }
  bool operator()(const xValKey& v) const {return (v.getPhys() != phys_id);}
  xValKey::ids_size_t phys_id;
};
struct xIsGeom : public key_filter_t {
  xIsGeom(const std::string& s)  { geom_id = xKeyInfo::getGeomId(s); }
  bool operator()(const xValKey& v) const {return (v.getGeom() == geom_id);}
  xValKey::ids_size_t geom_id;
};
struct xBelongsToEntity : public key_filter_t {
  xBelongsToEntity(AOMD::mEntity* e_) : e(e_) {}
  bool operator()(const xValKey& key) const {return (e == key.getEnti() || 
                             (e->find(key.getEnti()) && ( e->getLevel() >
                              key.getEnti()->getLevel()) )  );}
  AOMD::mEntity* e;
};



class xValKeyExtend 
{
public:
  xValKeyExtend(const std::string& e);
  virtual ~xValKeyExtend()= default;;
  void operator()(xValKey& key); 
  virtual void addToExtension(const std::string& a) {std::cout<<"No implementation of addToExtension for xValKeyExtend is done, use derived class"; throw; }
  virtual void resetExtension(void) {std::cout<<"No implementation of resetExtension for xValKeyExtend is done, use derived class"; throw; }
protected:
  std::string extension;
};

class xValKeyMultiExtend : public xValKeyExtend
{
    public:
        xValKeyMultiExtend(const std::string& e):xValKeyExtend(e),base_extension(e) {}
  ~xValKeyMultiExtend() override{
    //std::cout <<"  ~xValKeyMultiExtend " << std::endl;
  };
        void addToExtension(const std::string& a) override {extension=base_extension+"_"+a;} 
        void resetExtension(void) override {extension=base_extension;} 
    private:
        std::string base_extension;
};

typedef std::function<void (xValKey&)>  ValKeyModifier_t;
typedef xValKeyExtend*  ValKeyModifierPtr;

struct xAcceptAllKeys : public key_filter_t {
  xAcceptAllKeys()  { };
  bool operator()(const xValKey& v) const {return true;}
};

struct xAcceptOnPhysString : public key_filter_t {
  xAcceptOnPhysString(const std::string & physString_): physString(physString_)  { };
  bool operator()(const xValKey& v) const {return xKeyInfo::getPhysName(v.getPhys())== physString;}
  
private:
  std::string physString;
};


} // end of namespace

#endif









