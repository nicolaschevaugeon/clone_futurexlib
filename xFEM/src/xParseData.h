/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.
*/

#ifndef _xParseData_
#define _xParseData_
#include <map>
#include <string>
#include <list>
#include <vector>
#include "xVector.h"
#include "xTensor2.h"

//! xParseData is a small class to store data of vrious type together in a data structure and read them from an input file.
/*!
  when a xParseData is created, the user can then register data, with a name and a type, using the registers functions. Defaults values can be given at the registering stage.
 when this is done, function parse will scan an inputfile for the data, report if some data are missing in the input, or are not registered. The the user can retrive data from their name and type, using the get functions.
*/
struct xParseData{
    enum datatype {INT, DOUBLE, XVECTOR, STRING, LISTSTRING, UNDEFINED, LISTXVECTOR, XTENSOR2, MAPSTRINGDOUBLE};
  void registerInt(std::string name, int defaultval=0);
  void registerDouble(std::string name, double defaultval=0.);
  void registerString(std::string name, std::string defaultval="");
  void registerListString(std::string name, std::list < std::string> defaultval=std::list< std::string>() );
  void registerVector(std::string name, xtensor::xVector<> defaultval= xtensor::xVector<>(0.,0.,0.));
  void registerListVector(std::string name, std::list<xtensor::xVector<> > defaultval=std::list<xtensor::xVector<> >() );
  void registerTensor2(std::string name, xtensor::xTensor2<> defaultval= xtensor::xTensor2<>(0.));
  void registerMapStringDouble(std::string name, std::map<std::string, double> defaultval= std::map<std::string, double>());//GREG
  std::string getString(std::string name) const;
  void setString(std::string name,const std::string &s );
  std::list<std::string > getListString(std::string name) const;
  std::list<xtensor::xVector<> > getListVector(std::string name) const;
  std::map<std::string, double> getMapStringDouble(std::string name) const;//GREG
  double getDouble(std::string name) const;
  void setDouble(std::string name, double d);
  xtensor::xVector<> getVector(std::string name) const;
  int getInt(std::string name) const;
  void setInt(std::string name, int i);
  datatype getType(std::string name) const;
  xtensor::xTensor2<> getTensor2(std::string name) const;
  int parse(std::istream &);
private:
  std::map<std::string, int>  intdata;
  std::map<std::string, double>  doubledata;
  std::map<std::string, std::string>  stringdata;
  std::map<std::string, std::list<std::string> >  liststringdata;
  std::map<std::string, xtensor::xVector<> > xvectordata;
  std::map<std::string, std::list<xtensor::xVector<> > > listvectordata;
  std::map<std::string, int> isdefined;
  std::map<std::string, xtensor::xTensor2<> > xtensor2data;
  std::map<std::string, std::map<std::string, double> > mapstringdoubledata;//GREG
};
#endif
