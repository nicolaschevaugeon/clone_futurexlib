/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/

#ifndef __Variables_H
#define __Variables_H

#include <boost/mpl/identity.hpp>
#include <cassert>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

#include "xTable.h"
#include "xTensor2.h"
#include "xTensor4.h"
#include "xTensorsPtr.h"
#include "xValue.h"
#include "xVector.h"

namespace xfem
{
class xTensors;
class xPieceWiseLinear
{
  public:
   xPieceWiseLinear(){};
   xPieceWiseLinear(const std::map<double, double>& val);
   void read(std::istream&);
   double operator()(double) const;
   std::map<double, double> values;
};

class xTensorsSignature
{
  public:
   typedef std::unordered_map<std::string, int> container_t;
   typedef container_t::const_iterator const_iterator;

  public:
   void register_scalar(const std::string& s) { scalars.insert(make_pair(s, scalars.size())); }
   int size_scalar() const { return scalars.size(); }
   const_iterator begin_scalar() const { return scalars.begin(); }
   const_iterator end_scalar() const { return scalars.end(); }
   bool exist_scalar(const std::string& s) const
   {
      //      printf("exist scalar\n");
      if (scalars.find(s) == scalars.end()) return false;
      return true;
   }
   int scalar(const std::string& s) const
   {
      const_iterator it = scalars.find(s);
      if (it == scalars.end()) assert(0);
      return it->second;
   }

   void register_vector(const std::string& s) { vectors.insert(make_pair(s, vectors.size())); }
   int size_vector() const { return vectors.size(); }
   const_iterator begin_vector() const { return vectors.begin(); }
   const_iterator end_vector() const { return vectors.end(); }
   bool exist_vector(const std::string& s) const
   {
      //      std::cout<<"exist vector "<<s<<std::endl;
      if (vectors.find(s) == vectors.end()) return false;
      return true;
   }
   int vector(const std::string& s) const
   {
      const_iterator it = vectors.find(s);
      if (it == vectors.end()) assert(0);
      return it->second;
   }

   void register_tensor2(const std::string& s) { tensor2s.insert(make_pair(s, tensor2s.size())); }
   int size_tensor2() const { return tensor2s.size(); }
   const_iterator begin_tensor2() const { return tensor2s.begin(); }
   const_iterator end_tensor2() const { return tensor2s.end(); }
   int tensor2(const std::string& s) const
   {
      const_iterator it = tensor2s.find(s);
      if (it == tensor2s.end()) assert(0);
      return it->second;
   }

   void register_tensor4(const std::string& s) { tensor4s.insert(make_pair(s, tensor4s.size())); }
   int size_tensor4() const { return tensor4s.size(); }
   const_iterator begin_tensor4() const { return tensor4s.begin(); }
   const_iterator end_tensor4() const { return tensor4s.end(); }
   int tensor4(const std::string& s) const
   {
      const_iterator it = tensor4s.find(s);
      if (it == tensor4s.end()) assert(0);
      return it->second;
   }

   void register_string(const std::string& s) { strings.insert(make_pair(s, strings.size())); }
   int size_string() const { return strings.size(); }
   const_iterator begin_string() const { return strings.begin(); }
   const_iterator end_string() const { return strings.end(); }
   bool exist_string(const std::string& s) const
   {
      if (strings.find(s) == strings.end()) return false;
      return true;
   }
   int astring(const std::string& s) const
   {
      const_iterator it = strings.find(s);
      if (it == strings.end()) assert(0);
      return it->second;
   }

   void register_table(const std::string& s) { tables.insert(make_pair(s, tables.size())); }
   int size_table() const { return tables.size(); }
   const_iterator begin_table() const { return tables.begin(); }
   const_iterator end_table() const { return tables.end(); }
   bool exist_table(const std::string& s) const
   {
      if (tables.find(s) == tables.end()) return false;
      return true;
   }
   int table(const std::string& s) const
   {
      const_iterator it = tables.find(s);
      if (it == tables.end()) assert(0);
      return it->second;
   }

   void register_piecewiselinear(const std::string& s) { piecewiselinears.insert(make_pair(s, piecewiselinears.size())); }
   int size_piecewiselinear() const { return piecewiselinears.size(); }
   const_iterator begin_piecewiselinear() const { return piecewiselinears.begin(); }
   const_iterator end_piecewiselinear() const { return piecewiselinears.end(); }
   bool exist_piecewiselinear(const std::string& s) const
   {
      if (piecewiselinears.find(s) == piecewiselinears.end()) return false;
      return true;
   }
   int piecewiselinear(const std::string& s) const
   {
      const_iterator it = piecewiselinears.find(s);
      if (it == piecewiselinears.end()) assert(0);
      return it->second;
   }

  private:
   // optimisation
   // remplacer par vector ce sera plus rapide
   container_t scalars;
   container_t vectors;
   container_t tensor2s;
   container_t tensor4s;
   container_t strings;
   container_t tables;
   container_t piecewiselinears;
};

class xTensors
{
  public:
   xTensors() : sgn(nullptr) {}
   // attention : *s ne devrait pas etre modifie par la suite ( Tenssig * const pas pratique)
   void setSignature(xTensorsSignature* s)
   {
      if (s) sgn = s;

      scalars.resize(sgn->size_scalar());
      fill(scalars.begin(), scalars.end(), 0);
      vectors.resize(sgn->size_vector());
      tensor2s.resize(sgn->size_tensor2());
      tensor4s.resize(sgn->size_tensor4());
      strings.resize(sgn->size_string());
      tables.resize(sgn->size_table());
      piecewiselinears.resize(sgn->size_piecewiselinear());
   }

   const xTensorsSignature* getSignature() const { return sgn; }

   xTensors(xTensorsSignature* s) { setSignature(s); }

   double& scalar(const std::string& name) { return scalars[sgn->scalar(name)]; }

   const double& scalar(const std::string& name) const { return scalars[sgn->scalar(name)]; }

   xtensor::xVector<>& vector(const std::string& name) { return vectors[sgn->vector(name)]; }

   const xtensor::xVector<>& vector(const std::string& name) const { return vectors[sgn->vector(name)]; }

   xtensor::xTensor2<>& tensor2(const std::string& name) { return tensor2s[sgn->tensor2(name)]; }

   const xtensor::xTensor2<>& tensor2(const std::string& name) const { return tensor2s[sgn->tensor2(name)]; }

   xtensor::xTensor4<>& tensor4(const std::string& name) { return tensor4s[sgn->tensor4(name)]; }

   const xtensor::xTensor4<>& tensor4(const std::string& name) const { return tensor4s[sgn->tensor4(name)]; }

   std::string& astring(const std::string& name) { return strings[sgn->astring(name)]; }

   const std::string& astring(const std::string& name) const { return strings[sgn->astring(name)]; }

   xTable& table(const std::string& name) { return tables[sgn->table(name)]; }

   const xTable& table(const std::string& name) const { return tables[sgn->table(name)]; }

   xPieceWiseLinear& piecewiselinear(const std::string& name) { return piecewiselinears[sgn->piecewiselinear(name)]; }

   const xPieceWiseLinear& piecewiselinear(const std::string& name) const { return piecewiselinears[sgn->piecewiselinear(name)]; }

   void read(const std::string& filename);
   void read(FILE* fp, const std::string& name);
   void read_scalar(FILE* fp, double& val);
   void read_vector(FILE* fp, xtensor::xVector<>& val);
   void read_string(FILE* fp, std::string& val);
   void read_table(FILE* fp, xTable& val);
   void read_piecewiselinear(FILE* fp, xPieceWiseLinear& val);
   // helper function to code with templates
   double& get(const std::string& name, boost::mpl::identity<double>) { return scalar(name); }
   const double& get(const std::string& name, boost::mpl::identity<double>) const { return scalar(name); }
   xtensor::xVector<>& get(const std::string& name, boost::mpl::identity<xtensor::xVector<>>) { return vector(name); }
   const xtensor::xVector<>& get(const std::string& name, boost::mpl::identity<xtensor::xVector<>>) const { return vector(name); }
   xtensor::xTensor2<>& get(const std::string& name, boost::mpl::identity<xtensor::xTensor2<>>) { return tensor2(name); }
   const xtensor::xTensor2<>& get(const std::string& name, boost::mpl::identity<xtensor::xTensor2<>>) const
   {
      return tensor2(name);
   }
   xtensor::xTensor4<>& get(const std::string& name, boost::mpl::identity<xtensor::xTensor4<>>) { return tensor4(name); }
   const xtensor::xTensor4<>& get(const std::string& name, boost::mpl::identity<xtensor::xTensor4<>>) const
   {
      return tensor4(name);
   }
   std::string& get(const std::string& name, boost::mpl::identity<std::string>) { return astring(name); }
   const std::string& get(const std::string& name, boost::mpl::identity<std::string>) const { return astring(name); }
   xTable& get(const std::string& name, boost::mpl::identity<xTable>) { return table(name); }
   const xTable& get(const std::string& name, boost::mpl::identity<xTable>) const { return table(name); }
   const xPieceWiseLinear& get(const std::string& name, boost::mpl::identity<xPieceWiseLinear>) const
   {
      return piecewiselinear(name);
   }

   friend std::ostream& operator<<(std::ostream& s, const xfem::xTensors& t);

  private:
   xTensorsSignature* sgn;
   std::vector<double> scalars;
   std::vector<xtensor::xVector<>> vectors;
   std::vector<xtensor::xTensor2<>> tensor2s;
   std::vector<xtensor::xTensor4<>> tensor4s;
   std::vector<std::string> strings;
   std::vector<xTable> tables;
   std::vector<xPieceWiseLinear> piecewiselinears;
};

class xTensorsValuePtr : public xfem::xValue<tensorsPtr_t>
{
  public:
   xTensorsValuePtr();

   tensorsPtr_t getVal() const override;
   void setVal(tensorsPtr_t v) override;
   std::ostream& printVal(std::ostream& o) const override;

  private:
   tensorsPtr_t value;
};

}  // namespace xfem

#endif
