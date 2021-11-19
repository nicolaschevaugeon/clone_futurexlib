/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/

#include "xTensors.h"

#include <iostream>

namespace xfem
{
using std::cerr;
using std::cout;
using std::endl;
using std::string;

void xTensors::read(const string& filename)
{
   if (filename[0] == '*') return;
   FILE* fp = fopen(filename.c_str(), "r");
   if (fp == nullptr)
   {
      cerr << "The data file " << filename << " cannot be opened\n";
      assert(1 == 0);
   }
   char key_char[256];
   string key;
   int i;
   // int tmp;
   char c;

   while ((i = fscanf(fp, "%[1234567890A-Z_]", key_char)) != EOF)
   {
      if (i == 0)
      {
         c = fgetc(fp);
         // if ( c != '#' && c !='\n' && c != ' ')
         // {
         //   fprintf(stderr, "The following character is not known : %c\n", c);
         //   assert(1 == 0);
         // }
         if (c == '#') fscanf(fp, "%*[^\n] \n");
      }
      else
      {
         key = key_char;
         read(fp, key);
      }
   }
   fclose(fp);
}

void xTensors::read(FILE* fp, const string& name)
{
   // cout << "reading " << name << endl;
   if (sgn->exist_scalar(name))
      read_scalar(fp, scalar(name));
   else if (sgn->exist_vector(name))
      read_vector(fp, vector(name));
   // else if (exist_tensor2(name)) read_tensor2(fp, sgn->tensor2(name));
   else if (sgn->exist_table(name))
      read_table(fp, table(name));
   else if (sgn->exist_string(name))
      read_string(fp, astring(name));
   else if (sgn->exist_piecewiselinear(name))
      read_piecewiselinear(fp, piecewiselinear(name));

   // else {
   //   cerr << "the name " << name << " is not known by the signature of Tensors" << std::endl;
   //   assert(0);
   // }
   return;
}
void xTensors::read_scalar(FILE* fp, double& val)
{
   //  int tmp=fscanf(fp, "%*[\n= ] %lf", &val);
   fscanf(fp, "%*[\n= ] %lf", &val);
}

void xTensors::read_vector(FILE* fp, xtensor::xVector<>& val)
{
   int dim;
   // int tmp;
   // tmp=fscanf(fp, "%*[\n= ] %d", &dim);
   fscanf(fp, "%*[\n= ] %d", &dim);
   assert((dim <= 3) && (dim > 0));
   for (int i = 0; i < dim; i++)
   {
      //  tmp=fscanf(fp, "%*[\n= ] %lf", &val[i]);
      fscanf(fp, "%*[\n= ] %lf", &val[i]);
   }
}

void xTensors::read_string(FILE* fp, string& val)
{
   char name_char[256];
   // int tmp = fscanf(fp, "%*[\n= ] %s[a-z]", name_char );
   fscanf(fp, "%*[\n= ] %s[a-z]", name_char);
   val = name_char;
}
void xTensors::read_table(FILE* fp, xTable& val)
{
   int dim;
   std::vector<int> dimensions;
   //  int tmp=fscanf(fp, "%*[\n= ] %d", &dim);
   fscanf(fp, "%*[\n= ] %d", &dim);
   //  cout << dim << " ";
   dimensions.resize(dim);
   for (int i = 0; i < dim; i++)
   {
      //    tmp=fscanf(fp, "%*[\n= ] %d", &dimensions[i]);
      fscanf(fp, "%*[\n= ] %d", &dimensions[i]);
   }
   if (dim > 0)
   {
      double factor;
      // tmp=fscanf(fp, "%*[\n= ] %lf", &factor);
      fscanf(fp, "%*[\n= ] %lf", &factor);
      val.setdimnd(dimensions);
      int tot = 1;
      for (int i = 0; i < dim; ++i) tot *= dimensions[i];
      for (int i = 0; i < tot; ++i)
      {
         double v;
         // tmp=fscanf(fp, "%*[\n= ] %lf", &v);
         fscanf(fp, "%*[\n= ] %lf", &v);
         val.set(i, v * factor);
      }
   }
}

void xTensors::read_piecewiselinear(FILE* fp, xPieceWiseLinear& val)
{
   val.values.clear();
   int dim;
   // int tmp=fscanf(fp, "%*[\n= ] %d", &dim);
   fscanf(fp, "%*[\n= ] %d", &dim);
   for (int i = 0; i < dim; ++i)
   {
      double t, v;
      //     tmp=fscanf(fp, "%*[\n= ] %lf", &t);
      // tmp=fscanf(fp, "%*[\n= ] %lf", &v);
      fscanf(fp, "%*[\n= ] %lf", &t);
      fscanf(fp, "%*[\n= ] %lf", &v);
      val.values[t] = v;
   }
   return;
}

std::ostream& operator<<(std::ostream& s, const xTensors& t)
{
   if (int i = t.sgn->size_scalar())
   {
      s << i << " scalar values " << endl;
      xTensorsSignature::const_iterator it = t.sgn->begin_scalar();
      xTensorsSignature::const_iterator itend = t.sgn->end_scalar();
      for (; it != itend; ++it)
      {
         s << it->first << endl;
         s << (t.scalars)[it->second] << endl;
      }
   }
   if (int i = t.sgn->size_vector())
   {
      s << i << " vector values " << endl;
      xTensorsSignature::const_iterator it = t.sgn->begin_vector();
      xTensorsSignature::const_iterator itend = t.sgn->end_vector();
      for (; it != itend; ++it)
      {
         s << it->first << endl;
         s << (t.vectors)[it->second] << endl;
      }
   }
   if (int i = t.sgn->size_tensor2())
   {
      s << i << " tensor2 values " << endl;
      xTensorsSignature::const_iterator it = t.sgn->begin_tensor2();
      xTensorsSignature::const_iterator itend = t.sgn->end_tensor2();
      for (; it != itend; ++it)
      {
         s << it->first << endl;
         s << (t.tensor2s)[it->second] << endl;
      }
   }
   if (int i = t.sgn->size_string())
   {
      s << i << " string values " << endl;
      xTensorsSignature::const_iterator it = t.sgn->begin_string();
      xTensorsSignature::const_iterator itend = t.sgn->end_string();
      for (; it != itend; ++it)
      {
         s << it->first << endl;
         s << (t.strings)[it->second] << endl;
      }
   }
   return s;
}

xTensorsValuePtr::xTensorsValuePtr() : xfem::xValue<tensorsPtr_t>()
{
   tensorsPtr_t Foo(new xTensors);
   value = Foo;
}

tensorsPtr_t xTensorsValuePtr::getVal() const { return value; }

void xTensorsValuePtr::setVal(tensorsPtr_t v) { value = v; }

std::ostream& xTensorsValuePtr::printVal(std::ostream& o) const
{
   o << "xValueTensorsPtr::printVal not implemented" << endl;
   return o;
}

void xPieceWiseLinear::read(std::istream& f)
{
   values.clear();
   int n;
   f >> n;
   for (int i = 0; i < n; ++i)
   {
      double t, v;
      f >> t;
      f >> v;
      values[t] = v;
   }
   return;
}

xPieceWiseLinear::xPieceWiseLinear(const std::map<double, double>& val) : values(val)
{
   if (values.size() < 2) throw;
}

double xPieceWiseLinear::operator()(double T) const
{
   std::map<double, double>::const_iterator it = values.lower_bound(T);
   std::map<double, double>::const_iterator it2 = values.upper_bound(T);

   if (it == it2 && (it == values.begin() || it == values.end()))
   {
      // std::cout << "warning interpolating out of  bound : "<< T  << std::endl;
      if (it == values.begin())
      {
         double vb = (*it).second;
         if (fabs(it->first - T) > 1.e-5)
            std::cout << "warning interpolating out of  bound : " << T << " " << it->first << std::endl;
         return vb;
      }
      if (it == values.end())
      {
         double ve = (*(--it)).second;
         if (fabs(it->first - T) > 1.e-5)
            std::cout << "warning interpolating out of  bound : " << T << " " << it->first << std::endl;
         return ve;
      }
   }
   if (it == values.begin()) return (*it).second;

   std::pair<double, double> after = *it;
   std::pair<double, double> before = *(--it);

   double alpha = (T - before.first) / (after.first - before.first);
   return before.second * (1. - alpha) + after.second * alpha;
}

}  // namespace xfem
