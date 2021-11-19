/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/

#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <map>
#include <string>

// aomd
#include "mEdge.h"
#include "mEntityContainer.h"
#include "mFace.h"
#include "mIterator.h"
#include "mTet.h"
#include "mVertex.h"
#include "xRegion.h"

// xfem
#include "xDebug.h"
#include "xElement.h"
#include "xVectorField.h"

using std::cout;
using std::endl;
using std::ostream;
using namespace AOMD;

namespace xfem
{
xtensor::xVector<>& xVectorField::operator()(mEntity* e) { return f[e]; }

void xVectorField::setSupport(const xRegion& m, const xtensor::xVector<>& val)
{
   support = m;
   clear();
   for (xIter it = support.begin(0); it != support.end(0); ++it)
   {
      mVertex* v = (mVertex*)*it;
      f[v] = val;
   }
}

const xtensor::xVector<>& xVectorField::operator()(mEntity* e) const { return f.find(e)->second; }

// mMesh* xVectorField::getMesh() const {return mesh;}

// void xVectorField::setMesh(mMesh* m) {mesh = m; clear(); }

void xVectorField::clear() { f.clear(); }

xVectorField::xVectorField(const xRegion& s, const xtensor::xVector<>& val) : support(s)
{
   for (xIter it = support.begin(0); it != support.end(0); ++it)
   {
      // for(mIteratorPtr it(mesh->getLevelIterator(0)); !it->done(); it->next() ) {
      mVertex* v = (mVertex*)*it;
      f[v] = val;
   }
}

std::vector<xtensor::xVector<>> xVectorField::getVals(mEntity* e) const
{
   std::vector<xtensor::xVector<>> vals;
   for (int i = 0; i < e->size(0); i++)
   {
      mVertex* v = (mVertex*)e->get(0, i);
      vals.push_back(f.find(v)->second);
   }
   return vals;
}

// anthony
void xVectorField::exportGmsh(const std::string& name) { exportGmsh(name, support); }

void xVectorField::exportGmsh(const std::string& name, xRegion sub)
{
   printf("in xVectorField::exportGmsh\n");
   std::string filename = name + ".pos";
   FILE* fout = fopen(filename.c_str(), "w");
   fprintf(fout, "View \"%s\" { \n", name.c_str());
   for (mEntity* e : sub)
   {
      std::vector<xtensor::xPoint> p(e->size(0));
      std::vector<xtensor::xVector<>> lev = getVals(e);
      for (int j = 0; j < e->size(0); j++)
      {
         mVertex* v = static_cast<mVertex*>(e->get(0, j));
         p[j] = v->point();
      }
      switch (e->getLevel())
      {
         case 1:
            fprintf(fout, "VL(%12.5e, %12.5e, %12.5e,", p[0](0), p[0](1), p[0](2));
            fprintf(fout, "%12.5e, %12.5e, %12.5e)", p[1](0), p[1](1), p[1](2));
            fprintf(fout, "{%12.5e, %12.5e, %12.5e, %12.5e, %12.5e, %12.5e};\n", lev[0](0), lev[0](1), lev[0](2), lev[1](0),
                    lev[1](1), lev[1](2));
            break;
         case 2:
            fprintf(fout, "VT(%12.5e, %12.5e, %12.5e,", p[0](0), p[0](1), p[0](2));
            fprintf(fout, "%12.5e, %12.5e, %12.5e,", p[1](0), p[1](1), p[1](2));
            fprintf(fout, "%12.5e, %12.5e, %12.5e)", p[2](0), p[2](1), p[2](2));
            fprintf(fout, "{%12.5e, %12.5e, %12.5e, %12.5e, %12.5e, %12.5e, %12.5e, %12.5e, %12.5e};\n", lev[0](0), lev[0](1),
                    lev[0](2), lev[1](0), lev[1](1), lev[1](2), lev[2](0), lev[2](1), lev[2](2));
            break;
         case 3:
            fprintf(fout, "VS(%12.5e, %12.5e, %12.5e,", p[0](0), p[0](1), p[0](2));
            fprintf(fout, "%12.5e, %12.5e, %12.5e,", p[1](0), p[1](1), p[1](2));
            fprintf(fout, "%12.5e, %12.5e, %12.5e,", p[2](0), p[2](1), p[2](2));
            fprintf(fout, "%12.5e, %12.5e, %12.5e)", p[3](0), p[3](1), p[3](2));
            fprintf(fout,
                    "{%12.5e, %12.5e, %12.5e, %12.5e, %12.5e, %12.5e, %12.5e, "
                    "%12.5e, %12.5e, %12.5e, %12.5e, %12.5e};\n",
                    lev[0](0), lev[0](1), lev[0](2), lev[1](0), lev[1](1), lev[1](2), lev[2](0), lev[2](1), lev[2](2), lev[3](0),
                    lev[3](1), lev[3](2));
            break;
         default:
            assert(0);
            break;
      }
   }
   fprintf(fout, "};\n");
   fclose(fout);
   return;
}

void xVectorField::printDebug() const
{
   for (mEntity* e : support.range(0))
   {
      mVertex* v = static_cast<mVertex*>(e);
      xtensor::xPoint p = v->point();
      std::cout << p << f.find(v)->second << endl;
   }
}

void xVectorField::exportMatlab(ostream& fout, const std::string& field_name, int level)
{
   fout << field_name << "= [ ";
   for (mEntity* e : support.range(level)) fout << f.find(e)->second << ";\n";
   fout << "];\n";
   return;
}

void xVectorField::importGmsh(const std::string& filename)
{
   FILE* fin = fopen(filename.c_str(), "r");
   assert(fin != nullptr);
   // int ok;
   // char temp[256],temp1[256];
   //  fscanf(fin, "View %s[\"0-9a-z_] %s[{\n]", temp, temp1);
   //  ok = fscanf(fin, "%*[View \"] %*[abcdefghijklmnopqrstuvwxyz0123456789_] %*[\"\n={ ]");
   fscanf(fin, "%*[View \"] %*[abcdefghijklmnopqrstuvwxyz0123456789_] %*[\"\n={ ]");
   // printf("read view %s \n", temp);
   double a0, a1, a2, b0, b1, b2;
   xtensor::xVector<> a, b;
   for (mEntity* e : support)
   {
      mEntity* v0 = e->get(0, 0);
      mEntity* v1 = e->get(0, 1);
      switch (e->getLevel())
      {
         case 1:
            //	ok = fscanf(fin, "VL(%*f %*[,] %*f %*[,] %*f %*[,] %*f %*[,] %*f %*[,] %*f %*[){] %lf %*[\n=, ] %lf %*[\n=, ] %lf
            //%*[\n=, ] %lf %*[\n=, ] %lf %*[\n=, ] %lf %*[};\n=, ]",
            fscanf(fin,
                   "VL(%*f %*[,] %*f %*[,] %*f %*[,] %*f %*[,] %*f %*[,] %*f %*[){] %lf %*[\n=, ] %lf %*[\n=, ] %lf %*[\n=, ] "
                   "%lf %*[\n=, ] %lf %*[\n=, ] %lf %*[};\n=, ]",
                   &a0, &a1, &a2, &b0, &b1, &b2);
            a(0) = a0;
            a(1) = a1;
            a(2) = a2;
            b(0) = b0;
            b(1) = b1;
            b(2) = b2;
            cout << "speed read " << a << " " << b << endl;
            f[v0] = a;
            f[v1] = b;
            break;
         default:
            assert(0);
            break;
      }
   }
   // ok = fscanf(fin, "};\n");
   fscanf(fin, "};\n");
   fclose(fin);
   return;
}

xRegion xVectorField::getSupport() const { return support; }

xtensor::xVector<> xVectorField::getVal(mEntity* e, const xtensor::xPoint& uvw) const
{
   // const bool debug = xdebug_flag;
   xElement elem(e);
   elem.setUvw(uvw);
   std::vector<xtensor::xVector<>> vals = getVals(e);
   xtensor::xVector<> out = elem.getInterpoVec(vals);
   return out;
}

xtensor::xTensor2<double> xVectorField::getGrad(mEntity* e) const
{
   // const bool debug = xdebug_flag;
   xElement elem(e);
   elem.setUvw(xtensor::xPoint(0., 0., 0.));
   std::vector<xtensor::xVector<>> vals = getVals(e);
   xtensor::xTensor2<> out;
   elem.getGradInterpoVec(vals, out);
   return out;
}

}  // namespace xfem
