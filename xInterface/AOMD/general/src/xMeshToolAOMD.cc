/*
   This file is a part of eXlibris C++ Library
   under the GNU Lesser General Public License.
   See the NOTICE.md & LICENSE.md files for terms
   and conditions.
*/

#include "mEdge.h"
#include "mFace.h"
#include "mMesh.h"
#include "mRegion.h"
#include "mTet.h"
// #include "xDebug.h"
//#include "xMesh.h"

#include "xAOMDEntityUtil.h"
#include "xMesh.h"
#include "xMeshToolAOMD.h"

namespace xinterface
{
namespace aomd
{
using namespace AOMD;

void modifyAllState(AOMD::mMesh& mmesh)
{
   mmesh.modifyState(3, 2, true);
   mmesh.modifyState(3, 1, true);
   mmesh.modifyState(3, 0, true);
   mmesh.modifyState(2, 1, true);
   mmesh.modifyState(2, 0, true);
   mmesh.modifyState(1, 0, true);
   mmesh.modifyState(0, 1, true);
   mmesh.modifyState(0, 2, true);
   mmesh.modifyState(0, 3, true);
   mmesh.modifyState(1, 2, true);
   mmesh.modifyState(1, 3, true);
   mmesh.modifyState(2, 3, true);
}

void modifyAllStateFalse(AOMD::mMesh& mmesh)
{
   mmesh.modifyState(3, 2, false);
   mmesh.modifyState(3, 1, false);
   mmesh.modifyState(2, 1, false);
   mmesh.modifyState(2, 0, false);
   mmesh.modifyState(0, 1, false);
   mmesh.modifyState(0, 2, false);
   mmesh.modifyState(1, 2, false);
   mmesh.modifyState(1, 3, false);
   mmesh.modifyState(2, 3, false);
   mmesh.modifyState(2, 1, false);
   mmesh.modifyState(2, 0, false);
   mmesh.modifyState(1, 0, false);
   mmesh.modifyState(0, 1, false);
   mmesh.modifyState(0, 2, false);
   mmesh.modifyState(1, 2, false);
}

void xCopyMesh(const xfem::xMesh& in, xfem::xMesh& out)
{
   xAttachedDataManagerAOMD<AOMD::mEntity*> is_copied_to;
   xAttachedDataManagerAOMD<AOMD::mEntity*> is_the_copy_of;
   xCopyMesh(in, out, is_copied_to, is_the_copy_of);
}

void xCopyMesh(const xfem::xMesh& in, xfem::xMesh& out, xAttachedDataManagerAOMD<AOMD::mEntity*>& is_copied_to,
               xAttachedDataManagerAOMD<AOMD::mEntity*>& is_the_copy_of)
{
//    const bool debug = xfem::xdebug_flag;
    const bool debug = false;
   if (&in == &out)
   {
      std::cout << "warnning : copying a mesh onto itself" << std::endl;
      return;
   }
   AOMD::mMesh& mout = out.getMesh();

   for (mEntity* pvin : in.range(0))
   {
      mVertex& vin = static_cast<mVertex&>(*pvin);
      mVertex& vout = *mout.createVertex(vin.point(), vin.getClassification());
      is_copied_to.setData(vin) = &vout;
      is_the_copy_of.setData(vout) = &vin;
      if (debug)
      {
         cout << " vin is : " << &vin << endl;
         vin.print();
         cout << " vout is : " << &vout << endl;
         vout.print();
      }
   }
   for (mEntity* pein : in.range(1))
   {
      mEdge& ein = static_cast<mEdge&>(*pein);
      mVertex& v0in = static_cast<mVertex&>(*ein.get(0, 0));
      mVertex& v1in = static_cast<mVertex&>(*ein.get(0, 1));
      mVertex& v0out = static_cast<mVertex&>(*is_copied_to.at(v0in));
      mVertex& v1out = static_cast<mVertex&>(*is_copied_to.at(v1in));
      mEdge& eout = *mout.createEdge(&v0out, &v1out, ein.getClassification());
      is_copied_to.setData(ein) = &eout;
      is_the_copy_of.setData(eout) = &ein;
   }
   for (mEntity* pfin : in.range(2))
   {
      mFace& fin = static_cast<mFace&>(*pfin);
      switch (fin.size(0))
      {
         case 3:
         {
            mVertex& v0in = static_cast<mVertex&>(*fin.get(0, 0));
            mVertex& v1in = static_cast<mVertex&>(*fin.get(0, 1));
            mVertex& v2in = static_cast<mVertex&>(*fin.get(0, 2));
            mVertex& v0out = static_cast<mVertex&>(*is_copied_to.at(v0in));
            mVertex& v1out = static_cast<mVertex&>(*is_copied_to.at(v1in));
            mVertex& v2out = static_cast<mVertex&>(*is_copied_to.at(v2in));
            mFace& fout = *mout.createFaceWithVertices(&v0out, &v1out, &v2out, fin.getClassification());
            is_copied_to.setData(fin) = &fout;
            is_the_copy_of.setData(fout) = &fin;
            break;
         }
         case 4:
         {
            mVertex& v0in = static_cast<mVertex&>(*fin.get(0, 0));
            mVertex& v1in = static_cast<mVertex&>(*fin.get(0, 1));
            mVertex& v2in = static_cast<mVertex&>(*fin.get(0, 2));
            mVertex& v3in = static_cast<mVertex&>(*fin.get(0, 3));
            mVertex& v0out = static_cast<mVertex&>(*is_copied_to.at(v0in));
            mVertex& v1out = static_cast<mVertex&>(*is_copied_to.at(v1in));
            mVertex& v2out = static_cast<mVertex&>(*is_copied_to.at(v2in));
            mVertex& v3out = static_cast<mVertex&>(*is_copied_to.at(v3in));
            mFace& fout = *mout.createFaceWithVertices(&v0out, &v1out, &v2out, &v3out, fin.getClassification());
            is_copied_to.setData(fin) = &fout;
            is_the_copy_of.setData(fout) = &fin;
            break;
         }
         default:
         {
            cout << "xCopyMesh, 2D : copying face with size(0) != 3 or 4  not coded yet  " << __FILE__ << ":" << __LINE__ << endl;
            throw 1;
         }
      }
   }
   for (mEntity* prin : in.range(3))
   {
      mRegion& rin = static_cast<mRegion&>(*prin);
      switch (rin.size(0))
      {
         case 4:
         {
            mVertex& v0in = static_cast<mVertex&>(*rin.get(0, 0));
            mVertex& v1in = static_cast<mVertex&>(*rin.get(0, 1));
            mVertex& v2in = static_cast<mVertex&>(*rin.get(0, 2));
            mVertex& v3in = static_cast<mVertex&>(*rin.get(0, 3));
            mVertex& v0out = static_cast<mVertex&>(*is_copied_to.at(v0in));
            mVertex& v1out = static_cast<mVertex&>(*is_copied_to.at(v1in));
            mVertex& v2out = static_cast<mVertex&>(*is_copied_to.at(v2in));
            mVertex& v3out = static_cast<mVertex&>(*is_copied_to.at(v3in));
            mRegion& rout = *mout.createTetWithVertices(&v0out, &v1out, &v2out, &v3out, rin.getClassification());
            is_copied_to.setData(rin) = &rout;
            is_the_copy_of.setData(rout) = &rin;
            break;
         }
         default:
         {
            cout << "xCopyMesh, 3D : copying region with size(0) != 4 not coded yet  " << __FILE__ << ":" << __LINE__ << endl;
            throw 1;
         }
      }
   }
   return;
}

void ModifyAllState(AOMD::mMesh& mesh)
{
   mesh.modifyState(3, 2, true);
   mesh.modifyState(3, 1, true);
   mesh.modifyState(3, 0, true);
   mesh.modifyState(2, 1, true);
   mesh.modifyState(2, 0, true);
   mesh.modifyState(1, 0, true);
   mesh.modifyState(0, 1, true);
   mesh.modifyState(0, 2, true);
   mesh.modifyState(0, 3, true);
   mesh.modifyState(1, 2, true);
   mesh.modifyState(1, 3, true);
   mesh.modifyState(2, 3, true);
}

void Simplexify(AOMD::mMesh& mesh)
{
   std::cout << "=========Begin of Simplexify()=========" << std::endl;
   //  PrintMesh();

   list<mEntity*> toTriangularize;
   for (auto it = mesh.begin(2); it != mesh.end(2); ++it)
   {
      mFace* _F = (mFace*)(*it);
      if (_F->getType() == mEntity::QUAD)
      {
         toTriangularize.push_back(_F);
      }
   }
   for (auto it = toTriangularize.begin(); it != toTriangularize.end(); ++it)
   {
      mFace* _Q = (mFace*)(*it);
      QuadToTri(mesh, _Q);
      DeleteEntity(mesh, _Q);
   }

   list<mEntity*> toTetrize;
   for (auto it = mesh.begin(3); it != mesh.end(3); ++it)
   {
      mRegion* _R = (mRegion*)(*it);

      switch (_R->getType())
      {
         case mEntity::HEX:
         {
            toTetrize.push_back(_R);
            break;
         }
         case mEntity::PRISM:
         {
            toTetrize.push_back(_R);
            break;
         }
         case mEntity::PYRAMID:
         {
            toTetrize.push_back(_R);
            break;
         }
         default:
         {
         }
      }
   }
   for (auto it = toTetrize.begin(); it != toTetrize.end(); ++it)
   {
      mHex* _H = (mHex*)(*it);
      switch (_H->getType())
      {
         case mEntity::HEX:
         {
            HexToTet(mesh, _H);
            DeleteEntity(mesh, _H);
            break;
         }
         case mEntity::PRISM:
         {
            //	    PrismToTet( _H ) ;
            //	    deleteEntity( _H );
            try
            {
               throw logic_error("ERROR: PRISM element are not treated yet in  Xtest::SimplexifyMesh.");
            }
            catch (exception& e)
            {
               cerr << e.what() << endl;
               exit(0);
            }
            break;
         }
         case mEntity::PYRAMID:
         {
            //	    PyramidToTet( _H ) ;
            //	    deleteEntity( _H );
            try
            {
               throw logic_error("ERROR: PYRAMID element are not treated yet in  Xtest::SimplexifyMesh.");
            }
            catch (exception& e)
            {
               cerr << e.what() << endl;
               exit(0);
            }

            break;
         }
         default:
         {
         }
      }
   }

   std::cout << "=========End of Simplexify()=========" << std::endl;
   //      PrintMesh();
};

void DeleteEntity(AOMD::mMesh& mesh, mEntity* pent) { mesh.DEL(pent); };

void QuadToTri(AOMD::mMesh& mesh, mFace* pent)
{
   vector<int> connect;
   connect.clear();
   for (int i = 0; i < 4; ++i)
   {
      connect.push_back((pent->get(0, i))->getId());
   }
   int imin = distance(connect.begin(), min_element(connect.begin(), connect.end()));

   vector<int> new_connect;
   new_connect.reserve(4);
   new_connect.clear();
   for (int i = 0; i < 4; ++i)
   {
      new_connect[(i - imin + 4) % 4] = connect[i];
   }

   mesh.createFaceWithVertices(new_connect[0], new_connect[2], new_connect[3], pent->getClassification());

   mesh.createFaceWithVertices(new_connect[0], new_connect[1], new_connect[2], pent->getClassification());
};

void HexToTet(mMesh& mesh, mHex* pent)
{
   vector<int> connect;
   connect.reserve(8);
   connect.clear();
   for (int i = 0; i < 8; ++i)
   {
      connect.push_back((pent->get(0, i))->getId());
   }
   int imin = distance(connect.begin(), min_element(connect.begin(), connect.end()));
   vector<int> new_connect;
   new_connect.reserve(8);
   new_connect.clear();
   if (imin < 4)
   {
      for (int i = 0; i < 4; ++i)
      {
         new_connect[(i - imin + 4) % 4] = connect[i];
      }
      for (int i = 4; i < 8; ++i)
      {
         new_connect[(i - imin + 4) % 4 + 4] = connect[i];
      }
   }
   if (3 < imin)
   {
      for (int i = 4; i < 8; ++i)
      {
         new_connect[(-i + imin + 4) % 4] = connect[i];
      }
      for (int i = 0; i < 4; ++i)
      {
         new_connect[(-i + imin + 4) % 4 + 4] = connect[i];
      }
   }
   int N0 = new_connect[0];
   int N1 = new_connect[1];
   int N2 = new_connect[2];
   int N3 = new_connect[3];
   int N4 = new_connect[4];
   int N5 = new_connect[5];
   int N6 = new_connect[6];
   int N7 = new_connect[7];

   int selection[] = {2, 3, 6, 7};
   vector<int> index;
   index.assign(selection, selection + 4);
   imin = IndexOfMinAmong(new_connect, index);
   if (imin == 6 || imin == 3)
   {
      CreateTetWithPrismVertices(mesh, N0, N5, N1, N3, N6, N2, pent->getClassification());
      CreateTetWithPrismVertices(mesh, N5, N0, N4, N6, N3, N7, pent->getClassification());
   }
   else
   {
      int selection[] = {4, 5, 6, 7};
      vector<int> index;
      index.assign(selection, selection + 4);
      imin = IndexOfMinAmong(new_connect, index);
      if (imin == 4 || imin == 6)
      {
         CreateTetWithPrismVertices(mesh, N4, N6, N5, N0, N2, N1, pent->getClassification());
         CreateTetWithPrismVertices(mesh, N4, N7, N6, N0, N3, N2, pent->getClassification());
      }
      else
      {
         int selection[] = {5, 6, 2, 1};
         vector<int> index;
         index.assign(selection, selection + 4);
         imin = IndexOfMinAmong(new_connect, index);
         if (imin == 1 || imin == 6)
         {
            CreateTetWithPrismVertices(mesh, N5, N6, N1, N4, N7, N0, pent->getClassification());
            CreateTetWithPrismVertices(mesh, N1, N6, N2, N0, N7, N3, pent->getClassification());
         }
         else
         {
            // xEntity* tet;

            mesh.createTetWithVertices(N5, N4, N7, N0, pent->getClassification());
            mesh.createTetWithVertices(N0, N7, N5, N2, pent->getClassification());
            mesh.createTetWithVertices(N5, N0, N2, N1, pent->getClassification());
            mesh.createTetWithVertices(N0, N7, N2, N3, pent->getClassification());
            mesh.createTetWithVertices(N5, N2, N7, N6, pent->getClassification());
         }
      }
   }
};

void CreateTetWithPrismVertices(mMesh& mesh, int N0, int N1, int N2, int N3, int N4, int N5, pGEntity classif)
{
   vector<int> connect;
   connect.reserve(6);
   connect.clear();
   connect.push_back(N0);
   connect.push_back(N1);
   connect.push_back(N2);
   connect.push_back(N3);
   connect.push_back(N4);
   connect.push_back(N5);
   int imin = min_element(connect.begin(), connect.end()) - connect.begin();
   vector<int> new_connect;
   new_connect.reserve(6);
   new_connect.clear();
   if (imin < 3)
   {
      for (int i = 0; i < 3; ++i)
      {
         new_connect[(i - imin + 3) % 3] = connect[i];
      }
      for (int i = 3; i < 6; ++i)
      {
         new_connect[(i - imin + 3) % 3 + 3] = connect[i];
      }
   }
   if (2 < imin)
   {
      for (int i = 3; i < 6; ++i)
      {
         new_connect[(-i + imin + 3) % 3] = connect[i];
      }
      for (int i = 0; i < 3; ++i)
      {
         new_connect[(-i + imin + 3) % 3 + 3] = connect[i];
      }
   }
   N0 = new_connect[0];
   N1 = new_connect[1];
   N2 = new_connect[2];
   N3 = new_connect[3];
   N4 = new_connect[4];
   N5 = new_connect[5];

   mesh.createTetWithVertices(N3, N5, N4, N0, classif);

   CreateTetWithPyramidVertices(mesh, N4, N5, N2, N1, N0, classif);
};

void CreateTetWithPyramidVertices(mMesh& mesh, int N0, int N1, int N2, int N3, int N4, pGEntity classif)
{
   vector<int> connect;
   connect.reserve(5);
   connect.clear();
   connect.push_back(N0);
   connect.push_back(N1);
   connect.push_back(N2);
   connect.push_back(N3);
   connect.push_back(N4);
   int selection[] = {0, 1, 2, 3};
   vector<int> index;
   index.assign(selection, selection + 4);
   int imin = IndexOfMinAmong(connect, index);
   vector<int> new_connect;
   new_connect.reserve(5);
   new_connect.clear();
   for (int i = 0; i < 4; ++i)
   {
      new_connect[(i - imin + 4) % 4] = connect[i];
   }
   N0 = new_connect[0];
   N1 = new_connect[1];
   N2 = new_connect[2];
   N3 = new_connect[3];
   N4 = connect[4];

   mesh.createTetWithVertices(N0, N2, N3, N4, classif);
   mesh.createTetWithVertices(N0, N1, N2, N4, classif);
};

int IndexOfMinAmong(vector<int>& iValue, vector<int>& index)
{
   int imin = 0;
   int iValMin = iValue[index[imin]];
   for (size_t i = 1; i < index.size(); ++i)
   {
      if (iValue[index[i]] < iValMin)
      {
         imin = i;
         iValMin = iValue[index[imin]];
      }
   }
   return index[imin];
};
/*
void PrintMesh( mMesh& mesh )
{
    MeshQueryInterface query( mesh );
    query.printMesh();
};
*/
void printMesh(std::ostream& os, xfem::xMesh& mesh, bool deep)
{
   os << std::endl << "Export mesh" << std::endl;
   for (int i = 0; i < 4; ++i)
      for (AOMD::mEntity* pe : mesh.range(i)) printEntity(os, pe, deep);
}
/*
    void readMesh(const string & filename) {
        std::ifstream meshfile(filename.c_str());
        if (!meshfile.is_open())
{
    cout << "can't open meshfile "<< filename << " in readMesh " << __FILE__ << " " << __LINE__ << endl;
    throw;
}
        meshfile.close();
        AOMD::AOMD_Util::Instance()->import(filename.c_str(), &mesh);
        cout<<" begin modifyAllState()"<<endl;
        modifyAllState();
        cout<<" end modifyAllState()"<<endl;
    }

*/

/*
    void getPartition(const aomdMeshQueryInterface& in, const aomdMeshQueryInterface& out) {
        auto it = in.begin(in.dim());
        auto itend = in.end(in.dim());

        while (it != itend)
{

                //  void * pointerToData = (*it).getAttachedPointer(xfem::xMesh::get_partition_tag());
                if (! (*it).hasData(xfem::xMesh::get_partition_tag()) )
                {
                        std::cout << "No data attached to entity "; ( *it ).print(); std::cout << " with tag xMesh::partition_tag"
                        std::cout<< std::endl;
                        throw;
                }
                //  xfem::xMesh* subpat = static_cast<xfem::xMesh*>(pointerToData);
                xfem::xMesh* subpat = (*it).getData<xfem::xMesh*>( xfem::xMesh::get_partition_tag() );
                xCopyMesh(  dynamic_cast<const aomdMeshQueryInterface&>( subpat->getQuery() ) , out);

    ++it;
}
        return;
    }

*/

/*
    void  xCopyMesh(const aomdMeshQueryInterface& in, const aomdMeshQueryInterface& out, bool clean_tag_on_source_mesh)
    {
        const bool debug = false;
         if (  &(in.mesh) == &(out.mesh) )
{
    std::cout << "warnning : copying a mesh onto itself" << std::endl;
    return;
}
        //  should check first if in and out are not  0 ....
        // xMesh & in = *pin;
        // xMesh & out = *pout;

        auto it = in.beginVertex();
        auto ite = in.endVertex();
        while (it != ite)
{
    xVertex vin = *it ;
    if (debug) cout << " vin is : " ; vin.print(); cout << endl;
    xtensor::xPoint pt=vin.point();
    AOMD::pGEntity pge = (any_cast<AOMD::mEntity*>( vin.getEntityIdentifier()) )->getClassification() ;
    xVertex vout( out.mesh.createVertex(pt(0),pt(1),pt(2), pge), out );
    if (debug) cout << " vout is : " ; vout.print(); cout << endl;
        vin.attachData<xEntity>(xfem::xMesh::get_has_a_copy_tag(), vout);
        vout.attachData<xEntity>(xfem::xMesh::get_is_the_copy_of_tag(),vin);
    ++it;
}

        it = in.beginEdge();
        ite = in.endEdge();
        while (it != ite)
{
    xEdge ein =  *it ;
    xVertex v0in = ein.getVertex(0);
    xVertex v1in = ein.getVertex(1);
        xVertex v0out = v0in.getData<xEntity>(xfem::xMesh::get_has_a_copy_tag());
        xVertex v1out = v1in.getData<xEntity>(xfem::xMesh::get_was_created_by_tag());
    AOMD::pGEntity pge = (any_cast<AOMD::mEntity*>( ein.getEntityIdentifier()) )->getClassification() ;
    xEdge eout( out.mesh.createEdge( any_cast<AOMD::mVertex*>(v0out.getEntityIdentifier()),
                     any_cast<AOMD::mVertex*>(v1out.getEntityIdentifier()),
                     pge),
                out);
        ein.attachData<xEntity>(xfem::xMesh::get_was_created_by_tag(),eout);
        eout.attachData<xEntity>(xfem::xMesh::get_is_the_copy_of_tag(),ein);
    ++it;
}

        it = in.beginFace();
        ite = in.endFace();
        while (it != ite)
{
    xFace fin = *it ;
    int nnode = fin.size(0);
    switch (nnode)
        {
        case 3 : {
            xVertex v0in =  fin.getVertex(0);
            xVertex v1in =  fin.getVertex(1);
            xVertex v2in =  fin.getVertex(2);
                xVertex v0out = v0in.getData<xEntity>(xfem::xMesh::get_was_created_by_tag());
                xVertex v1out = v1in.getData<xEntity>(xfem::xMesh::get_was_created_by_tag());
                xVertex v2out = v2in.getData<xEntity>(xfem::xMesh::get_was_created_by_tag());
            AOMD::pGEntity pge = (any_cast<AOMD::mEntity*>( fin.getEntityIdentifier()) )->getClassification() ;
            xFace fout( out.mesh.createFaceWithVertices( any_cast<AOMD::mVertex*>(v0out.getEntityIdentifier()),
                             any_cast<AOMD::mVertex*>(v1out.getEntityIdentifier()),
                             any_cast<AOMD::mVertex*>(v2out.getEntityIdentifier()),
                             pge),
            out);
                fin.attachData<xEntity>(xfem::xMesh::get_was_created_by_tag(),fout);
                fout.attachData<xEntity>(xfem::xMesh::get_is_the_copy_of_tag(),fin);
            break;
        }
        case 4 : {
            cout << "xCopyMesh , quadrangles : PAS ENCORE PROGRAMME ! " << endl;
            throw 1;
        }
        default : {
            throw 1;
        }
        }
    ++it;
}

        it = in.beginSolid();
        ite = in.endSolid();
        while (it != ite)
{
    xSolid tetin =  *it ;
    int nnode = tetin.size(0);
    switch (nnode)
        {
        case 4 : {
            xVertex v0in = tetin.getVertex(0);
            xVertex v1in = tetin.getVertex(1);
            xVertex v2in = tetin.getVertex(2);
            xVertex v3in = tetin.getVertex(3);
                xVertex v0out = v0in.getData<xEntity>(xfem::xMesh::get_was_created_by_tag());
                xVertex v1out = v1in.getData<xEntity>(xfem::xMesh::get_was_created_by_tag());
                xVertex v2out = v2in.getData<xEntity>(xfem::xMesh::get_was_created_by_tag());
                xVertex v3out = v3in.getData<xEntity>(xfem::xMesh::get_was_created_by_tag());
            AOMD::pGEntity pge = (any_cast<AOMD::mEntity*>( tetin.getEntityIdentifier()) )->getClassification() ;
            xSolid tetout(out.mesh.createTetWithVertices( any_cast<AOMD::mVertex*>(v0out.getEntityIdentifier()),
                                any_cast<AOMD::mVertex*>(v1out.getEntityIdentifier()),
                                any_cast<AOMD::mVertex*>(v2out.getEntityIdentifier()),
                                any_cast<AOMD::mVertex*>(v3out.getEntityIdentifier()),
                                pge) ,
                out);
                tetin.attachData<xEntity>(xfem::xMesh::get_was_created_by_tag(),tetout);
                tetout.attachData<xEntity>(xfem::xMesh::get_is_the_copy_of_tag(),tetin);
            break;
        }
        case 5 : {
            cout << "xCopyMesh ,3D : PAS ENCORE PROGRAMME ! " << endl;
            throw 1;
        }
        default : {
            throw 1;
        }
        }
    ++it;
}
        if (clean_tag_on_source_mesh)
{
    for (int i = 0; i < 4; ++i)
        {
            auto it = in.begin(i);
            auto ite = in.end(i);
            while (it != ite)
    {
        xVertex vin = *it ;
                vin.deleteData<xEntity>(xfem::xMesh::get_was_created_by_tag());
        ++it;
    }
        }
}

        return;
    }
 */

}  // namespace aomd

}  // namespace xinterface
