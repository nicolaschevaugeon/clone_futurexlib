/*
  xfem : C++ Finite Element Library
  developed under the GNU Lesser General Public License
  See the NOTICE & LICENSE files for conditions.
*/
#include <iostream>

#include "xSplitAOMD.h"

// AOMD
#include "mEntity.h"
#include "mEdge.h"
#include "mFace.h"
#include "mVertex.h"

#include "xAspectRatioAOMD.h"

using namespace AOMD;
using namespace std;

namespace xinterface
{

  namespace xtemplaterefinemesh
  {

    //===== xMeshSplitWarperAOMD ===============================================================================================
    //
    xMeshSplitWarperAOMD::xMeshSplitWarperAOMD(AOMD::mMesh &m_) : mesh(m_),nb_pc(0){pc[0] = pc[1] = pc[2] = 0.; }

    AOMD::mEntity * xMeshSplitWarperAOMD::createCOGVertex (AOMD::mEntity &e)
    {
      Trellis_Util::mPoint pp (mesh.COG(&e));
      return mesh.createVertex(pp,mesh.getGEntity(0,0));
    }

    AOMD::mEntity * xMeshSplitWarperAOMD::createEdge(AOMD::mEntity *v1, AOMD::mEntity *v2)
    {
      mEntity * edge = (mEntity *) mesh.createEdge( (mVertex *) v1,(mVertex *) v2,mesh.getGEntity(0,1));
      v1->add(edge);
      v2->add(edge);
      return edge;
    }

    AOMD::mEntity * xMeshSplitWarperAOMD::createFaceWithEdge(AOMD::mEntity *e1, AOMD::mEntity *e2,  AOMD::mEntity *e3)
    {
      mEntity * face = (mEntity *) mesh.createFaceWithEdges ((mEdge *) e1,(mEdge *) e2,(mEdge *) e3,mesh.getGEntity(0,2));
      e1->add(face);
      e2->add(face);
      e3->add(face);
      for (int i = 0; i < 3; ++i) face->get(0,i)->add(face);
      return face;
    }
    AOMD::mEntity * xMeshSplitWarperAOMD::createTetWithFace(AOMD::mEntity *f1, AOMD::mEntity *f2, AOMD::mEntity *f3, AOMD::mEntity *f4)
    {
      mEntity * t = (mEntity *) mesh.createTetWithFaces ((mFace *) f1,(mFace *) f2,(mFace *) f3,(mFace *) f4,mesh.getGEntity(0,3));
      mFace *faces[4] = {nullptr,nullptr,nullptr,nullptr};
      for (int i = 0; i < 4; ++i)
	{
	  // adjancy
	  // 3->2 done by construction
	  // 2->3
	  faces[i] = (mFace *) t->get(2,i);
	  faces[i]->add(t);
	}
      // adjancy
      // 3->1 + 1->3
      for (int i = 1; i < 4; ++i)
	{
	  mEdge * e = faces[0]->commonEdge(faces[i]);
	  assert(e);
	  t->add(e);
	  e->add(t);
	}
      {
        mEdge * e = faces[3]->commonEdge(faces[1]);
        assert(e);
        t->add(e);
        e->add(t);
      }
      for (int i = 1; i < 3; ++i)
	{
	  mEdge * e = faces[i]->commonEdge(faces[i+1]);
	  assert(e);
	  t->add(e);
	  e->add(t);
	}
      // adjancy
      // 3->0 + 0->3
      for (int i = 0; i < 4; ++i)
	{
	  mVertex * v = faces[i/3]->commonVertex(faces[1+( 2*i )/3],faces[( i%2 ) ? 2 : 3]);
	  assert(v);
	  t->add(v);
	  v->add(t);
	}
      return t;
    }

    AOMD::mEntity * xMeshSplitWarperAOMD::getEdgeFromVertex(AOMD::mEntity *v1, AOMD::mEntity *v2)
    {
      mEntity * edge = mesh.getEdge((mVertex *) v1,(mVertex *) v2);
      // Note : here testing swap is rather specific to Split ...
      if (!edge) edge = mesh.getEdge((mVertex *) v2,(mVertex *) v1);
      return edge;
    }

    AOMD::mEntity * xMeshSplitWarperAOMD::getFaceFromVertex(AOMD::mEntity *v1, AOMD::mEntity *v2,  AOMD::mEntity *v3)
    {
      AOMD::mEntity *f = (AOMD::mEntity *) ( mesh.getTri((mVertex *) v1,(mVertex *) v2,(mVertex *) v3));
      // Note : here testing swap is rather specific to Split ...
      if (f) return f;
      f = (AOMD::mEntity *) ( mesh.getTri((mVertex *) v3,(mVertex *) v2,(mVertex *) v1));
      if (f) return f;
      f = (AOMD::mEntity *) ( mesh.getTri((mVertex *) v2,(mVertex *) v1,(mVertex *) v3));
      if (f) return f;
      f = (AOMD::mEntity *) ( mesh.getTri((mVertex *) v1,(mVertex *) v3,(mVertex *) v2));
      if (f) return f;
      f = (AOMD::mEntity *) ( mesh.getTri((mVertex *) v2,(mVertex *) v3,(mVertex *) v1));
      if (f) return f;
      f = (AOMD::mEntity *) ( mesh.getTri((mVertex *) v3,(mVertex *) v1,(mVertex *) v2));
      return f;
    }
    void xMeshSplitWarperAOMD::removeEntityWithAdj(AOMD::mEntity &e)
    {
      mesh.DEL_updateAdj(&e);
    }
    auto xMeshSplitWarperAOMD::begin(int what)->iter_type
    {
      return mesh.begin(what);
    }
    auto xMeshSplitWarperAOMD::end(int what)->iter_type
    {
      return mesh.end(what);
    }
    void xMeshSplitWarperAOMD::printEntity(const AOMD::mEntity *e, bool deep)
    {
      if (deep) cout<<endl;
      int n;
      switch (e->getLevel())
	{
        case 0 :
	  {
	    cout<<"node "<<e<<" "<<e->getId()<<endl;

	    for (int i = 1; deep && i < 4; ++i)
	      {
		if (e->isAdjacencyCreated(i))
		  {
		    cout <<" adjancy "<<i<<endl;
		    for (int k = 0,m = e->size(i); k < m; ++k) printEntity(e->get(i,k),false);
		  }
	      }
	    break;
	  }
        case 1 :
	  {
	    n = 2;
	    cout<<"edge "<<e<<" "<<e->getId();
	    if (e->isAdjacencyCreated(0))
	      {
		if (e->size(0) == n)
		  {
		    cout<<" ("<<e->get(0,0)->getId();
		    for (int k = 1; k < n; ++k) cout<<"-"<<e->get(0,k)->getId();
		    cout<<")"<<endl;
		  }
		else
		  cout<<" (wrong number of nodes !! sould be "<<n<<" put is "<<e->size(0)<<")"<<endl;
	      }

	    for (int i = 1; deep && i < 4; ++i)
	      {
		if (e->isAdjacencyCreated(i))
		  {
		    cout <<" adjancy "<<i<<endl;
		    for (int k = 0,m = e->size(i); k < m; ++k) printEntity(e->get(i,k),false);
		  }
	      }
	    break;
	  }
        case 2 :
	  {
	    n = 3;
	    cout<<"tria "<<e<<" "<<e->getId();
	    if (e->isAdjacencyCreated(0))
	      {
		if (e->size(0) == n)
		  {
		    cout<<" ("<<e->get(0,0)->getId();
		    for (int k = 1; k < n; ++k) cout<<"-"<<e->get(0,k)->getId();
		    cout<<")"<<endl;
		  }
		else
		  cout<<" (wrong number of nodes !! sould be "<<n<<" put is "<<e->size(0)<<")"<<endl;
	      }

	    for (int i = 1; deep && i < 4; ++i)
	      {
		if (e->isAdjacencyCreated(i))
		  {
		    cout <<" adjancy "<<i<<endl;
		    for (int k = 0,m = e->size(i); k < m; ++k) printEntity(e->get(i,k),false);
		  }
	      }
	    break;
	  }
        case 3 :
	  {
	    n = 4;
	    cout<<"tet "<<e<<" "<<e->getId();
	    if (e->isAdjacencyCreated(0))
	      {
		if (e->size(0) == n)
		  {
		    cout<<" ("<<e->get(0,0)->getId();
		    for (int k = 1; k < n; ++k) cout<<"-"<<e->get(0,k)->getId();
		    cout<<")"<<endl;
		  }
		else
		  cout<<" (wrong number of nodes !! sould be "<<n<<" put is "<<e->size(0)<<")"<<endl;
	      }

	    for (int i = 1; deep && i < 4; ++i)
	      {
		if (e->isAdjacencyCreated(i))
		  {
		    cout <<" adjancy "<<i<<endl;
		    for (int k = 0,m = e->size(i); k < m; ++k) printEntity(e->get(i,k),false);
		  }
	      }
	    break;
	  }
	}
    }
    void xMeshSplitWarperAOMD::printMesh(bool deep)
    {
      cout<<endl<<"Export mesh"<<endl;
      for (int i = 0; i < 4; ++i)
	{
	  iter_type b = mesh.begin(i);
	  iter_type e = mesh.end(i);
	  for (; b != e; ++b) printEntity(*b,deep);
	}
    }
    void xMeshSplitWarperAOMD::writeToFile(const std::string &file_name)
    {
      AOMD::AOMD_Util::Instance()->exportGmshFile(file_name.c_str(),&mesh);
    }
    double xMeshSplitWarperAOMD::tetAspectRatio (AOMD::mEntity *tet)
    {
      return xtool::xAspectRatio<AOMD::mTet,xtool::enum_xAspectRatio::MAXLENGTHMINHEIGHTRATIO >::aspectRatio(tet);
    }
    double xMeshSplitWarperAOMD::tetAspectRatio (AOMD::mEntity *v1, AOMD::mEntity *v2,AOMD::mEntity *v3, AOMD::mEntity *v4)
    {
      return xtool::xAspectRatio<AOMD::mTet,xtool::enum_xAspectRatio::MAXLENGTHMINHEIGHTRATIO >::aspectRatio(v1,v2,v3,v4);
    }

    void xMeshSplitWarperAOMD::accumulatePointCoord(const AOMD::mEntity * e)
    {
      assert(dynamic_cast < const mVertex * >( e ));
      Trellis_Util::mPoint pp (((mVertex *)(e))->point());
      pc[0]+=pp(0);
      pc[1]+=pp(1);
      pc[2]+=pp(2);
      ++nb_pc;
    }
    void xMeshSplitWarperAOMD::resetPointCoord()
    {
      nb_pc = 0;
      pc[0] = pc[1] = pc[2] = 0.;
    }
    AOMD::mEntity * xMeshSplitWarperAOMD::createMeanPointCoord ()
    {
      assert(nb_pc);
      return ((mEntity *) mesh.createVertex(pc[0]/nb_pc,pc[1]/nb_pc,pc[2]/nb_pc,mesh.getGEntity(0,0)));
    }
    //===== xSplitTransferAOMDGeom =============================================================================================
    void xSplitTransferAOMDGeom::collect(AOMD::mEntity *e)
    {
      g = e->getClassification();
    }
    /// A method which transfer classification previously collected
    void xSplitTransferAOMDGeom::transfer(AOMD::mEntity *se)
    {
      se->classify(g);
    }

  }         // end subnamespace

}         // end namespace
