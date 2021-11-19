/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.
*/

#include <cassert>
#include "lPhysSurf.h"
#include "lSimpleGeometry.h"
#include "mxIterator.h"
#include "Field.h"
#include "lOperator.h"
#include <fstream>

//for adaptation
#include "MeshAdapt.h"
#include "AdaptUtil.h"
#include "AOMD_Internals.h"
#include "MeshTools.h"



#if 1
int CB_count=0;
extern "C" void myCallback(pPList oldCavity, pPList newcavity,void *userdata)
{
  CB_count++;
  // PList_printx(oldCavity);
  // PList_printx(newcavity);
}
#endif


double taille(mEntity* edge)
{
  mVertex* v1 = (mVertex*) edge->get(0,0);
  mVertex* v2 = (mVertex*) edge->get(0,1);
  Trellis_Util::mPoint p1 = v1->point();
  Trellis_Util::mPoint p2 = v2->point();
  mVector vec(p1, p2);
  return std::sqrt(vec*vec);
} 


void refine_surface(const lPointToDouble& surf, int nb_pass, mxMesh* mesh) {

  const bool debug = true;
  for (int i = 0; i < nb_pass; ++i) 
    {
      std::hash_map<mEntity*,int,EntityHashKey,EntityEqualKey> cut;
      lField ls(mesh);
      ls.load(surf);
      //lFitToVertices fit(1.e-2);
      //ls.accept(fit);
      int nb = 0;
      double tmax = 0., tmin = 1e6;
      for(mxIter it = mesh->begin(1); it != mesh->end(1); ++it)
	{
	  mEntity* e = *it;
	  std::vector<double> lsv = ls.getVals(e);
	  cut[e] = 0;
	  //if (lsv[0] == 0. || lsv[1] == 0. || lsv[0] * lsv[1] < 0.0)  
	  if (fabs(lsv[0]) < 0.1 || fabs(lsv[1]) < 0.1)
	    { 
	      cut[e] = 1; 
	      nb++; 
	      tmax = std::max(taille(e), tmax);
	      tmin = std::min(taille(e), tmin);
	    }
	}
      cout << "nb edges to cut " << nb << endl;     
      cout << "tmax " << tmax << " tmin " << tmin << endl; 
 
      MeshAdapt rdr(mesh,0,0,1);    // snap off  
      rdr.setCallback(myCallback,0);
      ///////////
	EIter editer=M_edgeIter(mesh);
      pEdge edge;
      int nbmax = 400, nbe = 0;
      double frac = 0.2;
      double tcut = tmin + (1.-frac) * (tmax-tmin); 
      while( edge=EIter_next(editer) )
	{ 
	  if (debug && cut[edge] == 1 && taille(edge) > tcut) cout << "edge of taille " << 
                                 taille(edge) << " will be cut\n";
	  if (cut[edge] == 1 && taille(edge) > tcut)   
	    { 
	      rdr.setAdaptLevel(edge,1); 
	      nbe++;
	      //we also refine the edges linked to the nodes of 
	      std::list<mEntity*> mirrors;
	      mesh->lookupForMirrorEntities(edge, mirrors);
	      std::list<mEntity*>::const_iterator it;
	      for (it = mirrors.begin() ; it != mirrors.end(); ++it)
		{
		  rdr.setAdaptLevel(*it,1);
		  nbe++;
		}
	      if (nbe == nbmax) break;
	    }
          else                  rdr.setAdaptLevel(edge,0);
	}
      EIter_delete(editer);
      
      // do refinement
      cout<<"-------Begin refining and snapping-----"<<endl;
      rdr.run();
      cout << " refinement and snapping done !!!" << endl;
      cout << " call Callback function "<< CB_count << " times" <<endl;
      cout << "after step refine_surface mesh number of nodes " << mesh->size(0) << endl;
    }
  cout << "after total refine_surface mesh number of nodes " << mesh->size(0) << endl;

}

void derefine_surface(const lPointToDouble& surf, int nb_pass, mxMesh* mesh) {

  const bool debug = true;
  for (int i = 0; i < nb_pass; ++i) 
    {
      std::hash_map<mEntity*,int,EntityHashKey,EntityEqualKey> deref;
      lField ls(mesh);
      ls.load(surf);
      //lFitToVertices fit(1.e-2);
      //ls.accept(fit);
      int nb = 0;
      double tmax = 0., tmin = 1e6;
      for(mxIter it = mesh->begin(1); it != mesh->end(1); ++it)
	{
	  mEntity* e = *it;
	  std::vector<double> lsv = ls.getVals(e);
	  deref[e] = 0;
	  //if (lsv[0] == 0. || lsv[1] == 0. || lsv[0] * lsv[1] < 0.0)  
	  if (fabs(lsv[0]) > 0.1 && fabs(lsv[1]) > 0.1)
	    { 
	      deref[e] = 1; 
	      nb++; 
	      tmax = std::max(taille(e), tmax);
	      tmin = std::min(taille(e), tmin);
	    }
	}
      cout << "nb edges to cut " << nb << endl;     
      cout << "tmax " << tmax << " tmin " << tmin << endl; 
 
      MeshAdapt rdr(mesh,0,0,1);    // snap off  
      rdr.setCallback(myCallback,0);
      ///////////
	EIter editer=M_edgeIter(mesh);
      pEdge edge;
      int nbmax = 400, nbe = 0;
      double frac = 0.2;
      double tcut = tmin + (1.-frac) * (tmax-tmin); 
      while( edge=EIter_next(editer) )
	{ 
	  //if (debug && cut[edge] == 1 && taille(edge) > tcut) cout << "edge of taille " << 
            //                     taille(edge) << " will be cut\n";
	  if (deref[edge] == 1)   
	    { 
	      rdr.setAdaptLevel(edge,-1); 
	      nbe++;
	      //we also refine the edges linked to the nodes of 
	      std::list<mEntity*> mirrors;
	      mesh->lookupForMirrorEntities(edge, mirrors);
	      std::list<mEntity*>::const_iterator it;
	      for (it = mirrors.begin() ; it != mirrors.end(); ++it)
		{
		  rdr.setAdaptLevel(*it,-1);
		  nbe++;
		}
	      if (nbe == nbmax) break;
	    }
          else                  rdr.setAdaptLevel(edge,0);
	}
      EIter_delete(editer);
      
      // do refinement
      cout<<"-------Begin refining and snapping-----"<<endl;
      rdr.run();
      cout << " derefinement and snapping done !!!" << endl;
      cout << " call Callback function "<< CB_count << " times" <<endl;
      cout << "after step derefine_surface mesh number of nodes " << mesh->size(0) << endl;
    }
  cout << "after total derefine_surface mesh number of nodes " << mesh->size(0) << endl;

}

void refine_interaction(const std::vector<lPointToDouble>& surfs, int nb_pass, mxMesh* mesh) {

  const bool debug = true;
  std::vector<lPointToDouble>::const_iterator its;
  double tmax, tmin;
  for (int i = 0; i < nb_pass; ++i) 
    {
      std::hash_map<mEntity*,int,EntityHashKey,EntityEqualKey>  nbcuts;
      for (its = surfs.begin(); its != surfs.end(); ++its)
	{
	  lField ls(mesh);
	  ls.load(*its);
	  //lFitToVertices fit(1.e-2);
	  //ls.accept(fit);
	  tmax = 0.; tmin = 1e6;
	  for(mxIter it = mesh->begin(1); it != mesh->end(1); ++it)
	    {
	      mEntity* e = *it;
	      std::vector<double> lsv = ls.getVals(e);
	      //nbcuts[e] = 0.;
              if (lsv[0] == 0.)                  nbcuts[e]++;
              if (lsv[1] == 0.)                  nbcuts[e]++;
	      if (lsv[0] * lsv[1] < 0.0)         nbcuts[e]++;
	      tmax = std::max(taille(e), tmax);
	      tmin = std::min(taille(e), tmin);
	    }
	}
      
      MeshAdapt rdr(mesh,0,0,1);    // snap off  
      rdr.setCallback(myCallback,0);
      ///////////
#if 1
      if (debug) {
        int nbtocut=0;
	EIter editer_dbg=M_edgeIter(mesh);
	pEdge ed;
	while( ed=EIter_next(editer_dbg) )
	  { 
	    if (nbcuts[ed] > 1) nbtocut++;
	  }
        cout << "nb of edges to nbcut " << nbtocut << endl;
	EIter_delete(editer_dbg);
      }
#endif
      int nbmax = 50, nbe = 0;
      double frac = 0.5;
      double tcut = tmin + (1.-frac) * (tmax-tmin); 
      cout << "tmax " << tmax << " tmin " << tmin << endl;
      if (debug) cout << "tcut is " << tcut << endl;
      EIter editer=M_edgeIter(mesh);
      pEdge edge;
      while( edge=EIter_next(editer) )
	{ 
	  if (debug && nbcuts[edge] > 1) cout << "edge of taille " << 
                                 taille(edge) << " should be cut\n";
	  if (nbcuts[edge] > 1 && taille(edge) > tcut)   
	    {
	      rdr.setAdaptLevel(edge,1);
	      nbe++;
	      std::list<mEntity*> mirrors;
	      mesh->lookupForMirrorEntities(edge, mirrors);
	      std::list<mEntity*>::const_iterator it;
	      for (it = mirrors.begin() ; it != mirrors.end(); ++it)
		{
		  rdr.setAdaptLevel(*it,1);
		}
	      if (nbe == nbmax) break;
	    }
	}
      EIter_delete(editer);
      
      // do refinement
      cout<<"-------Begin refining and snapping-----"<<endl;
      rdr.run();
      cout << " refinement and snapping done !!!" << endl;
      cout << " call Callback function "<< CB_count << " times" <<endl;
      cout << "after step of refine_interaction mesh number of nodes " << mesh->size(0) << endl;
    }
  cout << "after total  refine_interaction mesh number of nodes " << mesh->size(0) << endl;

}


int main(int argc, char *argv[])  
{  
  const bool debug = true;
  std::ifstream in(argv[1]);
  std::string mesh_file; 
  in >> mesh_file;
  cout << mesh_file << endl;


#if 1 //orthognal non woven
  double radius = 0.45;
  lCylinder cylx(Trellis_Util::mPoint(0.,0.,0.),    mVector(1.,0.,0.), radius);
  lCylinder cylz1(Trellis_Util::mPoint(1.,1.,0.),   mVector(0.,0.,1.), radius);
  lCylinder cylz2(Trellis_Util::mPoint(-1.,-1.,0.), mVector(0.,0.,1.), radius);
  lCylinder cylz3(Trellis_Util::mPoint(1.,-1.,0.),  mVector(0.,0.,1.), radius);
  lCylinder cylz4(Trellis_Util::mPoint(-1.,1.,0.),  mVector(0.,0.,1.), radius);
  lCylinder cyly1(Trellis_Util::mPoint(0.,0.,-1.),  mVector(0.,1.,0.), radius);
  lCylinder cyly2(Trellis_Util::mPoint(0.,0.,1.),   mVector(0.,1.,0.), radius);
  lUnion u1(cylz1, cylz2);
  lUnion u2(cylz3, cylz4);
  lUnion u3(cyly1, cyly2);
  lUnion u4(cylx, u1);
  lUnion u5(u2, u3);
  lUnion total(u4, u5);
  std::vector<lPointToDouble> surfs;
  surfs.push_back(cylx);
  surfs.push_back(cyly1);
  surfs.push_back(cyly2);
  surfs.push_back(cylz1);
  surfs.push_back(cylz2);
  surfs.push_back(cylz3);
  surfs.push_back(cylz4);
#endif

#if 0 //tubular core
  double radius = 0.495;
  lCylinder cylz1(Trellis_Util::mPoint(0.,0.,0.),   mVector(0.,0.,1.), radius);
  lCylinder cylz2(Trellis_Util::mPoint(-0.5,-1.,0.), mVector(0.,0.,1.), radius);
  lCylinder cylz3(Trellis_Util::mPoint( 0.5,-1.,0.),  mVector(0.,0.,1.), radius);
  lCylinder cylz4(Trellis_Util::mPoint(-0.5, 1.,0.), mVector(0.,0.,1.), radius);
  lCylinder cylz5(Trellis_Util::mPoint( 0.5, 1.,0.),  mVector(0.,0.,1.), radius);
  lCylinder cylz6(Trellis_Util::mPoint(-1., 0.,0.),  mVector(0.,0.,1.), radius);
  lCylinder cylz7(Trellis_Util::mPoint( 1., 0.,0.),  mVector(0.,0.,1.), radius);
  lUnion u1(cylz1, cylz2);
  lUnion u2(cylz3, cylz4);
  lUnion u3(cylz5, cylz6);
  lUnion u4(u1, u2);
  lUnion u5(u3, cylz7);
  lUnion total(u4, u5);
  std::vector<lPointToDouble> surfs;
  surfs.push_back(cylz1);
  surfs.push_back(cylz2);
  surfs.push_back(cylz3);
  surfs.push_back(cylz4);
  surfs.push_back(cylz5);
  surfs.push_back(cylz6);
  surfs.push_back(cylz7);
#endif

#if 0 //honeycomb core
  double thickness = 0.2;
  double radius = 0.5-thickness/std::sqrt(3.);
  lHexaCylinder hex1(Trellis_Util::mPoint(0.,0.,0.), mVector(0.,0.,1.), mVector(1.,0.,0.), radius);
  lHexaCylinder hex2(Trellis_Util::mPoint(0.,-1.,0.), mVector(0.,0.,1.), mVector(1.,0.,0.), radius);
  lHexaCylinder hex3(Trellis_Util::mPoint(0.,1.,0.), mVector(0.,0.,1.), mVector(1.,0.,0.), radius);
  lHexaCylinder hex4(Trellis_Util::mPoint(-1.,-0.5,0.), mVector(0.,0.,1.), mVector(1.,0.,0.), radius);
  lHexaCylinder hex5(Trellis_Util::mPoint(-1., 0.5,0.), mVector(0.,0.,1.), mVector(1.,0.,0.), radius);
  lHexaCylinder hex6(Trellis_Util::mPoint(1.,-0.5,0.), mVector(0.,0.,1.), mVector(1.,0.,0.), radius);
  lHexaCylinder hex7(Trellis_Util::mPoint(1., 0.5,0.), mVector(0.,0.,1.), mVector(1.,0.,0.), radius);
  lUnion3 u1(hex1, hex2, hex3);
  lUnion3 u2(hex4, hex5, hex6);
  lUnion3 total(u1, u2, hex7);
#endif

#if 0 //square core
  double thickness = 0.01;
  //double thickness = 0.6;
  double radius = std::sqrt(2.)/2.-thickness/2.;
  double sqr2 = std::sqrt(2.)/2.;
  lSquareCylinder sq1(Trellis_Util::mPoint(0.,-1.,0.), mVector(0.,0.,1.), mVector(sqr2,sqr2,0.), radius);
  lSquareCylinder sq2(Trellis_Util::mPoint(1., 0.,0.), mVector(0.,0.,1.), mVector(sqr2,sqr2,0.), radius);
  lSquareCylinder sq3(Trellis_Util::mPoint(0., 1.,0.), mVector(0.,0.,1.), mVector(sqr2,sqr2,0.), radius);
  lSquareCylinder sq4(Trellis_Util::mPoint(-1.,0.,0.), mVector(0.,0.,1.), mVector(sqr2,sqr2,0.), radius);

  lUnion u1(sq1, sq2);
  lUnion u2(sq3, sq4);
  lUnion total(u1, u2);

  std::vector<lPointToDouble> surfs;
  surfs.push_back(sq1);
  surfs.push_back(sq2);
  surfs.push_back(sq3);
  surfs.push_back(sq4);
#endif


#if 0 // sinus
  double L = 2.;
  double l = L/2.;
  double per = 2.0;
  double radius = L/20.;
  double ampl = l/4.; 
  //sinus dans le plan xy
  lSinusCylinder sxy1(mVector(1.,0.,0.), mVector(0.,1.,0.), Trellis_Util::mPoint(0.5, 0.5, 1.), 
		       ampl, per, radius);
  lSinusCylinder sxy2(mVector(1.,0.,0.), mVector(0.,1.,0.), Trellis_Util::mPoint(-0.5, 0.5, 0.), 
		       ampl, per, radius);
  lSinusCylinder sxy3(mVector(1.,0.,0.), mVector(0.,1.,0.), Trellis_Util::mPoint(0.5, 0.5, -1.), 
		       ampl, per, radius);

  lSinusCylinder sxy4(mVector(1.,0.,0.), mVector(0.,1.,0.), Trellis_Util::mPoint(0.5, -0.5, 1.), 
		       ampl, per, radius);
  lSinusCylinder sxy5(mVector(1.,0.,0.), mVector(0.,1.,0.), Trellis_Util::mPoint(-0.5, -0.5, 0.), 
		       ampl, per, radius);
  lSinusCylinder sxy6(mVector(1.,0.,0.), mVector(0.,1.,0.), Trellis_Util::mPoint(0.5, -0.5, -1.), 
		       ampl, per, radius);

  lUnion3 sxya(sxy1, sxy2, sxy3);
  lUnion3 sxyb(sxy4, sxy5, sxy6);
  lUnion sxy(sxya, sxyb);

  //sinus dans le plan yz
  lSinusCylinder syz1(mVector(0.,0.,1.), mVector(0.,-1.,0.), Trellis_Util::mPoint(1.0, 0.5, 0.5), 
		       ampl, per, radius);
  lSinusCylinder syz2(mVector(0.,0.,1.), mVector(0.,-1.,0.), Trellis_Util::mPoint(0., 0.5, -0.5), 
		       ampl, per, radius);
  lSinusCylinder syz3(mVector(0.,0.,1.), mVector(0.,-1.,0.), Trellis_Util::mPoint(-1., 0.5, 0.5), 
		       ampl, per, radius);
  lSinusCylinder syz4(mVector(0.,0.,1.), mVector(0.,-1.,0.), Trellis_Util::mPoint(1.0, -0.5, 0.5), 
		       ampl, per, radius);
  lSinusCylinder syz5(mVector(0.,0.,1.), mVector(0.,-1.,0.), Trellis_Util::mPoint(0., -0.5, -0.5), 
		       ampl, per, radius);
  lSinusCylinder syz6(mVector(0.,0.,1.), mVector(0.,-1.,0.), Trellis_Util::mPoint(-1., -0.5, 0.5), 
		       ampl, per, radius);

  lUnion3 syza(syz1, syz2, syz3);
  lUnion3 syzb(syz4, syz5, syz6);
  lUnion syz(syza, syzb);


  lUnion total(sxy, syz);

#endif

#if 0 //random
  //volume rtio de particules 0.2672
  double frac = 0.2672;
  int    nb_sphere = 8;
  double min_dist = 0.05;
  lRandomPeriodicSpheres total(nb_sphere, frac, 2., min_dist);
  //lPhysSurf surface(data->mesh, 
		     // total,
	             // mxClassifyOn("inclusion", data->mesh), 
		     // mxClassifyOn("matrix", data->mesh));
  //std::vector<lPointToDouble> surfs(total.getSpheres());
#endif


  mxMesh *mesh = new mxMesh(mesh_file);
  cout << "done reading mesh_file" << endl;
  cout << "init  mesh number of nodes " << mesh->size(0) << endl;

  mesh->periodicAssociation();  

  if (debug) cout << "refine interactio" << endl;
  //refine_interaction(surfs, 1, mesh);
  //refine_surface(total, 1, mesh);
  derefine_surface(total, 1, mesh);

  if (debug) cout << "refine surface" << endl;
  //refine_surface(total, 7, mesh);

  cout << "final mesh number of nodes " << mesh->size(0) << endl;
  cout << "before write mesh " << endl;
  M_writeSMS(mesh,"mesh-refined",2);
  cout << "after  write mesh " << endl;

  mesh->periodicAssociation();  


  lPhysSurf lsfinal(mesh, total);
  


  return 0;
  
}









