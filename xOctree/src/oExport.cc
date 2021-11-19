/*
  octree is a subproject of  xfem : C++ Finite Element Library
  developed under the GNU Lesser General Public License
  See the NOTICE, CONTRIBUTORS & LICENSE files for conditions.
*/

#include <cstdio>
#include <cmath>
#include <iostream>
#include "oOctree.h"
#include "oExport.h"
#include "oLevelSet.h"
#include "oTopo.h"
#include "oKeyManager.h"

using namespace std;

namespace xoctree {



  void ExportGMSHAsciiNodes (const oOctree& octree, const oKeyManager& key_manager, const string & fn)
  {
    string floc = fn + ".pos";
    FILE *f = fopen (floc.c_str(), "w");
    
    oKeyManager::const_iterator it =  key_manager.begin();
    oKeyManager::const_iterator ite = key_manager.end();
    
    const oMapping& mapping = octree.getMapping();

    fprintf(f,"View \"nodes\" {\n");
    double xyz[3];
    for ( ;it != ite; ++it)
      {
        const oKey * key = *it;
        mapping.ijk2xyz(key->getIJK(), octree.getLevelMax(), xyz);
        int nodeID = key->getId();
        fprintf(f,"SP(%g,%g,%g) {%d};\n",xyz[0], xyz[1],xyz[2],nodeID);
      }
    fprintf(f,"};\n");
    fclose (f);
  }



  void ExportGMSHAsciiOctreeLevels (const oOctree& octree, const string& fn, bool simplex) 
  {
    string floc = fn + ".pos";
    FILE *f = fopen (floc.c_str(),"w");

    fprintf(f,"View \"element level\" {\n");
    const oMapping& mapping = octree.getMapping();

    for (int l=0;l<=octree.getLevelMax();l++)
      {
        oOctree::const_iterator it  = octree.begin(l);
        oOctree::const_iterator ite = octree.end  (l);
        double dl = (double) l;
        while (it != ite)
	  {
            if (*it == 1)
	      {
                double inf[3],sup[3];
                int ijk[3];
                octree.octree2cartesian(it, l, ijk);
                mapping.getBox(ijk, l,inf,sup);
                if (octree.getDim() == 1)
		  {
                    fprintf(f,"VL(%g,%g,%g,%g,%g,%g)  {%g,%g,%g,%g,%g,%g};\n",
                            inf[0],inf[1],0.,
                            sup[0],inf[1],0.,
                            0.,dl,0.,0.,dl,0.);
		  }
                else if (octree.getDim() == 2)
		  {
                    if (!simplex)
		      {
                        fprintf(f,"SQ(%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g)  {%g,%g,%g,%g};\n",
                                inf[0],inf[1],0.,
                                sup[0],inf[1],0.,
                                sup[0],sup[1],0.,
                                inf[0],sup[1],0.,
                                dl,dl,dl,dl);
		      }
                    else
		      {
                        if ( (ijk[0]+ijk[1]) % 2 == 0) // is i+j is even
			  {
                            fprintf(f,"ST(%g,%g,%g,%g,%g,%g,%g,%g,%g)  {%g,%g,%g};\n",
                                    inf[0],inf[1],0.,
                                    sup[0],inf[1],0.,
                                    sup[0],sup[1],0.,
                                    dl,dl,dl);

                            fprintf(f,"ST(%g,%g,%g,%g,%g,%g,%g,%g,%g)  {%g,%g,%g};\n",
                                    inf[0],inf[1],0.,
                                    sup[0],sup[1],0.,
                                    inf[0],sup[1],0.,
                                    dl,dl,dl);
			  }
                        else
			  {
                            fprintf(f,"ST(%g,%g,%g,%g,%g,%g,%g,%g,%g)  {%g,%g,%g};\n",
                                    inf[0],inf[1],0.,
                                    sup[0],inf[1],0.,
                                    inf[0],sup[1],0.,
                                    dl,dl,dl);

                            fprintf(f,"ST(%g,%g,%g,%g,%g,%g,%g,%g,%g)  {%g,%g,%g};\n",
                                    sup[0],inf[1],0.,
                                    sup[0],sup[1],0.,
                                    inf[0],sup[1],0.,
                                    dl,dl,dl);
			  }
		      }
		  }
                else if (octree.getDim() == 3)
		  {
                    // 			if (!simplex)
                    // 		  {
                    fprintf(f,"SH(%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g)  {%g,%g,%g,%g,%g,%g,%g,%g};\n",
                            inf[0],inf[1],inf[2],
                            sup[0],inf[1],inf[2],
                            sup[0],sup[1],inf[2],
                            inf[0],sup[1],inf[2],
                            inf[0],inf[1],sup[2],
                            sup[0],inf[1],sup[2],
                            sup[0],sup[1],sup[2],
                            inf[0],sup[1],sup[2],
                            dl,dl,dl,dl,dl,dl,dl,dl);
		  }
	      }
            ++it;
	  }
      }
    fprintf(f,"};\n");
    fclose(f);
  }

  void ExportGMSHAscii (const oLevelSet& ls, const oOctree& octree, const string &fn, bool simplex) 
  {
    string floc = fn + ".pos";
    FILE *f = fopen (floc.c_str(),"w");

    //fprintf(f,"View \"level set\" {\n");
    fprintf(f,"View \"%s\" {\n", fn.c_str());
    const oMapping& mapping = octree.getMapping();

    std::vector<double> lsv(octree.getTopo().pow_base2[octree.getDim()]);
    for (int l=0;l<=octree.getLevelMax();l++)
      {
        oOctree::const_iterator it  = octree.begin(l);
        oOctree::const_iterator ite = octree.end  (l);
        while (it != ite)
	  {
            if (*it == 1)
	      {
                double inf[3],sup[3];
                int ijk[3];
                octree.octree2cartesian(it, l, ijk);
                mapping.getBox(ijk, l,inf,sup);
                ls.getLevelSetValues (ijk, l, lsv);
                if (octree.getDim() == 1)
		  {
                    fprintf(f,"SL(%g,%g,%g,%g,%g,%g)  {%g,%g};\n",
                            inf[0],inf[1],0.,
                            sup[0],inf[1],0.,
                            lsv[0], lsv[1]);
		  }
                else if (octree.getDim() == 2)
		  {
                    if (!simplex)
		      {
                        fprintf(f,"SQ(%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g)  {%g,%g,%g,%g};\n",
                                inf[0],inf[1],0.,
                                sup[0],inf[1],0.,
                                sup[0],sup[1],0.,
                                inf[0],sup[1],0.,
                                lsv[0], lsv[1], lsv[2], lsv[3]);
		      }
                    else
		      {
                        if ( (ijk[0]+ijk[1]) % 2 == 0) // is i+j is even
			  {
                            fprintf(f,"ST(%g,%g,%g,%g,%g,%g,%g,%g,%g)  {%g,%g,%g};\n",
                                    inf[0],inf[1],0.,
                                    sup[0],inf[1],0.,
                                    sup[0],sup[1],0.,
                                    lsv[0], lsv[1], lsv[2]);

                            fprintf(f,"ST(%g,%g,%g,%g,%g,%g,%g,%g,%g)  {%g,%g,%g};\n",
                                    inf[0],inf[1],0.,
                                    sup[0],sup[1],0.,
                                    inf[0],sup[1],0.,
                                    lsv[0], lsv[2], lsv[3]);
			  }
                        else
			  {
                            fprintf(f,"ST(%g,%g,%g,%g,%g,%g,%g,%g,%g)  {%g,%g,%g};\n",
                                    inf[0],inf[1],0.,
                                    sup[0],inf[1],0.,
                                    inf[0],sup[1],0.,
                                    lsv[0], lsv[1], lsv[3]);

                            fprintf(f,"ST(%g,%g,%g,%g,%g,%g,%g,%g,%g)  {%g,%g,%g};\n",
                                    sup[0],inf[1],0.,
                                    sup[0],sup[1],0.,
                                    inf[0],sup[1],0.,
                                    lsv[1], lsv[2], lsv[3]);
			  }
		      }
		  }
                else if (octree.getDim() == 3)
		  {
                    assert(!simplex); //simplex not coded yet
                    fprintf(f,"SH(%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g)  {%g,%g,%g,%g,%g,%g,%g,%g};\n",
                            inf[0],inf[1],inf[2],
                            sup[0],inf[1],inf[2],
                            sup[0],sup[1],inf[2],
                            inf[0],sup[1],inf[2],
                            inf[0],inf[1],sup[2],
                            sup[0],inf[1],sup[2],
                            sup[0],sup[1],sup[2],
                            inf[0],sup[1],sup[2],
                            lsv[0],lsv[1],lsv[2],lsv[3],lsv[4],lsv[5],lsv[6],lsv[7]);
		  }
	      }
            ++it;
	  }
      }
    fprintf(f,"};\n");
    fclose(f);
  }


  void ExportGMSHBinary (const oLevelSet& ls,  const oOctree& octree, const string& fn)
  {
    const bool debug = false;
    string floc = fn + ".pos";
    FILE *file = fopen (floc.c_str(),"wb");

    fprintf(file, "$PostFormat\n");
    fprintf(file, "%g %d %d\n", 1.4, 1, (int)(sizeof(double)));
    fprintf(file, "$EndPostFormat\n");
    fprintf(file, "$View\n");

    string view_name = fn;

    int one = 1;
    int nb_time_steps = 1;
    int nb_scalar_points=0; int nb_vector_points=0; int nb_tensor_points=0;
    int nb_scalar_lines=0; int nb_vector_lines=0; int nb_tensor_lines=0;
    int nb_scalar_triangles=0; int nb_vector_triangles=0; int nb_tensor_triangles=0;
    int nb_scalar_quadrangles=0; int nb_vector_quadrangles=0; int nb_tensor_quadrangles=0;
    int nb_scalar_tetrahedra=0; int nb_vector_tetrahedra=0; int nb_tensor_tetrahedra=0;
    int nb_scalar_hexahedra=0; int nb_vector_hexahedra=0; int nb_tensor_hexahedra=0;
    int nb_scalar_prisms=0; int nb_vector_prisms=0; int nb_tensor_prisms=0;
    int nb_scalar_pyramids=0; int nb_vector_pyramids=0; int nb_tensor_pyramids=0;
    int nb_scalar_lines2=0; int nb_vector_lines2=0; int nb_tensor_lines2=0;
    int nb_scalar_triangles2=0; int nb_vector_triangles2=0; int nb_tensor_triangles2=0;
    int nb_scalar_quadrangles2=0; int nb_vector_quadrangles2=0; int nb_tensor_quadrangles2=0;
    int nb_scalar_tetrahedra2=0; int nb_vector_tetrahedra2=0; int nb_tensor_tetrahedra2=0;
    int nb_scalar_hexahedra2=0; int nb_vector_hexahedra2=0; int nb_tensor_hexahedra2=0;
    int nb_scalar_prisms2=0; int nb_vector_prisms2=0; int nb_tensor_prisms2=0;
    int nb_scalar_pyramids2=0; int nb_vector_pyramids2=0; int nb_tensor_pyramids2=0;
    int nb_text2d=0; int nb_text2d_chars=0; int nb_text3d=0; int nb_text3d_chars=0;




    std::vector<double> scalar_quadrangles, scalar_hexahedra;

    int count_active = 0;
    for (int l=0;l<=octree.getLevelMax();l++)
      {
        oOctree::const_iterator it  = octree.begin(l);
        oOctree::const_iterator ite = octree.end  (l);
        for (;it != ite; ++it) { if (*it == 1) count_active++;}
      }
    int nb_vals;
    if (octree.getDim() == 2)
      {
        nb_scalar_quadrangles=count_active;
        nb_vals = count_active*16;
        scalar_quadrangles.resize(nb_vals);
      }
    else
      {
        nb_scalar_hexahedra=count_active;
        nb_vals = count_active*32;
        scalar_hexahedra.resize(nb_vals);
      }


    double time_step_values[1];
    time_step_values[0]= 0.;
    int count = 0;
    std::vector<double> lsv(octree.getTopo().pow_base2[octree.getDim()]);
    const oMapping& mapping = octree.getMapping();

    for (int l=0;l<=octree.getLevelMax();l++)
      {
        oOctree::const_iterator it  = octree.begin(l);
        oOctree::const_iterator ite = octree.end  (l);
        for (;it != ite;++it)
	  {
            if (*it == 1)
	      {
                double inf[3],sup[3];
                int ijk[3];
                octree.octree2cartesian(it, l, ijk);
                mapping.getBox(ijk, l,inf,sup);
                if (octree.getDim() == 2) {
		  scalar_quadrangles[count++] = inf[0];
		  scalar_quadrangles[count++] = sup[0];
		  scalar_quadrangles[count++] = sup[0];
		  scalar_quadrangles[count++] = inf[0];
		  scalar_quadrangles[count++] = inf[1];
		  scalar_quadrangles[count++] = inf[1];
		  scalar_quadrangles[count++] = sup[1];
		  scalar_quadrangles[count++] = sup[1];
		  scalar_quadrangles[count++] = inf[2];
		  scalar_quadrangles[count++] = inf[2];
		  scalar_quadrangles[count++] = inf[2];
		  scalar_quadrangles[count++] = inf[2];

		  ls.getLevelSetValues (ijk, l, lsv);
		  for (int i =0; i < 4; ++i) scalar_quadrangles[count++] = lsv[i];
                }
                else if (octree.getDim() == 3)
		  {
                    scalar_hexahedra[count++] = inf[0];
                    scalar_hexahedra[count++] = sup[0];
                    scalar_hexahedra[count++] = sup[0];
                    scalar_hexahedra[count++] = inf[0];
                    scalar_hexahedra[count++] = inf[0];
                    scalar_hexahedra[count++] = sup[0];
                    scalar_hexahedra[count++] = sup[0];
                    scalar_hexahedra[count++] = inf[0];
                    scalar_hexahedra[count++] = inf[1];
                    scalar_hexahedra[count++] = inf[1];
                    scalar_hexahedra[count++] = sup[1];
                    scalar_hexahedra[count++] = sup[1];
                    scalar_hexahedra[count++] = inf[1];
                    scalar_hexahedra[count++] = inf[1];
                    scalar_hexahedra[count++] = sup[1];
                    scalar_hexahedra[count++] = sup[1];
                    scalar_hexahedra[count++] = inf[2];
                    scalar_hexahedra[count++] = inf[2];
                    scalar_hexahedra[count++] = inf[2];
                    scalar_hexahedra[count++] = inf[2];
                    scalar_hexahedra[count++] = sup[2];
                    scalar_hexahedra[count++] = sup[2];
                    scalar_hexahedra[count++] = sup[2];
                    scalar_hexahedra[count++] = sup[2];
                    ls.getLevelSetValues (ijk, l, lsv);
                    for (int i =0; i < 8; ++i) scalar_hexahedra[count++] = lsv[i];
		  }
	      }
	  }
      }

    if (debug) cout << " nb quads " << nb_scalar_quadrangles << " nb hex " << nb_scalar_hexahedra << endl;
    if (debug) cout << " size quad " << scalar_quadrangles.size() << " size hex " << scalar_hexahedra.size() << endl;


    fprintf(file, "%s %d "
	    "%d %d %d %d %d %d %d %d %d "
	    "%d %d %d %d %d %d %d %d %d "
	    "%d %d %d %d %d %d %d %d %d "
	    "%d %d %d %d %d %d %d %d %d "
	    "%d %d %d %d %d %d %d %d %d "
	    "%d %d %d %d\n",
            view_name.c_str(), nb_time_steps,
            nb_scalar_points, nb_vector_points, nb_tensor_points,
            nb_scalar_lines, nb_vector_lines, nb_tensor_lines,
            nb_scalar_triangles, nb_vector_triangles, nb_tensor_triangles,
            nb_scalar_quadrangles, nb_vector_quadrangles, nb_tensor_quadrangles,
            nb_scalar_tetrahedra, nb_vector_tetrahedra, nb_tensor_tetrahedra,
            nb_scalar_hexahedra, nb_vector_hexahedra, nb_tensor_hexahedra,
            nb_scalar_prisms, nb_vector_prisms, nb_tensor_prisms,
            nb_scalar_pyramids, nb_vector_pyramids, nb_tensor_pyramids,
            nb_scalar_lines2, nb_vector_lines2, nb_tensor_lines2,
            nb_scalar_triangles2, nb_vector_triangles2, nb_tensor_triangles2,
            nb_scalar_quadrangles2, nb_vector_quadrangles2, nb_tensor_quadrangles2,
            nb_scalar_tetrahedra2, nb_vector_tetrahedra2, nb_tensor_tetrahedra2,
            nb_scalar_hexahedra2, nb_vector_hexahedra2, nb_tensor_hexahedra2,
            nb_scalar_prisms2, nb_vector_prisms2, nb_tensor_prisms2,
            nb_scalar_pyramids2, nb_vector_pyramids2, nb_tensor_pyramids2,
            nb_text2d, nb_text2d_chars, nb_text3d, nb_text3d_chars);
    size_t fwriteflag;
    fwriteflag = fwrite(&one, sizeof(int), 1, file);
    fwriteflag = fwrite(time_step_values, sizeof(double), nb_time_steps, file);
    fwriteflag = fwrite(&(*scalar_quadrangles.begin()), sizeof(double), scalar_quadrangles.size(),file);
    fwriteflag = fwrite(&(*scalar_hexahedra.begin()), sizeof(double), scalar_hexahedra.size(),file);
    if (debug) cout << "Last fwriteflag: "<<fwriteflag; // to avoid -Wunused-but-set-variable, not realy usefful but leave fwriteflag available under debugger step by step
    fprintf(file, "\n$EndView\n");
    fclose(file);
  }



  void ExportENSIGHTAscii (const oLevelSet& ls, const oOctree& octree, const string& fn) 
  {
    const bool debug = false;

    string fgeo_name = fn + ".geo";
    string fcase_name = fn + ".case";
    string fscl_name = fn + ".scl";

    //Cr�ation fichier Case ensight
    FILE *fcase = fopen (fcase_name.c_str(),"w");
    fprintf(fcase, "#\n#fichier : %s\n#\n", fn.c_str());
    fprintf(fcase, "FORMAT\n");
    fprintf(fcase, "type:\t\t\tensight gold\n");
    fprintf(fcase, "\nGEOMETRY\n");
    fprintf(fcase, "model:\t\t\t%s\n\n", fgeo_name.c_str());
    fprintf(fcase, "VARIABLE\n");
    //aplication du champ scalaire.
    fprintf(fcase, "scalar per node:%15s\t%s\n", "champsca", fscl_name.c_str());
    fclose(fcase);


    //debut du fichier contenant le scalaire.
    FILE *fscl = fopen (fscl_name.c_str(),"w");
    fprintf(fscl, "Octree_scl variable scalar\n");
    fprintf(fscl, "part\n");

    //
    //

    FILE *fgeo = fopen (fgeo_name.c_str(),"w");

    fprintf(fgeo,"Octree_mesh geometry\n");
    fprintf(fgeo,"\n");
    fprintf(fgeo,"node id assign\n");
    fprintf(fgeo,"element id assign\n");
    fprintf(fgeo,"part\n");
    fprintf(fgeo,"%10d\n", 1);
    fprintf(fgeo,"my mesh\n");
    fprintf(fgeo,"coordinates\n");

    int NbElts = 0;
    for (int l=0;l<=octree.getLevelMax();l++)
      {
        oOctree::const_iterator it  = octree.begin(l);
        oOctree::const_iterator ite = octree.end  (l);
        for (;it != ite; ++it) { if (*it == 1) NbElts++;}
      }
    int NbNodes = NbElts * 4;
    if (octree.getDim() == 3) NbNodes = NbElts * 8;
    if (debug) cout << " NbElts " << NbElts << endl;


    vector<double> x_node;
    vector<double> y_node;
    vector<double> z_node;
    vector<double> ls_node;
    if (debug) cout << "  size x_node " << x_node.size() << endl;

    const oMapping& mapping = octree.getMapping();
    fprintf(fgeo,"%10d\n", NbNodes);
    //inscription du nombre de noeuds dans le fichier scalaire
    fprintf(fscl,"%10d\n", 1);
    fprintf(fscl,"coordinates\n");
    // export of the nodes
    std::vector<double> lsv(octree.getTopo().pow_base2[octree.getDim()]);
    for (int l=0;l<=octree.getLevelMax();l++)
      {
        oOctree::const_iterator it  = octree.begin(l);
        oOctree::const_iterator ite = octree.end  (l);
        for (;it != ite;++it)
	  {
            if (*it == 1)
	      {
                double inf[3],sup[3];
                int ijk[3];
                octree.octree2cartesian(it, l, ijk);
                mapping.getBox(ijk, l,inf,sup);
                if (octree.getDim() == 2)
		  {
                    x_node.push_back(inf[0]); y_node.push_back(inf[1]); z_node.push_back(inf[2]);
                    x_node.push_back(sup[0]); y_node.push_back(inf[1]); z_node.push_back(inf[2]);
                    x_node.push_back(sup[0]); y_node.push_back(sup[1]); z_node.push_back(inf[2]);
                    x_node.push_back(inf[0]); y_node.push_back(sup[1]); z_node.push_back(inf[2]);
                    //on met chaque valeur scalaire correspondant au noeud.
                    ls.getLevelSetValues (ijk, l, lsv);
                    fprintf(fscl, "%12.5e\n%12.5e\n%12.5e\n%12.5e\n", lsv[0], lsv[1], lsv[2], lsv[3]);
		  }
                else
		  {
                    x_node.push_back(inf[0]); y_node.push_back(inf[1]); z_node.push_back(inf[2]);
                    x_node.push_back(sup[0]); y_node.push_back(inf[1]); z_node.push_back(inf[2]);
                    x_node.push_back(sup[0]); y_node.push_back(sup[1]); z_node.push_back(inf[2]);
                    x_node.push_back(inf[0]); y_node.push_back(sup[1]); z_node.push_back(inf[2]);
                    x_node.push_back(inf[0]); y_node.push_back(inf[1]); z_node.push_back(sup[2]);
                    x_node.push_back(sup[0]); y_node.push_back(inf[1]); z_node.push_back(sup[2]);
                    x_node.push_back(sup[0]); y_node.push_back(sup[1]); z_node.push_back(sup[2]);
                    x_node.push_back(inf[0]); y_node.push_back(sup[1]); z_node.push_back(sup[2]);
                    ls.getLevelSetValues (ijk, l, lsv);
                    fprintf(fscl, "%12.5e\n%12.5e\n%12.5e\n%12.5e\n%12.5e\n%12.5e\n%12.5e\n%12.5e\n",
                            lsv[0],lsv[1],lsv[2],lsv[3],lsv[4],lsv[5],lsv[6],lsv[7]);

		  }
	      }
	  }
      }
    //export of the nodes coordinates
    vector<double>::const_iterator it;
    vector<double>::const_iterator itex = x_node.end(); //same of all vector<double>
    for (it = x_node.begin(); it != itex; ++it)
      {
        fprintf(fgeo, "%12.5e\n", *it);
      }
    vector<double>::const_iterator itey = y_node.end(); //same of all vector<double>
    for (it = y_node.begin(); it != itey; ++it)
      {
        fprintf(fgeo, "%12.5e\n", *it);
      }
    vector<double>::const_iterator itez = z_node.end(); //same of all vector<double>
    for (it = z_node.begin(); it != itez; ++it)
      {
        fprintf(fgeo, "%12.5e\n", *it);
      }



    //export of the elements
    if (octree.getDim() == 2) fprintf(fgeo,"quad4\n");
    else fprintf(fgeo,"hexa8\n");
    fprintf(fgeo,"%10d\n", NbElts);

    int count = 1;
    for (int l=0;l<=octree.getLevelMax();l++)
      {
        oOctree::const_iterator it  = octree.begin(l);
        oOctree::const_iterator ite = octree.end  (l);
        for (;it != ite;++it)
	  {
            if (*it == 1)
	      {
                double inf[3],sup[3];
                int ijk[3];
                octree.octree2cartesian(it, l, ijk);
                mapping.getBox(ijk, l,inf,sup);
                if (octree.getDim() == 2)
		  {
                    //		  fprintf(fgeo,"%10d%10d%10d%10d\n",
                    //  count++, count++, count++, count++);
                    fprintf(fgeo,"%10d%10d%10d%10d\n",
                            count, count+1, count+2, count+3);
                    count +=4;
		  }
                else
		  {
                    //  fprintf(fgeo,"%10d%10d%10d%10d%10d%10d%10d%10d\n",
                    //  count++, count++, count++, count++,
                    //  count++, count++, count++, count++);
                    fprintf(fgeo,"%10d%10d%10d%10d%10d%10d%10d%10d\n",
                            count, count+1, count+2, count+3,
                            count+4, count+5, count+6, count+7);
                    count +=8;
		  }
	      }
	  }
      }
    fclose(fgeo);
    fclose(fscl);

  }



  void ExportENSIGHTBinary (const oLevelSet& ls, const oOctree& octree, const string& fn) 
  {
    const bool debug = false;

    string fgeo_name = fn + ".geo";
    string fcase_name = fn + ".case";
    string fscl_name = fn + ".scl";

    //Cr�ation fichier Case ensight
    FILE *fcase = fopen (fcase_name.c_str(),"w");
    fprintf(fcase, "#\n#fichier : %s\n#\n", fn.c_str());
    fprintf(fcase, "FORMAT\n");
    fprintf(fcase, "type:\t\t\tensight gold\n");
    fprintf(fcase, "\nGEOMETRY\n");
    fprintf(fcase, "model:\t\t\t%s\n\n", fgeo_name.c_str());
    fprintf(fcase, "VARIABLE\n");
    //aplication du champ scalaire.
    fprintf(fcase, "scalar per node:%15s\t%s\n", "champsca", fscl_name.c_str());
    fclose(fcase);


    //debut du fichier contenant le scalaire.
    FILE *fscl = fopen (fscl_name.c_str(),"wb");
    size_t fwriteflag;
    fwriteflag = fwrite("Octree_scl variable scalar", sizeof(char), 80, fscl);
    fwriteflag = fwrite("part\n", sizeof(char), 80, fscl);

    //
    //

    FILE *fgeo = fopen (fgeo_name.c_str(),"wb");
    fwriteflag = fwrite("C binary", sizeof(char), 80, fgeo);
    fwriteflag = fwrite("Octree_mesh geometry", sizeof(char), 80, fgeo);
    fwriteflag = fwrite("", sizeof(char), 80, fgeo);
    fwriteflag = fwrite("node id assign", sizeof(char), 80, fgeo);
    fwriteflag = fwrite("element id assign", sizeof(char), 80, fgeo);
    int buffInt;
    buffInt=1;
    fwriteflag = fwrite("part", sizeof(char), 80, fgeo);
    fwriteflag = fwrite(&buffInt, sizeof(int),1,fgeo);
    fwriteflag = fwrite("my mesh", sizeof(char), 80, fgeo);
    fwriteflag = fwrite("coordinates", sizeof(char), 80, fgeo);

    int NbElts = 0;
    for (int l=0;l<=octree.getLevelMax();l++)
      {
        oOctree::const_iterator it  = octree.begin(l);
        oOctree::const_iterator ite = octree.end  (l);
        for (;it != ite; ++it) { if (*it == 1) NbElts++;}
      }
    int NbNodes = NbElts * 4;
    if (octree.getDim() == 3) NbNodes = NbElts * 8;
    if (debug) cout << " NbElts " << NbElts << endl;


    vector<double> x_node;
    vector<double> y_node;
    vector<double> z_node;
    vector<double> ls_node;
    if (debug) cout << "  size x_node " << x_node.size() << endl;

    fwriteflag = fwrite(&NbNodes, sizeof(int),1,fgeo);
    //inscription du nombre de noeuds dans le fichier scalaire
    fwriteflag = fwrite(&buffInt,sizeof(int),1,fscl);
    fwriteflag = fwrite("coordinates\n", sizeof(char), 80, fscl);
    // export of the nodes
    float buffFloat;
    std::vector<double> lsv(octree.getTopo().pow_base2[octree.getDim()]);

    const oMapping& mapping = octree.getMapping();

    for (int l=0;l<=octree.getLevelMax();l++)
      {
        oOctree::const_iterator it  = octree.begin(l);
        oOctree::const_iterator ite = octree.end  (l);
        for (;it != ite;++it)
	  {
            if (*it == 1)
	      {
                double inf[3],sup[3];
                int ijk[3];
                octree.octree2cartesian(it, l, ijk);
                mapping.getBox(ijk, l,inf,sup);
                if (octree.getDim() == 2)
		  {
                    x_node.push_back(inf[0]); y_node.push_back(inf[1]); z_node.push_back(inf[2]);
                    x_node.push_back(sup[0]); y_node.push_back(inf[1]); z_node.push_back(inf[2]);
                    x_node.push_back(sup[0]); y_node.push_back(sup[1]); z_node.push_back(inf[2]);
                    x_node.push_back(inf[0]); y_node.push_back(sup[1]); z_node.push_back(inf[2]);
                    //on met chaque valeur scalaire correspondant au noeud.
                    ls.getLevelSetValues (ijk, l, lsv);
                    buffFloat=(float)lsv[0];
                    fwriteflag = fwrite(&buffFloat, sizeof(float), 1, fscl);
                    buffFloat=(float)lsv[1];
                    fwriteflag = fwrite(&buffFloat, sizeof(float), 1, fscl);
                    buffFloat=(float)lsv[2];
                    fwriteflag =  fwrite(&buffFloat, sizeof(float), 1, fscl);
                    buffFloat=(float)lsv[3];
                    fwriteflag = fwrite(&buffFloat, sizeof(float), 1, fscl);
		  }
                else
		  {
                    x_node.push_back(inf[0]); y_node.push_back(inf[1]); z_node.push_back(inf[2]);
                    x_node.push_back(sup[0]); y_node.push_back(inf[1]); z_node.push_back(inf[2]);
                    x_node.push_back(sup[0]); y_node.push_back(sup[1]); z_node.push_back(inf[2]);
                    x_node.push_back(inf[0]); y_node.push_back(sup[1]); z_node.push_back(inf[2]);
                    x_node.push_back(inf[0]); y_node.push_back(inf[1]); z_node.push_back(sup[2]);
                    x_node.push_back(sup[0]); y_node.push_back(inf[1]); z_node.push_back(sup[2]);
                    x_node.push_back(sup[0]); y_node.push_back(sup[1]); z_node.push_back(sup[2]);
                    x_node.push_back(inf[0]); y_node.push_back(sup[1]); z_node.push_back(sup[2]);


                    ls.getLevelSetValues (ijk, l, lsv);

                    buffFloat=(float)lsv[0];
                    fwriteflag = fwrite(&buffFloat, sizeof(float), 1, fscl);
                    buffFloat=(float)lsv[1];
                    fwriteflag = fwrite(&buffFloat, sizeof(float), 1, fscl);
                    buffFloat=(float)lsv[2];
                    fwriteflag = fwrite(&buffFloat, sizeof(float), 1, fscl);
                    buffFloat=(float)lsv[3];
                    fwriteflag = fwrite(&buffFloat, sizeof(float), 1, fscl);
                    buffFloat=(float)lsv[4];
                    fwriteflag = fwrite(&buffFloat, sizeof(float), 1, fscl);
                    buffFloat=(float)lsv[5];
                    fwriteflag = fwrite(&buffFloat, sizeof(float), 1, fscl);
                    buffFloat=(float)lsv[6];
                    fwriteflag = fwrite(&buffFloat, sizeof(float), 1, fscl);
                    buffFloat=(float)lsv[7];
                    fwriteflag = fwrite(&buffFloat, sizeof(float), 1, fscl);
		  }
	      }
	  }
      }
    //export of the nodes coordinates
    vector<double>::const_iterator it;
    vector<double>::const_iterator itex = x_node.end(); //same of all vector<double>

    for (it = x_node.begin(); it != itex; ++it)
      {
        buffFloat=(float)*it;
        fwriteflag = fwrite(&buffFloat, sizeof(float), 1,fgeo);
      }
    vector<double>::const_iterator itey = y_node.end(); //same of all vector<double>
    for (it = y_node.begin(); it != itey; ++it)
      {
        buffFloat=(float)*it;
        fwriteflag = fwrite(&buffFloat, sizeof(float), 1,fgeo);
      }
    vector<double>::const_iterator itez = z_node.end(); //same of all vector<double>
    for (it = z_node.begin(); it != itez; ++it)
      {
        buffFloat=(float)*it;
        fwriteflag = fwrite(&buffFloat, sizeof(float), 1,fgeo);
      }



    //export of the elements
    if (octree.getDim() == 2) fwriteflag = fwrite("quad4", sizeof(char),80,fgeo);
    else fwriteflag = fwrite("hexa8\n", sizeof(char),80,fgeo);
    fwriteflag = fwrite(&NbElts, sizeof(int), 1, fgeo);

    int count = 1;
    for (int l=0;l<=octree.getLevelMax();l++)
      {
        oOctree::const_iterator it  = octree.begin(l);
        oOctree::const_iterator ite = octree.end  (l);
        for (;it != ite;++it)
	  {
            if (*it == 1)
	      {
                double inf[3],sup[3];
                int ijk[3];
                octree.octree2cartesian(it, l, ijk);
                mapping.getBox(ijk, l,inf,sup);
                if (octree.getDim() == 2)
		  {
                    fwriteflag = fwrite(&count, sizeof(int), 1, fgeo);
                    count++;
                    fwriteflag = fwrite(&count, sizeof(int), 1, fgeo);
                    count++;
                    fwriteflag = fwrite(&count, sizeof(int), 1, fgeo);
                    count++;
                    fwriteflag = fwrite(&count, sizeof(int), 1, fgeo);

		  }
                else
		  {
                    fwriteflag = fwrite(&count, sizeof(int), 1, fgeo);
                    count++;
                    fwriteflag = fwrite(&count, sizeof(int), 1, fgeo);
                    count++;
                    fwriteflag = fwrite(&count, sizeof(int), 1, fgeo);
                    count++;
                    fwriteflag = fwrite(&count, sizeof(int), 1, fgeo);
                    count++;
                    fwriteflag = fwrite(&count, sizeof(int), 1, fgeo);
                    count++;
                    fwriteflag = fwrite(&count, sizeof(int), 1, fgeo);
                    count++;
                    fwriteflag = fwrite(&count, sizeof(int), 1, fgeo);
                    count++;
                    fwriteflag = fwrite(&count, sizeof(int), 1, fgeo);
                    count++;

		  }
	      }
	  }
      }
    fclose(fgeo);
    fclose(fscl);

    if (debug) cout << "Last fwriteflag: "<<fwriteflag; // to avoid -Wunused-but-set-variable, not realy usefful but leave fwriteflag available under debugger step by step
  }





  void ExportGMSHMesh(const oOctree& octree,  const oKeyManager& key_manager, const string& fn)
  {
    string floc = fn + ".msh";
    FILE *fcase = fopen (floc.c_str(),"w");
    fprintf(fcase, "$MeshFormat\n");
    fprintf(fcase, "2 0 8\n");
    fprintf(fcase, "$EndMeshFormat\n");
    fprintf(fcase, "$Nodes\n");
    fprintf(fcase, "%d\n", (int)(key_manager.size()));
    oKeyManager::const_iterator it  = key_manager.begin();
    oKeyManager::const_iterator ite = key_manager.end();
    double xyz[3];
    const oMapping& mapping = octree.getMapping();
    for (;it != ite;++it)
      {
        const oKey * key = *it;
        mapping.ijk2xyz(key->getIJK(), octree.getLevelMax(), xyz);
        int nodeID = key->getId();
        fprintf(fcase,"%d %g %g %g\n", nodeID,xyz[0], xyz[1],xyz[2]);
      }
    int NbElts = octree.getNbActiveCells();

    fprintf(fcase, "$EndNodes\n");
    fprintf(fcase, "$Elements\n");
    fprintf(fcase, "%d\n",NbElts);

    int count = 1;

    for (int l=0;l<=octree.getLevelMax();l++)
      {
        oOctree::const_iterator it  = octree.begin(l);
        oOctree::const_iterator ite = octree.end  (l);
        for (;it != ite;++it)
	  {
            if (*it == 1)
	      {
                int ijk[3]; int node_ids[8];
                octree.octree2cartesian(it, l, ijk);
                key_manager.getNodeIdsOnElt (ijk,  l, node_ids);
                if (octree.getDim() == 2)
		  fprintf(fcase, "%d 3 3 1 1 1 %d %d %d %d\n",
			  count++,node_ids[0],node_ids[1],node_ids[2],node_ids[3]);
                else if (octree.getDim() == 3) fprintf(fcase, "%d 5 3 1 1 1 %d %d %d %d %d %d %d %d\n",
                                                       count++,node_ids[0],node_ids[1],node_ids[2],
						       node_ids[3],node_ids[4],node_ids[5],node_ids[6],node_ids[7]);
                else{
		  fclose(fcase);
		  throw;
                }
	      }
	  }
      }

    fprintf(fcase, "$EndElements\n");

    fclose(fcase);
  }

  void ExportGMSHAsciiCurvature (const oLevelSet& ls, const oOctree& octree, const string &fn, bool simplex)
  {
    string floc = fn + ".pos";
    FILE *f = fopen (floc.c_str(),"w");

    //fprintf(f,"View \"level set\" {\n");
    fprintf(f,"View \"%s\" {\n", fn.c_str());
    const oMapping& mapping = octree.getMapping();

    std::vector<double> lsv(octree.getTopo().pow_base2[octree.getDim()]);
    std::vector<double> lsc(octree.getTopo().pow_base2[octree.getDim()]);
    for (int l=0;l<=octree.getLevelMax();l++)
      {
        oOctree::const_iterator it  = octree.begin(l);
        oOctree::const_iterator ite = octree.end  (l);
        while (it != ite)
	  {
            if (*it == 1)
	      {
                double inf[3],sup[3];
                int ijk[3];
                octree.octree2cartesian(it, l, ijk);
                mapping.getBox(ijk, l,inf,sup);
                ls.getLevelSetValues (ijk, l, lsv);
                ls.getLevelSetCurvatures(ijk, l, lsc);
                if (octree.getDim() == 1)
		  {
                    fprintf(f,"SL(%g,%g,%g,%g,%g,%g)  {%g,%g};\n",
                            inf[0],inf[1],0.,
                            sup[0],inf[1],0.,
                            lsc[0], lsc[1]);
		  }
                else if (octree.getDim() == 2)
		  {
                    if (!simplex)
		      {
                        fprintf(f,"SQ(%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g)  {%g,%g,%g,%g};\n",
                                inf[0],inf[1],0.,
                                sup[0],inf[1],0.,
                                sup[0],sup[1],0.,
                                inf[0],sup[1],0.,
                                lsc[0], lsc[1], lsc[2], lsc[3]);
		      }
                    else
		      {
                        if ( (ijk[0]+ijk[1]) % 2 == 0) // is i+j is even
			  {
                            fprintf(f,"ST(%g,%g,%g,%g,%g,%g,%g,%g,%g)  {%g,%g,%g};\n",
                                    inf[0],inf[1],0.,
                                    sup[0],inf[1],0.,
                                    sup[0],sup[1],0.,
                                    lsc[0], lsc[1], lsc[2]);

                            fprintf(f,"ST(%g,%g,%g,%g,%g,%g,%g,%g,%g)  {%g,%g,%g};\n",
                                    inf[0],inf[1],0.,
                                    sup[0],sup[1],0.,
                                    inf[0],sup[1],0.,
                                    lsc[0], lsc[2], lsc[3]);
			  }
                        else
			  {
                            fprintf(f,"ST(%g,%g,%g,%g,%g,%g,%g,%g,%g)  {%g,%g,%g};\n",
                                    inf[0],inf[1],0.,
                                    sup[0],inf[1],0.,
                                    inf[0],sup[1],0.,
                                    lsc[0], lsc[1], lsc[3]);

                            fprintf(f,"ST(%g,%g,%g,%g,%g,%g,%g,%g,%g)  {%g,%g,%g};\n",
                                    sup[0],inf[1],0.,
                                    sup[0],sup[1],0.,
                                    inf[0],sup[1],0.,
                                    lsc[1], lsc[2], lsc[3]);
			  }
		      }
		  }
                else if (octree.getDim() == 3)
		  {
                    assert(!simplex); //simplex not coded yet
                    fprintf(f,"SH(%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g)  {%g,%g,%g,%g,%g,%g,%g,%g};\n",
                            inf[0],inf[1],inf[2],
                            sup[0],inf[1],inf[2],
                            sup[0],sup[1],inf[2],
                            inf[0],sup[1],inf[2],
                            inf[0],inf[1],sup[2],
                            sup[0],inf[1],sup[2],
                            sup[0],sup[1],sup[2],
                            inf[0],sup[1],sup[2],
                            lsc[0],lsc[1],lsc[2],lsc[3],lsc[4],lsc[5],lsc[6],lsc[7]);
		  }
	      }
            ++it;
	  }
      }
    fprintf(f,"};\n");
    fclose(f);
  }


  void ExportGMSHAsciiOctreeOnFinestLevel (const oOctree& octree, const std::vector<double> &valVec, const string& fn, bool simplex ){
    string floc = fn + ".pos";
    FILE *f = fopen (floc.c_str(),"w");

    string nameView = "View \"" + fn + "\" {\n";
    fprintf(f,"%s",nameView.c_str());
    //    fprintf(f,"View \"element level\" {\n");
    const oMapping& mapping = octree.getMapping();
    int l = octree.getLevelMax();
    const int powp = octree.getTopo().pow_base2[l];

    oOctree::const_iterator it  = octree.begin(l);
    oOctree::const_iterator ite = octree.end  (l);
    while (it != ite)
      {
        double inf[3],sup[3];
        int ijk[3];
        octree.octree2cartesian(it, l, ijk);
        mapping.getBox(ijk, l,inf,sup);

        int idx = (ijk[0] + powp * ijk[1]+ powp*powp*ijk[2]);
        double val = valVec[idx];

        if (octree.getDim() == 1)
	  {
            fprintf(f,"VL(%g,%g,%g,%g,%g,%g)  {%g,%g,%g,%g,%g,%g};\n",
                    inf[0],inf[1],0.,
                    sup[0],inf[1],0.,
                    0.,val,0.,0.,val,0.);
	  }
        else if (octree.getDim() == 2)
	  {
            if (!simplex)
	      {
                fprintf(f,"SQ(%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g)  {%g,%g,%g,%g};\n",
                        inf[0],inf[1],0.,
                        sup[0],inf[1],0.,
                        sup[0],sup[1],0.,
                        inf[0],sup[1],0.,
                        val,val,val,val);
	      }
            else
	      {
                if ( (ijk[0]+ijk[1]) % 2 == 0) // is i+j is even
		  {
                    fprintf(f,"ST(%g,%g,%g,%g,%g,%g,%g,%g,%g)  {%g,%g,%g};\n",
                            inf[0],inf[1],0.,
                            sup[0],inf[1],0.,
                            sup[0],sup[1],0.,
                            val,val,val);

                    fprintf(f,"ST(%g,%g,%g,%g,%g,%g,%g,%g,%g)  {%g,%g,%g};\n",
                            inf[0],inf[1],0.,
                            sup[0],sup[1],0.,
                            inf[0],sup[1],0.,
                            val,val,val);
		  }
                else
		  {
                    fprintf(f,"ST(%g,%g,%g,%g,%g,%g,%g,%g,%g)  {%g,%g,%g};\n",
                            inf[0],inf[1],0.,
                            sup[0],inf[1],0.,
                            inf[0],sup[1],0.,
                            val,val,val);

                    fprintf(f,"ST(%g,%g,%g,%g,%g,%g,%g,%g,%g)  {%g,%g,%g};\n",
                            sup[0],inf[1],0.,
                            sup[0],sup[1],0.,
                            inf[0],sup[1],0.,
                            val,val,val);
		  }
	      }
	  }
        else if (octree.getDim() == 3)
	  {
            // 			if (!simplex)
            // 		  {
            fprintf(f,"SH(%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g)  {%g,%g,%g,%g,%g,%g,%g,%g};\n",
                    inf[0],inf[1],inf[2],
                    sup[0],inf[1],inf[2],
                    sup[0],sup[1],inf[2],
                    inf[0],sup[1],inf[2],
                    inf[0],inf[1],sup[2],
                    sup[0],inf[1],sup[2],
                    sup[0],sup[1],sup[2],
                    inf[0],sup[1],sup[2],
                    val,val,val,val,val,val,val,val);
	  }

        ++it;
      }

    fprintf(f,"};\n");
    fclose(f);
  }



  void ExportGMSHAscii (const oCellFieldOnFine& field,  const oOctree& octree, const string& fn)
  {
    const int levelMax = octree.getLevelMax();
    //Only work for "isotropic" octree
    //  ^
    // /!\  3D assert not OK !!!
    // ---
    assert(field.size() == octree.getTopo().pow_base2[levelMax] * octree.getTopo().pow_base2[levelMax]);

    string floc = fn + ".pos";
    FILE *f = fopen (floc.c_str(),"w");

    fprintf(f,"View \"element level\" {\n");
    const oMapping& mapping = octree.getMapping();

    oOctree::const_iterator it  = octree.begin(levelMax);
    oOctree::const_iterator ite = octree.end  (levelMax);
    //        double dl = (double) l;
    while (it != ite)
      {
        double inf[3],sup[3];
        int ijk[3];
        octree.octree2cartesian(it, levelMax, ijk);
        mapping.getBox(ijk, levelMax,inf,sup);
        double val = field.getVal(ijk);
        if (octree.getDim() == 1)
	  {
            throw;
	  }
        else if (octree.getDim() == 2)
	  {
            fprintf(f,"SQ(%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g)  {%g,%g,%g,%g};\n",
                    inf[0],inf[1],0.,
                    sup[0],inf[1],0.,
                    sup[0],sup[1],0.,
                    inf[0],sup[1],0.,
                    val,val,val,val);
	  }
        else if (octree.getDim() == 3)
	  {
            fprintf(f,"SH(%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g)  {%g,%g,%g,%g,%g,%g,%g,%g};\n",
                    inf[0],inf[1],inf[2],
                    sup[0],inf[1],inf[2],
                    sup[0],sup[1],inf[2],
                    inf[0],sup[1],inf[2],
                    inf[0],inf[1],sup[2],
                    sup[0],inf[1],sup[2],
                    sup[0],sup[1],sup[2],
                    inf[0],sup[1],sup[2],
                    val,val,val,val,val,val,val,val);
	  }

        ++it;
      }

    fprintf(f,"};\n");
    fclose(f);
  }

}//end namespace



