/* 
   This file is a part of eXlibris C++ Library
   under the GNU Lesser General Public License.
   See the NOTICE.md & LICENSE.md files for terms 
   and conditions.
*/


#include "xRefCut.h"
#include "xRefMesh.h"


namespace xcut
{


  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////// function implementation ////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //
  ///////////////
  /// 1D Edge
  ///////////////
  int cutEdgeRefByLevelSet ( const std::vector < double > & ls_value, const std::vector < int > & node_label, xRefMesh & cutmesh)
  {
    //
    //     0------1    --> u
    //         0
    //
    // mapping u : [-1,1] <=> [x0,x1]
    //
    double lsv0 = ls_value[0];
    double lsv1 = ls_value[1];
    const double one = 1.0E+0;
    const double zero = 0.0E+0;

    // the iso-zero cut the edge ////////////////////////////////////////////////////////////////////////////
    if ( ( lsv0*lsv1 ) < zero )
      {
        xPoint UVW;

        //1D => v,w=0
        UVW.v = zero;
        UVW.w = zero;

        xRefMesh::elemdef_t conectivity(2);
        xRefMesh::elemdef_t iso_zero_conectivity(1);

        //store the constitutive nodes
        // node 0
        UVW.u = zero;
        cutmesh.addPoint(0,UVW);
        // node 1
        UVW.u = one;
        cutmesh.addPoint(1,UVW);
        // add the cutting node
        UVW.u = -2.0*lsv0/( lsv1 - lsv0 )-1.0; // lineare interpolation in [-1,1]
        cutmesh.addPoint(2,UVW);


        // sign of ls give the "in" "out"
        if (lsv0 < zero)
	  {
            // add element "in"
            conectivity[0] = 2;
            conectivity[1] = 0;
            cutmesh.addElem(0,conectivity);
            // add element "out"
            conectivity[1] = 1;
            cutmesh.addElem(1,conectivity);
	  }
        else
	  {
            // add element "in"
            conectivity[0] = 2;
            conectivity[1] = 1;
            cutmesh.addElem(0,conectivity);
            // add element "out"
            conectivity[1] = 0;
            cutmesh.addElem(1,conectivity);
	  }

        // add iso zero "element"
        iso_zero_conectivity[0] = 2;
        cutmesh.addElem(2,iso_zero_conectivity);

        // separation betwen "in" an "out" element
        cutmesh.setLimite(1);

        // iso zero element is at the end
        cutmesh.setIsozero(2);


      }
    // the iso-zero pass thru 2 nodes  /////////////////////////////////////////////////////////////////////
    // or this function shouldn't have been called as level set don't cut the element /////////////////////
    else
      {
        // no cut
        return ( 0 );
      }


    // normal endding
    return ( 1 );

  }

  //
  //
  //
  //
  //
  ///////////////
  /// 2D triangle
  ///////////////
  int cutTriRefByLevelSet ( const std::vector < double > & ls_value, const std::vector < int > &node_label, xRefMesh & cutmesh)
  {
    /********************************************************************************************************************
     * DEFINITION
     * ==========
     *
     *
     *        v
     *        ^
     *        |
     *
     *        2
     *        | \
     *        |  \
     *       2|   \ 1
     *        |    \
     *        |     \
     *        0------1    --> u
     *            0
     *
     *
     *    mapping u : [0,1] <=> [x0,x1]
     *            v : [0,1] <=> [y0,y2]
     *
     *
     *    edge conectivity :  edge 0 is defined by vertice 0 and 1
     *                        edge 1 is defined by vertice 1 and 2
     *                        edge 2 is defined by vertice 0 and 2
     *
     * CUT SUMARY
     * ==========
     *
     *    There is 6 known possibles cut of a triangle by a iso-zero surface by linear interpolation
     *
     *    All of them can be reduce to 2 references cases by rotations/permutations
     *
     *    The 2 cases differs by there number of edge cut :
     *
     *      case 1 (C1) : 1 edge only is cut, the iso-zero surface passe thru the opposite coin of the triangle
     *      case 2 (C2) : 2 edges are cut, the iso-zero surface passe thru them
     *
     * CONVENTIONS
     * ===========
     *
     *
     *    construction of subentities of subentities is folowing these conventions :
     *                       - edges conectivity start from cut node of the original edge to nodes of the triangle
     *                       - numbering of edges use original edge number + 10 * logical order from "in" to "out" (starting at 1)
     *
     *    from iso-zero edges number "i", original edge number is given by i%10
     */
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////
    // permutation/rotation of a triangle  : first tree termes for C1, all termes for C2
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////
    const unsigned char perm_rot[4][5] = {
      // node  0  1  2  3  4
      {2, 0, 1, 4, 3},                                // C1 : edge 0 cut pass thru 2; C2 : edge 1,2 cut
      {0, 1, 2, 3, 4},                                // reference case no permutation C1 : edge 1 cut pass thru 0; C2 : edge 0,2 cut 
      {1, 2, 0, 4, 3},                                // C1 : edge 2 cut pass thru 1; C2 : edge 0,1 cut
      {0, 0, 0, 0, 0}                                 // Error case
    };
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////
    // indirection to permutation/rotation of a triangle for case 2
    //  table have 2 index of 2 terms long (no comon index possible = 2 same edge stored in id_cut is imposible)
    //  index2 is always greater then index1
    //  out of this 2X2=4 terms only 3 are relevant
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////
    const unsigned char indirect[2][2] = {
      // edge    1   2     cut
      {  2,  1},                                  // edge 0 cut
      {  3,  0}                                   // edge 1 cut
    };
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////
    // sign of cut node give "in" and "out" for C1
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////
    const bool sign_c1[4][2] = {
      {false,true},                            // 0
      {false,true},                            // 1
      {true,false},                            // 2
      {false,true}                             // 3
    };
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////
    // sign of first cut node give "in" and "out" for C2
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////
    const bool sign_c2[4][2] = {
      {true,false},                            // 0
      {false,true},                            // 1
      {true,false},                            // 2
      {false,true}                             // 3
    };
#ifdef XREFMESH_WITH_SUB
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////
    // vertex to edge table
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////
    const unsigned char v_to_e[3][3] = {
      // vertex     0  1  2
      {6, 0, 2},                                // vertex 0
      {0, 6, 1},                                // vertex 1
      {2, 1, 6}                                 // vertex 2
    };
#endif
    unsigned char nb_cut = 0;
    unsigned char id_cut[2];
    unsigned char sign_edge[2];
    xPoint cut_point[2];
    xPoint UVW;
    double lsv0 = ls_value[0];
    double lsv1 = ls_value[1];
    double lsv2 = ls_value[2];
    const double one = 1.0E+0;
    const double zero = 0.0E+0;

    xRefMesh::elemdef_t conectivity(3);
    xRefMesh::elemdef_t edge_conectivity(2);

    // the iso-zero cut the edge 0 ////////////////////////////////////////////////////////////////////////////
    if ( ( lsv0*lsv1 ) < zero )
      {

        // lineare interpolation
        cut_point[nb_cut].u = -lsv0/( lsv1 - lsv0 );
        cut_point[nb_cut].v = zero;
        //2D => w=0
        cut_point[nb_cut].w = zero;

        sign_edge[nb_cut] = ( lsv0 < zero ) ? 1 : 0;
        id_cut[nb_cut] = 0;
        ++nb_cut;
      }

    // the iso-zero cut the edge 1 ////////////////////////////////////////////////////////////////////////////
    if ( ( lsv1*lsv2 ) < 0 )
      {
        // lineare interpolation
        cut_point[nb_cut].v = -lsv1/( lsv2 - lsv1 );
        cut_point[nb_cut].u = one -  cut_point[nb_cut].v;
        //2D => w=0
        cut_point[nb_cut].w = zero;

        sign_edge[nb_cut] = ( lsv1 < zero ) ? 1 : 0;
        id_cut[nb_cut] = 1;
        ++nb_cut;
      }

    // the iso-zero cut the edge 2 ////////////////////////////////////////////////////////////////////////////
    if ( ( lsv0*lsv2 ) < 0 )
      {
        // lineare interpolation
        cut_point[nb_cut].u = zero;
        cut_point[nb_cut].v = -lsv0/( lsv2 - lsv0 );
        //2D => w=0
        cut_point[nb_cut].w = zero;

        sign_edge[nb_cut] = ( lsv0 < zero ) ? 1 : 0;
        id_cut[nb_cut] = 2;
        ++nb_cut;
      }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////
    // creation of the  cut mesh //////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////
    switch (nb_cut)
      {
        //  1 edge only is cut, the iso-zero surface passe thru  the opposites coins of the triangle
        //  submesh => 2 triangles; iso-zero => 1 edge
        //  subedge => 1 edge cut => 2 edges
        // based on the the reference case : edge 1 cut, pass thru node 0
      case 1 :
	{
	  // local
	  const unsigned char i_edge = id_cut[0];
	  const unsigned char * perm_node = &( perm_rot[i_edge][0] );
	  const unsigned char i0 = perm_node[0];
	  const unsigned char i1 = perm_node[1];
	  const unsigned char i2 = perm_node[2];

	  // adding the new nodes
	  cutmesh.addPoint(3,cut_point[0]);

	  // add iso zero element
	  edge_conectivity[0] = i0;
	  edge_conectivity[1] = 3;
	  cutmesh.addElem(2,edge_conectivity);
	  // iso zero element is at the end
	  cutmesh.setIsozero(2);
	  // separation betwen "in" an "out" element
	  cutmesh.setLimite(1);
#ifdef XREFMESH_WITH_SUB
	  // separation betwen "in" an "out" edge
	  cutmesh.setEdgeLimite(20+i_edge);
#endif
	  // conectivity of generated element all start by node 0
	  conectivity[0] = i0;

	  // if node i1 is "in"
	  if (sign_c1[i_edge][sign_edge[0]])
	    {
	      // add element "in"
	      conectivity[1] = i1;
	      conectivity[2] = 3;
	      cutmesh.addElem(0,conectivity);
	      // add element "out"
	      conectivity[1] = 3;
	      conectivity[2] = i2;
	      cutmesh.addElem(1,conectivity);
#ifdef XREFMESH_WITH_SUB
	      // add edge "in"
	      edge_conectivity[0] = 3;
	      edge_conectivity[1] = i1;
	      cutmesh.addEdge(10+i_edge,edge_conectivity);
	      // add edge "out"
	      edge_conectivity[1] = i2;
	      cutmesh.addEdge(20+i_edge,edge_conectivity);
#endif
	    }
	  else
	    {
	      // add element "in"
	      conectivity[1] = 3;
	      conectivity[2] = i2;
	      cutmesh.addElem(0,conectivity);
	      // add element "out"
	      conectivity[1] = i1;
	      conectivity[2] = 3;
	      cutmesh.addElem(1,conectivity);
#ifdef XREFMESH_WITH_SUB
	      // add edge "in"
	      edge_conectivity[0] = 3;
	      edge_conectivity[1] = i2;
	      cutmesh.addEdge(10+i_edge,edge_conectivity);
	      // add edge "out"
	      edge_conectivity[1] = i1;
	      cutmesh.addEdge(20+i_edge,edge_conectivity);
#endif
	    }
	  break;
	}
        //  2 edges are cut
        //  submesh => 3 triangles; iso-zero => 1 edge
        //  subedge => 2 edge cut => 4 edges
        // based on the the reference case : edge 0,2 cut
      case 2 :
	{
	  // local
	  const unsigned char i = indirect[id_cut[0]][id_cut[1]-1];
	  const unsigned char * perm_node = &( perm_rot[i][0] );
	  const unsigned char i0 = perm_node[0];
	  const unsigned char i1 = perm_node[1];
	  const unsigned char i2 = perm_node[2];
	  const unsigned char i3 = perm_node[3];
	  const unsigned char i4 = perm_node[4];
	  unsigned char n1,n2,n3;
	  //	     unsigned char n4;
#ifdef XREFMESH_WITH_SUB
	  // set edge number
	  const unsigned char i_edge0 = v_to_e[i0][i1];
	  const unsigned char i_edge2 = v_to_e[i0][i2];
#endif

	  // adding the new nodes
	  cutmesh.addPoint(3,cut_point[0]);
	  cutmesh.addPoint(4,cut_point[1]);

	  // add iso zero element
	  edge_conectivity[0] = 3;
	  edge_conectivity[1] = 4;
	  cutmesh.addElem(3,edge_conectivity);
	  // iso zero element is at the end
	  cutmesh.setIsozero(3);

	  // for conform mesh use node label to follow a uniforme cut strategie among triangle and theire supported tetraedron if any
	  // use of the lowest label as starting point of the cut
	  if (node_label[i1] < node_label[i2])
	    {
	      // i1 have the lowest label
	      n1 = i4;
	      n2 = i3;
	      n3 = i1;
	    }
	  else
	    {
	      // i2 have the lowest label
	      n1 = i3;
	      n2 = i2;
	      n3 = i4;
	    }

	  // sign of first cut edge give  "in"/"out" domaine
	  // reference : i0 is "in"
	  if (sign_c2[i][sign_edge[0]])
	    {
	      // add element "in"
	      conectivity[0] = i0;
	      conectivity[1] = i3;
	      conectivity[2] = i4;
	      cutmesh.addElem(0,conectivity);
	      // separation betwen "in" an "out" element
	      cutmesh.setLimite(1);
	      // add element "out"
	      conectivity[0] = n1;
	      conectivity[1] = i1;
	      conectivity[2] = i2;
	      cutmesh.addElem(1,conectivity);
	      conectivity[1] = n2;
	      conectivity[2] = n3;
	      cutmesh.addElem(2,conectivity);
#ifdef XREFMESH_WITH_SUB
	      // add edge "in"
	      edge_conectivity[0] = i3;
	      edge_conectivity[1] = i0;
	      cutmesh.addEdge(10+i_edge0,edge_conectivity);
	      edge_conectivity[0] = i4;
	      cutmesh.addEdge(20+i_edge2,edge_conectivity);
	      // add edge "out"
	      edge_conectivity[1] = i2;
	      cutmesh.addEdge(30+i_edge2,edge_conectivity);
	      edge_conectivity[0] = i3;
	      edge_conectivity[1] = i1;
	      cutmesh.addEdge(40+i_edge0,edge_conectivity);
	      // separation betwen "in" an "out" edge
	      cutmesh.setEdgeLimite(30+i_edge2);
#endif
	    }
	  else
	    {
	      // add element "in"
	      conectivity[0] = n1;
	      conectivity[1] = i1;
	      conectivity[2] = i2;
	      cutmesh.addElem(0,conectivity);
	      conectivity[1] = n2;
	      conectivity[2] = n3;
	      cutmesh.addElem(1,conectivity);
	      // separation betwen "in" an "out" element
	      cutmesh.setLimite(2);
	      // add element "out"
	      conectivity[0] = i0;
	      conectivity[1] = i3;
	      conectivity[2] = i4;
	      cutmesh.addElem(2,conectivity);
#ifdef XREFMESH_WITH_SUB
	      // add edge "out"
	      edge_conectivity[0] = i3;
	      edge_conectivity[1] = i0;
	      cutmesh.addEdge(30+i_edge0,edge_conectivity);
	      edge_conectivity[0] = i4;
	      cutmesh.addEdge(40+i_edge2,edge_conectivity);
	      // add edge "in"
	      edge_conectivity[1] = i2;
	      cutmesh.addEdge(10+i_edge2,edge_conectivity);
	      edge_conectivity[0] = i3;
	      edge_conectivity[1] = i1;
	      cutmesh.addEdge(20+i_edge0,edge_conectivity);
	      // separation betwen "in" an "out" edge
	      cutmesh.setEdgeLimite(30+i_edge0);
#endif
	    }
	  break;
	}

        // no cut
      default :
	{
	  return( 0 );
	  break;
	}
      }

    //2D => w=0
    UVW.w = zero;

    //store the constitutive nodes
    // node 2
    UVW.u = zero;
    UVW.v = one;
    cutmesh.addPoint(2,UVW);
    // node 1
    UVW.u = one;
    UVW.v = zero;
    cutmesh.addPoint(1,UVW);
    // node 0
    UVW.u = zero;
    cutmesh.addPoint(0,UVW);

    // normal endding
    return ( 1 );

  }


  ///////////////////
  /// 3D tetrahedron
  ///////////////////
  /********************************************************************************************************************
   * DEFINITION
   * ==========
   *
   *
   *        v
   *        ^
   *        |
   *
   *        2
   *        | \
   *        |  \
   *       2|   \ 1                             face 0
   *        |    \
   *        |     \
   *        0------1    --> u
   *            0
   *
   *        w
   *        ^
   *        |
   *
   *        3
   *        | \
   *        |  \
   *       5|   \ 3                             face 1
   *        |    \
   *        |     \
   *        0------1    --> u
   *            0
   *
   *
   *              w
   *              ^
   *              |
   *
   *              3
   *             / \
   *            /   \
   *          3/     \4                         face 2
   *          /       \
   *         /         \
   *        1-----------2
   *              1
   *
   *
   *
   *        w
   *        ^
   *        |
   *
   *        3
   *        | \
   *        |  \
   *       5|   \ 4                             face 3
   *        |    \
   *        |     \
   *        0------2    --> v
   *            2
   *
   *
   *
   *    mapping u : [0,1] <=> [x0,x1]
   *            v : [0,1] <=> [y0,y2]
   *            w : [0,1] <=> [z0,z3]
   *
   *
   *
   *
   *    edge conectivity :  edge 0 is defined by vertice 0 and 1
   *                        edge 1 is defined by vertice 1 and 2
   *                        edge 2 is defined by vertice 0 and 2
   *                        edge 3 is defined by vertice 1 and 3
   *                        edge 4 is defined by vertice 2 and 3
   *                        edge 5 is defined by vertice 0 and 3
   *
   *
   *    face conectivity : (normal going out of tetaedron )
   *                        face 0 is defined by vertice 0,2 and 1
   *                        face 1 is defined by vertice 0,1 and 3
   *                        face 2 is defined by vertice 1,2 and 3
   *                        face 3 is defined by vertice 2,0 and 3
   *
   *    tetra conectivity : tetra defined by vertice 0, 1, 2 and 3
   *
   *
   * CUT SUMARY
   * ==========
   *
   *    There is 25 known possibles cut of a tetraedron by a iso-zero surface by linear interpolation
   *
   *    All of them can be reduce to 4 references cases by rotations/permutations
   *
   *    The 4 cases differs by there number of edge cut :
   *
   *      case 1 (C1) : 1 edge only is cut, the iso-zero surface passe thru the 2 opposites coins of the tetrahedron
   *      case 2 (C2) : 2 edges are cut, the iso-zero surface passe thru the opposite coin of the tetrahedron
   *      case 3 (C3) : 3 edges are cut, the iso-zero surface passe thru them
   *      case 4 (C4) : 4 edges are cut, the iso-zero surface passe thru them
   *
   * CONVENTIONS
   * ===========
   *
   *  construction of iso-zero is folowing these conventions :
   *      - normal of faces are always from "in" to "out"
   *      - in the case C4 the first 2 nodes of the last iso-zero triangle define intersection edge of those 2 iso-zero triangles
   *  construction of subentities of subentities is folowing these conventions :
   *      - edges conectivity start from cut node of the original edge to nodes of the tetraedron
   *      - faces conectivity start from cut node of the original face/edge to nodes of the tetraedron (nodes of the tetraedron are ending)
   *      - numbering of edges use original edge number + 10 * logical order from "in" to "out" (starting at 1)
   *      - numbering of faces use original face number + 10 * logical order from "in" to "out" (starting at 1)
   *
   *    from iso-zero edges number "i", original edge number is given by i%10
   *    from iso-zero face number "i", original face number is given by i%10
   *
   ********************************************************************************************************************/
  int cutTetRefByLevelSet ( const std::vector < double > & ls_value, const std::vector < int > & node_label, xRefMesh & cutmesh)
  {
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////
    // permutation/rotation for all case
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////
    const unsigned char perm_rot[13][4] = {
      {0, 2, 3, 1},                                // C1 edge 0 cut, C2 edge 0,1 cut, C3 edge 0,1,3 cut
      {1, 0, 3, 2},                                // C1 edge 1 cut, C2 edge 1,2 cut
      {0, 3, 1, 2},                                // C1 edge 2 cut, C2 edge 2,4 cut, C3 edge 1,2,4 cut, C4 edge 0,2,3,4 cut
      {1, 2, 0, 3},                                // C1 edge 3 cut, C2 edge 3,4 cut
      {2, 0, 1, 3},                                // C1 edge 4 cut, C2 edge 4,5 cut, C4 edge 0,1,4,5 cut
      {0, 1, 2, 3},                                // reference case no permutation
      {0, 0, 0, 0},                                // Error case
      {2, 1, 3, 0},                                // C2 edge 0,2 cut, C3 edge 0,2,5 cut
      {3, 0, 2, 1},                                // C2 edge 0,3 cut
      {1, 3, 2, 0},                                // C2 edge 0,5 cut
      {2, 3, 0, 1},                                // C2 edge 1,3 cut
      {3, 1, 0, 2},                                // C2 edge 1,4 cut
      {3, 2, 1, 0}                                 // C2 edge 2,5 cut
    };
#ifdef XREFMESH_WITH_SUB
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////
    // vertex to edge table
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////
    const unsigned char v_to_e[4][4] = {
      // vertex     0  1  2  3
      {6, 0, 2, 5},                                // vertex 0
      {0, 6, 1, 3},                                // vertex 1
      {2, 1, 6, 4},                                // vertex 2
      {5, 3, 4, 6}                                 // vertex 3
    };
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////
    // face permutation/rotation for all case
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////
    const unsigned char perm_face[13][4] = {
      {3, 0, 2, 1},                                // C1 edge 0 cut, C2 edge 0,1 cut, C3 edge 0,1,3 cut
      {1, 0, 3, 2},                                // C1 edge 1 cut, C2 edge 1,2 cut
      {1, 3, 2, 0},                                // C1 edge 2 cut, C2 edge 2,4 cut, C3 edge 1,2,4 cut, C4 edge 0,2,3,4 cut
      {0, 2, 3, 1},                                // C1 edge 3 cut, C2 edge 3,4 cut
      {0, 3, 1, 2},                                // C1 edge 4 cut, C2 edge 4,5 cut, C4 edge 0,1,4,5 cut
      {0, 1, 2, 3},                                // reference case no permutation
      {6, 6, 6, 6},                                // Error case
      {2, 0, 1, 3},                                // C2 edge 0,2 cut, C3 edge 0,2,5 cut
      {3, 1, 0, 2},                                // C2 edge 0,3 cut
      {2, 1, 3, 0},                                // C2 edge 0,5 cut
      {3, 2, 1, 0},                                // C2 edge 1,3 cut
      {1, 2, 0, 3},                                // C2 edge 1,4 cut
      {2, 3, 0, 1}                                 // C2 edge 2,5 cut
    };
#endif
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////
    //  Acces table to permutation table
    //  case 1 no table as  permutation table can be accessed directly
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////
    //  Acces table to permutation table
    //  case 2 : table have 2 index of 5 terms long (no comon index possible = 2 same edge stored in id_cut)
    //           out of this 5X5=25 terms only 12 are relevant
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////
    const unsigned char indirect_c2[5][5] = {
      // edge    1   2   3   4   5  cut
      {  0,  7,  8,  6,  9},                                  // edge 0 cut
      {  6,  1, 10, 11,  6},                                  // edge 1 cut
      {  6,  6,  6,  2, 12},                                  // edge 2 cut
      {  6,  6,  6,  3,  5},                                  // edge 3 cut
      {  6,  6,  6,  6,  4}                                   // edge 4 cut
    };
    const bool sign_c2[13][2] = {
      {true,false},                             // 0
      {false,true},                             // 1
      {true,false},                             // 2
      {false,true},                             // 3
      {false,true},                             // 4
      {false,true},                             // 5
      {false,true},                             // 6
      {true,false},                             // 7
      {true,false},                             // 8
      {true,false},                             // 9
      {true,false},                             // 10
      {true,false},                             // 11
      {true,false}                              // 12
    };
    const unsigned char perm_c2[13][2] = {
      { 5, 4},                         // 0
      { 5, 4},                         // 1
      { 5, 4},                         // 2
      { 5, 4},                         // 3
      { 5, 4},                         // 4
      { 4, 5},                         // 5
      { 4, 5},                         // 6
      { 4, 5},                         // 7
      { 4, 5},                         // 8
      { 5, 4},                         // 9
      { 5, 4},                         // 10
      { 4, 5},                         // 11
      { 4, 5}                          // 12
    };
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////
    //  Acces table to permutation table
    //  case 3 : table have 3 index of 4 terms long (no comon index possible = 2 same edge stored in id_cut)
    //           out of this 4X4x4=64 terms only 4 are relevant
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////
    const unsigned char indirect_c3[4][4][4] = {
      // edge        1                  2                 3                4                       cut
      // edge  2   3   4   5      2   3   4   5     2   3   4   5     2   3   4   5                cut
      { {  6,  0,  6,  6}, {  6,  6,  6,  7},{  6,  6,  6,  6},{  6,  6,  6,  6}, },                                  // edge 0 cut
      { {  6,  6,  6,  6}, {  6,  6,  2,  6},{  6,  6,  6,  6},{  6,  6,  6,  6}, },                                  // edge 1 cut
      { {  6,  6,  6,  6}, {  6,  6,  6,  6},{  6,  6,  6,  6},{  6,  6,  6,  6}, },                                  // edge 2 cut
      { {  6,  6,  6,  6}, {  6,  6,  6,  6},{  6,  6,  6,  6},{  6,  6,  6,  5}, }                                   // edge 3 cut
    };
    const bool sign_c3[8][2] = {
      {true,false},                            // 0
      {false,true},                            // 1
      {true,false},                            // 2
      {false,true},                            // 3
      {false,true},                            // 4
      {false,true},                            // 5
      {false,true},                            // 6
      {true,false}                             // 7
    };
    const unsigned char perm_c3[8][3] = {
      { 5, 6, 4},                         // 0
      { 4, 5, 6},                         // 1
      { 6, 4, 5},                         // 2
      { 4, 5, 6},                         // 3
      { 4, 5, 6},                         // 4
      { 4, 5, 6},                         // 5
      { 4, 5, 6},                         // 6
      { 4, 6, 5}                          // 7
    };
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////
    //  Acces table to permutation table
    //  case 4 : table have 4 index of 3 terms long (no comon index possible = 2 same edge stored in id_cut)
    //           out of this 3X3x3x3=81 terms only 3 are relevant
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////
    const unsigned char indirect_c4[3][3][3][3] = {
      // edge                     1                                           2                                           3                                cut
      // edge       2             3            4                2             3            4                2             3            4                   cut
      // edge    3  4  5       3  4  5      3  4  5          3  4  5       3  4  5      3  4  5          3  4  5       3  4  5      3  4  5                cut
      {  { {  6, 6, 6},  {  6, 6, 6}, {  6, 6, 4} }, { {  6, 6, 6},  {  6, 2, 6}, {  6, 6, 6} }, { {  6, 6, 6},  {  6, 6, 6}, {  6, 6, 6} } }, // edge 0 cut
      {  { {  6, 6, 6},  {  6, 6, 6}, {  6, 6, 6} }, { {  6, 6, 6},  {  6, 6, 5}, {  6, 6, 6} }, { {  6, 6, 6},  {  6, 6, 6}, {  6, 6, 6} } }, // edge 1 cut
      {  { {  6, 6, 6},  {  6, 6, 6}, {  6, 6, 6} }, { {  6, 6, 6},  {  6, 6, 6}, {  6, 6, 6} }, { {  6, 6, 6},  {  6, 6, 6}, {  6, 6, 6} } } // edge 2 cut
    };
    const bool sign_c4[6][2] = {
      {false,true},                            // not used
      {false,true},                            // not used
      {true,false},                            // C4 edge 0,2,3,4 cut
      {false,true},                            // not used
      {false,true},                            // C4 edge 0,1,4,5 cut
      {false,true}                             // reference case no permutation
    };
    const unsigned char perm_c4[6][4] = {
      { 4, 5, 6, 7},               // not used
      { 4, 5, 6, 7},               // not used
      { 6, 4, 7, 5},               // C4 edge 0,2,3,4 cut
      { 4, 5, 6, 7},               // not used
      { 4, 5, 7, 6},               // C4 edge 0,1,4,5 cut
      { 4, 5, 6, 7}                // reference case no permutation
    };
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////
    //
    // local
    unsigned char nb_cut = 0;
    unsigned char id_cut[4];
    unsigned char sign_edge[4];
    xPoint cut_point[4];
    double lsv0 = ls_value[0];
    double lsv1 = ls_value[1];
    double lsv2 = ls_value[2];
    double lsv3 = ls_value[3];
    const double one = 1.0E+0;
    const double zero = 0.0E+0;
    double r;
    xRefMesh::elemdef_t conectivity(4);
    xRefMesh::elemdef_t face_conectivity(3);
    xRefMesh::elemdef_t edge_conectivity(2);


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Determination of the cutting points ////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////

    // the iso-zero cut the edge 0 ////////////////////////////////////////////////////////////////////////////
    if ( ( lsv0*lsv1 ) < 0 )
      {

        // lineare interpolation
        cut_point[nb_cut].u = -lsv0/( lsv1 - lsv0 );
        cut_point[nb_cut].v = zero;
        cut_point[nb_cut].w = zero;

        sign_edge[nb_cut] = ( lsv0 < zero ) ? 1 : 0;
        id_cut[nb_cut] = 0;
        ++nb_cut;
      }

    // the iso-zero cut the edge 1 ////////////////////////////////////////////////////////////////////////////
    if ( ( lsv1*lsv2 ) < 0 )
      {
        // lineare interpolation
        cut_point[nb_cut].v = r = -lsv1/( lsv2 - lsv1 );
        cut_point[nb_cut].u = one - r;
        cut_point[nb_cut].w = zero;

        sign_edge[nb_cut] = ( lsv1 < zero ) ? 1 : 0;
        id_cut[nb_cut] = 1;
        ++nb_cut;
      }

    // the iso-zero cut the edge 2 ////////////////////////////////////////////////////////////////////////////
    if ( ( lsv0*lsv2 ) < 0 )
      {
        // lineare interpolation
        cut_point[nb_cut].u = zero;
        cut_point[nb_cut].v = -lsv0/( lsv2 - lsv0 );
        cut_point[nb_cut].w = zero;

        sign_edge[nb_cut] = ( lsv0 < zero ) ? 1 : 0;
        id_cut[nb_cut] = 2;
        ++nb_cut;
      }

    // the iso-zero cut the edge 3 ////////////////////////////////////////////////////////////////////////////
    if ( ( lsv1*lsv3 ) < 0 )
      {
        // lineare interpolation
        cut_point[nb_cut].w = r = -lsv1/( lsv3 - lsv1 );
        cut_point[nb_cut].u = one - r;
        cut_point[nb_cut].v = zero;

        sign_edge[nb_cut] = ( lsv1 < zero ) ? 1 : 0;
        id_cut[nb_cut] = 3;
        ++nb_cut;
      }

    // the iso-zero cut the edge 4 ////////////////////////////////////////////////////////////////////////////
    if ( ( lsv2*lsv3 ) < 0 )
      {
        // lineare interpolation
        cut_point[nb_cut].w = r = -lsv2/( lsv3 - lsv2 );
        cut_point[nb_cut].v = one - r;
        cut_point[nb_cut].u = zero;

        sign_edge[nb_cut] = ( lsv2 < zero ) ? 1 : 0;
        id_cut[nb_cut] = 4;
        ++nb_cut;
      }

    // the iso-zero cut the edge 5 ////////////////////////////////////////////////////////////////////////////
    if ( ( lsv0*lsv3 ) < 0 )
      {
        // lineare interpolation
        cut_point[nb_cut].w = -lsv0/( lsv3 - lsv0 );
        cut_point[nb_cut].u = zero;
        cut_point[nb_cut].v = zero;

        sign_edge[nb_cut] = ( lsv0 < zero ) ? 1 : 0;
        id_cut[nb_cut] = 5;
        ++nb_cut;
      }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////
    // creation of the  cut mesh //////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////
    switch (nb_cut)
      {
        //  1 edge only is cut, the iso-zero surface passe thru the 2 opposites coins of the tetrahedron
        //  submesh => 2 tetrahedrons; iso-zero => 1 triangle
        //  subedge => 1 edge cut => 2 edges
        //  subface => 2 faces cut => 4 faces
        // based on the the reference case : edge 5 cut, pass thru node 1,2
      case 1 :
	{
	  // local
	  // the first 6 row of perm_rot can be accessed directely
	  const unsigned char i_edge = id_cut[0];
	  const unsigned char * perm_node = &( perm_rot[i_edge][0] );
	  const unsigned char i0 = perm_node[0];
	  const unsigned char i1 = perm_node[1];
	  const unsigned char i2 = perm_node[2];
	  const unsigned char i3 = perm_node[3];


	  // adding the new node
	  cutmesh.addPoint(4,cut_point[0]);

	  // face start in all case by 4
	  face_conectivity[0] = 4;

#ifdef XREFMESH_WITH_SUB


	  // set face number
	  const unsigned char i_face1 = perm_face[i_edge][1];
	  const unsigned char i_face3 = perm_face[i_edge][3];

	  // separation betwen "in" an "out" edge
	  cutmesh.setEdgeLimite(20+i_edge);

	  // separation betwen "in" an "out" face
	  cutmesh.setFaceLimite(30+i_face3);

	  // iso zero face begin (edges)
	  cutmesh.setIsozeroFace(50+i_face1);

	  // add iso zero edge to faces container
	  // first add edge not cutting a face : the iso-zero edge corresponding to edge i1,i2
	  // is id is k*10 + 4 + id edge 
	  // 4 is a shift to avoid to use this as a id of a face
	  edge_conectivity[0] = i1;
	  edge_conectivity[1] = i2;
	  cutmesh.addFace(70+v_to_e[i1][i2]+4,edge_conectivity);
	  // edge start in all case by 4 now
	  edge_conectivity[0] = 4;
	  // second add edge cutting a face
	  // i4,i2
	  cutmesh.addFace(60+i_face3,edge_conectivity);
	  edge_conectivity[1] = i1;
	  cutmesh.addFace(50+i_face1,edge_conectivity);


#endif

	  //
	  // node 0 is "in"
	  if (sign_edge[0])
	    {
	      // add element "in"
	      conectivity[0] = i0;
	      conectivity[1] = i1;
	      conectivity[2] = i2;
	      conectivity[3] = 4;
	      cutmesh.addElem(0,conectivity);
	      // separation betwen "in" an "out" element
	      cutmesh.setLimite(1);
	      // add element "out"
	      conectivity[0] = 4;
	      conectivity[3] = i3;
	      cutmesh.addElem(1,conectivity);
	      // iso zero element is at the end
	      cutmesh.setIsozero(2);
	      // add iso zero element
	      face_conectivity[1] = i1;
	      face_conectivity[2] = i2;
	      cutmesh.addElem(2,face_conectivity);

#ifdef XREFMESH_WITH_SUB
	      // add edge "in"
	      edge_conectivity[1] = i0;
	      cutmesh.addEdge(10+i_edge,edge_conectivity);
	      // add edge "out"
	      edge_conectivity[1] = i3;
	      cutmesh.addEdge(20+i_edge,edge_conectivity);


	      // add face "in"
	      face_conectivity[1] = i0;
	      face_conectivity[2] = i1;
	      cutmesh.addFace(10+i_face1,face_conectivity);
	      face_conectivity[1] = i2;
	      face_conectivity[2] = i0;
	      cutmesh.addFace(20+i_face3,face_conectivity);
	      // add face "out"
	      face_conectivity[1] = i3;
	      face_conectivity[2] = i2;
	      cutmesh.addFace(30+i_face3,face_conectivity);
	      face_conectivity[1] = i1;
	      face_conectivity[2] = i3;
	      cutmesh.addFace(40+i_face1,face_conectivity);
#endif
	    }
	  else
	    {
	      // add element "in"
	      conectivity[0] = 4;
	      conectivity[1] = i1;
	      conectivity[2] = i2;
	      conectivity[3] = i3;
	      cutmesh.addElem(0,conectivity);
	      // separation betwen "in" an "out" element
	      cutmesh.setLimite(1);
	      // add element "out"
	      conectivity[0] = i0;
	      conectivity[3] = 4;
	      cutmesh.addElem(1,conectivity);
	      // iso zero element is at the end
	      cutmesh.setIsozero(2);
	      // add iso zero element
	      face_conectivity[1] = i2;
	      face_conectivity[2] = i1;
	      cutmesh.addElem(2,face_conectivity);

#ifdef XREFMESH_WITH_SUB
	      // add edge "in"
	      edge_conectivity[1] = i3;
	      cutmesh.addEdge(10+i_edge,edge_conectivity);
	      // add edge "out"
	      edge_conectivity[1] = i0;
	      cutmesh.addEdge(20+i_edge,edge_conectivity);

	      // add face "in"
	      face_conectivity[1] = i1;
	      face_conectivity[2] = i3;
	      cutmesh.addFace(10+i_face1,face_conectivity);
	      face_conectivity[1] = i3;
	      face_conectivity[2] = i2;
	      cutmesh.addFace(20+i_face3,face_conectivity);
	      // add face "out"
	      face_conectivity[1] = i2;
	      face_conectivity[2] = i0;
	      cutmesh.addFace(30+i_face3,face_conectivity);
	      face_conectivity[1] = i0;
	      face_conectivity[2] = i1;
	      cutmesh.addFace(40+i_face1,face_conectivity);
#endif
	    }


	  break;
	}
        //  2 edges are cut, the iso-zero surface passe thru the opposite coin of the tetrahedron
        //  submesh => 3 tetrahedrons; iso-zero => 1 triangle
        //  subedge => 2 edge cut => 4 edges
        //  subface => 3 faces cut => 7 faces
        // based on the the reference case : edge 3,5 cut, pass thru node 2
      case 2 :
	{
	  // local
	  const unsigned char i = indirect_c2[id_cut[0]][id_cut[1]-1];
	  const unsigned char * perm_node = &( perm_rot[i][0] );
	  const unsigned char i0 = perm_node[0];
	  const unsigned char i1 = perm_node[1];
	  const unsigned char i2 = perm_node[2];
	  const unsigned char i3 = perm_node[3];
	  const unsigned char i4 = perm_c2[i][0];
	  const unsigned char i5 = perm_c2[i][1];
	  unsigned char n1,n2,n3,n4;

	  // adding the new nodes
	  cutmesh.addPoint(4,cut_point[0]);
	  cutmesh.addPoint(5,cut_point[1]);

#ifdef XREFMESH_WITH_SUB
	  // set edge number
	  const unsigned char i_edge3 = v_to_e[i1][i3];
	  const unsigned char i_edge5 = v_to_e[i0][i3];

	  // set face number
	  const unsigned char i_face1 = perm_face[i][1];
	  const unsigned char i_face2 = perm_face[i][2];
	  const unsigned char i_face3 = perm_face[i][3];

	  //conform
	  unsigned char n5;

	  // separation betwen "in" an "out" edge
	  cutmesh.setEdgeLimite(30+i_edge5);

	  // iso zero face begin (edges)
	  cutmesh.setIsozeroFace(80+i_face2);

	  // add iso zero edge to faces container
	  edge_conectivity[0] = i4;
	  edge_conectivity[1] = i2;
	  cutmesh.addFace(80+i_face2,edge_conectivity);
	  edge_conectivity[0] = i5;
	  cutmesh.addFace(90+i_face3,edge_conectivity);
	  edge_conectivity[1] = i4;
	  cutmesh.addFace(100+i_face1,edge_conectivity);
#endif

	  // for conform mesh use node label to follow a uniforme cut strategie among 2 tetraedron relate by the same face
	  // use of the lowest label as starting point of the cut
	  if (node_label[i0] < node_label[i1])
	    {
	      // i0 have the lowest label
	      n1 = i4;
	      n2 = i0;
	      n3 = i4;
	      n4 = i5;
#ifdef XREFMESH_WITH_SUB
	      n5 = i0;
#endif
	    }
	  else
	    {
	      // i1 have the lowest label
	      n1 = i5;
	      n2 = i5;
	      n3 = i1;
	      n4 = i4;
#ifdef XREFMESH_WITH_SUB
	      n5 = i1;
#endif
	    }

	  //
	  // node 0 is "in"
	  if (sign_c2[i][sign_edge[1]])
	    {
	      // add elements "in"
	      conectivity[0] = i0;
	      conectivity[1] = i1;
	      conectivity[2] = i2;
	      conectivity[3] = n1;
	      cutmesh.addElem(0,conectivity);
	      conectivity[0] = n2;
	      conectivity[1] = n3;
	      conectivity[3] = n4;
	      cutmesh.addElem(1,conectivity);
	      // separation betwen "in" an "out" element
	      cutmesh.setLimite(2);
	      // add element "out"
	      conectivity[0] = i5;
	      conectivity[1] = i4;
	      conectivity[3] = i3;
	      cutmesh.addElem(2,conectivity);
	      // iso zero element is at the end
	      cutmesh.setIsozero(3);
	      // add iso zero element
	      face_conectivity[0] = i5;
	      face_conectivity[1] = i4;
	      face_conectivity[2] = i2;
	      cutmesh.addElem(3,face_conectivity);

#ifdef XREFMESH_WITH_SUB
	      // add edge "in"
	      edge_conectivity[0] = i5;
	      edge_conectivity[1] = i0;
	      cutmesh.addEdge(10+i_edge5,edge_conectivity);
	      // add edge "out"
	      edge_conectivity[1] = i3;
	      cutmesh.addEdge(30+i_edge5,edge_conectivity);
	      edge_conectivity[0] = i4;
	      cutmesh.addEdge(40+i_edge3,edge_conectivity);
	      // add edge "in"
	      edge_conectivity[1] = i1;
	      cutmesh.addEdge(20+i_edge3,edge_conectivity);

	      // separation betwen "in" an "out" face
	      cutmesh.setFaceLimite(50+i_face3);

	      // starting from i5
	      // add face "in"
	      face_conectivity[0] = i5;
	      face_conectivity[1] = i2;
	      face_conectivity[2] = i0;
	      cutmesh.addFace(10+i_face3,face_conectivity);
	      // add face "out"
	      face_conectivity[1] = i3;
	      face_conectivity[2] = i2;
	      cutmesh.addFace(50+i_face3,face_conectivity);
	      face_conectivity[1] = i4;
	      face_conectivity[2] = i3;
	      cutmesh.addFace(60+i_face1,face_conectivity);

	      // starting from i5 or i4 depending on cut choise
	      // add face "in"
	      face_conectivity[0] = n1;
	      face_conectivity[1] = i0;
	      face_conectivity[2] = i1;
	      cutmesh.addFace(20+i_face1,face_conectivity);

	      // starting from i4
	      // add face "in"
	      face_conectivity[0] = i4;
	      face_conectivity[1] = i5;
	      face_conectivity[2] = n5;
	      cutmesh.addFace(30+i_face1,face_conectivity);
	      face_conectivity[1] = i1;
	      face_conectivity[2] = i2;
	      cutmesh.addFace(40+i_face2,face_conectivity);
	      // add face "out"
	      face_conectivity[1] = i2;
	      face_conectivity[2] = i3;
	      cutmesh.addFace(70+i_face2,face_conectivity);
#endif
	    }
	  else
	    {
	      // add elements "in"
	      conectivity[0] = i5;
	      conectivity[1] = i4;
	      conectivity[2] = i2;
	      conectivity[3] = i3;
	      cutmesh.addElem(0,conectivity);
	      // separation betwen "in" an "out" element
	      cutmesh.setLimite(1);
	      // add element "out"
	      conectivity[0] = i0;
	      conectivity[1] = i1;
	      conectivity[3] = n1;
	      cutmesh.addElem(1,conectivity);
	      conectivity[0] = n2;
	      conectivity[1] = n3;
	      conectivity[3] = n4;
	      cutmesh.addElem(2,conectivity);
	      // iso zero element is at the end
	      cutmesh.setIsozero(3);
	      // add iso zero element
	      face_conectivity[0] = i4;
	      face_conectivity[1] = i5;
	      face_conectivity[2] = i2;
	      cutmesh.addElem(3,face_conectivity);

#ifdef XREFMESH_WITH_SUB
	      // add edge "out"
	      edge_conectivity[0] = i5;
	      edge_conectivity[1] = i0;
	      cutmesh.addEdge(30+i_edge5,edge_conectivity);
	      // add edge "in"
	      edge_conectivity[1] = i3;
	      cutmesh.addEdge(10+i_edge5,edge_conectivity);
	      edge_conectivity[0] = i4;
	      cutmesh.addEdge(20+i_edge3,edge_conectivity);
	      // add edge "out"
	      edge_conectivity[1] = i1;
	      cutmesh.addEdge(40+i_edge3,edge_conectivity);

	      // separation betwen "in" an "out" face
	      cutmesh.setFaceLimite(40+i_face3);

	      // starting from i5
	      // add face "out"
	      face_conectivity[0] = i5;
	      face_conectivity[1] = i2;
	      face_conectivity[2] = i0;
	      cutmesh.addFace(40+i_face3,face_conectivity);
	      // add face "in"
	      face_conectivity[1] = i3;
	      face_conectivity[2] = i2;
	      cutmesh.addFace(10+i_face3,face_conectivity);
	      face_conectivity[1] = i4;
	      face_conectivity[2] = i3;
	      cutmesh.addFace(20+i_face1,face_conectivity);

	      // starting from i5 or i4 depending on cut choise
	      // add face "out"
	      face_conectivity[0] = n1;
	      face_conectivity[1] = i0;
	      face_conectivity[2] = i1;
	      cutmesh.addFace(50+i_face1,face_conectivity);

	      // starting from i4
	      // add face "out"
	      face_conectivity[0] = i4;
	      face_conectivity[1] = i5;
	      face_conectivity[2] = n5;
	      cutmesh.addFace(60+i_face1,face_conectivity);
	      face_conectivity[1] = i1;
	      face_conectivity[2] = i2;
	      cutmesh.addFace(70+i_face2,face_conectivity);
	      // add face "in"
	      face_conectivity[1] = i2;
	      face_conectivity[2] = i3;
	      cutmesh.addFace(30+i_face2,face_conectivity);

#endif
	    }
	  break;
	}
        //  3 edges are cut, the iso-zero surface passe thru them
        //  submesh => 4 tetrahedrons; iso-zero => 1 triangle
        //  subedge => 3 edge cut => 6 edges
        //  subface => 3 faces cut => 9 faces
        // based on the the reference case : edge 3,4,5 cut
      case 3 :
	{

	  // local
	  const unsigned char i = indirect_c3[id_cut[0]][id_cut[1]-1][id_cut[2]-2];
	  const unsigned char * perm_node = &( perm_rot[i][0] );
	  const unsigned char i0 = perm_node[0];
	  const unsigned char i1 = perm_node[1];
	  const unsigned char i2 = perm_node[2];
	  const unsigned char i3 = perm_node[3];
	  const unsigned char i4 = perm_c3[i][0];
	  const unsigned char i5 = perm_c3[i][1];
	  const unsigned char i6 = perm_c3[i][2];
	  unsigned char I0,I1,I2,I4,I5,I6;
	  unsigned char n1,n2,n3,n4;

	  // adding the new nodes
	  cutmesh.addPoint(4,cut_point[0]);
	  cutmesh.addPoint(5,cut_point[1]);
	  cutmesh.addPoint(6,cut_point[2]);

#ifdef XREFMESH_WITH_SUB
	  // set edge number
	  const unsigned char i_edge3 = v_to_e[i1][i3];
	  const unsigned char i_edge4 = v_to_e[i2][i3];
	  const unsigned char i_edge5 = v_to_e[i0][i3];

	  // set face number
	  const unsigned char i_face1 = perm_face[i][1];
	  const unsigned char i_face2 = perm_face[i][2];
	  const unsigned char i_face3 = perm_face[i][3];

	  unsigned char n5,I_face1,I_face2,I_face3;

	  // separation betwen "in" an "out" edge
	  cutmesh.setEdgeLimite(40+i_edge5);

	  // iso zero face begin (edges)
	  cutmesh.setIsozeroFace(100+i_face1);

	  // add iso zero edge to faces container
	  edge_conectivity[0] = i6;
	  edge_conectivity[1] = i4;
	  cutmesh.addFace(100+i_face1,edge_conectivity);
	  edge_conectivity[0] = i5;
	  cutmesh.addFace(110+i_face2,edge_conectivity);
	  edge_conectivity[1] = i6;
	  cutmesh.addFace(130+i_face3,edge_conectivity);
#endif

	  // for conform mesh use node label to follow a uniforme cut strategie among 2 tetraedron relate by the same face
	  // use of the lowest label as starting point of the cut
	  // permutation from actual to reference which is i1 having the minimum label
	  if (node_label[i0] < node_label[i1])
	    {
	      if (node_label[i0] < node_label[i2])
		{
		  // i0 is the minimum
		  I0 = i2;
		  I1 = i0;
		  I2 = i1;
		  I4 = i6;
		  I5 = i4;
		  I6 = i5;
#ifdef XREFMESH_WITH_SUB
		  I_face1 = i_face3;
		  I_face2 = i_face1;
		  I_face3 = i_face2;
#endif
		  if (node_label[I0] < node_label[I2])
		    {
		      n1 = I5;
		      n2 = I0;
		      n3 = I5;
		      n4 = I6;
#ifdef XREFMESH_WITH_SUB
		      n5 = I0;
#endif
		    }
		  else
		    {
		      n1 = I6;
		      n2 = I6;
		      n3 = I2;
		      n4 = I5;
#ifdef XREFMESH_WITH_SUB
		      n5 = I2;
#endif
		    }
		}
	      else
		{
		  // i2 is the minimum and i0<i1
		  I0 = i1;
		  I1 = i2;
		  I2 = i0;
		  I4 = i5;
		  I5 = i6;
		  I6 = i4;
#ifdef XREFMESH_WITH_SUB
		  I_face1 = i_face2;
		  I_face2 = i_face3;
		  I_face3 = i_face1;
		  n5 = I2;
#endif
		  n1 = I6;
		  n2 = I6;
		  n3 = I2;
		  n4 = I5;
		}
	    }
	  else if (node_label[i1] < node_label[i2])
	    {
	      // i1 is the minimum
	      I0 = i0;
	      I1 = i1;
	      I2 = i2;
	      I4 = i4;
	      I5 = i5;
	      I6 = i6;
#ifdef XREFMESH_WITH_SUB
	      I_face1 = i_face1;
	      I_face2 = i_face2;
	      I_face3 = i_face3;
#endif
	      if (node_label[i0] < node_label[i2])
		{
		  n1 = I5;
		  n2 = I0;
		  n3 = I5;
		  n4 = I6;
#ifdef XREFMESH_WITH_SUB
		  n5 = I0;
#endif
		}
	      else
		{
		  n1 = I6;
		  n2 = I6;
		  n3 = I2;
		  n4 = I5;
#ifdef XREFMESH_WITH_SUB
		  n5 = I2;
#endif
		}
	    }
	  else
	    {
	      // i2 is the minimum and i1<i0
	      I0 = i1;
	      I1 = i2;
	      I2 = i0;
	      I4 = i5;
	      I5 = i6;
	      I6 = i4;
#ifdef XREFMESH_WITH_SUB
	      I_face1 = i_face2;
	      I_face2 = i_face3;
	      I_face3 = i_face1;
	      n5 = I0;
#endif
	      n1 = I5;
	      n2 = I0;
	      n3 = I5;
	      n4 = I6;
	    }

	  //
	  // node 0 is "in"
	  if (sign_c3[i][sign_edge[2]])
	    {
	      // add elements "in"
	      conectivity[0] = I0;
	      conectivity[1] = I1;
	      conectivity[2] = I2;
	      conectivity[3] = n1;
	      cutmesh.addElem(0,conectivity);
	      conectivity[0] = n2;
	      conectivity[2] = n3;
	      conectivity[3] = n4;
	      cutmesh.addElem(1,conectivity);
	      conectivity[0] = I6;
	      conectivity[2] = I5;
	      conectivity[3] = I4;
	      cutmesh.addElem(2,conectivity);
	      // separation betwen "in" an "out" element
	      cutmesh.setLimite(3);
	      // add element "out"
	      conectivity[1] = I4;
	      conectivity[3] = i3;
	      cutmesh.addElem(3,conectivity);
	      // iso zero element is at the end
	      cutmesh.setIsozero(4);
	      // add iso zero element
	      face_conectivity[0] = i6;
	      face_conectivity[1] = i4;
	      face_conectivity[2] = i5;
	      cutmesh.addElem(4,face_conectivity);

#ifdef XREFMESH_WITH_SUB
	      // add edge "in"
	      edge_conectivity[0] = i6;
	      edge_conectivity[1] = i0;
	      cutmesh.addEdge(10+i_edge5,edge_conectivity);
	      // add edge "out"
	      edge_conectivity[1] = i3;
	      cutmesh.addEdge(40+i_edge5,edge_conectivity);
	      edge_conectivity[0] = i4;
	      cutmesh.addEdge(50+i_edge3,edge_conectivity);
	      edge_conectivity[0] = i5;
	      cutmesh.addEdge(60+i_edge4,edge_conectivity);
	      // add edge "in"
	      edge_conectivity[1] = i2;
	      cutmesh.addEdge(20+i_edge4,edge_conectivity);
	      edge_conectivity[0] = i4;
	      edge_conectivity[1] = i1;
	      cutmesh.addEdge(30+i_edge3,edge_conectivity);

	      // separation betwen "in" an "out" face
	      cutmesh.setFaceLimite(70+I_face1);

	      // starting from i6
	      // add face "in"
	      face_conectivity[0] = I6;
	      face_conectivity[1] = I5;
	      face_conectivity[2] = n5;
	      cutmesh.addFace(10+I_face3,face_conectivity);
	      face_conectivity[1] = I0;
	      face_conectivity[2] = I1;
	      cutmesh.addFace(20+I_face1,face_conectivity);
	      // add face "out"
	      face_conectivity[1] = I4;
	      face_conectivity[2] = i3;
	      cutmesh.addFace(70+I_face1,face_conectivity);

	      // starting from i5 or i6 depending on cut choise
	      // add face "in"
	      face_conectivity[0] = n1;
	      face_conectivity[1] = I2;
	      face_conectivity[2] = I0;
	      cutmesh.addFace(30+I_face3,face_conectivity);

	      // starting from i5
	      // add face "in"
	      face_conectivity[0] = I5;
	      face_conectivity[1] = I1;
	      face_conectivity[2] = I2;
	      cutmesh.addFace(40+I_face2,face_conectivity);
	      face_conectivity[1] = I4;
	      face_conectivity[2] = I1;
	      cutmesh.addFace(50+I_face2,face_conectivity);
	      // add face "out"
	      face_conectivity[1] = I6;
	      face_conectivity[2] = i3;
	      cutmesh.addFace(80+I_face3,face_conectivity);

	      // starting from i4
	      // add face "in"
	      face_conectivity[0] = I4;
	      face_conectivity[1] = I6;
	      face_conectivity[2] = I1;
	      cutmesh.addFace(60+I_face1,face_conectivity);
	      // add face "out"
	      face_conectivity[1] = I5;
	      face_conectivity[2] = i3;
	      cutmesh.addFace(90+I_face2,face_conectivity);
#endif
	    }
	  else
	    {
	      // add elements "in"
	      conectivity[0] = I6;
	      conectivity[1] = I4;
	      conectivity[2] = I5;
	      conectivity[3] = i3;
	      cutmesh.addElem(0,conectivity);
	      // separation betwen "in" an "out" element
	      cutmesh.setLimite(1);
	      // add element "out"
	      conectivity[1] = I1;
	      conectivity[3] = I4;
	      cutmesh.addElem(1,conectivity);
	      conectivity[0] = I0;
	      conectivity[2] = I2;
	      conectivity[3] = n1;
	      cutmesh.addElem(2,conectivity);
	      conectivity[0] = n2;
	      conectivity[2] = n3;
	      conectivity[3] = n4;
	      cutmesh.addElem(3,conectivity);
	      // iso zero element is at the end
	      cutmesh.setIsozero(4);
	      // add iso zero element
	      face_conectivity[0] = i6;
	      face_conectivity[1] = i5;
	      face_conectivity[2] = i4;
	      cutmesh.addElem(4,face_conectivity);

#ifdef XREFMESH_WITH_SUB
	      // add edge "out"
	      edge_conectivity[0] = i6;
	      edge_conectivity[1] = i0;
	      cutmesh.addEdge(40+i_edge5,edge_conectivity);
	      // add edge "in"
	      edge_conectivity[1] = i3;
	      cutmesh.addEdge(10+i_edge5,edge_conectivity);
	      edge_conectivity[0] = i4;
	      cutmesh.addEdge(20+i_edge3,edge_conectivity);
	      edge_conectivity[0] = i5;
	      cutmesh.addEdge(30+i_edge4,edge_conectivity);
	      // add edge "out"
	      edge_conectivity[1] = i2;
	      cutmesh.addEdge(50+i_edge4,edge_conectivity);
	      edge_conectivity[0] = i4;
	      edge_conectivity[1] = i1;
	      cutmesh.addEdge(60+i_edge3,edge_conectivity);

	      // separation betwen "in" an "out" face
	      cutmesh.setFaceLimite(40+I_face3);

	      // starting from i6
	      // add face "out"
	      face_conectivity[0] = I6;
	      face_conectivity[1] = I5;
	      face_conectivity[2] = n5;
	      cutmesh.addFace(40+I_face3,face_conectivity);
	      face_conectivity[1] = I0;
	      face_conectivity[2] = I1;
	      cutmesh.addFace(50+I_face1,face_conectivity);
	      // add face "in"
	      face_conectivity[1] = I4;
	      face_conectivity[2] = i3;
	      cutmesh.addFace(10+I_face1,face_conectivity);

	      // starting from i5 or i6 depending on cut choise
	      // add face "out"
	      face_conectivity[0] = n1;
	      face_conectivity[1] = I2;
	      face_conectivity[2] = I0;
	      cutmesh.addFace(60+I_face3,face_conectivity);

	      // starting from i5
	      // add face "out"
	      face_conectivity[0] = I5;
	      face_conectivity[1] = I1;
	      face_conectivity[2] = I2;
	      cutmesh.addFace(70+I_face2,face_conectivity);
	      face_conectivity[1] = I4;
	      face_conectivity[2] = I1;
	      cutmesh.addFace(80+I_face2,face_conectivity);
	      // add face "in"
	      face_conectivity[1] = I6;
	      face_conectivity[2] = i3;
	      cutmesh.addFace(20+I_face3,face_conectivity);

	      // starting from i4
	      // add face "out"
	      face_conectivity[0] = I4;
	      face_conectivity[1] = I6;
	      face_conectivity[2] = I1;
	      cutmesh.addFace(90+I_face1,face_conectivity);
	      // add face "in"
	      face_conectivity[1] = I5;
	      face_conectivity[2] = i3;
	      cutmesh.addFace(30+I_face2,face_conectivity);
#endif
	    }
	  break;
	}
        //   4 edges are cut, the iso-zero surface passe thru them
        //  submesh => 6 tetrahedrons; iso-zero => 2 triangle
        //  subedge => 4 edge cut => 8 edges
        //  subface => 4 faces cut => 9 faces
        // based on the the reference case : edge 1,2,3,5 cut
      case 4 :
	{

	  // local
	  const unsigned char i = indirect_c4[id_cut[0]][id_cut[1]-1][id_cut[2]-2][id_cut[3]-3];
	  const unsigned char * perm_node = &( perm_rot[i][0] );
	  const unsigned char i0 = perm_node[0];
	  const unsigned char i1 = perm_node[1];
	  const unsigned char i2 = perm_node[2];
	  const unsigned char i3 = perm_node[3];
	  const unsigned char i4 = perm_c4[i][0];
	  const unsigned char i5 = perm_c4[i][1];
	  const unsigned char i6 = perm_c4[i][2];
	  const unsigned char i7 = perm_c4[i][3];
	  unsigned char n1,n2,n3,n4,n5,n6;
	  unsigned char m1,m2,m3,m4,m5;


	  // adding the new nodes
	  cutmesh.addPoint(4,cut_point[0]);
	  cutmesh.addPoint(5,cut_point[1]);
	  cutmesh.addPoint(6,cut_point[2]);
	  cutmesh.addPoint(7,cut_point[3]);

#ifdef XREFMESH_WITH_SUB
	  // set edge number
	  const unsigned char i_edge1 = v_to_e[i1][i2];
	  const unsigned char i_edge2 = v_to_e[i0][i2];
	  const unsigned char i_edge3 = v_to_e[i1][i3];
	  const unsigned char i_edge5 = v_to_e[i0][i3];

	  // set face number
	  const unsigned char i_face0 = perm_face[i][0];
	  const unsigned char i_face1 = perm_face[i][1];
	  const unsigned char i_face2 = perm_face[i][2];
	  const unsigned char i_face3 = perm_face[i][3];

	  // conform topology
	  unsigned char n7,m6;
	  // separation betwen "in" an "out" edge
	  cutmesh.setEdgeLimite(50+i_edge1);

	  // separation betwen "in" an "out" face
	  cutmesh.setFaceLimite(70+i_face1);

	  // iso zero face begin (edges)
	  cutmesh.setIsozeroFace(130+i_face1);

	  // add iso zero edge to faces container
	  edge_conectivity[0] = i7;
	  edge_conectivity[1] = i6;
	  cutmesh.addFace(130+i_face1,edge_conectivity);
	  edge_conectivity[0] = i4;
	  cutmesh.addFace(140+i_face2,edge_conectivity);
	  edge_conectivity[1] = i5;
	  cutmesh.addFace(150+i_face0,edge_conectivity);
	  edge_conectivity[0] = i7;
	  cutmesh.addFace(160+i_face3,edge_conectivity);

#endif
	  // for conform mesh use node label to follow a uniforme cut strategie among 2 tetraedron relate by the same face
	  // use of the lowest label as starting point of the cut
	  if (node_label[i0] < node_label[i1])
	    {
	      // i0 have the lowest label
	      n1 = i0;
	      n2 = i1;
	      n3 = i4;
	      n4 = i0;
	      n5 = i4;
	      n6 = i6;
#ifdef XREFMESH_WITH_SUB
	      n7 = i6;
#endif
	    }
	  else
	    {
	      // i1 have the lowest label
	      n1 = i1;
	      n2 = i4;
	      n3 = i5;
	      n4 = i7;
	      n5 = i1;
	      n6 = i1;
#ifdef XREFMESH_WITH_SUB
	      n7 = i7;
#endif
	    }
	  if (node_label[i2] < node_label[i3])
	    {
	      // i2 have the lowest label
	      m1 = i2;
	      m2 = i7;
	      m3 = i6;
	      m4 = i2;
	      m5 = i3;
#ifdef XREFMESH_WITH_SUB
	      m6 = i6;
#endif
	    }
	  else
	    {
	      // i3 have the lowest label
	      m1 = i3;
	      m2 = i5;
	      m3 = i3;
	      m4 = i4;
	      m5 = i2;
#ifdef XREFMESH_WITH_SUB
	      m6 = i4;
#endif
	    }

	  //
	  // node 0 is "in"
	  if (sign_c4[i][sign_edge[3]])
	    {
	      // add elements "in"
	      conectivity[0] = n1;
	      conectivity[1] = n2;
	      conectivity[2] = n3;
	      conectivity[3] = i6;
	      cutmesh.addElem(0,conectivity);
	      conectivity[0] = n4;
	      conectivity[1] = n5;
	      conectivity[2] = i5;
	      cutmesh.addElem(1,conectivity);
	      conectivity[0] = i0;
	      conectivity[1] = n6;
	      conectivity[3] = i7;
	      cutmesh.addElem(2,conectivity);
	      // separation betwen "in" an "out" element
	      cutmesh.setLimite(3);
	      // add element "out"
	      conectivity[0] = i7;
	      conectivity[1] = i6;
	      conectivity[3] = m1;
	      cutmesh.addElem(3,conectivity);
	      conectivity[0] = i5;
	      conectivity[2] = i4;
	      cutmesh.addElem(4,conectivity);
	      conectivity[0] = m2;
	      conectivity[1] = m3;
	      conectivity[2] = m4;
	      conectivity[3] = m5;
	      cutmesh.addElem(5,conectivity);
	      // iso zero element is at the end
	      cutmesh.setIsozero(6);
	      // add iso zero element
	      face_conectivity[0] = i7;
	      face_conectivity[1] = i6;
	      face_conectivity[2] = i5;
	      cutmesh.addElem(6,face_conectivity);
	      face_conectivity[0] = i5;
	      face_conectivity[2] = i4;
	      cutmesh.addElem(7,face_conectivity);

#ifdef XREFMESH_WITH_SUB
	      // add edge "in"
	      edge_conectivity[0] = i6;
	      edge_conectivity[1] = i1;
	      cutmesh.addEdge(10+i_edge3,edge_conectivity);
	      edge_conectivity[0] = i4;
	      cutmesh.addEdge(20+i_edge1,edge_conectivity);
	      // add edge "out"
	      edge_conectivity[1] = i2;
	      cutmesh.addEdge(50+i_edge1,edge_conectivity);
	      edge_conectivity[0] = i5;
	      cutmesh.addEdge(60+i_edge2,edge_conectivity);
	      // add edge "in"
	      edge_conectivity[1] = i0;
	      cutmesh.addEdge(30+i_edge2,edge_conectivity);
	      edge_conectivity[0] = i7;
	      cutmesh.addEdge(40+i_edge5,edge_conectivity);
	      // add edge "out"
	      edge_conectivity[1] = i3;
	      cutmesh.addEdge(70+i_edge5,edge_conectivity);
	      edge_conectivity[0] = i6;
	      cutmesh.addEdge(80+i_edge3,edge_conectivity);

	      // starting from i7
	      // add face "in"
	      face_conectivity[0] = i7;
	      face_conectivity[1] = i5;
	      face_conectivity[2] = i0;
	      cutmesh.addFace(10+i_face3,face_conectivity);
	      // add face "out"
	      face_conectivity[1] = i6;
	      face_conectivity[2] = i3;
	      cutmesh.addFace(70+i_face1,face_conectivity);

	      // starting from i7 or i6 depending on cut choise
	      // add face "in"
	      face_conectivity[0] = n7;
	      face_conectivity[1] = i0;
	      face_conectivity[2] = i1;
	      cutmesh.addFace(20+i_face1,face_conectivity);

	      // starting from i6
	      // add face "in"
	      face_conectivity[0] = i6;
	      face_conectivity[1] = i7;
	      face_conectivity[2] = n1;
	      cutmesh.addFace(30+i_face1,face_conectivity);
	      // add face "out"
	      face_conectivity[1] = i4;
	      face_conectivity[2] = m1;
	      cutmesh.addFace(80+i_face2,face_conectivity);

	      // starting from i6 or i4 depending on cut choise
	      // add face "out"
	      face_conectivity[0] = m6;
	      face_conectivity[1] = i2;
	      face_conectivity[2] = i3;
	      cutmesh.addFace(90+i_face2,face_conectivity);

	      // starting from i4
	      // add face "in"
	      face_conectivity[0] = i4;
	      face_conectivity[1] = i6;
	      face_conectivity[2] = i1;
	      cutmesh.addFace(40+i_face2,face_conectivity);
	      // add face "out"
	      face_conectivity[1] = i5;
	      face_conectivity[2] = i2;
	      cutmesh.addFace(100+i_face0,face_conectivity);

	      // starting from i5 or i4 depending on cut choise
	      // add face "in"
	      face_conectivity[0] = n3;
	      face_conectivity[1] = i1;
	      face_conectivity[2] = i0;
	      cutmesh.addFace(50+i_face0,face_conectivity);

	      // starting from i5
	      // add face "in"
	      face_conectivity[0] = i5;
	      face_conectivity[1] = i4;
	      face_conectivity[2] = n1;
	      cutmesh.addFace(60+i_face0,face_conectivity);
	      // add face "out"
	      face_conectivity[1] = i7;
	      face_conectivity[2] = m1;
	      cutmesh.addFace(110+i_face3,face_conectivity);

	      // starting from i5 or i7 depending on cut choise
	      // add face "out"
	      face_conectivity[0] = m2;
	      face_conectivity[1] = i3;
	      face_conectivity[2] = i2;
	      cutmesh.addFace(120+i_face3,face_conectivity);

#endif
	    }
	  else
	    {
	      // add elements "in"
	      conectivity[0] = i7;
	      conectivity[1] = i6;
	      conectivity[2] = i5;
	      conectivity[3] = m1;
	      cutmesh.addElem(0,conectivity);
	      conectivity[0] = i5;
	      conectivity[2] = i4;
	      cutmesh.addElem(1,conectivity);
	      conectivity[0] = m2;
	      conectivity[1] = m3;
	      conectivity[2] = m4;
	      conectivity[3] = m5;
	      cutmesh.addElem(2,conectivity);
	      // separation betwen "in" an "out" element
	      cutmesh.setLimite(3);
	      // add element "out"
	      conectivity[0] = n1;
	      conectivity[1] = n2;
	      conectivity[2] = n3;
	      conectivity[3] = i6;
	      cutmesh.addElem(3,conectivity);
	      conectivity[0] = n4;
	      conectivity[1] = n5;
	      conectivity[2] = i5;
	      cutmesh.addElem(4,conectivity);
	      conectivity[0] = i0;
	      conectivity[1] = n6;
	      conectivity[3] = i7;
	      cutmesh.addElem(5,conectivity);
	      // iso zero element is at the end
	      cutmesh.setIsozero(6);
	      // add iso zero element
	      face_conectivity[0] = i6;
	      face_conectivity[1] = i7;
	      face_conectivity[2] = i5;
	      cutmesh.addElem(6,face_conectivity);
	      face_conectivity[1] = i5;
	      face_conectivity[2] = i4;
	      cutmesh.addElem(7,face_conectivity);

#ifdef XREFMESH_WITH_SUB
	      // add edge "out"
	      edge_conectivity[0] = i6;
	      edge_conectivity[1] = i1;
	      cutmesh.addEdge(60+i_edge3,edge_conectivity);
	      edge_conectivity[0] = i4;
	      cutmesh.addEdge(50+i_edge1,edge_conectivity);
	      // add edge "in"
	      edge_conectivity[1] = i2;
	      cutmesh.addEdge(10+i_edge1,edge_conectivity);
	      edge_conectivity[0] = i5;
	      cutmesh.addEdge(20+i_edge2,edge_conectivity);
	      // add edge "out"
	      edge_conectivity[1] = i0;
	      cutmesh.addEdge(70+i_edge2,edge_conectivity);
	      edge_conectivity[0] = i7;
	      cutmesh.addEdge(80+i_edge5,edge_conectivity);
	      // add edge "in"
	      edge_conectivity[1] = i3;
	      cutmesh.addEdge(30+i_edge5,edge_conectivity);
	      edge_conectivity[0] = i6;
	      cutmesh.addEdge(40+i_edge3,edge_conectivity);

	      // starting from i7
	      // add face "out"
	      face_conectivity[0] = i7;
	      face_conectivity[1] = i5;
	      face_conectivity[2] = i0;
	      cutmesh.addFace(80+i_face3,face_conectivity);
	      // add face "in"
	      face_conectivity[1] = i6;
	      face_conectivity[2] = i3;
	      cutmesh.addFace(10+i_face1,face_conectivity);

	      // starting from i7 or i6 depending on cut choise
	      // add face "out"
	      face_conectivity[0] = n7;
	      face_conectivity[1] = i0;
	      face_conectivity[2] = i1;
	      cutmesh.addFace(70+i_face1,face_conectivity);

	      // starting from i6
	      // add face "out"
	      face_conectivity[0] = i6;
	      face_conectivity[1] = i7;
	      face_conectivity[2] = n1;
	      cutmesh.addFace(90+i_face1,face_conectivity);
	      // add face "in"
	      face_conectivity[1] = i4;
	      face_conectivity[2] = m1;
	      cutmesh.addFace(20+i_face2,face_conectivity);

	      // starting from i6 or i4 depending on cut choise
	      // add face "in"
	      face_conectivity[0] = m6;
	      face_conectivity[1] = i2;
	      face_conectivity[2] = i3;
	      cutmesh.addFace(30+i_face2,face_conectivity);

	      // starting from i4
	      // add face "out"
	      face_conectivity[0] = i4;
	      face_conectivity[1] = i6;
	      face_conectivity[2] = i1;
	      cutmesh.addFace(100+i_face2,face_conectivity);
	      // add face "in"
	      face_conectivity[1] = i5;
	      face_conectivity[2] = i2;
	      cutmesh.addFace(40+i_face0,face_conectivity);

	      // starting from i5 or i4 depending on cut choise
	      // add face "out"
	      face_conectivity[0] = n3;
	      face_conectivity[1] = i1;
	      face_conectivity[2] = i0;
	      cutmesh.addFace(110+i_face0,face_conectivity);

	      // starting from i5
	      // add face "out"
	      face_conectivity[0] = i5;
	      face_conectivity[1] = i4;
	      face_conectivity[2] = n1;
	      cutmesh.addFace(120+i_face0,face_conectivity);
	      // add face "in"
	      face_conectivity[1] = i7;
	      face_conectivity[2] = m1;
	      cutmesh.addFace(50+i_face3,face_conectivity);

	      // starting from i5 or i7 depending on cut choise
	      // add face "in"
	      face_conectivity[0] = m2;
	      face_conectivity[1] = i3;
	      face_conectivity[2] = i2;
	      cutmesh.addFace(60+i_face3,face_conectivity);
#endif
	    }
	  break;
	}
        // no cut
      default :
	{
	  return( 0 );
	  break;
	}
      }

    //store the constitutive nodes
    // node 2
    xPoint UVW;
    UVW.u = zero;
    UVW.v = one;
    UVW.w = zero;
    cutmesh.addPoint(2,UVW);
    // node 1
    UVW.u = one;
    UVW.v = zero;
    cutmesh.addPoint(1,UVW);
    // node 0
    UVW.u = zero;
    cutmesh.addPoint(0,UVW);
    // node 3
    UVW.w = one;
    cutmesh.addPoint(3,UVW);

    // normal endding
    return ( 1 );

  }


} // end of namespace
