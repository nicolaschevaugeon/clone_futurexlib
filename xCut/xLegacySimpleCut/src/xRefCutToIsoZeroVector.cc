/* 
   This file is a part of eXlibris C++ Library
   under the GNU Lesser General Public License.
   See the NOTICE.md & LICENSE.md files for terms 
   and conditions.
*/

#include <iostream>
#include <sstream>
#include <vector>


//xfem
#include "xRefCutToIsoZeroVector.h"
#include "xMesh.h"
#include "xLevelSet.h"

//AOMD
#include "mEntity.h"

using namespace xfem;

namespace xcut
{

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////// xRefCutToIsoZeroVector class implementation ///////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //
  /////////////////////////////////////// Static member //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // case C1 : giving id of the edge cut this table gives id of the 
  // 2 tetrahedron nodes in iso zero.
  unsigned char xRefCutToIsoZeroVector::c1_node[6][2]={
    {2,3},     // edge 0 cut
    {0,3},     // edge 1 cut
    {1,3},     // edge 2 cut
    {0,1},     // edge 3 cut
    {1,2},     // edge 4 cut
    {0,2}      // edge 5 cut
  };

  //  case 2 : giving id of the 2 edges cut this table gives the id
  //           of the tetrahedron node in iso zero
  //  there is 6 edges so from start it is a 6x6 matrix table.
  //  As  id_cut[0]<id_cut[1] by construction 
  //   => only upper termes are relevant
  //  As no comon index are possible = 2 same edge stored in id_cut is imposible by construction
  //   => first colonne and last line may be remove => 5X5 table where second index have to be shifted by 1
  //
  //   out of this (5)(5+1)/2=15 terms only 12 are relevant as somme conbinaition are impossible :
  //      opposite edge cut will never give a C2 case
  //   6 is a failure code for security
  unsigned char xRefCutToIsoZeroVector::c2_node[5][5]={
    // edge    1   2   3   4   5  cut
    {  3,  3,  6,  2,  2},     // edge 0 cut
    {  6,  3,  0,  6,  0},     // edge 1 cut
    {  6,  6,  1,  1,  6},     // edge 2 cut
    {  6,  6,  6,  1,  0},     // edge 3 cut
    {  6,  6,  6,  6,  2}      // edge 4 cut
  };
  // Node of the tetrahedron in reference coordinate
  double xRefCutToIsoZeroVector::tet_node[4][3]={
    {0.,0.,0.},    // node 0
    {1.,0.,0.},    // node 1
    {0.,1.,0.},    // node 2
    {0.,0.,1.}     // node 3
  };

  // size of 3 doubles (portable up to computer using 86 bytes for representation of double considering that at the same times this 
  // same computer keep a 1 byte representation for unsigned char ......) 
  unsigned char xRefCutToIsoZeroVector::size3d=(unsigned char)(3*sizeof(double));
  /////////////////////////////////////// End Static member //////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////// Constructor ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////// End constructor ////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////// Destructor /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////// End Destructor /////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////// Private methode ////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////// End Private methode ////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////// Public methode /////////////////////////////////////////////////////////////////////////////////////////////////////////////

  ///////////////
  /// 1D Edge
  ///////////////
  int xRefCutToIsoZeroVector::cutEdgeRefByLevelSet ( const std::vector < double > & ls_value, std::vector < double > & p ) const
  {
    //
    //     0------1    --> u
    //         0
    //
    // mapping u : [-1,1] <=> [x0,x1]
    //
    double lsv0 = ls_value[0];
    double lsv1 = ls_value[1];
    double prod = lsv0*lsv1;
    const double one = 1.0E+0;
    const double zero = 0.0E+0;
    assert(p.size()==3);

    // the iso-zero cut the edge ////////////////////////////////////////////////////////////////////////////
    if ( prod <= zero )
      {
        p[0]=-2.0*lsv0/( lsv1 - lsv0 )-one;
        p[1]=zero;
        p[2]=zero;

        if ( prod < zero )
	  return 1;
        else if ( lsv0 == zero )
	  return -1;
        else
	  return -2;

      }

    // the iso-zero dont cut the edge ////////////////////////////////////////////////////////////////////////////
    return 0;

  }

  //
  //
  //
  //
  //
  ///////////////
  /// 2D triangle
  ///////////////
  int xRefCutToIsoZeroVector::cutTriRefByLevelSet ( const std::vector < double > & ls_value, std::vector < double > & p ) const
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
     */
    unsigned char nb_cut = 0,dp=0;
    unsigned char id_opposite_node[2];
    const double lsv0 = ls_value[0];
    const double lsv1 = ls_value[1];
    const double lsv2 = ls_value[2];
    const double one = 1.0E+0;
    const double zero = 0.0E+0;
    assert(p.size()==6);

    // the iso-zero cut the edge 0 ////////////////////////////////////////////////////////////////////////////
    if ( ( lsv0*lsv1 ) < zero )
      {
        p[0]=-lsv0/( lsv1 - lsv0 );
        p[1]=zero;
        p[2]=zero;

        id_opposite_node[0]=2;
        ++nb_cut;
        dp+=3;
      }

    // the iso-zero cut the edge 1 ////////////////////////////////////////////////////////////////////////////
    if ( ( lsv1*lsv2 ) < 0 )
      {
        const double v = -lsv1/( lsv2 - lsv1 );
        p[dp++]=one-v;
        p[dp++]=v;
        p[dp++]=zero;

        id_opposite_node[nb_cut]=0;
        ++nb_cut;
      }

    // the iso-zero cut the edge 2 ////////////////////////////////////////////////////////////////////////////
    if ( ( lsv0*lsv2 ) < 0 )
      {
        p[dp++]=zero;
        p[dp++]=-lsv0/( lsv2 - lsv0 );
        p[dp++]=zero;

        id_opposite_node[nb_cut]=1;
        ++nb_cut;
      }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////
    // if no cut check if iso zero is passing on one node/edge ////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////
    if (!nb_cut)
      {
        if (lsv0 == zero)
	  {
            p[0]=zero;
            p[1]=zero;
            p[2]=zero;

            if (lsv1 == zero)
	      {
                p[3]=one;
                p[4]=zero;
                p[5]=zero;
                return -1;
	      }
            else if (lsv2 == zero)
	      {
                p[3]=zero;
                p[4]=one;
                p[5]=zero;
                return -3;
	      }
            else
	      {
                // only node 0 is in isozero
                return -4;
	      }
	  }
        // node 0 is not in isozero
        else if (lsv1 == zero )
	  {
            p[0]=one;
            p[1]=zero;
            p[2]=zero;
            if (lsv2 == zero)
	      {
                p[3]=zero;
                p[4]=one;
                p[5]=zero;
                return -2;
	      }
            else
	      {
                // only node 1 is in isozero
                return -5;
	      }
	  }
        // node 0 and 1 are not in isozero
        else if (lsv2 == zero )
	  {
            // only node 2 is in isozero
            p[0]=zero;
            p[1]=one;
            p[2]=zero;
            return -6;
	  }
        //  node 0,1 and 2 are not in isozero and no cut
        else
	  return 0;

      }
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////
    // if only one cut adding  node passing true iso zero  ////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////
    else if (nb_cut<2)
      {
        p[3]=zero;
        p[4]=zero;
        p[5]=zero;

        const int iso_node=id_opposite_node[0];
        if ( iso_node > 1)
	  p[4]=one;
        else if ( iso_node > 0)
	  p[3]=one;
      }

    // normal endding
    return 1;

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
   *       4|   \ 5                            face 1
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
   *          5/     \3                        face 2
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
   *       4|   \ 3                            face 3
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
   *                        edge 3 is defined by vertice 2 and 3
   *                        edge 4 is defined by vertice 0 and 3
   *                        edge 5 is defined by vertice 1 and 3
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
   ********************************************************************************************************************/
  int xRefCutToIsoZeroVector::cutTetRefByLevelSet ( const std::vector < double > & ls_value, std::vector < double > & p ) const
  {
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////
    //
    // local
    unsigned char nb_cut = 0,dp=0;
    unsigned char id_cut[4];
    double lsv0 = ls_value[0];
    double lsv1 = ls_value[1];
    double lsv2 = ls_value[2];
    double lsv3 = ls_value[3];
    const double one = 1.0E+0;
    const double zero = 0.0E+0;
    double r;
    assert(p.size()==12);


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Determination of the cutting points ////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////

    // the iso-zero cut the edge 0 ////////////////////////////////////////////////////////////////////////////
    if ( ( lsv0*lsv1 ) < zero )
      {

        p[0]=-lsv0/( lsv1 - lsv0 );
        p[1]=zero;
        p[2]=zero;

        dp+=3;
        id_cut[0]=0;
        ++nb_cut;
      }

    // the iso-zero cut the edge 1 ////////////////////////////////////////////////////////////////////////////
    if ( ( lsv1*lsv2 ) < zero )
      {
        // lineare interpolation
        r = -lsv1/( lsv2 - lsv1 );
        p[dp++]=one-r;
        p[dp++]=r;
        p[dp++]=zero;

        id_cut[nb_cut]=1;
        ++nb_cut;
      }

    // the iso-zero cut the edge 2 ////////////////////////////////////////////////////////////////////////////
    if ( ( lsv0*lsv2 ) < zero )
      {
        p[dp++]=zero;
        p[dp++]=-lsv0/( lsv2 - lsv0 );
        p[dp++]=zero;

        id_cut[nb_cut]=2;
        ++nb_cut;
      }

    // the iso-zero cut the edge 3 ////////////////////////////////////////////////////////////////////////////
    if ( ( lsv2*lsv3 ) < zero )
      {
        r = -lsv2/( lsv3 - lsv2 );
        p[dp++]=zero;
        p[dp++]=one-r;
        p[dp++]=r;

        id_cut[nb_cut]=3;
        ++nb_cut;
      }

    // the iso-zero cut the edge 4 ////////////////////////////////////////////////////////////////////////////
    if ( ( lsv0*lsv3 ) < 0 )
      {
        // lineare interpolation
        p[dp++]=zero;
        p[dp++]=zero;
        p[dp++]=-lsv0/( lsv3 - lsv0 );

        id_cut[nb_cut]=4;
        ++nb_cut;
      }

    // the iso-zero cut the edge 5 ////////////////////////////////////////////////////////////////////////////
    if ( ( lsv1*lsv3 ) < zero )
      {
        r = -lsv1/( lsv3 - lsv1 );
        p[dp++]=one-r;
        p[dp++]=zero;
        p[dp++]=r;

        id_cut[nb_cut]=5;
        ++nb_cut;
      }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////
    // setting return value for CASE4 and /////////////////////////////////////////////////////////////////////
    // if no cut check if iso zero is passing on one face /////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////
    switch (nb_cut)
      {
        //  1 edge only is cut, the iso-zero surface passe thru the 2 opposites coins of the tetrahedron
        //  iso-zero => 1 triangle
      case 1 :
	{
	  const unsigned char i=id_cut[0];
	  memcpy(&p[3],&tet_node[c1_node[i][0]][0],size3d);
	  memcpy(&p[6],&tet_node[c1_node[i][1]][0],size3d);
	  return 1;
	  break;
	}
        //  2 edges are cut, the iso-zero surface passe thru the opposite coin of the tetrahedron
        //  iso-zero => 1 triangle
      case 2 :
	{
	  memcpy(&p[6],&tet_node[c2_node[id_cut[0]][id_cut[1]-1]][0],size3d);
	  return 1;
	  break;
	}
        //  3 edges are cut, the iso-zero surface passe thru them
        //  iso-zero => 1 triangle
      case 3 :
	{
	  return 1;
	  break;
	}
        //  4 edges are cut, the iso-zero surface passe thru them
        //  iso-zero => 2 triangle
        //  3 cases possible :
        //   edge 1,2,4,5 cut
        //   edge 0,1,3,4 cut
        //   edge 0,2,3,5 cut
        //     ? ? ? look like with this edge ordering, whatever the cut is, we ave the same sequence
        //     of cutting node gently turning around the polygone !!!
        //     so return point may in all case be taken as 0 1 2 and 0 2 3 to form 2 correct triangle ... 
        //     TODO : check and if problem use intermediate point container and order testing to obtaine the espected 
        //     0 1 2 and 0 2 3 sequence.
      case 4 :
	{
	  return 2;
	  break;
	}
        // no cut
      default :
	{
	  if (lsv0 == zero)
	    {
	      p[0]=zero;
	      p[1]=zero;
	      p[2]=zero;
	      if (lsv1 == zero)
		{
		  p[3]=one;
		  p[4]=zero;
		  p[5]=zero;
		  if (lsv2 == zero)
		    {
		      p[6]=zero;
		      p[7]=one;
		      p[8]=zero;
		      return -1;
		    }
		  else if (lsv3 == zero)
		    {
		      p[6]=zero;
		      p[7]=zero;
		      p[8]=one;
		      return -2;
		    }
		  else // only edge 0 is in isozero
		    return -5;
		}
	      else if (lsv2 == zero)
		{
		  p[3]=zero;
		  p[4]=one;
		  p[5]=zero;
		  if (lsv3 == zero)
		    {
		      p[6]=zero;
		      p[7]=zero;
		      p[8]=one;
		      return -4;
		    }
		  else // only edge 2 is in isozero
		    return -7;
		}
	      else if (lsv3 == zero)
		{
		  // only edge 4 is in isozero
		  p[3]=zero;
		  p[4]=zero;
		  p[5]=one;
		  return -9;
		}
	      else
		{
		  // only node 0 is in isozero
		  return -11;
		}
	    }
	  // node 0 is not in isozero
	  else if (lsv1 == zero)
	    {
	      p[0]=one;
	      p[1]=zero;
	      p[2]=zero;
	      if (lsv2 == zero)
		{
		  p[3]=zero;
		  p[4]=one;
		  p[5]=zero;
		  if (lsv3 == zero)
		    {
		      p[6]=zero;
		      p[7]=zero;
		      p[8]=one;
		      return -3;
		    }
		  else
		    {
		      // only edge 1 is in isozero
		      return -6;
		    }
		}
	      else if (lsv3 == zero)
		{
		  // only edge 5 is in isozero
		  p[3]=zero;
		  p[4]=zero;
		  p[5]=one;
		  return -10;
		}
	      else
		{
		  // only node 1 is in isozero
		  return -12;
		}
	    }
	  //  node 0 and 1 are not in isozero
	  else if (lsv2 == zero)
	    {
	      p[0]=zero;
	      p[1]=one;
	      p[2]=zero;
	      if (lsv3 == zero)
		{
		  // only edge 3 is in isozero
		  p[3]=zero;
		  p[4]=zero;
		  p[5]=one;
		  return -8;
		}
	      else
		{
		  // only node 2 is in isozero
		  return -13;
		}
	    }
	  //  node 0,1 and 2 are not in isozero
	  else if (lsv3 == zero)
	    {
	      p[0]=zero;
	      p[1]=zero;
	      p[2]=one;
	      // only node 3 is in isozero
	      return -14;
	    }
	  //  node 0,1,2 and 3 are not in isozero and no cut
	  else
	    return 0;
	  break;
	}
      }

  }
  /////////////////////////////////////// End Public methode /////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

} // end of namespace
