/*
    octree is a subproject of  xfem : C++ Finite Element Library
    developed under the GNU Lesser General Public License
    See the NOTICE, CONTRIBUTORS & LICENSE files for conditions.
*/

#include "oTopo.h"

namespace xoctree{


    const int oTopo::base2_ijk[8][3] = { {0,0,0}, 
					 {1,0,0}, 
					 {0,1,0}, 
					 {1,1,0}, 
					 {0,0,1}, 
					 {1,0,1}, 
					 {0,1,1}, 
					 {1,1,1} };


    const int oTopo::nb_entities[4][3] = { {1,0,0},   //dim 0
                                           {2,1,0},   //dim 1
					   {4,4,1},   //dim 2
					   {8,12,6} }; //dim 3
    

    const int oTopo::ijk_node_shift[3][12][3] = { { {0,0,0},   //nodes
						    {1,0,0}, 
						    {1,1,0}, 
						    {0,1,0}, 
						    {0,0,1}, 
						    {1,0,1}, 
						    {1,1,1}, 
						    {0,1,1} },  //end nodes
						  { {1,0,0},    //edges
						    {2,1,0},    
						    {1,2,0},    
						    {0,1,0},    
						    {1,0,2},    
						    {2,1,2},    
						    {1,2,2},    
						    {0,1,2},    
						    {0,0,1},    
						    {2,0,1},    
						    {2,2,1},    
						    {0,2,1} }, //end edges
						  { {1,0,1},   // faces   
						    {2,1,1},    
						    {1,2,1},    
						    {0,1,1}, 
						    {1,1,0},   
						    {1,1,2} } }; //end faces


    const int oTopo::ijk_center_shift[4][3] = {
                                                  {},          // dim 0   
                                                  {1, 0, 0},   // dim 1   
                                                  {1, 1, 0},   // dim 2   
                                                  {1, 1, 1} }; // dim 3   


    const int oTopo::edge_connect[12][2] = { {0,1}, 
					     {1,2}, 
					     {2,3}, 
					     {3,0}, 
					     {4,5}, 
					     {5,6}, 
					     {6,7}, 
					     {7,4}, 
					     {0,4}, 
					     {1,5}, 
					     {2,6}, 
					     {3,7} };






    const int oTopo::face_connect[6][4] = { {0,1,5,4}, 
					    {1,2,6,5}, 
					    {2,3,7,6}, 
					    {3,0,4,7}, 
					    {1,0,3,2},
					    {4,5,6,7} };


    // Pay attention to the fact that description of the face is made consistent in between face_connect and face_edge_connect.
    // Around a face, edge and node are traversed in the same way and start by the same topological location :
    //               for face 0 node description is 0, 1, 5, 4 and
    //               in between 0 and 1 there is edge 0, then in between 1 and 5 there is edge 9,
    //               then in between 5 and 4 there is edge 4, then in between 4 and 0 there is edge 8.
    const int oTopo::face_edge_connect[6][4] = { {0, 9, 4, 8},
						 {1,10, 5, 9}, 
						 {2,11, 6,10}, 
						 {3, 8, 7,11}, 
						 {0, 3, 2, 1},
						 {4, 5, 6, 7} };
    
  
      const int oTopo::ijk_edge_neighbor_shift_2D[4][3] = { 
	                               { 0, -1, 0}, //Edge = 0
							       { 1,  0, 0}, //Edge = 1
							       { 0,  1, 0}, //Edge = 2
							       {-1,  0, 0}  //Edge = 3
							        };


      // for Dim = 3
      const int oTopo::ijk_edge_neighbor_shift_3D[12][3][3] =    { 
      //  id neighbor       0             1             2
	                { { 0, -1,  0}, { 0, -1, -1}, { 0,  0, -1}, },    //Edge = 0
	                { { 0,  0, -1}, { 1,  0, -1}, { 1,  0,  0}, },    //Edge = 1
	                { { 0,  0, -1}, { 0,  1, -1}, { 0,  1,  0}, },    //Edge = 2
	                { {-1,  0,  0}, {-1,  0, -1}, { 0,  0, -1}, },    //Edge = 3
	                { { 0, -1,  0}, { 0, -1,  1}, { 0,  0,  1}, },    //Edge = 4
	                { { 1,  0,  0}, { 1,  0,  1}, { 0,  0,  1}, },    //Edge = 5
	                { { 0,  1,  0}, { 0,  1,  1}, { 0,  0,  1}, },    //Edge = 6
	                { {-1,  0,  0}, {-1,  0,  1}, { 0,  0,  1}, },    //Edge = 7
	                { {-1, -1,  0}, {-1,  0,  0}, { 0, -1,  0}, },    //Edge = 8
	                { { 0, -1,  0}, { 1, -1,  0}, { 1,  0,  0}, },    //Edge = 9
	                { { 1,  0,  0}, { 1,  1,  0}, { 0,  1,  0}, },    //Edge = 10
	                { {-1,  0,  0}, {-1,  1,  0}, { 0,  1,  0}, }     //Edge = 11

                                                                    };


      const int oTopo::ijk_face_neighbor_shift[6][3] =  { 
                             { 0, -1,  0},  //Face 0 
							 { 1,  0,  0},  //Face 1
							 { 0,  1,  0},  //Face 2
							 {-1,  0,  0},  //Face 3
							 { 0,  0, -1},  //Face 4
							 { 0,  0,  1}   //Face 5
                                                        };

    const int oTopo::nb_nodes_connected_to_node[4] = { 0, 1, 2, 3};
    const int oTopo::connected_nodes[4][8][3] = {
	                                                     {  },              //Dim = 0            
	                                                     { {1}, 
							       {0} },           //Dim = 1
	                                                     { {1,3}, 
							       {0,2}, 
							       {1,3}, 
							       {0,2} },           //Dim = 2
	                                                     { {1,3,4}, 
							       {0,2,5}, 
							       {1,3,6}, 
							       {0,2,7}, 
							       {0,5,7}, 
							       {1,4,6}, 
							       {2,5,7}, 
							       {3,4,6}  } };           //Dim = 3
      

    const int oTopo::stencil_location[4][8][3] = {
	                                                     {  },              //Dim = 0                   
	                                                     { {1}, 
							       {0} },           //Dim = 1
	                                                     { {1,3}, 
							       {0,3}, 
							       {2,0}, 
							       {2,1} },           //Dim = 2
	                                                     { {1,3,5}, 
							       {0,3,5}, 
							       {1,3,6}, 
							       {2,0,5}, 
							       {4,1,3}, 
							       {4,0,3}, 
							       {4,2,0}, 
							       {4,2,1}  } };           //Dim = 3


//for (int l=0; l <= topo::MAX_DEPTH ; ++l) {  index_max_for_level[l] = powint(2,l)-1;}
    //modified as following : 

    
  //first line is power of 2, second line power of 4, third line power of 8

    const int oTopo::generation_size[3][MAX_DEPTH+1] = { {1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768/*, 65536, 131072*/},
  //Dim = 1
  {1, 4, 16, 64, 256, 1024, 4096, 16384, 65536, 262144, 1048576, 4194304, 16777216, 67108864, 268435456, 1073741824 /*,4294967296, 17179869184*/},
  //Dim = 2
  {1, 8, 64, 512, 4096, 32768, 262144, 2097152, 16777216, 134217728, 1073741824, -1, -1, -1, -1, -1/*, -1, -1*/} };
// ATTENTION : topo::MAX_DEPTH in 3D = 10 !!!!

  
//   const int oTopo::pow_base2[topo::MAX_DEPTH+1] ={ 1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768} ;

  const int oTopo::pow_base2[MAX_ABSOLUTE_DEPTH+1] ={ 1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536, 131072, 262144, 524288, 1048576, 2097152, 4194304, 8388608, 16777216, 33554432, 67108864, 134217728/*, 268435456, 536870912*/} ;

  const int oTopo::index_max_for_level[MAX_ABSOLUTE_DEPTH+1] = { 0, 1, 3, 7, 15, 31, 63, 127, 255, 511, 1023, 2047, 4095, 8191, 16383, 32767, 65535, 131071, 262143, 524287, 1048575, 2097151, 4194303, 8388607, 16777215, 33554431, 67108863, 134217727/*, 268435455, 536870911 */};

}
