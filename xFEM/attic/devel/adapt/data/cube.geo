x0 = -1;
y0 = -1;
z0 = -1;

L1 = 2;
L2 = 2;
L3 = 2;

h = 2; 

nelt    = 32;
nbelt_x = nelt;
nbelt_y = nelt;
nbelt_z = nelt;


nbpt_x = nbelt_x + 1;
nbpt_y = nbelt_y + 1;
nbpt_z = nbelt_z + 1;

Point(1) = {x0,y0,z0,h};
Point(2) = {x0+L1,y0,z0,h};
Point(3) = {x0+L1,y0+L2,z0,h};
Point(4) = {x0,y0+L2,z0,h};
Point(5) = {x0,y0,z0+L3,h};
Point(6) = {x0+L1,y0,z0+L3,h};
Point(7) = {x0+L1,y0+L2,z0+L3,h};
Point(8) = {x0,y0+L2,z0+L3,h};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};
Line(5) = {1,5};
Line(6) = {2,6};
Line(7) = {3,7};
Line(8) = {4,8};
Line(9) = {5,6};
Line(10) = {6,7};
Line(11) = {7,8};
Line(12) = {8,5};

Line Loop(13) = {1,2,3,4};    Plane Surface(19) = {13};
Line Loop(14) = {-9,-5,1,6};  Plane Surface(20) = {14};
Line Loop(15) = {-10,-6,2,7}; Plane Surface(21) = {15};
Line Loop(16) = {11,7,-3,-8}; Plane Surface(22) = {16};
Line Loop(17) = {12,-5,-4,8}; Plane Surface(23) = {17};
Line Loop(18) = {11,12,9,10}; Plane Surface(24) = {18};

Surface Loop(25) = {19,-22,21,24,-20,-23}; Complex Volume(26) = {25};

/* Transfinite Line{5,1,6,9} = nbpt Using Power 1.0; */
/* Transfinite Line{8,3,7,11} = nbpt Using Power 1.0; */
/* Transfinite Line{4,2,10,12} = nbpt Using Power 1.0; */

Transfinite Line{1,9,3,11}  = nbpt_x Using Power 1.0; 
Transfinite Line{10,2,4,12} = nbpt_y Using Power 1.0; 
Transfinite Line{5,8,7,6}   = nbpt_z Using Power 1.0; 

Transfinite Surface{20} = {1,2,6,5};
Transfinite Surface{19} = {1,2,3,4};
Transfinite Surface{22} = {4,3,7,8};
Transfinite Surface{24} = {7,6,5,8};
Transfinite Surface{23} = {4,8,5,1};
Transfinite Surface{21} = {3,7,6,2};

/*Recombine Surface {20,19,22,24,23,21};*/

Transfinite Volume{26} = {5,6,2,1,8,7,3,4}; 


Physical Point   (101)  = {1,2,3,4,5,6,7,8} ;

/* Physical Point   (108)  = {1,2,3,4,5,6,7,8} ; */

/* Physical Line(100) = {1,2,3}; */

Physical Surface(200) = {19,20,21,22,23,24}; 

/* Physical Surface(219) = {19}; */
/* Physical Surface(224) = {24}; */

Physical Volume(300) = {26}; 









