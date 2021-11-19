nbele = 4;
x0 = -1.;
y0 = -1.;
z0 = -1.;
L1 = 2.;
L2 = 2.;
L3 = 2.;
h = 2/nbele;

Point(1) = {x0,y0,z0,h};
Point(2) = {x0+L1,y0,z0,h};
Point(3) = {x0+L1,y0+L2,z0,h};
Point(4) = {x0,y0+L2,z0,h};
Point(5) = {x0,y0,z0+L3,h};
Point(6) = {x0+L1,y0,z0+L3,h};
Point(7) = {x0+L1,y0+L2,z0+L3,h};
Point(8) = {x0,y0+L2,z0+L3,h};
Point(9) = {x0-L1,y0,z0,h};
Point(10) = {x0-L1,y0+L2,z0,h};
	 


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

Line(13) = {1, 9};
Line(14) = {9, 10};
Line(15) = {10, 4};



Line Loop(16) = {1,2,3,4};
Plane Surface(201) = {16};

Line Loop(17) = {-9,-5,1,6};
Plane Surface(202) = {17};

Line Loop(18) = {-10,-6,2,7};
Plane Surface(203) = {18};

Line Loop(19) = {-11,-7,3,8};
Plane Surface(204) = {19};

Line Loop(20) = {12,-5,-4,8};
Plane Surface(205) = {20};

Line Loop(21) = {11,12,9,10};
Plane Surface(206) = {21};

Line Loop(22) = {13,14,15,4};
Plane Surface(207) = {22};

Surface Loop(31) = {201,202,203,204,205,206};
Volume(301) = {31};

Physical Point   (101)  = {1} ;
Physical Point   (102)  = {2} ;
Physical Point   (103)  = {3} ;
Physical Point   (104)  = {4} ;
Physical Point   (105)  = {5} ;
Physical Point   (106)  = {6} ;
Physical Point   (107)  = {7} ;
Physical Point   (108)  = {8} ;
Physical Point   (109)  = {9} ;
Physical Point   (110)  = {10};

Physical Line   (201)  = {1} ;
Physical Line   (202)  = {2} ;
Physical Line   (203)  = {3} ;
Physical Line   (204)  = {4} ;
Physical Line   (205)  = {5} ;
Physical Line   (206)  = {6} ;
Physical Line   (207)  = {7} ;
Physical Line   (208)  = {8} ;
Physical Line   (209)  = {9} ;
Physical Line   (210)  = {10} ;
Physical Line   (211)  = {11} ;
Physical Line   (212)  = {12} ;	 	 
Physical Line   (213)  = {13} ;	 	 

Physical Line   (214)  = {14} ;	 	 
Physical Line   (215)  = {15} ;	 	 

Physical Surface   (414) = {201};
Physical Surface   (416) = {202};
Physical Surface   (418) = {203};
Physical Surface   (420) = {204};
Physical Surface   (422) = {205};
Physical Surface   (424) = {206};

Physical Surface   (426) = {207};

Physical Volume    (801) = {301};









