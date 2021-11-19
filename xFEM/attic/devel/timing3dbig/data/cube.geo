nbele = 10;
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

Line Loop(13) = {1,2,3,4};
Plane Surface(20) = {13};
Line Loop(15) = {-9,-5,1,6};
Plane Surface(16) = {15};
Line Loop(17) = {-10,-6,2,7};
Plane Surface(18) = {17};
Line Loop(19) = {-11,-7,3,8};
Plane Surface(22) = {19};
Line Loop(21) = {12,-5,-4,8};
Plane Surface(14) = {21};
Line Loop(23) = {11,12,9,10};
Plane Surface(24) = {23};

Surface Loop(25) = {14,16,18,20,22,24};
Complex Volume(26) = {25};

Physical Point   (101)  = {1} ;
Physical Point   (102)  = {2} ;
Physical Point   (103)  = {3} ;

Physical Surface   (214) = {14};
Physical Surface   (109) = {16};
Physical Surface   (218) = {18};
Physical Surface   (220) = {20};
Physical Surface   (111) = {22};
Physical Surface   (224) = {24};


Physical Volume    (121) = {26};
/*Physical Line      (400) = {1,2,3};*/










