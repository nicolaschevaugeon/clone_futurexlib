dens=0.5;

Point(1) = {0, 0, 0, dens};
Point(2) = {0, 4.5, 0, dens};
Point(3) = {5, 4.2, 0, dens};
Point(4) = {5, 6.2, 0, dens};
Line (1) = {1, 3};
Line (2) = {3, 4};
Line (3) = {4, 2};
Line (4) = {2, 1};


Physical Point(101) = {1};
Physical Point(103) = {3};
Physical Point(104) = {4};
Physical Point(102) = {2};

Physical Line(111) = {1};
Physical Line(112) = {2};
Physical Line(113) = {3};
Physical Line(114) = {4};
Line Loop(6) = {2,3,4,1};
Plane Surface(7) = {6};

Physical Surface(121) = {7};
