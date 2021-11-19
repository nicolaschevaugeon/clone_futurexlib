lc = 0.1;
Point(1) = {0., 0., 0, lc};
Point(2) = {1., 0., 0, lc};
Point(3) = {1., 1., 0, lc};
Point(4) = {0., 1., 0, lc};	

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line Loop(5) = {3, 4, 1, 2};
Plane Surface(6) = {5};
Physical Point(101) = {1};
Physical Point(102) = {2};
Physical Point(103) = {3};
Physical Point(104) = {4};

Physical Line(201) = {1};
Physical Line(202) = {2};
Physical Line(203) = {3};
Physical Line(204) = {4};

Physical Surface(301) = {6};
