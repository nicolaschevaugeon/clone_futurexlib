lc = 0.1;
Point(1) = {-1, -1, 0, lc};
Point(2) = {1, -1, 0, lc};
Point(3) = {1, -0, 0, lc};
Point(4) = {0, 0, 0, lc};
Point(5) = {0, 1, 0, lc};
Point(6) = {-1, 1, 0, lc};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 1};
Point(7) = {0.1, -0.3, 0, 1.0};
Point(8) = {-0.2, -0.5, 0, 1.0};
Point(9) = {-0.6, -0.4, 0, 1.0};
Point(10) = {-0.6, -0.7, 0, 1.0};
Point(11) = {0, -0.8, 0, 1.0};
Point(12) = {0.4, -0.7, 0, 1.0};
Point(13) = {0.7, -0.3, 0, 1.0};
Point(14) = {0.4, -0.2, 0, 1.0};
Point(15) = {-0.3, -0.7, 0, 1.0};

Line Loop(7) = {6, 1, 2, 3, 4, 5};
Spline(11) = {15, 8, 7, 14, 13, 12, 11, 15};
Line Loop(12) = {11};
Plane Surface(8) = {7, 12};
Physical Line(9) = {1, 2, 3, 4, 5, 6};
Physical Surface(10) = {8};


