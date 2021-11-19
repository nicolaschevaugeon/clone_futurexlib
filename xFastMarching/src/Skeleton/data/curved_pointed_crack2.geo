ta = 0.1;
r = 1;
Point(1) = {0, 0, 0, ta};
Point(2) = {-r, 0, 0, ta};
Point(3) = {0, r, 0, ta};
Point(4) = {0, -r, 0, ta};
Point(5) = {10, -r, 0, ta};
Point(6) = {10, 1, 0, ta};
Point(7) = {4.5, 1.2, 0, 1.0};
Point(8) = {4.5, -0.8, 0, 1.0};
Point(9) = {-1.2, 0.7, 0, 1.0};
Point(10) = {-1.8, -0.2, 0, 1.0};
Point(11) = {-1.1, -0.9, 0, 1.0};
BSpline(1) = {3, 9, 10, 11, 4};
Spline(3) = {4, 8, 5};
Spline(4) = {3, 7, 6};
Line(5) = {6, 5};
Line Loop(6) = {4, 5, -3, -1};
Plane Surface(7) = {6};
Physical Line(8) = {4, 2, 1, 3};
Physical Line(12) = {5};
Physical Surface(9) = {7};




