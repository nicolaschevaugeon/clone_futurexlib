ta = 0.1;
r = 1;
Point(1) = {0, 0, 0, ta};
Point(2) = {-r, 0, 0, ta};
Point(3) = {0, r, 0, ta};
Point(4) = {0, -r, 0, ta};
Point(5) = {10, -r, 0, ta};
Point(6) = {10, 1, 0, ta};
Circle(1) = {4, 1, 2};
Circle(2) = {2, 1, 3};
Line(3) = {5, 4};
Line(4) = {3, 6};
Line(5) = {6, 5};
Line Loop(6) = {4, 5, 3, 1, 2};
Plane Surface(7) = {6};
Physical Line(8) = {4, 2, 1, 3};
Physical Line(12) = {5};
Physical Surface(9) = {7};
