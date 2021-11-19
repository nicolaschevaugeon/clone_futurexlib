nbele = 50;
nbpt = nbele+1;
ta = 2/nbele;
r = 0.1;
/* Point      1 */
Point(newp) = {-1.0,-1.0,0.0,ta};
/* Point      2 */
Point(newp) = {1.0,-1.0,0.0,ta};
/* Point      3 */
Point(newp) = {1.0,1.0,0.0,ta};
/* Point      4 */
Point(newp) = {-1.0,1.0,0.0,ta};

Line(9)  = {1,2};
Line(10) = {2,3};
Line(11) = {4,3};
Line(12) = {1,4};
  


Physical Line    (109)  = {9};
Physical Line    (110)  = {10};
Physical Line    (111)  = {11};
Physical Line    (112)  = {12};












Point(5) = {0, 0, 0, 1.0};
Point(6) = {0, r, 0, ta};
Point(7) = {r, 0, 0, ta};
Point(8) = {0, -r, 0, ta};
Point(9) = {-r, 0, 0, ta};
Circle(113) = {6, 5, 7};
Circle(114) = {7, 5, 8};
Circle(115) = {8, 5, 9};
Circle(116) = {9, 5, 6};
Line Loop(117) = {12, 11, -10, -9};
Line Loop(118) = {116, 113, 114, 115};
Plane Surface(119) = {117, 118};
Physical Line    (200)  = {113,114,115,116};
Physical Surface(400) = 119;