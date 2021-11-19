nbele = 19;
nbpt = nbele+1;
ta = 2/nbele;
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
  
Line Loop(18) = {-12,9,10,-11};
Ruled Surface(21) = {18};
/*Plane Surface(21) = {18};*/ 

/*
Transfinite Line {9,11,10,12} = nbpt Using Power 1.0; 
Transfinite Surface {21} = {1,2,3,4}; 
*/

/*Recombine Surface {21} = 90;*/

Physical Point   (101)  = {1} ;
Physical Point   (102)  = {2} ;

Physical Line    (109)  = {9};
Physical Line    (110)  = {10};
Physical Line    (111)  = {11};
Physical Line    (112)  = {12};

Physical Surface (121)  = {21} ;











