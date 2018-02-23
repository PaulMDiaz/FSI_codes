// ---- 2D Circular Cylinder Gmsh Tutorial ----
// 2D_cylinder_tutorial.geo
// Creates a mesh with an inner structured-quad region and
// an radius unstructured tri region
//
// Created 11/26/2014 by Jacob Crabill
// Aerospace Computing Lab, Stanford University
// --------------------------------------------

// Gmsh allows variables; these will be used to set desired
// element sizes at various Points.
// Adjust as necessary

cl1 = 0.025; //entry
cl2 = .01; // inner circles
cl3 = 0.05;  //exit
cl4 = 0.05;  //structure
nc1 = 20;  //number of nodes on 1 quarter of cylinder


// number of nodes on structure, vertical and horizontal
N_s_v = 4;
N_s_h = 40;

//cl1 = 0.0125; //entry
//cl2 = .005; // inner circles
//cl3 = 0.025;  //exit
//cl4 = 0.005;  //structure

//nc1 = 40;

radius = 0.05;

struct_l = 0.35;
struct_h = 0.02;

//numinner = 0;
//extr = -1;

// Exterior (bounding box) of mesh
Point(1) = {0, 0, 0, cl1};
Point(2) = { 2.2, 0, 0, cl3};
Point(3) = { 2.2,  0.41, 0, cl3};
Point(4) = {0,  0.41, 0, cl1};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

// Cylinder
Point(5) = {0.2,   0.2, 0, cl2};

Point(8) = {0.2 + (radius*radius-0.01*0.01)^(0.5), 0.2 + struct_h/2, 0, cl4};
Point(9) = {0.2 + (radius*radius-0.01*0.01)^(0.5), 0.2 - struct_h/2, 0, cl4};



Point(11) = {0.2, 0.2-radius, 0, cl2};
Point(12) = {0.2-radius,  0.2, 0, cl2};
Point(13) = {0.2, 0.2 + radius, 0, cl2};

Circle(7) = {11, 5, 9};
Circle(8) = {12, 5, 11};
Circle(9) = {13, 5, 12};
Circle(10) = {8, 5, 13};


Circle(11) = {8, 5, 9};

//Circle(15) = {13, 5, 9};
//Circle(16) = {12, 5, 8};

Transfinite Line {7,8,9,10} = nc1;

// Structure use points 8 and 9
//Point(24) = {0.2 + (radius*radius-0.01*0.01)^(0.5), 0.2 - struct_h/2, 0, cl4};
Point(25) = {0.2+radius+struct_l, 0.2-  struct_h/2, 0, cl4};
Point(26) = {0.2+radius+struct_l, 0.2+struct_h/2, 0, cl4};
//Point(27) = {0.2 + (radius*radius-0.01*0.01)^(0.5), 0.2+struct_h/2, 0, cl4};
Line(24) = {9, 25};
Line(25) = {25, 26};
Line(26) = {26, 8};
//Line(27) = {27, 24}; //have to revisit this, replaced by circle 8.

//Circle(17) = {27, 5, 8};

Transfinite Line {24,26} = N_s_h;
//Transfinite Line {25,27} = N_s_v;
Transfinite Line {25, 11} = N_s_v;

//Transfinite Line {5,6,7,8,13,14,15,16} = N; // We want 40 points along each of these lines
//Transfinite Line {9,10,11,12} = numinner Using Progression 1.1;    // And 10 points along each of these lines

//Using Progression 1.1

// Each region which to be independently meshed must have a line loop
// Regions which will be meshed with Transfinite Surface must have 4 lines
// and be labeled in CCW order, with the correct orientation of each edge
//Line Loop(1) = {1, 2, 3, 4, 7, 16, 8, 15}; // Exterior
//Line Loop(2) = {10, 8, -11, -5}; // RH side of quad region - note ordering
//Line Loop(3) = {7, -12, -6, 9}; // LH side of quad region - note ordering
//Line Loop(4) = {-10, -14, 12, 15}; // RH side of quad region - note ordering
//Line Loop(5) = {16, -9, -13, 11}; // LH side of quad region - note ordering

//I think these are the loops that define regions.
Line Loop(1) = {1, 2, 3, 4, 10, 9, 8, 7, 24, 25, 26}; // Exterior
Line Loop(2) = {24, 25, 26, 11}; // Structure boundary

//Line Loop(2) = {10, 8, -11, -5}; // RH side of quad region - note ordering
//Line Loop(3) = {7, -12, -6, 9}; // LH side of quad region - note ordering
//Line Loop(4) = {-10, -14, 12, 15}; // RH side of quad region - note ordering
//Line Loop(5) = {16, -9, -13, 11}; // LH side of quad region - note ordering

Plane Surface(1) = {1}; // radius unstructured region
Plane Surface(2) = {2}; // structure perhaps

Transfinite Surface{2};

//Plane Surface(2) = {2}; // RH inner structured region
//Plane Surface(3) = {3}; // LH inner structured region
//Plane Surface(4) = {4}; // RH inner structured region
//Plane Surface(5) = {5}; // LH inner structured region

// Mesh these surfaces in a structured manner
//Transfinite Surface{2,3,4,5};

// Turn into quads (optional, but Transfinite Surface looks best with quads)
//Recombine Surface {2,3,4,5};
// Turn radius region into unstructured quads (optional)
//Recombine Surface {1};

// Change layer to inc1rease z subdivision
//Extrude {0, 0, extr} { Surface{1,2,3,4,5}; Layers{1}; Recombine;}


// Apply boundary conditions
// Note: Can change names later at top of .msh file
// Each boundary in gmsh must be labeled differently
// rename the boundaries manually in the resulting .msh file
//Physical Line("Bottom") = {1};
//Physical Line("Right")  = {2};
//Physical Line("Top")    = {3};
//Physical Line("Left")   = {4};
//Physical Line("Circle") = {5,6};
// Alternate version - make all 4 radius bounds part of the same B.C.:
//Physical Line("Char") = {1,2,3,4};

// IMPORTANT: "FLUID" MUST contain all fluid surfaces(2D)/volumes(3D)


//Physical Surface("wall") = {115, 79, 97,141};
//Physical Surface("inflow") = {41, 37, 29};
//Physical Surface("outflow") = {33};
//Physical Surface("periodic_0_r") = {1,2,3,4,5};
//Physical Surface("periodic_0_l") = {58,124,80,102,146};
//Physical Volume("fluid") = {1, 2, 4, 3, 5};


//Physical Surface("fluid") = {1};
//+
//Physical Surface("structure") = {1};
//+
Physical Line("Walls") = {3, 1};
//+
Physical Line("inlet") = {4};
//+
Physical Line("outlet") = {2};
//+
Physical Line("cylinder") = {10, 9, 8, 7};
//+
Physical Line("end") = {11};
//+
Physical Line("FSI") = {26, 25, 24};
//+
Physical Surface("structure") = {2};
//+
Physical Surface("fluid") = {1};
