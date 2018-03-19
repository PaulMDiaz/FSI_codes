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

//fine
//cl1 = 0.0125; //entry
//cl2 = .0025; // circle

//cl3 = 0.025; //exit
//cl4 = 0.0025; //structure
//cl5 = 0.00125; //structure fluid end (larger velocity)

//medium
cl1 = 0.0125; //entry
cl2 = .005; // circle

cl3 = 0.025; //exit
cl4 = 0.005; //structure
cl5 = 0.0025; //structure fluid end (larger velocity)


//course
//cl1 = 0.05; //entry
//cl2 = .01; // circle

//cl3 = 0.1; //exit
//cl4 = 0.01; //structure

// can reintroduce this node counting if I want.
//nc1 = 20;  //number of nodes on 1 quarter of cylinder

// number of nodes on structure, vertical and horizontal
//N_s_v = 4;
//N_s_h = 40;

radius = 0.05;    //radius of circle

struct_l = 0.35;  //length of structure
struct_h = 0.02;  //height of structure

// Exterior (bounding box) of mesh
Point(1) = {0, 0, 0, cl1};
Point(2) = { 2.5, 0, 0, cl3};
Point(3) = { 2.5,  0.41, 0, cl3};
Point(4) = {0,  0.41, 0, cl1};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

// Cylinder
Point(5) = {0.2,   0.2, 0, cl2};

// points on edge of cylinder and end of structure
Point(6) = {0.2 + (radius*radius-0.01*0.01)^(0.5), 0.2 + struct_h/2, 0, cl4}; //bar end top
Point(7) = {0.2 + (radius*radius-0.01*0.01)^(0.5), 0.2 - struct_h/2, 0, cl4}; //bar end lower

// points spaced around cylinder
Point(8) = {0.2, 0.2-radius, 0, cl2};    //bottom
Point(9) = {0.2-radius,  0.2, 0, cl2};   //left
Point(10) = {0.2, 0.2 + radius, 0, cl2};  //top

// Define curves
Circle(5) = {8, 5, 7};   // lower right
Circle(6) = {9, 5, 8};  // lower left
Circle(7) = {10, 5, 9};  // top left
Circle(8) = {6, 5, 10};  // top right

// Structure use points 8 and 9
Point(11) = {0.2+radius+struct_l, 0.2-  struct_h/2, 0, cl5}; //lower right
Point(12) = {0.2+radius+struct_l, 0.2+struct_h/2, 0, cl5};  //lower left

Line(9) = {7, 11};   //top
Line(10) = {11, 12};  //end
Line(11) = {12, 6};   //lower

//I think these are the loops that define regions.
Line Loop(12) = {1, 2, 3, 4}; //Exterior
Line Loop(13) = {8, 7, 6, 5, 9, 10, 11}; // cylinder and structure
Plane Surface(14) = {12, 13}; // plane surface on fluid region

//Plane Surface(15) = {13, 14};
//Line Loop(1) = {1, 2, 3, 4}
//Line Loop(2) = {8, 7, 6, 5, 9, 10, 11}; // Circle and structure.

// Structure
//Line Loop(2) = {24, 25, 26, 11}; // Structure boundary

// Introduce refinement at end of bar

Point(15) = {0.65,  0.2, 0, cl5};
Point{15} In Surface{14};


// Apply boundary conditions
// IMPORTANT: "FLUID" MUST contain all fluid surfaces(2D)/volumes(3D)


//"fluid"
Physical Surface(1) = {14};
//+"structure"
//Physical Surface(1) = {2};
//+ walls
Physical Line(15) = {3, 1};
//+inlet
Physical Line(16) = {4};
//+ outlet
Physical Line(17) = {2};
//+cylinder
Physical Line(18) = {8, 7, 6, 5};
//+bar end
//Physical Line(19) = {11};
//+FSI
Physical Line(20) = {11, 10, 9};
