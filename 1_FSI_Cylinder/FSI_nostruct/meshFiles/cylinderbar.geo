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
cl1 = 0.0125; //entry
cl2 = .0025; // circle

cl3 = 0.025; //exit
cl4 = 0.0025; //structure

cl5 = .00125; // bar end refinement

//medium
//cl1 = 0.025; //entry
//cl2 = .005; // circle

//cl3 = 0.05; //exit
//cl4 = 0.005; //structure cylinder end (smaller velocity)
//cl5 = 0.0025; //structure fluid end (larger velocity)


//course
//cl1 = 0.05; //entry
//cl2 = .01; // circle

//cl3 = 0.1; //exit
//cl4 = 0.01; //structure cylinder end (smaller velocity)


// number of nodes on structure, vertical and horizontal
//N_s_v = 4;
//N_s_h = 40;

radius = 0.05;    //radius of circle

struct_l = 0.35;  //length of structure
struct_h = 0.02;  //height of structure

// Cylinder
Point(5) = {0.2,   0.2, 0, cl2};

// points on edge of cylinder and end of structure
Point(6) = {0.2 + (radius*radius-0.01*0.01)^(0.5), 0.2 + struct_h/2, 0, cl4}; //bar end top
Point(7) = {0.2 + (radius*radius-0.01*0.01)^(0.5), 0.2 - struct_h/2, 0, cl4}; //bar end lower

Circle(8) = {6, 5, 7};

// Structure use points 8 and 9
Point(11) = {0.2+radius+struct_l, 0.2-  struct_h/2, 0, cl5}; //lower right
Point(12) = {0.2+radius+struct_l, 0.2+struct_h/2, 0, cl5};  //lower left

Line(9) = {7, 11};   //top
Line(10) = {11, 12};  //end
Line(11) = {12, 6};   //lower

Line Loop(12) = {9, 10, 11, 8}; // Structure boundary
Plane Surface(13) = {12}; //

// Introduce refinement at end of bar

//Point(14) = {0.35,  0.2, 0, cl5};
//Point(14) In Surface{13};


// Apply boundary conditions
// IMPORTANT: "FLUID" MUST contain all fluid surfaces(2D)/volumes(3D)

//"fluid"
//Physical Surface(1) = {14};
//+"structure"
Physical Surface(1) = {13};
//+ walls
//Physical Line(15) = {3, 1};
//+inlet
//Physical Line(16) = {4};
//+ outlet
//Physical Line(17) = {2};
//+cylinder
//Physical Line(18) = {8, 7, 6, 5};
//+bar end
Physical Line(19) = {8};
//+FSI
Physical Line(20) = {11, 10, 9};
