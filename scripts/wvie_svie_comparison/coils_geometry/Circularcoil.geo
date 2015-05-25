
 Mesh.Format      = 1; // 1=.msh format by default
 Mesh.Algorithm   = 1; // 2D mesh algorithm (1=MeshAdapt, 5=Delaunay, 6=Frontal) Default value: 1 

 Zpos = 0.12;
 RAD = 0.04;
 W = 0.001;
 Rext = RAD + W/2;
 Rint = RAD - W/2;


 // --------------------------------------------------------------- 
 // COIL 1

 // center of coil
 Center1 = 10001;
 Point(Center1) = {0 , 0, Zpos};

 // port 1
 P1ext = 10002;
 Point(P1ext) = {Rext*Cos(1e-6), Rext*Sin(1e-6), Zpos};
 P1int = 10003;
 Point(P1int) = {Rint*Cos(1e-6), Rint*Sin(1e-6), Zpos};
 P4ext = 10004;
 Point(P4ext) = {Rext*Cos(2*Pi-1e-6), Rext*Sin(2*Pi-1e-6), Zpos};
 P4int = 10005;
 Point(P4int) = {Rint*Cos(2*Pi-1e-6), Rint*Sin(2*Pi-1e-6), Zpos};

 // port 2
 P2ext = 10006;
 Point(P2ext) = {Rext*Cos(Pi-1e-6), Rext*Sin(Pi-1e-6), Zpos};
 P2int = 10007;
 Point(P2int) = {Rint*Cos(Pi-1e-6), Rint*Sin(Pi-1e-6), Zpos};
 P3ext = 10008;
 Point(P3ext) = {Rext*Cos(Pi+1e-6), Rext*Sin(Pi+1e-6), Zpos};
 P3int = 10009;
 Point(P3int) = {Rint*Cos(Pi+1e-6), Rint*Sin(Pi+1e-6), Zpos};

 L1 = 101;
 Line(L1) = {P1int, P1ext};
 C1 = 201;
 Circle(C1) = { P1ext, Center1, P2ext};
 L2 = 102;
 Line(L2) = {P2ext, P2int};
 C2 = 202;
 Circle(C2) = { P2int, Center1, P1int}; 

 L3 = 103;
 Line(L3) = {P3ext, P3int};
 C3 = 203;
 Circle(C3) = { P3int, Center1, P4int};
 L4 = 104;
 Line(L4) = {P4int, P4ext};
 C4 = 204;
 Circle(C4) = { P4ext, Center1, P3ext}; 

 S1 = 1001;
 Line Loop(S1) = { L1, C1, L2, C2};

 S2 = 1002;
 Line Loop(S2) = { L3, C3, L4, C4};

 CoilS1 = 1003;
 Plane Surface(CoilS1) = {S1};

 CoilS2 = 1004;
 Plane Surface(CoilS2) = {S2};

 Port1 = 1;
 Physical Line(Port1) = {L1};
 Port2 = 2;
 Physical Line(Port2) = {L4};
 Port3 = 3;
 Physical Line(Port3) = {L3};
 Port4 = 4;
 Physical Line(Port4) = {L2};

 External = 0;
 Physical Line(External) = { C1, C2, C3, C4};

 PCoil1 = 1;
 Physical Surface(PCoil1) = {CoilS1, CoilS2};
