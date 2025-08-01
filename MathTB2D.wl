(* ::Package:: *)

(* MathTB2D package *)
BeginPackage["MathTB2D`"];

MathTBVersion = "1.0";

SetModel::usage = 
  "SetModel[latticeVectors, sitePositions] validates and stores the lattice basis and site positions for subsequent operations.";

CheckLatticeStructure::usage = 
  "CheckLatticeStructure[] checks that the stored lattice structure is valid. Returns True if valid, False otherwise.";

PlotLatticeStructure::usage = 
  "PlotLatticeStructure[{n1, n2}] plots the 2D lattice repeated n1\[Times]n2 times based on the stored lattice structure.";

GetModel::usage = 
  "Returns system's parameters.";
  
AddHoppings::usage=
	"Add hopping terms";
	
ResetHoppings::usage=
	"Reset all the hoppings";

ResetOnsite::usage=
	"Reset all the onsite energy";
	
AddOnsite::usage=
	"Add onsite energy";
	
GetHK::usage=
	"Return the Hamiltonian";
	
PlotBands::usage=
	"Plot band spectrum";
	
kx::usage=
	"Crystal momentum";
	
ky::usage=
	"Crystal momentum";
  
SetModel::badBasis = "Lattice vectors `1` must be two non-collinear 2D numeric vectors.";
SetModel::collinear = "Lattice vectors `1` are collinear (determinant ~ 0).";
SetModel::badRelative = "Site positions `1` must each be a pair {a,b} with 0 \[LessEqual] a < 1 and 0 \[LessEqual] b < 1.";
SetModel::duplicateSites = "Duplicate site positions found in `1`.";
CheckLatticeStructure::noData = "No lattice structure has been set. Please call SetModel first.";

Begin["`Private`"];

(* Internal variables to store lattice data *)
(*Input form: $LatticeVectors= {v1, v2}, $SitePositions={s1, s2, ...} *)
(*All v1, v2, s1, s2,... are 2D vector represented by a list*)
ClearAll[$LatticeVectors, $SitePositions, $Orb, $Ht, $Hd, kx, ky];

(*Construct lattice*)
SetModel[latticeVectors : {_, _}, sitePositions_List, orb_List] := Module[{},
  $LatticeVectors = latticeVectors;
  $SitePositions = sitePositions;
  If[Length[$SitePositions]==Length[orb] && AllTrue[orb, (IntegerQ[#] && # > 0) &],
  ($Orb = orb;$Ht=ConstantArray[0, {Total[$Orb], Total[$Orb]}];$Hd=ConstantArray[0, {Total[$Orb], Total[$Orb]}]),
  Print["Error! The orbital input must be positive interger and match the site! "]];
];

(*Prints the lattice structure*)
GetModel[] := Module[{},
						 Print["Lattice Vectors=", $LatticeVectors];
						 Print["Sites=", $SitePositions];
						 Print["Obital number=", $Orb]];

(* Internal validation function *)
validateLattice[latticeVectors : {v1_List, v2_List}, sitePositions_List] := Module[{det},
  (* Check that v1 and v2 are 2D numeric vectors *)
  If[!(VectorQ[v1, NumericQ] && VectorQ[v2, NumericQ] && Length[v1] == 2 && Length[v2] == 2),
    Message[SetModel::badBasis, latticeVectors];
    Return[False];
  ];
  (* Check that v1 and v2 are not collinear *)
  det = v1[[1]] v2[[2]] - v1[[2]] v2[[1]];
  If[Abs[det] < 10^-8,
    Message[SetModel::collinear, latticeVectors];
    Return[False];
  ];
  (* Check that site positions are relative coords in [0,1] x [0,1] *)
  If[! And @@ (MatchQ[#, {_?NumericQ, _?NumericQ}] && 0 <= #[[1]] < 1 && 0 <= #[[2]] < 1 & /@ sitePositions),
    Message[SetModel::badRelative, sitePositions];
    Return[False];
  ];
  (* Check for duplicate sites *)
  If[Length[sitePositions] != Length[DeleteDuplicates[sitePositions]],
    Message[SetModel::duplicateSites, sitePositions];
    Return[False];
  ];
  True
];

(*Validate the structure*)
CheckLatticeStructure[] := Module[{},
  If[! ValueQ[$LatticeVectors] || ! ValueQ[$SitePositions],
    Message[CheckLatticeStructure::noData];
    Return[False];
  ];
  validateLattice[$LatticeVectors, $SitePositions]
];


(* Interactive plotting with neighbor bonds *)
PlotLatticeStructure[] := Module[{},
  (* Ensure valid data *)
  If[! CheckLatticeStructure[], Return[$Failed]];
  Manipulate[
    Module[{v1, v2, sp, realPos, pts, distList, shellDists, pairs, bondsList, bondPrims, colors, pointPrims},
      {v1, v2} = $LatticeVectors;
      sp = $SitePositions;
      (* Map relative coords to real-space positions *)
      realPos[{a_, b_}] := a v1 + b v2;
      (* Prepare fixed colors for each site index *)
      colors = Table[ColorData["BrightBands"][k/(Length[sp] + 1)], {k, Length[sp]}];
      (* Point primitives for sites *)
      pointPrims = Flatten[
        Table[
          {colors[[k]], PointSize[Large], Point[realPos[sp[[k]]] + i v1 + j v2]},
          {i, 0, n1 - 1}, {j, 0, n2 - 1}, {k, 1, Length[sp]}
        ], 2
      ];
      (* Collect all point coordinates for bond calculations *)
      pts = Flatten[Table[realPos[sp[[k]]] + i v1 + j v2, {i, 0, n1 - 1}, {j, 0, n2 - 1}, {k, 1, Length[sp]}], 2];
      (* Generate all unique point pairs *)
      pairs = Subsets[pts, {2}];
      (* Compute sorted unique distances *)
      distList = Sort[DeleteDuplicates[N[Norm /@ (Subtract @@@ pairs)]]];
      (* Take the first NN shells, allowing NN=0 *)
      shellDists = If[NN > 0, Take[distList, Min[NN, Length[distList]]], {}];
      (* For each shell, select bonds matching the shell distance *)
      bondsList = Table[Select[pairs, Abs[Norm[#[[1]] - #[[2]]] - shellDists[[k]]] < (10^-6)*Min[shellDists] &], {k, Length[shellDists]}];
      (* Line primitives: darker for lower shell index *)
      bondPrims = If[NN > 0, Flatten[Table[{GrayLevel[(k - 1)/(Length[shellDists] + 1)], Line[bondsList[[k]]]},{k, Length[shellDists]}],1],{}(* no bonds when NN==0 *)];
      (* Render bonds and points together *)
      Graphics[{bondPrims, pointPrims}, Axes -> False, AspectRatio -> Automatic, PlotRange -> All]
    ],
    {{n1, 3, "Cells in v1"}, 1, 20, 1},
    {{n2, 3, "Cells in v2"}, 1, 20, 1},
    {{NN, 1, "Neighbor bonds"}, 0, 5, 1},
    SaveDefinitions -> True
  ]
];



(* Add hoppings to the model*)
AddHoppings[t_,i_,j_,R_] := Module[{v1,v2,N,Ns,actualPos,\[Delta],hop,hopd,row,col},
								v1 = $LatticeVectors[[1]];
								v2 = $LatticeVectors[[2]];
								N = Total[$Orb];
								Ns = Length[$SitePositions];
								If[!(0<i<=Ns && 0<j<=Ns), (Print["Error! designated orbitals exceeds range."];Return[False])];
								If[!(IntegerQ[i] && IntegerQ[j]), (Print["Error! designated orbitals index must be integer!"];Return[False])];
								If[!AllTrue[R, IntegerQ], (Print["Error! Hopping R must be integer vector!"];Return[False])];
								If[i==j && Norm[R]==0,(Print["Error! To set onsite energy, use AddOnsite[] instead!"]; Return[False])];
								If[Dimensions[t][[1]]!=$Orb[[i]] || Dimensions[t][[2]]!=$Orb[[j]],(Print["Error! Wrong size of hopping matrix!"]; Return[False]) ];
								
								actualPos = Table[$SitePositions[[s, 1]]*v1 + $SitePositions[[s, 2]]*v2, {s, Ns}];
								\[Delta] = actualPos[[j]]-actualPos[[i]]+R[[1]]*v1+R[[2]]*v2;
								row = Sum[$Orb[[os]],{os,1,i}];
								col = Sum[$Orb[[os]],{os,1,j}];
								hop = Exp[I*\[Delta] . {kx,ky}]*t;
								hopd = Exp[-I*\[Delta] . {kx,ky}]*ConjugateTranspose[t];
								$Ht[[row-$Orb[[i]]+1;;row, col-$Orb[[j]]+1;;col]] += hop;
								$Ht[[col-$Orb[[j]]+1;;col, row-$Orb[[i]]+1;;row]] += hopd;
]

(*Set onsite energy*)
AddOnsite[onsite_]:= Module[{N},
					 If[!ValueQ[$Orb], (Print["Error! Orbitals and site position must be specified first!"];Return[False]), N=Total[$Orb]];
					 If[N!=Length[onsite], (Print["Error! Input must match the number of orbitals"];Return[False]) ];
					 If[AllTrue[onsite, Element[#, Reals]&], $Hd = DiagonalMatrix[onsite], (Print["Error! onsite energy must be real!"];Return[False])];
					 ]

(*Return total Hamiltonian in K space*)
GetHK[] := Module[{TotalH}, TotalH=$Ht+$Hd; Return[TotalH]];

(*Reset the hoppoings*)
ResetHoppings[] := Module[{N}, If[!ValueQ[$Orb],(Print["Error! Orbitals and site position must be specified first!"];Return[False]), N=Total[$Orb]];
						$Ht=ConstantArray[0, {N, N}];]
						
(*Reset the onsite energy*)
ResetOnsite[] := Module[{N}, If[!ValueQ[$Orb],(Print["Error! Orbitals and site position must be specified first!"];Return[False]), N=Total[$Orb]];
						$Hd=ConstantArray[0, {N, N}];]
						
						
						

(*Plot Band structure*)
PlotBands[xrange_,yrange_,Nx_,Ny_] :=  Module[{Ham,qx,qy,Ek,raw,eigen},
	If[Dimensions[xrange]!=2||xrange[[1]]>=xrange[2], (Print["Error! Invalid kx range!"];Return[False])];
	If[Dimensions[yrange]!=2||yrange[[1]]>=yrange[2], (Print["Error! Invalid ky range!"];Return[False])];
	If[!(IntegerQ[Nx]&& Nx>0), (Print["Error! Invalid Nx!"];Return[False])];
	If[!(IntegerQ[Ny]&& Ny>0), (Print["Error! Invalid Nx!"];Return[False])];
	
	
	qx = Subdivide[xrange[[1]],xrange[[2]],Nx-1]; qy = Subdivide[yrange[[1]],yrange[[2]],Ny-1];
	If[!MatrixQ[{qx}], (Print["Warning: qx is not a list:", qx];Return[False])];
	If[!MatrixQ[{qy}], (Print["Warning: qy is not a list:", qx];Return[False])];
	
	Ham = GetHK[]; eigen = Eigenvalues[Ham];
	DistributeDefinitions[Ham, eigen, qx, qy];
	raw = ParallelTable[Sort[Re[N[eigen/.{kx->qx[[i]],ky->qy[[j]]}]]],{j,1,Ny},{i,1,Nx}];
	
	Ek = Table[raw[[All,All,i]],{i,1,Length[Ham]}];
	ListPlot3D[Ek, DataRange->{{xrange[[1]],xrange[[2]]},{yrange[[1]],yrange[[2]]}},
	Mesh->None,ColorFunction->"Rainbow",ColorFunctionScaling->True,
	AxesLabel->{"kx","ky","E"},PlotRange->All,ImageSize->Large,PlotLegends->Automatic]
]


End[]; (* `Private` *)
Print["AFMVisual loaded successfully." <> " Version: " <> MathTBVersion];
Print["Developed by Junyu Tang (UCR). Licensed under MIT License."]
Print["Type ?MathTB2D`* to see all variables and available functions."];
Print["For more infos, visit: https://github.com/Junyu-Tang/MathTB2D"];
EndPackage[];
