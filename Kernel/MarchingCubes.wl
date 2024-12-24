BeginPackage["JerryI`MarchingCubes`"]

CMarchingCubes::usage = "CMarchingCubes[list3D_List, isoSurface_, opts___] returns a list of vertices corresponding to the sequence of triangles (non indexed)"

Begin["`Private`"]

getLibraryLinkVersion[] := 
Which[
    $VersionNumber >= 14.1, 
        With[{n = LibraryVersionInformation[FindLibrary["demo"] ]["WolframLibraryVersion"]},
            If[!NumberQ[n], 7, n]
        ], 
    $VersionNumber >= 13.1, 
        7, 
    $VersionNumber >= 12.1, 
        6, 
    $VersionNumber >= 12.0, 
        5, 
    $VersionNumber >= 11.2, 
        4, 
    $VersionNumber >= 10.0, 
        3, 
    $VersionNumber >= 9.0, 
        2, 
    True, 
        1
]; 


$directory = DirectoryName[If[$InputFileName == "", 
        NotebookFileName[], 
        $InputFileName
    ], 2]

lib = FileNameJoin[{
        $directory, 
        "LibraryResources", 
        $SystemID <> "-v" <> ToString[getLibraryLinkVersion[] ],
        "main" <> "." <> Internal`DynamicLibraryExtension[]
}]; 



compute = LibraryFunctionLoad[
 lib,  "process", 
 {  
   {Real, 3, "Constant"}, Real, Real 
 }, {Real, 2, Automatic}]; 

computeWithNormals = LibraryFunctionLoad[
 lib,  "processWithNormals", 
 {  
   {Real, 3, "Constant"}, Real, Real 
 }, {Real, 2, Automatic}]; 

getNormals = LibraryFunctionLoad[
 lib,  "getNormals", 
 {  
   Real 
 }, {Real, 2, Automatic}];


CMarchingCubes[img_List, iso_, OptionsPattern[] ] := With[{

},
    If[OptionValue["CalculateNormals"],
        {computeWithNormals[img, iso, 1.0], getNormals[1.0]}
    ,
        {compute[img, iso, 1.0], None}
    ]
] 

Options[CMarchingCubes] = {"CalculateNormals" -> True}

End[]
EndPackage[]
