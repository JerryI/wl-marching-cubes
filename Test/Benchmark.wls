Get[FileNameJoin[{$InputFileName // DirectoryName // ParentDirectory, "Kernel", "MarchingCubes.wl"}] ];


data3D = CellularAutomaton[{224, {2, {{2, 2, 2}, {2, 1, 2}, {2, 2, 2}}}, {1, 1}}, {SeedRandom[1]; RandomInteger[1, {15, 15}], 0}, 50];

version1 = ImageMesh[Image3D[data3D], Method -> "MarchingCubes"] // AbsoluteTiming;

version2 = CMarchingCubes[data3D, 0.5] // AbsoluteTiming;

version1//First//Print;
version2//First//Print;