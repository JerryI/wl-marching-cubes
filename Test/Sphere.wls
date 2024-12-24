Get[FileNameJoin[{$InputFileName // DirectoryName // ParentDirectory, "Kernel", "MarchingCubes.wl"}] ];

img = Table[Exp[-(i^2 + j^2 + k^2)], {i,-5,5,0.25}, {j, -5,5,0.25}, {k,-5,5,0.25}] // N;

{vertices, normals} = CMarchingCubes[img, 0.2];

graphics = 
  GraphicsComplex[vertices, {
    Polygon[ Partition[Range[vertices//Length], 3] ]
  }, "VertexNormals"->normals] // Graphics3D;


Export["Sphere.png", graphics]