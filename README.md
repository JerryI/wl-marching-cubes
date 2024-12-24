# Realtime Marching cubes for Wolfram Language
*Implemented in pure C and LibraryLink, aimed to be crossplatform*

**Alpha version**

![](./imgs/blob.gif)

## Interface
There is single function to calculate vertices and normals for the given 3D scalar field

```mathematica
{vertices, normals} = CMarchingCubes[scalar3D_List, isoSurface_Real, opts___]
```

where following options are supported

- `"CalculateNormals"` uses more advanced procedure to calculate normals per vertices as well. By the default is `True`. If `False` it will return `None` for normals.

## Benchmark

```mathematica
native = ImageMesh[Image3D[data3D], Method -> "MarchingCubes"] // AbsoluteTiming;

cversion = CMarchingCubes[data3D, 0.5] // AbsoluteTiming;
```

the results for 15x15x50 array on Mac Air M1 are following

```
{0.2767, ...} (* native *)
{0.0033, ...} (* cversion *)
```

## Platforms supported
- [x] OSX ARM64
- [ ] OSX x86
- [ ] Windows x86
- [ ] GNU/Linux x86

## Installation
### Option 1
Clone this repository to a folder, then load this folder using

```mathematica
PacletDirectoryLoad["path to cloned repo..."];

<<<JerryI`MarchingCubes`
```

### Option 2
Using [LPM](https://github.com/JerryI/wl-localpackages)

```mathematica
PacletRepositories[{
    Github -> "https://github.com/JerryI/wl-marching-cubes" -> "master"
}]

<<<JerryI`MarchingCubes`
```