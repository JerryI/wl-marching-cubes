#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "WolframLibrary.h"
#include "WolframNumericArrayLibrary.h"


#include <math.h>
#include <stdlib.h>

#define LINEAR_TABLE_SIZE 12
#define TRIANGULATION_TABLE_SIZE 256
#define MAX_VERTICES 3000000

#include "tables.h"

double *vertices;
double *normals;

mint triCount = 0;

mint normalBufferSize = 0;
mint vertexBufferSize = 1024*3;

// Linear interpolation
double Lerp(double start, double end, double amt) {
    return (1 - amt) * start + amt * end;
}

// Determine state index
int getState(double v0, double v1, double v2, double v3,
             double v4, double v5, double v6, double v7,
             double isovalue) {
    int state = 0;
    if (v0 >= isovalue) state |= 1;
    if (v1 >= isovalue) state |= 2;
    if (v2 >= isovalue) state |= 4;
    if (v3 >= isovalue) state |= 8;
    if (v4 >= isovalue) state |= 16;
    if (v5 >= isovalue) state |= 32;
    if (v6 >= isovalue) state |= 64;
    if (v7 >= isovalue) state |= 128;
    return state;
}

mint get3DIndex(int depth, int height, int width, mint* dims) {
    return (depth * dims[1] * dims[2]) + (height * dims[2]) + width;
}


double interpolate(double isovalue, double val1, double val2, double p1, double p2) {
    if (fabs(isovalue - val1) < 0.00001f)
        return p1;
    if (fabs(isovalue - val2) < 0.00001f)
        return p2;
    if (fabs(val1 - val2) < 0.00001f)
        return p1;
    return p1 + (isovalue - val1) * (p2 - p1) / (val2 - val1);
}


// Define cubeVertices at the top of your code
double cubeVertices[8][3] = {
    {0.0, 0.0, 0.0}, // v0
    {1.0, 0.0, 0.0}, // v1
    {1.0, 1.0, 0.0}, // v2
    {0.0, 1.0, 0.0}, // v3
    {0.0, 0.0, 1.0}, // v4
    {1.0, 0.0, 1.0}, // v5
    {1.0, 1.0, 1.0}, // v6
    {0.0, 1.0, 1.0}  // v7
};

// Edge to vertex mapping
int edgeVertices[12][2] = {
    {0, 1}, // Edge 0
    {1, 2}, // Edge 1
    {2, 3}, // Edge 2
    {3, 0}, // Edge 3
    {4, 5}, // Edge 4
    {5, 6}, // Edge 5
    {6, 7}, // Edge 6
    {7, 4}, // Edge 7
    {0, 4}, // Edge 8
    {1, 5}, // Edge 9 
    {2, 6}, // Edge 10
    {3, 7}  // Edge 11
};

void resizeBuffer(mint size) {

  double* newvertices = (double*) malloc(size * 3 * sizeof(double));
  memcpy((void*)newvertices, (void*)vertices, vertexBufferSize * 3 * sizeof(double));
  vertexBufferSize = size;
  free(vertices);
  vertices = newvertices;

  if (normalBufferSize > 0) {
    double* newnormals = (double*) malloc(size * 3 * sizeof(double));
    memcpy((void*)newnormals, (void*)normals, normalBufferSize * 3 * sizeof(double));
    normalBufferSize = size;
    free(normals);
    normals = newnormals;
  }
}

double *gradients = NULL;
mint gradientBufferSize = 0;

void resizeGradient(mint size) {
  if (gradients != NULL) free(gradients);

  gradientBufferSize = size;

  gradients = (double *)malloc(size * 3 * sizeof(double));
}

void calcGradients(double *img, mint* dims, double res) {

if (gradientBufferSize <= dims[0] * dims[1] * dims[2]) {
  resizeGradient(dims[0] * dims[1] * dims[2]);
}


for (int z = 1; z < dims[0] - 1; ++z) {
    for (int y = 1; y < dims[1] - 1; ++y) {
        for (int x = 1; x < dims[2] - 1; ++x) {
            int index = get3DIndex(z, y, x, dims);

            // Compute partial derivatives using central differences
            double dx = (img[get3DIndex(z, y, x + 1, dims)] - img[get3DIndex(z, y, x - 1, dims)]) / (2 * res);
            double dy = (img[get3DIndex(z, y + 1, x, dims)] - img[get3DIndex(z, y - 1, x, dims)]) / (2 * res);
            double dz = (img[get3DIndex(z + 1, y, x, dims)] - img[get3DIndex(z - 1, y, x, dims)]) / (2 * res);

            // Store the gradient
            gradients[index * 3 + 0] = dx;
            gradients[index * 3 + 1] = dy;
            gradients[index * 3 + 2] = dz;
        }
    }
}
}

void marchingCubesWithNormals(double *img, mint* dims, double res, double isovalue) {
    mint depth, height, width, ii;
    mint vertexCnt = 0;
    triCount = 0;

    calcGradients(img, dims, res);

    for (depth = 0; depth < dims[0] - 1; ++depth) {
        for (height = 0; height < dims[1] - 1; ++height) {
            for (width = 0; width < dims[2] - 1; ++width) {
                double x = width * res;
                double y = height * res;
                double z = depth * res;

                // Get scalar values at cube corners
                double val[8];
                
                val[0] = img[get3DIndex(depth,     height,     width,     dims)]; // v0
                val[1] = img[get3DIndex(depth,     height,     width + 1, dims)]; // v1
                val[2] = img[get3DIndex(depth,     height + 1, width + 1, dims)]; // v2
                val[3] = img[get3DIndex(depth,     height + 1, width,     dims)]; // v3
                val[4] = img[get3DIndex(depth + 1, height,     width,     dims)]; // v4
                val[5] = img[get3DIndex(depth + 1, height,     width + 1, dims)]; // v5
                val[6] = img[get3DIndex(depth + 1, height + 1, width + 1, dims)]; // v6
                val[7] = img[get3DIndex(depth + 1, height + 1, width,     dims)]; // v7

                // Compute cube state
                int state = getState(val[0], val[1], val[2], val[3], val[4], val[5], val[6], val[7], isovalue);

                // Skip if cube is entirely inside or outside
                if (state == 0 || state == 255)
                    continue;

                // Process triangles
                for (ii = 0; (triangulationTable[state][ii] != -1); ii += 3) {
                    if ((vertexCnt +3) >= vertexBufferSize) {
                      resizeBuffer((vertexCnt +3)*2);
                    }

                    int edgeIndices[3];
                    
                    edgeIndices[0] = triangulationTable[state][ii];
                    edgeIndices[1] = triangulationTable[state][ii + 1];
                    edgeIndices[2] = triangulationTable[state][ii + 2];

                    for (int e = 0; e < 3; ++e) {
                    
                        int edge = edgeIndices[e];
                        int v1 = edgeVertices[edge][0];
                        int v2 = edgeVertices[edge][1];

                        // Positions of the two vertices
                        double x1 = x + res * cubeVertices[v1][0];
                        double y1 = y + res * cubeVertices[v1][1];
                        double z1 = z + res * cubeVertices[v1][2];

                        double x2 = x + res * cubeVertices[v2][0];
                        double y2 = y + res * cubeVertices[v2][1];
                        double z2 = z + res * cubeVertices[v2][2];

                        // Scalar values at the vertices
                        double val1 = val[v1];
                        double val2 = val[v2];

                        // Interpolation factor
                        double t = (isovalue - val1) / (val2 - val1);

                        // Interpolate position
                        double vertPos[3];
                        vertPos[0] = x1 + t * (x2 - x1);
                        vertPos[1] = y1 + t * (y2 - y1);
                        vertPos[2] = z1 + t * (z2 - z1);

                        // Interpolate gradients (normals)
                        int ix1 = width + (int)cubeVertices[v1][0];
                        int iy1 = height + (int)cubeVertices[v1][1];
                        int iz1 = depth + (int)cubeVertices[v1][2];
                        int index1 = get3DIndex(iz1, iy1, ix1, dims);

                        int ix2 = width + (int)cubeVertices[v2][0];
                        int iy2 = height + (int)cubeVertices[v2][1];
                        int iz2 = depth + (int)cubeVertices[v2][2];
                        int index2 = get3DIndex(iz2, iy2, ix2, dims);

                        double grad1[3] = {gradients[index1 * 3 + 0], gradients[index1 * 3 + 1], gradients[index1 * 3 + 2]};
                        double grad2[3] = {gradients[index2 * 3 + 0], gradients[index2 * 3 + 1], gradients[index2 * 3 + 2]};

                        double normal[3];
                        normal[0] = grad1[0] + t * (grad2[0] - grad1[0]);
                        normal[1] = grad1[1] + t * (grad2[1] - grad1[1]);
                        normal[2] = grad1[2] + t * (grad2[2] - grad1[2]);

                        // Normalize the normal vector
                        double length = sqrt(normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]);
                        if (length != 0.0) {
                            normal[0] /= length;
                            normal[1] /= length;
                            normal[2] /= length;
                        }

                        // Store vertex position and normal
                        vertices[vertexCnt * 3 + 0] = vertPos[0];
                        vertices[vertexCnt * 3 + 1] = vertPos[1];
                        vertices[vertexCnt * 3 + 2] = vertPos[2];

                        normals[vertexCnt * 3 + 0] = normal[0];
                        normals[vertexCnt * 3 + 1] = normal[1];
                        normals[vertexCnt * 3 + 2] = normal[2];
                        
                        vertexCnt++;
                        
                    }
                    
                    triCount++;   
                }
            }
        }
    }

}

void marchingCubes(double *img, mint* dims, double res, double isovalue) {
    mint depth, height, width, ii;
    mint vertexCnt = 0;
    triCount = 0;

    for (depth = 0; depth < dims[0] - 1; ++depth) {
        for (height = 0; height < dims[1] - 1; ++height) {
            for (width = 0; width < dims[2] - 1; ++width) {
                double x = width * res;
                double y = height * res;
                double z = depth * res;

                // Get scalar values at cube corners
                double val[8];
                
                val[0] = img[get3DIndex(depth,     height,     width,     dims)]; // v0
                val[1] = img[get3DIndex(depth,     height,     width + 1, dims)]; // v1
                val[2] = img[get3DIndex(depth,     height + 1, width + 1, dims)]; // v2
                val[3] = img[get3DIndex(depth,     height + 1, width,     dims)]; // v3
                val[4] = img[get3DIndex(depth + 1, height,     width,     dims)]; // v4
                val[5] = img[get3DIndex(depth + 1, height,     width + 1, dims)]; // v5
                val[6] = img[get3DIndex(depth + 1, height + 1, width + 1, dims)]; // v6
                val[7] = img[get3DIndex(depth + 1, height + 1, width,     dims)]; // v7

                // Compute cube state
                int state = getState(val[0], val[1], val[2], val[3], val[4], val[5], val[6], val[7], isovalue);

                // Skip if cube is entirely inside or outside
                if (state == 0 || state == 255)
                    continue;

                // Process triangles
                for (ii = 0; (triangulationTable[state][ii] != -1); ii += 3) {
                
                    if ((vertexCnt +3) >= vertexBufferSize) {
                      resizeBuffer((vertexCnt +3)*2);
                    }

                    int edgeIndices[3];
                    
                    edgeIndices[0] = triangulationTable[state][ii];
                    edgeIndices[1] = triangulationTable[state][ii + 1];
                    edgeIndices[2] = triangulationTable[state][ii + 2];

                    for (int e = 0; e < 3; ++e) {
                    
                        int edge = edgeIndices[e];
                        int v1 = edgeVertices[edge][0];
                        int v2 = edgeVertices[edge][1];

                        // Positions of the two vertices
                        double x1 = x + res * cubeVertices[v1][0];
                        double y1 = y + res * cubeVertices[v1][1];
                        double z1 = z + res * cubeVertices[v1][2];

                        double x2 = x + res * cubeVertices[v2][0];
                        double y2 = y + res * cubeVertices[v2][1];
                        double z2 = z + res * cubeVertices[v2][2];

                        // Interpolate the position along the edge
                        double vertPos[3];
                        vertPos[0] = interpolate(isovalue, val[v1], val[v2], x1, x2);
                        vertPos[1] = interpolate(isovalue, val[v1], val[v2], y1, y2);
                        vertPos[2] = interpolate(isovalue, val[v1], val[v2], z1, z2);

                        // Store vertex
                        vertices[vertexCnt * 3 + 0] = vertPos[0];
                        vertices[vertexCnt * 3 + 1] = vertPos[1];
                        vertices[vertexCnt * 3 + 2] = vertPos[2];
                        vertexCnt++;
                        
                    }
                    
                    triCount++;   
                }
            }
        }
    }

}

DLLEXPORT mint WolframLibrary_getVersion() {
    return WolframLibraryVersion;
}

DLLEXPORT int WolframLibrary_initialize(WolframLibraryData libData) {
    vertices = (double *)malloc(vertexBufferSize * 3 * sizeof(double));

    return LIBRARY_NO_ERROR;
}

DLLEXPORT void WolframLibrary_uninitialize(WolframLibraryData libData) {
    return;
}



DLLEXPORT int process(WolframLibraryData libData, mint Argc, MArgument *Args, MArgument Res) {    

 int err = LIBRARY_NO_ERROR;

 MTensor in;
 mint in_type;
 mint in_rank;
 mint *in_dims;
 double *in_data;

 MTensor out;
 mint out_type;
 mint out_rank;
 mint out_dims[2];
 double *out_data;


 mint *tmp_dims;

 mint rows, row, cols, col;

 /* Code section */

 /* Get the input 3d array */
 in = MArgument_getMTensor(Args[0]);
 in_type = libData->MTensor_getType(in);
 in_rank = libData->MTensor_getRank(in);
 in_dims = libData->MTensor_getDimensions(in);
 in_data = libData->MTensor_getRealData(in);

 double res = MArgument_getReal(Args[2]);
 double th = MArgument_getReal(Args[1]);

 out_type = in_type;

 marchingCubes(in_data, in_dims, res, th);

 out_dims[0] = triCount*3;
 out_dims[1] = 3;

 /* Create the output array */
 err = libData->MTensor_new(out_type, 2, out_dims, &out);
 out_data = libData->MTensor_getRealData(out);

 memcpy(out_data, vertices, triCount*3*3*sizeof(double));

 MArgument_setMTensor(Res,out);
 return err;


}

mint global_type;

DLLEXPORT int getNormals(WolframLibraryData libData, mint Argc, MArgument *Args, MArgument Res) {    

 int err = LIBRARY_NO_ERROR;

 MTensor out;
 mint out_type;
 mint out_rank;
 mint out_dims[2];
 double *out_data;

 out_type = global_type;

 out_dims[0] = triCount*3;
 out_dims[1] = 3;

 /* Create the output array */
 err = libData->MTensor_new(out_type, 2, out_dims, &out);
 out_data = libData->MTensor_getRealData(out);

 memcpy(out_data, normals, triCount*3*3*sizeof(double));

 MArgument_setMTensor(Res,out);
 return err;


}

DLLEXPORT int processWithNormals(WolframLibraryData libData, mint Argc, MArgument *Args, MArgument Res) {    

 int err = LIBRARY_NO_ERROR;

 MTensor in;
 mint in_type;
 mint in_rank;
 mint *in_dims;
 double *in_data;

 MTensor out;
 mint out_type;
 mint out_rank;
 mint out_dims[2];
 double *out_data;


 mint *tmp_dims;

 mint rows, row, cols, col;

 if (normalBufferSize == 0) {
    normalBufferSize = vertexBufferSize;
    normals = (double*)malloc(vertexBufferSize * 3 * sizeof(double));
 }

 /* Code section */

 /* Get the input 3d array */
 in = MArgument_getMTensor(Args[0]);
 in_type = libData->MTensor_getType(in);
 in_rank = libData->MTensor_getRank(in);
 in_dims = libData->MTensor_getDimensions(in);
 in_data = libData->MTensor_getRealData(in);
 

 double res = MArgument_getReal(Args[2]);
 double th = MArgument_getReal(Args[1]);

 out_type = in_type;
 global_type = in_type;

 marchingCubesWithNormals(in_data, in_dims, res, th);

 out_dims[0] = triCount*3;
 out_dims[1] = 3;

 /* Create the output array */
 err = libData->MTensor_new(out_type, 2, out_dims, &out);
 out_data = libData->MTensor_getRealData(out);

 memcpy(out_data, vertices, triCount*3*3*sizeof(double));

 MArgument_setMTensor(Res,out);
 return err;


}